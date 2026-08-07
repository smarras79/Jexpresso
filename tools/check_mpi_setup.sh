#!/usr/bin/env bash
#
# check_mpi_setup.sh — find out WHY an MPI (or coupled Alya+Jexpresso) launch
# hangs, without hanging yourself.
#
# Every stage runs under a timeout, so a stage that would block forever is
# reported as HANG and the script moves on instead of freezing your terminal.
#
# Usage:
#   bash tools/check_mpi_setup.sh              # full ladder
#   NP=4 bash tools/check_mpi_setup.sh         # use 4 ranks for the parallel stages
#   TIMEOUT=120 bash tools/check_mpi_setup.sh  # allow longer per stage (cold Julia)
#
# Run it from the Jexpresso project root.
#
# Background: Jexpresso (via MPI.jl) and Alya.x are two separately linked
# binaries.  Every rank must be started by the launcher belonging to the SAME
# MPI library those binaries are linked against.  When they disagree, the ranks
# never form a world: MPI_Init either blocks forever waiting for a process
# manager that speaks a different protocol, or each rank silently comes up as
# its own world of size 1.  Both look like "the code hangs and never starts".
# See RUN-COUPLED.md §1-§3 and INSTALL.md §5.

set -u

NP="${NP:-2}"
TIMEOUT="${TIMEOUT:-90}"
JULIA="${JULIA:-julia}"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMPDIR_LOCAL="$(mktemp -d 2>/dev/null || mktemp -d -t jexpmpi)"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

if [ -t 1 ]; then
    R=$'\033[31m'; G=$'\033[32m'; Y=$'\033[33m'; B=$'\033[1m'; N=$'\033[0m'
else
    R=""; G=""; Y=""; B=""; N=""
fi

FAILURES=0
HANGS=0

hr()      { printf '%s\n' "----------------------------------------------------------------------"; }
stage()   { hr; printf '%s==> %s%s\n' "$B" "$*" "$N"; }
pass()    { printf '%s  PASS%s  %s\n' "$G" "$N" "$*"; }
fail()    { printf '%s  FAIL%s  %s\n' "$R" "$N" "$*"; FAILURES=$((FAILURES+1)); }
hang()    { printf '%s  HANG%s  %s\n' "$R" "$N" "$*"; HANGS=$((HANGS+1)); FAILURES=$((FAILURES+1)); }
warn()    { printf '%s  WARN%s  %s\n' "$Y" "$N" "$*"; }
info()    { printf '        %s\n' "$*"; }

# Run "$@" with a timeout. Output goes to $LAST_OUT, exit status in $LAST_RC
# (124 == timed out). Portable: macOS has no coreutils `timeout` by default.
LAST_OUT=""
LAST_RC=0
run_timeout() {
    local secs="$1"; shift
    local out="$TMPDIR_LOCAL/out.$$"
    : > "$out"
    "$@" > "$out" 2>&1 &
    local pid=$!
    local i=0
    while kill -0 "$pid" 2>/dev/null; do
        if [ "$i" -ge "$secs" ]; then
            # TERM first: mpirun/prterun forwards it and reaps its children.
            kill -TERM "$pid" 2>/dev/null
            sleep 3
            kill -KILL "$pid" 2>/dev/null
            wait "$pid" 2>/dev/null
            LAST_OUT="$(cat "$out")"
            LAST_RC=124
            return 124
        fi
        sleep 1
        i=$((i+1))
    done
    wait "$pid"
    LAST_RC=$?
    LAST_OUT="$(cat "$out")"
    return "$LAST_RC"
}

# Normalise an MPI banner to a family name: openmpi | mpich | unknown.
# MPICH derivatives (Intel MPI, MVAPICH, Cray MPICH) are ABI-compatible with
# MPICH and are reported as mpich.
mpi_family() {
    local s
    s="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
    case "$s" in
        *"open mpi"*|*openmpi*|*open-mpi*|*prterun*|*orterun*) echo openmpi ;;
        *mvapich*|*"intel(r) mpi"*|*"intel mpi"*|*"cray mpich"*|*mpich*) echo mpich ;;
        # Fall back to the library soname when the path carries no vendor name.
        *libmpi.40*|*libmpi_mpifh.40*|*libmpi.so.40*) echo openmpi ;;
        *libmpi.12*|*libmpifort.12*|*libmpi.so.12*)   echo mpich ;;
        *) echo unknown ;;
    esac
}

cd "$PROJECT_ROOT" || exit 1

printf '%sJexpresso MPI setup check%s   (project: %s)\n' "$B" "$N" "$PROJECT_ROOT"
printf 'ranks for parallel stages: %s     per-stage timeout: %ss\n' "$NP" "$TIMEOUT"

# ---------------------------------------------------------------------------
stage "Stage 0 — what is on PATH"
# ---------------------------------------------------------------------------
for tool in "$JULIA" mpirun mpiexec mpif90 mpicc; do
    p="$(command -v "$tool" 2>/dev/null)"
    if [ -n "$p" ]; then
        info "$(printf '%-8s %s' "$tool" "$p")"
    else
        info "$(printf '%-8s %s' "$tool" "(not found)")"
    fi
done

if ! command -v "$JULIA" >/dev/null 2>&1; then
    fail "julia not found on PATH — set JULIA=/path/to/julia and re-run."
    exit 1
fi

JULIA_BIN="$("$JULIA" -e 'print(joinpath(Sys.BINDIR, "julia"))' 2>/dev/null)"
if [ -z "$JULIA_BIN" ] || [ ! -x "$JULIA_BIN" ]; then
    warn "could not resolve the real julia binary; falling back to '$JULIA'"
    JULIA_BIN="$JULIA"
fi
info "real julia binary: $JULIA_BIN"
info "(this is the path to hand to mpirun — a juliaup shim on PATH can confuse"
info " OpenMPI 5's prterun, see RUN-COUPLED.md §6)"

LAUNCHER="$(command -v mpirun 2>/dev/null || command -v mpiexec 2>/dev/null || true)"
if [ -z "$LAUNCHER" ]; then
    fail "no mpirun/mpiexec on PATH — nothing can launch a parallel job."
    LAUNCHER_FAMILY=unknown
    LAUNCHER_BANNER=""
else
    LAUNCHER_BANNER="$("$LAUNCHER" --version 2>&1 | head -3)"
    LAUNCHER_FAMILY="$(mpi_family "$LAUNCHER_BANNER")"
    info "launcher: $LAUNCHER"
    info "launcher banner: $(printf '%s' "$LAUNCHER_BANNER" | head -1)"
    info "launcher family: $LAUNCHER_FAMILY"
fi

# Warn about several MPIs fighting over PATH.
_mpirun="$(command -v mpirun 2>/dev/null || true)"
_mpiexec="$(command -v mpiexec 2>/dev/null || true)"
_mpif90="$(command -v mpif90 2>/dev/null || true)"
if [ -n "$_mpirun" ] && [ -n "$_mpif90" ]; then
    if [ "$(dirname "$_mpirun")" != "$(dirname "$_mpif90")" ]; then
        warn "mpirun and mpif90 come from DIFFERENT installations:"
        info "  mpirun -> $_mpirun"
        info "  mpif90 -> $_mpif90"
        info "Fix PATH (or use absolute paths) before trusting anything below."
    fi
fi
if [ -n "$_mpirun" ] && [ -n "$_mpiexec" ] && \
   [ "$(dirname "$_mpirun")" != "$(dirname "$_mpiexec")" ]; then
    warn "mpirun and mpiexec come from different installations:"
    info "  mpirun  -> $_mpirun"
    info "  mpiexec -> $_mpiexec"
fi

# ---------------------------------------------------------------------------
stage "Stage 1 — what MPI.jl is bound to (serial, no launcher involved)"
# ---------------------------------------------------------------------------
if [ ! -f LocalPreferences.toml ]; then
    warn "no LocalPreferences.toml in the project root — MPI.jl is on its default"
    info "binding (the bundled JLL), which is almost never what a coupled run wants."
    info "See RUN-COUPLED.md §3 / INSTALL.md §5.2."
fi

run_timeout "$TIMEOUT" "$JULIA_BIN" --project=. -e '
    using MPIPreferences
    println("binary=", MPIPreferences.binary)
    using MPI
    impl, ver = MPI.identify_implementation()
    println("impl=", impl)
    println("version=", ver)
    println("libmpi=", MPI.API.libmpi)
    println("banner=", first(split(MPI.MPI_LIBRARY_VERSION_STRING, "\n")))
'
case "$LAST_RC" in
  0)
    printf '%s\n' "$LAST_OUT" | sed 's/^/        /'
    JL_BINARY="$(printf '%s' "$LAST_OUT" | sed -n 's/^binary=//p')"
    JL_IMPL="$(printf '%s'   "$LAST_OUT" | sed -n 's/^impl=//p')"
    JL_LIBMPI="$(printf '%s' "$LAST_OUT" | sed -n 's/^libmpi=//p')"
    JL_BANNER="$(printf '%s' "$LAST_OUT" | sed -n 's/^banner=//p')"
    JL_FAMILY="$(mpi_family "$JL_IMPL $JL_BANNER $JL_LIBMPI")"
    pass "MPI.jl loads; family = $JL_FAMILY"
    ;;
  124)
    hang "loading MPI.jl timed out after ${TIMEOUT}s."
    info "This is a SERIAL load — no launcher involved. If it hangs here the"
    info "problem is Julia-side: first-time precompilation of the project can"
    info "take many minutes. Precompile once, serially, then re-run:"
    info "    julia --project=. -e 'using Pkg; Pkg.precompile()'"
    JL_FAMILY=unknown; JL_LIBMPI=""
    ;;
  *)
    fail "MPI.jl failed to load:"
    printf '%s\n' "$LAST_OUT" | tail -20 | sed 's/^/        /'
    JL_FAMILY=unknown; JL_LIBMPI=""
    ;;
esac

# ---------------------------------------------------------------------------
stage "Stage 2 — do the launcher and MPI.jl belong to the same MPI?"
# ---------------------------------------------------------------------------
if [ "$JL_FAMILY" = unknown ] || [ "$LAUNCHER_FAMILY" = unknown ]; then
    warn "cannot compare (one side unknown): launcher=$LAUNCHER_FAMILY MPI.jl=$JL_FAMILY"
elif [ "$JL_FAMILY" = "$LAUNCHER_FAMILY" ]; then
    pass "both are $JL_FAMILY"
    info "Same family is necessary but not sufficient — two different builds of"
    info "the same family still break. Stage 4 is the real test."
else
    fail "MISMATCH: launcher is $LAUNCHER_FAMILY but MPI.jl is bound to $JL_FAMILY"
    info "This alone explains a launch that hangs forever with no output:"
    info "$LAUNCHER hands the ranks a process-manager protocol that"
    info "$JL_FAMILY does not speak, so MPI_Init never completes."
    info "Fix: rebind MPI.jl to the launcher's MPI (RUN-COUPLED.md §3, Option A)"
    info "     julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=[\"<libdir>\"])'"
    info "     julia --project=. -e 'using Pkg; Pkg.build(\"MPI\"; verbose=true)'"
    info "or launch with the launcher belonging to $JL_FAMILY."
fi

# Advice printed whenever a launch stalls. On macOS, Open MPI 5 (PRRTE) has
# two well-known ways to hang before a single line of user code runs.
prrte_hang_hints() {
    info "If you are on macOS with Open MPI 5 (PRRTE), try these in order:"
    info "  0) Stale daemons/session dirs from earlier killed runs. Every job you"
    info "     interrupted with ^C can leave a prte daemon and a session directory"
    info "     behind, and the next launch waits on them. Clear them first:"
    info "         pkill -9 prte; pkill -9 prted; pkill -9 prterun"
    info "         rm -rf /tmp/prte.* /tmp/ompi.* \"\${TMPDIR:-/tmp}\"/prte.* \"\${TMPDIR:-/tmp}\"/ompi.*"
    info "  a) Session-directory path too long. macOS TMPDIR is"
    info "     /var/folders/<...>/T/, and PMIx's rendezvous socket must fit in"
    info "     a 104-byte sun_path. When it does not, the server never binds and"
    info "     the ranks wait forever:"
    info "         export TMPDIR=/tmp"
    info "         mpirun --prtemca prte_tmpdir_base /tmp -np $NP ..."
    info "  b) Hostname does not resolve. PRRTE's out-of-band channel resolves"
    info "     the local hostname; on macOS the bare name often is not resolvable"
    info "     (VPN, mDNS off). Add it to /etc/hosts:"
    info "         echo \"127.0.0.1   \$(hostname -s)\" | sudo tee -a /etc/hosts"
    info "         echo \"127.0.0.1   \$(hostname)\"    | sudo tee -a /etc/hosts"
    info "     then: ping -c1 \$(hostname -s)"
    info "  c) Pin the loopback interface and shared-memory transport:"
    info "         mpirun --mca oob_tcp_if_include lo0 --mca btl self,sm -np $NP ..."
}

# ---------------------------------------------------------------------------
stage "Stage 2b — does the LAUNCHER work at all? (no MPI in the picture)"
# ---------------------------------------------------------------------------
# `mpirun -np N hostname` starts N processes that never call MPI_Init. If this
# hangs, the failure is entirely in the launcher/daemon layer (PRRTE, Hydra) —
# neither Julia, nor MPI.jl, nor Jexpresso is involved, and every later stage
# would hang for the same reason.
if [ -n "$LAUNCHER" ]; then
    run_timeout "$TIMEOUT" "$LAUNCHER" -np "$NP" /bin/hostname
    case "$LAST_RC" in
      0)   pass "launcher started $NP non-MPI processes"
           printf '%s\n' "$LAST_OUT" | head -"$NP" | sed 's/^/        /' ;;
      124) hang "the launcher could not even start $NP copies of /bin/hostname."
           info "Nothing above this layer can work. This is NOT a Jexpresso,"
           info "MPI.jl or Alya problem — it is the process manager itself."
           prrte_hang_hints ;;
      *)   fail "launcher exited $LAST_RC:"
           printf '%s\n' "$LAST_OUT" | tail -20 | sed 's/^/        /' ;;
    esac
fi

# ---------------------------------------------------------------------------
stage "Stage 2c — MPI_Init works from a NATIVE binary (C), 1 and $NP ranks"
# ---------------------------------------------------------------------------
# This is the test that splits "OpenMPI is broken on this machine" from
# "OpenMPI is fine, but MPI_Init misbehaves inside a Julia process".
#
# Julia ships its own copies of libraries OpenMPI also depends on (libevent,
# libhwloc, and on some builds pmix). MPI_Get_library_version — which Stage 1
# calls — touches none of them, so MPI.jl can load and report the right
# implementation while MPI_Init, which does the full PMIx handshake with the
# launcher's daemon, hangs. A C binary has no such bundled libraries, so if C
# passes here and Julia hangs in Stage 3/4, the fault is on the Julia side.
MPICC="$(command -v mpicc 2>/dev/null || true)"
if [ -z "$MPICC" ]; then
    info "no mpicc — skipping the native-binary comparison"
else
    CSRC="$TMPDIR_LOCAL/mpi_hello.c"
    CBIN="$TMPDIR_LOCAL/mpi_hello_c"
    cat > "$CSRC" <<'EOF'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("CHELLO rank=%d size=%d\n", rank, size);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
EOF
    if ! "$MPICC" "$CSRC" -o "$CBIN" >/dev/null 2>&1; then
        warn "mpicc could not build the C test — skipping"
    elif [ -n "$LAUNCHER" ]; then
        for n in 1 "$NP"; do
            run_timeout "$TIMEOUT" "$LAUNCHER" -np "$n" "$CBIN"
            case "$LAST_RC" in
              0)   # Exit 0 alone is not proof: a launcher that swallows its
                   # children's status, or ranks that came up as separate
                   # worlds, both look "successful". Require the expected
                   # number of CHELLO lines, each reporting a world of size n.
                   CN="$(printf '%s\n' "$LAST_OUT" | grep -c '^CHELLO' || true)"
                   CS="$(printf '%s\n' "$LAST_OUT" | sed -n 's/^CHELLO.*size=//p' | sort -u | tr -d '\n')"
                   if [ "$CN" = "$n" ] && [ "$CS" = "$n" ]; then
                       pass "C binary, $n rank(s): one world of size $n"
                   else
                       fail "C binary, $n rank(s): expected $n lines reporting size=$n, got $CN line(s) with size(s) '${CS:-none}'"
                       printf '%s\n' "$LAST_OUT" | tail -15 | sed 's/^/        /'
                       CSTAGE_BAD=1
                       break
                   fi ;;
              124) hang "C binary hung with $n rank(s) — the MPI install is broken here."
                   info "No Julia involved. Fix the MPI installation first."
                   CSTAGE_BAD=1
                   break ;;
              *)   fail "C binary exited $LAST_RC with $n rank(s):"
                   printf '%s\n' "$LAST_OUT" | tail -15 | sed 's/^/        /'
                   info "No Julia involved — the MPI installation itself cannot start"
                   info "a job. Stage 2d will try the known fixes."
                   CSTAGE_BAD=1
                   break ;;
            esac
        done
    fi
fi

# ---------------------------------------------------------------------------
# Stage 2d — only when the C test hung. At that point the fault is entirely in
# the MPI installation, so stop bisecting and start trying the known macOS /
# Open MPI 5 workarounds. Each is a single-rank C run under a timeout; the
# first one that completes is the fix.
# ---------------------------------------------------------------------------
if [ "${CSTAGE_BAD:-0}" = "1" ] && [ -n "${CBIN:-}" ] && [ -x "${CBIN:-/nonexistent}" ]; then
    stage "Stage 2d — trying known workarounds"

    info "environment facts first:"
    _tmpdir="${TMPDIR:-}"
    info "  TMPDIR   = ${_tmpdir:-<unset>}  (length ${#_tmpdir}; PMIx sockets must fit in 104 bytes)"
    _hshort="$(hostname -s 2>/dev/null || echo '')"
    info "  hostname = $(hostname 2>/dev/null || echo '?')  /  short: ${_hshort:-?}"
    if command -v ping >/dev/null 2>&1 && [ -n "$_hshort" ]; then
        if run_timeout 5 ping -c1 "$_hshort" >/dev/null 2>&1; then
            info "  short hostname resolves: yes"
        else
            warn "  short hostname does NOT resolve — workaround (b) below is likely your fix"
        fi
    fi
    if [ -x /usr/libexec/ApplicationFirewall/socketfilterfw ]; then
        info "  macOS firewall: $(/usr/libexec/ApplicationFirewall/socketfilterfw --getglobalstate 2>/dev/null | tail -1)"
    fi
    # libfabric is what MPICH's OFI netmod talks to. Which providers this build
    # actually has is the single most useful fact when MPI_Init fails inside
    # MPIDI_OFI_*, and it is not guessable.
    FI_INFO="$(command -v fi_info 2>/dev/null || true)"
    if [ -z "$FI_INFO" ] && command -v brew >/dev/null 2>&1; then
        _p="$(brew --prefix libfabric 2>/dev/null || true)"
        [ -n "$_p" ] && [ -x "$_p/bin/fi_info" ] && FI_INFO="$_p/bin/fi_info"
    fi
    if [ -n "$FI_INFO" ]; then
        info "  libfabric providers ($FI_INFO -l):"
        "$FI_INFO" -l 2>/dev/null | sed 's/^/            /' | head -12 || true
    else
        info "  fi_info not found — install it with 'brew install libfabric' to see"
        info "  which providers your MPICH can actually use."
    fi
    info "  FI_PROVIDER  = ${FI_PROVIDER:-<unset>}"
    info "  FI_TCP_IFACE = ${FI_TCP_IFACE:-<unset>}"

    STALE="$(pgrep -l 'prte|prted|prterun' 2>/dev/null | head -5 || true)"
    if [ -n "$STALE" ]; then
        warn "  stale launcher daemons are running right now:"
        printf '%s\n' "$STALE" | sed 's/^/            /'
        info "  Clear them (pkill -9 prte; pkill -9 prted) and re-run this script."
    fi

    # Short timeout: a working launch takes ~1 s, so there is no reason to wait
    # the full per-stage budget on each of six attempts.
    WT=$(( TIMEOUT < 25 ? TIMEOUT : 25 ))
    WORKAROUND_FOUND=""

    try_workaround() {
        local label="$1"; shift
        run_timeout "$WT" "$@"
        if [ "$LAST_RC" = 0 ] && printf '%s' "$LAST_OUT" | grep -q '^CHELLO'; then
            pass "$label  →  WORKS"
            [ -z "$WORKAROUND_FOUND" ] && WORKAROUND_FOUND="$label"
            return 0
        elif [ "$LAST_RC" = 124 ]; then
            info "  $label  →  still hangs"
        else
            info "  $label  →  exit $LAST_RC"
        fi
        return 1
    }

    # Singleton: no launcher at all. Isolates "MPI is broken" from
    # "the launcher/PMIx handshake is broken".
    try_workaround "singleton (no mpirun)" env -u PMIX_NAMESPACE "$CBIN" || true

    try_workaround "TMPDIR=/tmp" env TMPDIR=/tmp "$LAUNCHER" -np 1 "$CBIN" || true

    try_workaround "--prtemca prte_tmpdir_base /tmp" \
        "$LAUNCHER" --prtemca prte_tmpdir_base /tmp -np 1 "$CBIN" || true

    try_workaround "--mca oob_tcp_if_include lo0" \
        "$LAUNCHER" --mca oob_tcp_if_include lo0 -np 1 "$CBIN" || true

    try_workaround "TMPDIR=/tmp + oob lo0" \
        env TMPDIR=/tmp "$LAUNCHER" --mca oob_tcp_if_include lo0 -np 1 "$CBIN" || true

    try_workaround "--mca btl self,sm" \
        "$LAUNCHER" --mca btl self,sm -np 1 "$CBIN" || true

    # MPICH-flavoured: its OFI (libfabric) netmod picks a provider and a NIC at
    # MPI_Init, and there are three distinct ways that goes wrong:
    #
    #   nic=utun7 ...          it chose a VPN tunnel      -> pin the interface
    #   nic=(n/a) getinfo      no provider matched        -> the filter is wrong
    #                                                        (e.g. shm, which is
    #                                                        Linux-only)
    #   nic=tcp ep_enable      provider chosen, endpoint  -> the pinned interface
    #     "Bad file descriptor"  would not open              is itself unusable
    #
    # The last one is why FI_TCP_IFACE is NOT assumed to be an improvement:
    # pinning to lo0 fixes the first failure and can cause the third. So sweep
    # the combinations rather than betting on one.
    try_workaround "FI_PROVIDER=tcp (no interface pinned)" \
        env -u FI_TCP_IFACE FI_PROVIDER=tcp "$LAUNCHER" -np 1 "$CBIN" || true

    try_workaround "FI_PROVIDER=tcp + FI_TCP_IFACE=lo0" \
        env FI_PROVIDER=tcp FI_TCP_IFACE=lo0 "$LAUNCHER" -np 1 "$CBIN" || true

    # Whatever the machine's primary interface actually is.
    PRIMARY_IF=""
    if command -v route >/dev/null 2>&1; then
        PRIMARY_IF="$(route get default 2>/dev/null | awk '/interface:/{print $2}' | head -1)"
    fi
    [ -z "$PRIMARY_IF" ] && command -v ip >/dev/null 2>&1 && \
        PRIMARY_IF="$(ip route show default 2>/dev/null | awk '/dev/{for(i=1;i<=NF;i++) if($i=="dev") print $(i+1)}' | head -1)"
    if [ -n "$PRIMARY_IF" ]; then
        try_workaround "FI_PROVIDER=tcp + FI_TCP_IFACE=$PRIMARY_IF" \
            env FI_PROVIDER=tcp FI_TCP_IFACE="$PRIMARY_IF" "$LAUNCHER" -np 1 "$CBIN" || true
    fi

    try_workaround "FI_PROVIDER=sockets" \
        env -u FI_TCP_IFACE FI_PROVIDER=sockets "$LAUNCHER" -np 1 "$CBIN" || true

    try_workaround "FI_PROVIDER=udp" \
        env -u FI_TCP_IFACE FI_PROVIDER=udp "$LAUNCHER" -np 1 "$CBIN" || true

    if [ "$(uname -s)" = "Linux" ]; then
        try_workaround "FI_PROVIDER=shm (Linux only)" \
            env FI_PROVIDER=shm "$LAUNCHER" -np 1 "$CBIN" || true
    fi

    # A stale FI_PROVIDER already exported in the user's shell is itself a
    # common cause of the nic=(n/a) failure — check whether simply clearing it
    # fixes things.
    if [ -n "${FI_PROVIDER:-}" ]; then
        warn "FI_PROVIDER is already set in this shell to '${FI_PROVIDER}'"
        try_workaround "unset FI_PROVIDER" \
            env -u FI_PROVIDER -u FI_TCP_IFACE "$LAUNCHER" -np 1 "$CBIN" || true
    fi

    hr
    if [ -n "$WORKAROUND_FOUND" ]; then
        pass "WORKAROUND FOUND: $WORKAROUND_FOUND"
        case "$WORKAROUND_FOUND" in
          "singleton"*)
            info "MPI itself works; only the LAUNCHER path is broken. Reinstalling"
            info "Open MPI so its PRRTE/PMIx are rebuilt against the current"
            info "hwloc/libevent/pmix usually fixes it:"
            info "    brew reinstall open-mpi"
            info "If it persists, switch the whole stack to MPICH — RUN-COUPLED.md §4"
            info "already recommends MPICH on macOS:"
            info "    brew install mpich"
            info "    julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=[\"/opt/homebrew/lib\"])'"
            info "    julia --project=. -e 'using Pkg; Pkg.build(\"MPI\"; verbose=true)'"
            info "    cd AlyaProxy && MPIF90=/opt/homebrew/bin/mpif90 bash compilef90.sh && cd .."
            ;;
          TMPDIR*|--prtemca*)
            info "Make it permanent by adding to ~/.zshrc:"
            info "    export TMPDIR=/tmp"
            info "Then re-run this script to confirm every stage passes."
            ;;
          FI_PROVIDER*|"unset FI_PROVIDER")
            info "This was MPICH's libfabric layer. Export exactly the setting"
            info "named above in ~/.zshrc — and ONLY that one; the combinations"
            info "interact, so do not add FI_TCP_IFACE unless it is in the winner."
            info "Reference: shm is Linux-only (asking for it gives nic=(n/a));"
            info "pinning FI_TCP_IFACE avoids a VPN tunnel being chosen but can"
            info "itself fail with 'ep_enable ... Bad file descriptor' if that"
            info "interface cannot host an endpoint."
            ;;
          "unset FI_PROVIDER"*)
            info "The FI_PROVIDER already exported in your shell was the problem."
            info "Remove it from ~/.zshrc (and FI_TCP_IFACE with it)."
            ;;
          *)
            info "Add the same flags to every mpirun you launch, or set the"
            info "equivalent OMPI_MCA_* / PRTE_MCA_* variables in ~/.zshrc."
            ;;
        esac
    else
        fail "none of the workarounds helped — the Open MPI install is broken."
        info "Rebuild it against the current dependencies:"
        info "    brew reinstall open-mpi"
        info "and if that still hangs, move to MPICH (RUN-COUPLED.md §4 recommends"
        info "it on macOS anyway). Both Jexpresso and Alya.x must then be rebound"
        info "and rebuilt against MPICH — see RUN-COUPLED.md §3 and §5."
        info "To capture a stack for an upstream bug report, in another terminal"
        info "while a run is hung:"
        info "    sudo sample \$(pgrep -n mpi_hello_c) 3 -f /tmp/ompi.sample"
    fi
    prrte_hang_hints
fi

# ---------------------------------------------------------------------------
stage "Stage 3 — MPI_Init works with 1 rank"
# ---------------------------------------------------------------------------
HELLO="$TMPDIR_LOCAL/mpi_hello.jl"
cat > "$HELLO" <<'EOF'
using MPI
MPI.Init()
println("HELLO rank=", MPI.Comm_rank(MPI.COMM_WORLD),
        " size=", MPI.Comm_size(MPI.COMM_WORLD))
flush(stdout)
MPI.Barrier(MPI.COMM_WORLD)
MPI.Finalize()
EOF

if [ -n "$LAUNCHER" ]; then
    run_timeout "$TIMEOUT" "$LAUNCHER" -np 1 "$JULIA_BIN" --project=. "$HELLO"
    case "$LAST_RC" in
      0)   pass "1 rank OK"; printf '%s\n' "$LAST_OUT" | grep '^HELLO' | sed 's/^/        /' ;;
      124) hang "even 1 rank never finished MPI_Init (${TIMEOUT}s)."
           STAGE3_HUNG=1
           info "A single rank needs no rank-to-rank connection at all — it only"
           info "has to complete the PMIx handshake with the launcher's daemon."
           info "Stage 3b below will show exactly where it is stuck."
           prrte_hang_hints ;;
      *)   fail "1 rank exited $LAST_RC:"; printf '%s\n' "$LAST_OUT" | tail -20 | sed 's/^/        /' ;;
    esac
fi

# ---------------------------------------------------------------------------
# Stage 3b — only runs when Stage 3 hung. Re-launches the stuck case, lets it
# wedge, then samples it: a backtrace says which library it is blocked in, and
# the loaded-image list says WHICH COPY of that library it got. Julia bundles
# its own libevent / libhwloc / libz, and when MPI_Init picks one of those up
# instead of the ones OpenMPI was built against, it blocks in the PMIx
# handshake — exactly the "1 rank hangs" signature.
# ---------------------------------------------------------------------------
if [ "${STAGE3_HUNG:-0}" = "1" ] && [ -n "$LAUNCHER" ]; then
    stage "Stage 3b — where is it stuck?"
    SETTLE=25
    info "re-launching and letting it wedge for ${SETTLE}s …"
    "$LAUNCHER" -np 1 "$JULIA_BIN" --project=. "$HELLO" >/dev/null 2>&1 &
    WATCH_PID=$!
    sleep "$SETTLE"

    JPID="$(pgrep -f "$HELLO" 2>/dev/null | head -1 || true)"
    if [ -z "$JPID" ]; then
        warn "could not find the stuck julia process (it may have exited)"
    else
        info "stuck pid = $JPID"

        STACK="$TMPDIR_LOCAL/stack.txt"
        : > "$STACK"
        if command -v sample >/dev/null 2>&1; then                 # macOS
            sample "$JPID" 3 -f "$STACK" >/dev/null 2>&1 || true
        elif command -v eu-stack >/dev/null 2>&1; then             # Linux
            eu-stack -p "$JPID" > "$STACK" 2>&1 || true
        elif command -v gdb >/dev/null 2>&1; then
            gdb -p "$JPID" -batch -ex 'thread apply all bt' > "$STACK" 2>&1 || true
        fi

        if [ -s "$STACK" ]; then
            info "frames mentioning MPI / PMIx / event / hwloc:"
            grep -iE 'pmix|opal|ompi|libmpi|open-pal|event_base|hwloc|orte|prte' "$STACK" \
                | head -25 | sed 's/^/          /' || info "  (none — see $STACK)"
            cp "$STACK" ./jexpresso_mpi_stack.txt 2>/dev/null && \
                info "full sample written to ./jexpresso_mpi_stack.txt"
        else
            warn "could not sample the process."
            info "macOS may require: sudo sample $JPID 3"
            info "Run that by hand in another terminal while the job is hung."
        fi

        # Which copies of the shared libraries did it actually load? A path
        # under the Julia bundle for libevent/libhwloc/libpmix is the smoking
        # gun for a bundled-library clash.
        if command -v lsof >/dev/null 2>&1; then
            LIBS="$(lsof -p "$JPID" 2>/dev/null | grep -iE 'pmix|libevent|hwloc|open-pal|libmpi' | awk '{print $NF}' | sort -u)"
            if [ -n "$LIBS" ]; then
                info "MPI-related libraries loaded by the stuck process:"
                printf '%s\n' "$LIBS" | sed 's/^/          /'
                if printf '%s' "$LIBS" | grep -qiE 'julia.*(libevent|hwloc|pmix)'; then
                    fail "at least one of these came from the JULIA bundle, not from OpenMPI"
                    info "That is the hang: OpenMPI's PMIx client is running against"
                    info "Julia's bundled support library. Work around it by making the"
                    info "OpenMPI copies win at load time:"
                    info "    DYLD_INSERT_LIBRARIES=/opt/homebrew/lib/libpmix.dylib mpirun -np 1 ..."
                    info "or switch MPI.jl to the MPICH_jll route for standalone runs"
                    info "(INSTALL.md §5, Route C) — note that route is not usable for"
                    info "coupling unless Alya.x is rebuilt against the same MPICH."
                fi
            fi
        fi
    fi
    kill -TERM "$WATCH_PID" 2>/dev/null || true
    sleep 2
    kill -KILL "$WATCH_PID" 2>/dev/null || true
    wait "$WATCH_PID" 2>/dev/null || true
fi

# ---------------------------------------------------------------------------
stage "Stage 4 — MPI_Init works with $NP ranks  (THE decisive test)"
# ---------------------------------------------------------------------------
if [ -n "$LAUNCHER" ]; then
    run_timeout "$TIMEOUT" "$LAUNCHER" -np "$NP" "$JULIA_BIN" --project=. "$HELLO"
    NLINES="$(printf '%s\n' "$LAST_OUT" | grep -c '^HELLO' || true)"
    SIZES="$(printf '%s\n' "$LAST_OUT" | sed -n 's/^HELLO.*size=//p' | sort -u | tr '\n' ' ')"
    case "$LAST_RC" in
      0)
        printf '%s\n' "$LAST_OUT" | grep '^HELLO' | sed 's/^/        /'
        if [ "$NLINES" -eq "$NP" ] && [ "$(printf '%s' "$SIZES" | tr -d ' ')" = "$NP" ]; then
            pass "$NP ranks formed one world of size $NP"
            info "MPI itself is healthy. If your run still hangs, the cause is"
            info "inside Jexpresso, not the MPI setup — see 'Next steps' below."
        else
            fail "ranks did NOT form one world (saw $NLINES lines, sizes: $SIZES)"
            info "Each rank came up as its own world — the classic symptom of a"
            info "binary linked against one MPI being started by another MPI's"
            info "launcher. Jexpresso will then block on the first collective."
            info "Fix as described in Stage 2."
        fi
        ;;
      124)
        hang "$NP ranks never completed MPI_Init (${TIMEOUT}s)."
        info "THIS is your hang, and it is entirely below Jexpresso: the test"
        info "program is 6 lines of MPI.jl. Nothing in Jexpresso runs yet."
        if [ "${CSTAGE_BAD:-0}" = "1" ]; then
            info "Expected: Stage 2c already showed a plain C program hanging the"
            info "same way. This is the SAME failure, not a second one — nothing"
            info "here is Julia's or Jexpresso's doing. Act on Stage 2d."
        else
            info "Causes, in order of likelihood:"
            info "  1. launcher / MPI.jl mismatch (see Stage 2)"
            info "  2. the launcher daemon layer itself (see Stage 2b)"
            info "  3. hostname resolution or a firewall blocking the oob channel"
            prrte_hang_hints
            info "MPICH-specific: export MPICH_INTERFACE_HOSTNAME=127.0.0.1 (INSTALL.md §5.5)"
        fi
        ;;
      *)
        fail "$NP ranks exited $LAST_RC:"
        printf '%s\n' "$LAST_OUT" | tail -25 | sed 's/^/        /'
        ;;
    esac
fi

# ---------------------------------------------------------------------------
stage "Stage 5 — Alya.x (only relevant for coupled runs)"
# ---------------------------------------------------------------------------
if [ ! -x AlyaProxy/Alya.x ]; then
    info "AlyaProxy/Alya.x not built — skipping (build it with"
    info "  cd AlyaProxy && bash compilef90.sh)"
else
    if command -v otool >/dev/null 2>&1; then
        LINKED="$(otool -L AlyaProxy/Alya.x 2>/dev/null | grep -i mpi || true)"
    elif command -v ldd >/dev/null 2>&1; then
        LINKED="$(ldd AlyaProxy/Alya.x 2>/dev/null | grep -i mpi || true)"
    else
        LINKED=""
    fi
    if [ -n "$LINKED" ]; then
        info "Alya.x links against:"
        printf '%s\n' "$LINKED" | sed 's/^/          /'
        AX_FAMILY="$(mpi_family "$LINKED")"
        if [ -n "${JL_FAMILY:-}" ] && [ "$JL_FAMILY" != unknown ] && [ "$AX_FAMILY" != unknown ]; then
            if [ "$AX_FAMILY" = "$JL_FAMILY" ]; then
                pass "Alya.x and MPI.jl are both $AX_FAMILY"
                info "MPI.jl libmpi: ${JL_LIBMPI:-unknown}"
                info "Compare the paths above by eye — same family but different"
                info "installs still deadlocks. They must be the SAME library file."
            else
                fail "Alya.x is $AX_FAMILY but MPI.jl is $JL_FAMILY — coupled runs cannot work."
                info "Rebuild Alya.x with the matching mpif90 (RUN-COUPLED.md §5):"
                info "  cd AlyaProxy && MPIF90=<wrapper-of-$JL_FAMILY> bash compilef90.sh"
                FAMILY_MISMATCH=1
            fi
        fi
    else
        warn "could not read Alya.x's dynamic libraries (no otool/ldd, or static link)"
    fi

    if [ "${FAMILY_MISMATCH:-0}" = "1" ]; then
        info "Skipping the Alya run tests — Alya.x and the launcher are different"
        info "MPIs, so anything they produce is meaningless. Alya.x would come up"
        info "as N independent worlds of size 1 and still 'succeed'. Rebuild it"
        info "first, then re-run this script."
    elif [ -n "$LAUNCHER" ]; then
        # Alya alone: with no Julia ranks in the world it exchanges zero points
        # (npoin_recv is all-zero, so it posts no receives) and runs its whole
        # time loop to completion. That makes it a clean pass/fail test of the
        # SAME launcher driving a Fortran MPI binary — which is what separates
        # "MPI is broken here" from "MPI.jl / Julia is broken here".
        info "Alya alone with $NP ranks (no Julia — should run to completion)…"
        run_timeout "$TIMEOUT" "$LAUNCHER" -np "$NP" ./AlyaProxy/Alya.x
        case "$LAST_RC" in
          0)   # Exiting 0 is not enough: ranks that never joined a shared world
               # each run happily as their own world of size 1. Alya's rank 0
               # prints "=== Coupling labels (world size= N )", so a correct run
               # emits exactly ONE such line reporting N == NP.
               AW="$(printf '%s\n' "$LAST_OUT" | sed -n 's/.*world size= *\([0-9][0-9]*\).*/\1/p')"
               AW_N="$(printf '%s\n' "$AW" | grep -c '[0-9]' || true)"
               if [ "$AW_N" = 1 ] && [ "$(printf '%s' "$AW" | tr -d ' \n')" = "$NP" ]; then
                   pass "Alya alone completed in one world of size $NP"
                   info "If Stage 4 (Julia, $NP ranks) failed or hung while this passed,"
                   info "the fault is on the Julia side: MPI.jl's binding, or the Julia"
                   info "binary's architecture vs the MPI library's."
               else
                   fail "Alya alone exited 0 but did NOT form one world of size $NP"
                   info "world sizes reported: $(printf '%s' "$AW" | tr '\n' ' ')"
                   info "Each rank came up on its own — the binary and the launcher"
                   info "are not the same MPI. Rebuild Alya.x (RUN-COUPLED.md §5)."
               fi
               ;;
          124) hang "Alya alone never completed either."
               info "Both languages hang through this launcher, so the problem is"
               info "the MPI installation / launcher, not Jexpresso."
               [ "${CSTAGE_BAD:-0}" = "1" ] || prrte_hang_hints ;;
          *)   warn "Alya alone exited $LAST_RC:"
               printf '%s\n' "$LAST_OUT" | tail -15 | sed 's/^/        /' ;;
        esac
        # Clean up the .vts files that run just produced.
        rm -f alya_grid_*.vts 2>/dev/null || true

        info "MPMD smoke test: 1 Alya rank + 1 Julia rank sharing one world…"
        run_timeout "$TIMEOUT" "$LAUNCHER" -np 1 ./AlyaProxy/Alya.x \
            : -np 1 "$JULIA_BIN" --project=. "$HELLO"
        case "$LAST_RC" in
          124) hang "MPMD launch never got both sides through MPI_Init."
               if [ "${CSTAGE_BAD:-0}" = "1" ]; then
                   info "Expected — every launch hangs on this machine. Fix the MPI"
                   info "installation (Stage 2d) before reading anything into this."
               else
                   info "If Stage 4 passed, the Julia side is fine and Alya.x is linked"
                   info "against a different MPI — recompile it (RUN-COUPLED.md §5)."
               fi ;;
          *)   if printf '%s\n' "$LAST_OUT" | grep -q 'size=2'; then
                   pass "MPMD world formed correctly (size=2)"
               else
                   warn "MPMD launch returned $LAST_RC without a size=2 world:"
                   printf '%s\n' "$LAST_OUT" | tail -20 | sed 's/^/        /'
               fi ;;
        esac
    fi
fi

# ---------------------------------------------------------------------------
stage "Stage 6 — does JEXPRESSO_COUPLED survive the launch?"
# ---------------------------------------------------------------------------
# Jexpresso must see JEXPRESSO_COUPLED=1 or it takes its standalone path, runs
# its collectives on the full MPI_COMM_WORLD (Alya's ranks included) and
# deadlocks with no diagnostic at all.
#
# The documented way is to `export` it for the whole job, which works with any
# launcher — so that is what gets tested first. The per-app-context flag is
# launcher-specific (OpenMPI `-x VAR=VAL`, MPICH/Hydra `-env VAR VAL`) and is
# only probed as a convenience.
if [ -n "$LAUNCHER" ]; then
    ENVTEST="$TMPDIR_LOCAL/envtest.jl"
    cat > "$ENVTEST" <<'EOF'
println("COUPLED_SEEN=", get(ENV, "JEXPRESSO_COUPLED", "<unset>"))
EOF
    run_timeout "$TIMEOUT" env JEXPRESSO_COUPLED=1 "$LAUNCHER" \
        -np 1 /bin/echo alya-placeholder \
        : -np 1 "$JULIA_BIN" "$ENVTEST"
    if printf '%s\n' "$LAST_OUT" | grep -q 'COUPLED_SEEN=1'; then
        pass "exported JEXPRESSO_COUPLED reaches the Julia ranks (the documented way)"
    else
        fail "an EXPORTED JEXPRESSO_COUPLED did not reach the Julia ranks"
        printf '%s\n' "$LAST_OUT" | tail -10 | sed 's/^/        /'
        info "Your launcher is not forwarding the environment. Check its docs for"
        info "an explicit export flag (OpenMPI: -x VAR ; MPICH: -genv VAR VAL)."
    fi

    # Per-app-context flag, in the syntax this launcher actually speaks.
    case "$LAUNCHER_FAMILY" in
      mpich) ENVFLAG_DESC="-env JEXPRESSO_COUPLED 1 (MPICH/Hydra)"
             set -- -np 1 /bin/echo alya-placeholder : -np 1 -env JEXPRESSO_COUPLED 1 "$JULIA_BIN" "$ENVTEST" ;;
      *)     ENVFLAG_DESC="-x JEXPRESSO_COUPLED=1 (OpenMPI)"
             set -- -np 1 /bin/echo alya-placeholder : -np 1 -x JEXPRESSO_COUPLED=1 "$JULIA_BIN" "$ENVTEST" ;;
    esac
    run_timeout "$TIMEOUT" "$LAUNCHER" "$@"
    if printf '%s\n' "$LAST_OUT" | grep -q 'COUPLED_SEEN=1'; then
        pass "$ENVFLAG_DESC also works in the second app context"
    else
        warn "$ENVFLAG_DESC did NOT reach the Julia ranks"
        info "Harmless — RUN-COUPLED.md §6 tells you to export the variable"
        info "instead, precisely because this flag is not portable:"
        info "    export JEXPRESSO_COUPLED=1"
        info "    mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 \"\$JULIA_BIN\" --project=. ./src/Jexpresso.jl CompEuler 3dAlya"
    fi
fi

# ---------------------------------------------------------------------------
hr
if [ "$FAILURES" -eq 0 ]; then
    printf '%sAll stages passed.%s\n' "$G" "$N"
    printf 'MPI is healthy, so a hang in a real run is above this layer.\n'
    printf 'Next steps:\n'
    printf '  * Precompile serially first (a cold project can spend many minutes\n'
    printf '    precompiling, and under mpirun that output is buffered and invisible):\n'
    printf "      julia --project=. -e 'using Pkg; Pkg.precompile()'\n"
    printf '  * Re-launch with per-rank tagged, unbuffered output so you can see how\n'
    printf '    far it gets — OpenMPI: --output tag   MPICH: -prepend-rank\n'
    printf '  * Confirm the case runs on ONE rank first:\n'
    printf "      %s --project=. ./src/Jexpresso.jl CompEuler 3dAlya\n" "$JULIA_BIN"
else
    printf '%s%d stage(s) failed (%d hung).%s\n' "$R" "$FAILURES" "$HANGS" "$N"
    printf 'Fix the earliest failure first — later stages depend on it.\n'
fi
hr
exit 0

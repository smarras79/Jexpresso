#!/bin/bash
#==============================================================================
# compilef90.sh — build the Alya proxy (Alya.x) against the right MPI.
#
# Alya.x is built from myAlya.f90 (a symlink to alya_all2all_time_loop.f90).
#
# THE ONE RULE
# ------------
# A coupled run is a single MPMD job sharing one MPI_COMM_WORLD, so Alya.x and
# Jexpresso must use the SAME MPI installation — same implementation, same
# version, same ABI — and be launched by THAT installation's mpirun/mpiexec.
# Mismatches do not produce a clean error: the job hangs in MPI_Init, or every
# rank comes up in its own world of size 1. See ../RUN-COUPLED.md §1-§3.
#
# So by default this script does not guess. It asks Julia what MPI.jl is bound
# to and builds against that same implementation.
#
# USAGE
# -----
#   bash compilef90.sh                # match whatever MPI.jl is bound to
#   bash compilef90.sh mpich          # force MPICH
#   bash compilef90.sh openmpi        # force OpenMPI
#   bash compilef90.sh --help
#
# ENVIRONMENT (all optional)
# --------------------------
#   MPIF90=/path/to/mpif90   Use this exact wrapper. Highest precedence;
#                            skips all detection.
#   MPI=openmpi|mpich|auto   Same as the positional argument. Default: auto.
#   JULIA=/path/to/julia     Julia used to query MPI.jl's binding (default:
#                            julia). Only used when MPI=auto.
#   SKIP_JULIA_CHECK=1       Do not consult Julia at all, and do not verify the
#                            finished binary against it. Use when building on a
#                            machine that has no Jexpresso environment.
#   USE_MPIF08=0|1           Force the Fortran binding instead of probing:
#                              1 -> `use mpi_f08`      (typed handles)
#                              0 -> `include 'mpif.h'` (integer handles)
#   FFLAGS_EXTRA="..."       Extra flags appended to the compile line.
#
# WHAT IT HANDLES FOR YOU
# -----------------------
#   * Finds the wrapper belonging to the chosen implementation by PREFIX, not
#     by PATH — Homebrew can only link one MPI at a time, so `which mpif90`
#     may well be the other one's.
#   * Verifies the wrapper really is the implementation you asked for.
#   * Probes for the mpi_f08 module and falls back to mpif.h when a build does
#     not ship one (common with MPICH).
#   * On the mpif.h path, adds -fallow-argument-mismatch, without which
#     gfortran >= 10 rejects legal MPI usage.
#   * Compiles in a scratch directory so the vendored Open MPI mpif*.h in this
#     folder cannot shadow your MPI's own headers.
#   * Checks the finished binary links the same MPI family as MPI.jl.
#==============================================================================

set -e

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SRC_DIR/.." && pwd)"
OUT="$SRC_DIR/Alya.x"
JULIA="${JULIA:-julia}"

usage() { sed -n '2,/^#===/p' "${BASH_SOURCE[0]}" | sed 's/^#//;s/^ //'; exit 0; }

MPI_WANT="${MPI:-auto}"
case "${1:-}" in
    -h|--help|help) usage ;;
    "")             ;;
    openmpi|open-mpi|ompi) MPI_WANT=openmpi ;;
    mpich)                 MPI_WANT=mpich ;;
    auto)                  MPI_WANT=auto ;;
    *) echo "ERROR: unknown argument '$1'. Try: bash compilef90.sh --help" >&2
       exit 1 ;;
esac

# Classify an arbitrary string (a -show line, a library path, a banner) as
# openmpi | mpich | unknown. MPICH derivatives are ABI-compatible with MPICH.
mpi_family() {
    local s
    s="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
    case "$s" in
        *"open mpi"*|*openmpi*|*open-mpi*|*ompi*)                 echo openmpi ;;
        *mvapich*|*"intel"*|*"cray mpich"*|*mpich*)               echo mpich ;;
        *libmpi.40*|*libmpi_mpifh.40*|*libmpi.so.40*)             echo openmpi ;;
        *libmpi.12*|*libmpifort.12*|*libmpi.so.12*)               echo mpich ;;
        *)                                                        echo unknown ;;
    esac
}

# Print candidate bin directories for an implementation, most specific first.
candidate_bins() {
    local impl="$1" formula p
    case "$impl" in
        openmpi) formula=open-mpi ;;
        mpich)   formula=mpich ;;
        *)       return 0 ;;
    esac
    if command -v brew >/dev/null 2>&1; then
        p="$(brew --prefix "$formula" 2>/dev/null || true)"
        [ -n "$p" ] && echo "$p/bin"
    fi
    for p in "/usr/lib/$(uname -m)-linux-gnu/$impl/bin" \
             "/usr/lib64/$impl/bin" \
             "/usr/lib/$impl/bin" \
             "/opt/$impl/bin" \
             "/usr/local/$impl/bin"; do
        [ -d "$p" ] && echo "$p"
    done
    # MPI_HOME is what environment-module systems set.
    [ -n "${MPI_HOME:-}" ] && [ -d "$MPI_HOME/bin" ] && echo "$MPI_HOME/bin"
    return 0
}

# Find the mpif90 belonging to $1 (openmpi|mpich), verifying it really is that
# implementation before accepting it.
find_wrapper() {
    local impl="$1" bin cand fam
    while read -r bin; do
        [ -n "$bin" ] || continue
        cand="$bin/mpif90"
        [ -x "$cand" ] || cand="$bin/mpifort"
        [ -x "$cand" ] || continue
        fam="$(mpi_family "$("$cand" -show 2>/dev/null || true) $cand")"
        if [ "$fam" = "$impl" ]; then echo "$cand"; return 0; fi
    done <<< "$(candidate_bins "$impl")"

    # Last resort: whatever is on PATH, but only if it IS the right family.
    cand="$(command -v mpif90 2>/dev/null || command -v mpifort 2>/dev/null || true)"
    if [ -n "$cand" ]; then
        fam="$(mpi_family "$("$cand" -show 2>/dev/null || true) $cand")"
        [ "$fam" = "$impl" ] && { echo "$cand"; return 0; }
    fi
    return 1
}

#------------------------------------------------------------------------------
# Step 1 — decide which MPI to build against.
#------------------------------------------------------------------------------
JL_FAMILY=unknown
JL_LIBMPI=""
if [ "${SKIP_JULIA_CHECK:-0}" != "1" ] && command -v "$JULIA" >/dev/null 2>&1; then
    printf '==> Asking MPI.jl what Jexpresso is bound to ... '
    JL_OUT="$("$JULIA" --project="$PROJECT_DIR" -e '
        using MPI
        impl, _ = MPI.identify_implementation()
        println(impl); println(MPI.API.libmpi)' 2>/dev/null || true)"
    if [ -n "$JL_OUT" ]; then
        JL_IMPL="$(printf '%s\n' "$JL_OUT" | sed -n 1p)"
        JL_LIBMPI="$(printf '%s\n' "$JL_OUT" | sed -n 2p)"
        JL_FAMILY="$(mpi_family "$JL_IMPL $JL_LIBMPI")"
        echo "$JL_IMPL  ($JL_LIBMPI)"
    else
        echo "could not ask"
        echo "    (Julia present but MPI.jl would not load. Continuing without"
        echo "     the cross-check — pass 'openmpi' or 'mpich' to be explicit.)"
    fi
fi

if [ -n "${MPIF90:-}" ]; then
    echo "==> Using MPIF90 from the environment: $MPIF90"
elif [ "$MPI_WANT" = auto ]; then
    if [ "$JL_FAMILY" = unknown ]; then
        echo "ERROR: cannot determine which MPI to build against." >&2
        echo "  MPI.jl's binding could not be read, so there is nothing to match." >&2
        echo "  Say which one you want:" >&2
        echo "      bash compilef90.sh mpich" >&2
        echo "      bash compilef90.sh openmpi" >&2
        echo "  or point at a wrapper directly:" >&2
        echo "      MPIF90=/path/to/mpif90 bash compilef90.sh" >&2
        exit 1
    fi
    echo "==> Target: $JL_FAMILY (matching MPI.jl)"
    if ! MPIF90="$(find_wrapper "$JL_FAMILY")"; then
        echo "ERROR: MPI.jl is bound to $JL_FAMILY but no $JL_FAMILY mpif90 was found." >&2
        echo "  Looked in:" >&2
        candidate_bins "$JL_FAMILY" | sed 's/^/    /' >&2
        echo "  Install its Fortran wrapper, or pass MPIF90=/path/to/mpif90." >&2
        exit 1
    fi
else
    echo "==> Target: $MPI_WANT (requested)"
    if ! MPIF90="$(find_wrapper "$MPI_WANT")"; then
        echo "ERROR: no $MPI_WANT mpif90 found." >&2
        echo "  Looked in:" >&2
        candidate_bins "$MPI_WANT" | sed 's/^/    /' >&2
        echo "  Install it (brew install ${MPI_WANT/openmpi/open-mpi}), or pass" >&2
        echo "  MPIF90=/path/to/mpif90." >&2
        exit 1
    fi
    if [ "$JL_FAMILY" != unknown ] && [ "$JL_FAMILY" != "$MPI_WANT" ]; then
        echo
        echo "WARNING: you asked for $MPI_WANT but MPI.jl is bound to $JL_FAMILY."
        echo "         A coupled run needs both sides on the SAME MPI. Rebind the"
        echo "         Julia side (INSTALL.md §5.2) or build the other one here."
        echo
    fi
fi

echo "==> Fortran MPI wrapper: $MPIF90"
echo "==> It links against:"
"$MPIF90" -show
echo

WRAPPER_FAMILY="$(mpi_family "$("$MPIF90" -show 2>/dev/null || true) $MPIF90")"

#------------------------------------------------------------------------------
# Step 2 — build in a scratch directory, NOT in AlyaProxy/.
#
# This folder carries vendored copies of Open MPI's mpif.h / mpif-common.h /
# mpif-config.h. Fortran resolves `include 'mpif.h'` against the directory
# holding the source file FIRST, ahead of any -I path — so building the mpif.h
# variant in place would silently compile Open MPI's constants
# (MPI_COMM_WORLD=0) into a binary linked against a different MPI (MPICH uses
# 1140850688). That builds cleanly and then hands garbage handles to every MPI
# call. Compiling in a clean directory lets the wrapper's own -I win.
#------------------------------------------------------------------------------
BUILD_DIR="$(mktemp -d 2>/dev/null || mktemp -d -t alyaproxy)"
trap 'rm -rf "$BUILD_DIR"' EXIT

# macOS: gfortran hands the linker `-syslibroot $SDKROOT`. When SDKROOT is
# unset and the compiler's baked-in SDK path is stale — which any Xcode or
# Command Line Tools update can cause — the link dies with
#   ld: library not found for -lSystem
# even though compilation succeeded. Point it at the SDK that is actually
# installed. Only set what the user has not already chosen.
if [ "$(uname -s)" = "Darwin" ] && [ -z "${SDKROOT:-}" ] \
   && command -v xcrun >/dev/null 2>&1; then
    _sdk="$(xcrun --show-sdk-path 2>/dev/null || true)"
    if [ -n "$_sdk" ] && [ -d "$_sdk" ]; then
        export SDKROOT="$_sdk"
        echo "==> macOS: SDKROOT=$SDKROOT"
        echo
    fi
fi

# `cp` follows the myAlya.f90 -> alya_all2all_time_loop.f90 symlink.
cp "$SRC_DIR/myAlya.f90" "$BUILD_DIR/myAlya.f90"

if ls "$SRC_DIR"/mpif*.h >/dev/null 2>&1; then
    echo "==> Note: $SRC_DIR carries vendored mpif*.h (Open MPI's)."
    echo "    Building in a clean directory so your MPI's own headers are used."
    echo
fi

#------------------------------------------------------------------------------
# Step 3 — pick the Fortran binding.
#------------------------------------------------------------------------------
probe_mpif08() {
    printf 'program p\n  use mpi_f08\n  implicit none\nend program p\n' \
        > "$BUILD_DIR/probe.f90"
    ( cd "$BUILD_DIR" && "$MPIF90" -c probe.f90 -o probe.o ) >/dev/null 2>&1
}

# Does the compiler accept this flag? gfortran errors on unknown -f options,
# so a trivial compile is a reliable capability test.
probe_flag() {
    printf 'program p\n  implicit none\nend program p\n' > "$BUILD_DIR/flag.f90"
    ( cd "$BUILD_DIR" && "$MPIF90" "$1" -c flag.f90 -o flag.o ) >/dev/null 2>&1
}

if [ -n "${USE_MPIF08:-}" ]; then
    case "$USE_MPIF08" in
        1|true|yes|on)  USE_F08=1 ;;
        *)              USE_F08=0 ;;
    esac
    echo "==> Fortran binding: forced by USE_MPIF08=$USE_MPIF08"
else
    printf '==> Probing for the mpi_f08 module ... '
    if probe_mpif08; then
        USE_F08=1
        echo "available"
    else
        USE_F08=0
        echo "NOT available"
        echo "    This MPI ships no usable mpi_f08.mod (common with some MPICH"
        echo "    builds, and whenever gfortran is newer than the one that built"
        echo "    the MPI). Falling back to the legacy 'include mpif.h' bindings,"
        echo "    which the proxy supports natively. Nothing is lost — the MPI"
        echo "    calls are identical, only the handle types differ."
    fi
fi
rm -f "$BUILD_DIR"/probe.* 2>/dev/null || true

if [ "$USE_F08" = 1 ]; then
    FFLAGS=(-cpp -DUSEMPIF08)
    echo "==> Building with: use mpi_f08"
else
    FFLAGS=(-cpp)
    echo "==> Building with: include 'mpif.h'"

    # With mpif.h there are no interfaces, so every MPI routine is an external
    # procedure. gfortran >= 10 rejects a file that calls the same external
    # with different argument types — which is exactly what the MPI API is for
    # (MPI_Bcast on an INTEGER buffer and on a REAL(8) buffer):
    #
    #   Error: Type mismatch between actual argument at (1) and actual
    #   argument at (2) (REAL(8)/INTEGER(4)).
    #
    # Legal MPI usage, not a bug in the proxy; the mpi_f08 path avoids it via
    # generic interfaces. Downgrade the diagnostic where supported.
    printf '==> Probing for -fallow-argument-mismatch ... '
    if probe_flag -fallow-argument-mismatch; then
        FFLAGS+=(-fallow-argument-mismatch)
        echo "supported (enabled)"
    elif probe_flag -Wno-argument-mismatch; then
        FFLAGS+=(-Wno-argument-mismatch)
        echo "no, using -Wno-argument-mismatch (older gfortran)"
    else
        echo "not supported"
        echo "    If the build now fails with 'Type mismatch between actual"
        echo "    argument', your compiler needs its own equivalent flag —"
        echo "    add it via FFLAGS_EXTRA, e.g."
        echo "      FFLAGS_EXTRA=-some-flag bash compilef90.sh"
    fi
fi
rm -f "$BUILD_DIR"/flag.* 2>/dev/null || true

if [ -n "${FFLAGS_EXTRA:-}" ]; then
    # shellcheck disable=SC2206
    FFLAGS+=($FFLAGS_EXTRA)
    echo "==> Extra flags from FFLAGS_EXTRA: $FFLAGS_EXTRA"
fi

#------------------------------------------------------------------------------
# Step 4 — compile.
#------------------------------------------------------------------------------
echo "==> Compiling Alya.x ..."
COMPILE_LOG="$BUILD_DIR/compile.log"
set +e
( cd "$BUILD_DIR" && "$MPIF90" "${FFLAGS[@]}" myAlya.f90 -o "$OUT" ) 2>&1 | tee "$COMPILE_LOG"
COMPILE_RC=${PIPESTATUS[0]}
set -e

if [ "$COMPILE_RC" != 0 ]; then
    echo
    echo "==> BUILD FAILED (exit $COMPILE_RC). Reading the errors:"
    if grep -q 'library not found for -lSystem\|cannot find -lSystem' "$COMPILE_LOG"; then
        echo
        echo "  ld: library not found for -lSystem"
        echo
        echo "  This is a macOS toolchain problem, NOT an MPI one — compilation"
        echo "  succeeded and only the link failed. gfortran passes the linker"
        echo "  '-syslibroot \$SDKROOT'; the SDK it was built against is gone or"
        echo "  moved, which is what an Xcode / Command Line Tools update does."
        echo
        echo "  Fixes, in order:"
        echo "    1. Point at the SDK you actually have:"
        echo "         export SDKROOT=\$(xcrun --show-sdk-path)"
        echo "         bash compilef90.sh"
        echo "       (this script already tries that; if it still failed, xcrun"
        echo "        itself is returning a bad path — check: xcrun --show-sdk-path)"
        echo "    2. Repair the Command Line Tools:"
        echo "         sudo xcode-select --reset"
        echo "         xcode-select --install"
        echo "    3. Rebuild gfortran against the current SDK:"
        echo "         brew reinstall gcc"
    elif grep -q "Cannot open module file 'mpi_f08.mod'" "$COMPILE_LOG"; then
        echo "  The mpi_f08 probe passed but the real build could not use the"
        echo "  module. Force the legacy bindings:"
        echo "         USE_MPIF08=0 bash compilef90.sh"
    elif grep -q 'Error: Type mismatch between actual argument' "$COMPILE_LOG"; then
        echo "  Your compiler still rejects the mpif.h argument mismatches and"
        echo "  needs its own equivalent of -fallow-argument-mismatch:"
        echo "         FFLAGS_EXTRA=-some-flag bash compilef90.sh"
    else
        echo "  No known pattern matched — the compiler output above is the"
        echo "  authoritative message."
    fi
    exit "$COMPILE_RC"
fi

echo "==> Built $OUT"

#------------------------------------------------------------------------------
# Step 5 — verify the finished binary.
#------------------------------------------------------------------------------
echo "==> It links against:"
LINKED=""
if command -v otool >/dev/null 2>&1; then
    LINKED="$(otool -L "$OUT" 2>/dev/null | grep -i mpi || true)"
elif command -v ldd >/dev/null 2>&1; then
    LINKED="$(ldd "$OUT" 2>/dev/null | grep -i mpi || true)"
fi
if [ -n "$LINKED" ]; then
    printf '%s\n' "$LINKED" | sed 's/^[[:space:]]*/    /'
    BIN_FAMILY="$(mpi_family "$LINKED")"
else
    echo "    (could not read the binary's libraries)"
    BIN_FAMILY="$WRAPPER_FAMILY"
fi

echo
if [ "$JL_FAMILY" != unknown ] && [ "$BIN_FAMILY" != unknown ]; then
    if [ "$BIN_FAMILY" = "$JL_FAMILY" ]; then
        echo "==> OK: Alya.x and MPI.jl are both $BIN_FAMILY."
        echo "    MPI.jl libmpi: $JL_LIBMPI"
        echo "    Same family is necessary but not sufficient — check the paths"
        echo "    above are the SAME installation, then launch with that"
        echo "    installation's mpirun (put its bin first on PATH)."
    else
        echo "==> MISMATCH: Alya.x is $BIN_FAMILY but MPI.jl is $JL_FAMILY."
        echo "    A coupled run cannot work — it will hang in MPI_Init."
        echo "    Either rebuild here with '$JL_FAMILY', or rebind the Julia side"
        echo "    (INSTALL.md §5.2, case 2)."
        exit 1
    fi
else
    echo "    Compare these against MPI.jl by hand:"
    echo "      julia --project=$PROJECT_DIR -e 'using MPI; println(MPI.API.libmpi)'"
fi
echo
echo "==> Next steps. Run these from the project root:"
echo "      cd $PROJECT_DIR"
echo
echo "    Smoke-test Alya on its own — with no Jexpresso ranks in the world it"
echo "    exchanges nothing, runs its whole time loop and exits:"
echo "      mpirun -np 2 ./AlyaProxy/Alya.x"
echo
echo "    Full end-to-end check of both sides:"
echo "      bash tools/check_mpi_setup.sh"
echo
echo "    NOTE the leading './'. mpirun execs the program name directly and"
echo "    the current directory is not on PATH, so a bare 'mpirun -np 2 Alya.x'"
echo "    fails with:"
echo "      HYDU_create_process: execvp error on file Alya.x (No such file or directory)"

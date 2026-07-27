#!/usr/bin/env bash
#
# run_coupled.sh — launch Jexpresso (Julia) coupled with the Alya proxy
# (Fortran) as a single MPMD job, with preflight checks that turn the two
# classic "hangs forever and never starts" failures into an error message.
#
# Usage:
#   ./run_coupled.sh                       # 1 Alya rank + 2 Julia ranks, CompEuler 3dAlya
#   ./run_coupled.sh 2 2                   # 2 Alya ranks + 2 Julia ranks
#   ./run_coupled.sh 2 2 CompEuler thetaAlya
#
# Environment:
#   REBUILD_SYSIMAGE=1   build jexpresso.so first (removes Julia's cold-start JIT)
#   SYSIMAGE=path.so     use an existing sysimage (default: jexpresso.so if present)
#   SKIP_CHECKS=1        skip the preflight checks
#   JULIA=/path/to/julia override the Julia used
#   MPIRUN=/path/to/mpirun
#
# See RUN-COUPLED.md for the full setup guide, and
# `bash tools/check_mpi_setup.sh` for a deeper diagnosis when this script
# reports a problem.

set -euo pipefail

NALYA="${1:-1}"
NJULIA="${2:-2}"
EQS="${3:-CompEuler}"
CASE="${4:-3dAlya}"

cd "$(dirname "${BASH_SOURCE[0]}")"

JULIA="${JULIA:-julia}"
MPIRUN="${MPIRUN:-$(command -v mpirun || command -v mpiexec || true)}"

if [ -z "$MPIRUN" ]; then
    echo "ERROR: no mpirun/mpiexec on PATH. Set MPIRUN=/path/to/mpirun." >&2
    exit 1
fi

# Resolve the REAL julia binary. A juliaup shim on PATH is a shell wrapper that
# OpenMPI 5's prterun cannot always launch (RUN-COUPLED.md §6).
JULIA_BIN="$("$JULIA" -e 'print(joinpath(Sys.BINDIR, "julia"))')"

if [ ! -x AlyaProxy/Alya.x ]; then
    echo "ERROR: AlyaProxy/Alya.x not built. Build it with the mpif90 that matches" >&2
    echo "       MPI.jl's binding (RUN-COUPLED.md §5):" >&2
    echo "         cd AlyaProxy && bash compilef90.sh && cd .." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
if [ "${SKIP_CHECKS:-0}" != "1" ]; then
    echo "==> Preflight"

    # 1. Alya.x and MPI.jl must be the same MPI. A mismatch does not error at
    #    launch — it hangs inside MPI_Init, or silently gives each side its own
    #    world. This is the single most common coupled-run failure.
    if command -v otool >/dev/null 2>&1; then
        ALYA_LIBS="$(otool -L AlyaProxy/Alya.x 2>/dev/null | grep -i mpi || true)"
    elif command -v ldd >/dev/null 2>&1; then
        ALYA_LIBS="$(ldd AlyaProxy/Alya.x 2>/dev/null | grep -i mpi || true)"
    else
        ALYA_LIBS=""
    fi

    JL_INFO="$("$JULIA_BIN" --project=. -e '
        using MPIPreferences
        using MPI
        impl, _ = MPI.identify_implementation()
        println(impl, "|", MPI.API.libmpi, "|", MPIPreferences.binary)
    ' 2>/dev/null | tail -1)"

    JL_IMPL="${JL_INFO%%|*}"
    JL_REST="${JL_INFO#*|}"
    JL_LIBMPI="${JL_REST%%|*}"
    JL_BINARY="${JL_REST#*|}"

    echo "    MPI.jl : impl=${JL_IMPL:-?}  binary=${JL_BINARY:-?}"
    echo "             libmpi=${JL_LIBMPI:-?}"
    if [ -n "$ALYA_LIBS" ]; then
        echo "    Alya.x : $(printf '%s' "$ALYA_LIBS" | head -1 | sed 's/^[[:space:]]*//')"
    fi

    fam() {
        printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | \
        awk '/open ?mpi|open-mpi|libmpi\.40|libmpi\.so\.40/ {print "openmpi"; exit}
             /mpich|mvapich|intel|libmpi\.12|libmpi\.so\.12/ {print "mpich";   exit}
             {next} END {}'
    }
    A_FAM="$(fam "$ALYA_LIBS")"
    J_FAM="$(fam "$JL_IMPL $JL_LIBMPI")"
    L_FAM="$(fam "$("$MPIRUN" --version 2>&1 | head -2)")"

    if [ -n "$A_FAM" ] && [ -n "$J_FAM" ] && [ "$A_FAM" != "$J_FAM" ]; then
        echo "ERROR: Alya.x is linked against $A_FAM but MPI.jl is bound to $J_FAM." >&2
        echo "       A coupled run cannot work — it would hang in MPI_Init." >&2
        echo "       Reconcile them first (RUN-COUPLED.md §3), then rebuild Alya.x." >&2
        exit 1
    fi
    if [ -n "$L_FAM" ] && [ -n "$J_FAM" ] && [ "$L_FAM" != "$J_FAM" ]; then
        echo "ERROR: the launcher ($MPIRUN) is $L_FAM but MPI.jl is bound to $J_FAM." >&2
        echo "       The ranks would never form a shared world." >&2
        echo "       Run 'bash tools/check_mpi_setup.sh' for the fix." >&2
        exit 1
    fi
    echo "    OK"
fi

# ---------------------------------------------------------------------------
# Optional sysimage
# ---------------------------------------------------------------------------
SYSIMAGE="${SYSIMAGE:-jexpresso.so}"
if [ "${REBUILD_SYSIMAGE:-0}" = "1" ]; then
    echo "==> Building sysimage ($SYSIMAGE) — this takes a while, once."
    "$JULIA_BIN" --project=. create_Jexpresso_sysimage.jl
fi

# --startup-file=no skips ~/.julia/config/startup.jl on every rank; nothing in
# a batch run needs it and it is pure launch latency.
JULIA_ARGS=(--project=. --startup-file=no)
if [ -f "$SYSIMAGE" ]; then
    echo "==> Using sysimage $SYSIMAGE"
    JULIA_ARGS+=(--sysimage "$SYSIMAGE")
else
    echo "==> No sysimage ($SYSIMAGE not found). Every rank will JIT-compile"
    echo "    Jexpresso from scratch — tens of seconds per launch, every launch."
    echo "    Build it once with:  REBUILD_SYSIMAGE=1 ./run_coupled.sh"
fi

# ---------------------------------------------------------------------------
# Per-rank output tagging.
#
# Under a launcher, stdout is a pipe rather than a terminal, so it is
# block-buffered and interleaved: a working run can show nothing for a long
# time and then dump everything at once, which is indistinguishable from a
# hang. Tagging each line with its rank makes progress visible and shows which
# side of the coupling produced what. The flag is launcher-specific.
# Set TAG_OUTPUT=0 to turn it off.
LAUNCH_FLAGS=()
if [ "${TAG_OUTPUT:-1}" != "0" ]; then
    case "$("$MPIRUN" --version 2>&1 | head -2 | tr '[:upper:]' '[:lower:]')" in
        *"open mpi"*|*openmpi*|*prterun*) LAUNCH_FLAGS+=(--output tag) ;;
        *mpich*|*hydra*)                  LAUNCH_FLAGS+=(-prepend-rank) ;;
    esac
fi

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
# JEXPRESSO_COUPLED is exported for the WHOLE job rather than passed with `-x`
# inside the second app context: not every launcher honours `-x` after a `:`,
# and when it is dropped Jexpresso silently takes its standalone path and
# deadlocks against Alya with no diagnostic. Alya ignores the variable.
export JEXPRESSO_COUPLED=1

echo "==> $MPIRUN ${LAUNCH_FLAGS[*]-} -np $NALYA ./AlyaProxy/Alya.x : -np $NJULIA julia ... $EQS $CASE"
exec "$MPIRUN" ${LAUNCH_FLAGS[@]+"${LAUNCH_FLAGS[@]}"} \
    -np "$NALYA"  ./AlyaProxy/Alya.x \
  : -np "$NJULIA" "$JULIA_BIN" "${JULIA_ARGS[@]}" ./src/Jexpresso.jl "$EQS" "$CASE"

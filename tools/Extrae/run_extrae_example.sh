#!/bin/bash
#==============================================================================
# run_extrae_example.sh  --  convenience launcher for the MPI Extrae example
#
# Usage (from the Jexpresso project root):
#     ./tools/Extrae/run_extrae_example.sh [NRANKS]
#
# NRANKS defaults to 4 (one rank per core on a 4-core MacBook Air).
#
# This mirrors jexp_mpich.sh: it asks MPI.jl for its bundled `mpiexec` and
# launches the example under it, so it uses whatever MPI your MPIPreferences
# point at (MPICH_jll is the recommended route on macOS — see INSTALL.md).
#
# To capture a REAL Paraver trace, run this on Linux and point EXTRAE_LIB at
# the Extrae MPI tracing library (libmpitrace.so). The preload is attached to
# the *rank processes only* (via `env`), so Extrae is loaded after Julia but
# before the MPI library — the arrangement recommended by the Extrae.jl paper
# (this avoids the libstdc++ version clashes you get from preloading on the
# `julia` binary itself). EXTRAE_CONFIG_FILE should point at extrae.xml.
#
#   export EXTRAE_LIB=$EXTRAE_HOME/lib/libmpitrace.so
#   export EXTRAE_CONFIG_FILE=$PWD/tools/Extrae/extrae.xml
#   ./tools/Extrae/run_extrae_example.sh 32
#
# Leave EXTRAE_LIB unset to run without tracing (e.g. on macOS).
#==============================================================================
set -euo pipefail

# USER: change to your julia path if `julia` is not on your PATH ------------
JULIA="${JULIA:-julia}"
# ---------------------------------------------------------------------------

NRANKS="${1:-4}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
EXAMPLE="${SCRIPT_DIR}/extrae_mpi_jexpresso_pattern.jl"

# If EXTRAE_LIB is set, wrap each rank in `env LD_PRELOAD=<lib>` so the Extrae
# MPI interception library is loaded per rank. Otherwise launch plainly.
#
# EXTRAE_LIB must be the *C* MPI tracing library  libmpitrace.so
#   (NOT the Fortran one, libmpitracef.so — MPI.jl uses the C MPI ABI).
#
# EXTRAE_LIBPATH (optional) is prepended to LD_LIBRARY_PATH for the ranks so
# the loader can resolve Extrae's own dependencies (libunwind, PAPI, ...).
# With the Julia artifact, get it from the jll:
#   export EXTRAE_LIBPATH=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')
# With a system module it is usually  $EXTRAE_HOME/lib.
if [[ -n "${EXTRAE_LIB:-}" ]]; then
    case "${EXTRAE_LIB}" in
        *libmpitracef.so)
            echo "WARNING: EXTRAE_LIB points at the Fortran lib (libmpitracef.so);" >&2
            echo "         MPI.jl needs the C lib libmpitrace.so — MPI calls won't be traced." >&2
            ;;
    esac
    echo "Tracing ON: LD_PRELOAD=${EXTRAE_LIB}"
    echo "  EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE:-<unset!>}"
    ENV_ASSIGN="LD_PRELOAD=${EXTRAE_LIB}"
    if [[ -n "${EXTRAE_LIBPATH:-}" ]]; then
        echo "  EXTRAE_LIBPATH=${EXTRAE_LIBPATH}"
        ENV_ASSIGN="LD_LIBRARY_PATH=${EXTRAE_LIBPATH}:${LD_LIBRARY_PATH:-} ${ENV_ASSIGN}"
    fi
    PRELOAD_PREFIX="env ${ENV_ASSIGN}"
else
    echo "Tracing OFF (EXTRAE_LIB unset): running without an Extrae trace."
    PRELOAD_PREFIX=""
fi

echo "Launching Extrae MPI example with ${NRANKS} ranks ..."
cd "${PROJECT_ROOT}"
"${JULIA}" --project=. -e "
  using MPI
  run(\`\$(mpiexec()) -n ${NRANKS} ${PRELOAD_PREFIX} \$(Base.julia_cmd()) --project=. ${EXAMPLE}\`)
"

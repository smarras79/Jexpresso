#!/bin/bash -l
# ==============================================================================
# run_jexpresso_traced.sh  --  run a REAL Jexpresso case under Extrae tracing
#
# One command, no fragile quoting. It sets up the Extrae preload, turns the
# in-solver instrumentation on (JEXPRESSO_EXTRAE=1), and launches the solver
# with the Extrae library LD_PRELOAD-ed onto the rank processes only (loaded
# after Julia but before MPI — the arrangement the Extrae.jl paper recommends).
#
# Usage (from the Jexpresso project root):
#     ./tools/Extrae/run_jexpresso_traced.sh [NRANKS] [EQ] [CASE]
#
#   NRANKS  number of MPI ranks         (default 4)
#   EQ      equations dir under problems/ (default CompEuler)
#   CASE    case dir under problems/EQ/   (default theta)
#
# Examples:
#     ./tools/Extrae/run_jexpresso_traced.sh 4  CompEuler theta
#     ./tools/Extrae/run_jexpresso_traced.sh 16 CompEuler 3d
#
# EXTRAE_LIB / EXTRAE_LIBPATH / EXTRAE_CONFIG_FILE are auto-derived from
# Extrae_jll if you have not exported them; export your own to override (e.g.
# to use a system Extrae module instead of the Julia artifact).
# ==============================================================================
set -euo pipefail

JULIA="${JULIA:-julia}"

NRANKS="${1:-4}"
EQ="${2:-CompEuler}"
CASE="${3:-theta}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${PROJECT_ROOT}"

# --- auto-derive the Extrae preload from Extrae_jll (unless already set) ----
if [[ -z "${EXTRAE_LIB:-}" || -z "${EXTRAE_LIBPATH:-}" ]]; then
    ART="$("${JULIA}" --project=. -e 'using Extrae_jll; print(Extrae_jll.artifact_dir)')"
    : "${EXTRAE_LIB:=${ART}/lib/libmpitrace.so}"
    : "${EXTRAE_LIBPATH:=$("${JULIA}" --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')}"
fi
: "${EXTRAE_CONFIG_FILE:=${PROJECT_ROOT}/tools/Extrae/extrae.xml}"

if [[ ! -f "${EXTRAE_LIB}" ]]; then
    echo "ERROR: EXTRAE_LIB not found: ${EXTRAE_LIB}" >&2
    echo "       ls \"\$(dirname \"${EXTRAE_LIB}\")\"/libmpitrace*.so   to see what exists." >&2
    exit 1
fi
case "${EXTRAE_LIB}" in
    *libmpitracef.so) echo "WARNING: that is the Fortran lib; MPI.jl needs the C lib libmpitrace.so." >&2 ;;
esac

# Turn the in-solver instrumentation on (read by src/kernel/infrastructure/Profiling.jl).
export JEXPRESSO_EXTRAE=1
export EXTRAE_CONFIG_FILE

echo "EXTRAE_LIB=${EXTRAE_LIB}"
echo "EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE}"
echo "Tracing ${EQ}/${CASE} on ${NRANKS} ranks ..."

# The env-assignment string is built HERE in bash (no spaces in the paths, only
# colons), then embedded literally into the Julia backtick command — so neither
# bash nor Julia has to interpolate $LD_LIBRARY_PATH at the wrong time.
ENV_PREFIX="env -u OMP_NUM_THREADS LD_LIBRARY_PATH=${EXTRAE_LIBPATH}:${LD_LIBRARY_PATH:-} LD_PRELOAD=${EXTRAE_LIB}"

"${JULIA}" --project=. -e "
  using MPI
  run(\`\$(mpiexec()) -n ${NRANKS} ${ENV_PREFIX} \$(Base.julia_cmd()) --project=. src/Jexpresso.jl ${EQ} ${CASE}\`)
"

echo "--- Done. Paraver trace files in ${PROJECT_ROOT}: ---"
ls -la jexpresso-extrae.prv jexpresso-extrae.pcf jexpresso-extrae.row 2>/dev/null || \
    echo "  (no .prv yet — check output above; see tools/Extrae/README.md)"

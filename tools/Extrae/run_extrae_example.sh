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
# To capture a REAL Paraver trace you must run this on Linux with Extrae
# installed and LD_PRELOAD / MPIPreferences preloads configured (see README).
#==============================================================================
set -euo pipefail

# USER: change to your julia path if `julia` is not on your PATH ------------
JULIA="${JULIA:-julia}"
# ---------------------------------------------------------------------------

NRANKS="${1:-4}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
EXAMPLE="${SCRIPT_DIR}/extrae_mpi_jexpresso_pattern.jl"

echo "Launching Extrae MPI example with ${NRANKS} ranks ..."
cd "${PROJECT_ROOT}"
"${JULIA}" --project=. -e "
  using MPI
  run(\`\$(mpiexec()) -n ${NRANKS} \$(Base.julia_cmd()) --project=. ${EXAMPLE}\`)
"

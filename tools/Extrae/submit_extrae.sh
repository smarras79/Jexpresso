#!/bin/bash -l
# ==============================================================================
# submit_extrae.sh  --  SLURM batch job for the Extrae.jl MPI example
#
# Profiles tools/Extrae/extrae_mpi_jexpresso_pattern.jl on a SINGLE node and
# produces a Paraver trace (jexpresso-extrae.{prv,pcf,row}).
#
# It allocates one CPU and enough RAM PER RANK (Julia + Extrae need ~1 GB
# each — under-provisioning is what causes the oom_kill / EXIT CODE 9), sets
# up the Extrae preload automatically from Extrae_jll, and launches the
# example through tools/Extrae/run_extrae_example.sh.
#
# Submit from the Jexpresso project root:
#     sbatch tools/Extrae/submit_extrae.sh
#
# Change the rank count by editing --ntasks below (keep --nodes=1 for this
# single-node profiling example). Trace output lands in the submit directory.
# ==============================================================================
#SBATCH --job-name=jexp-extrae
#SBATCH --output=%x.%j.out          # %x=JobName  %j=JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=high_smarras
#SBATCH --account=smarras           # <-- set to your PI's allocation/account
#SBATCH --nodes=1                   # single node for this profiling example
#SBATCH --ntasks=32                 # = number of MPI ranks (one per core)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M         # ~4 GB/rank; lower to 2000M if your node is small
#SBATCH --time=01:00:00             # HH:MM:SS

set -euo pipefail

# --- locate the project (the dir you submitted from) -----------------------
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${PROJECT_ROOT}"

# --- environment -----------------------------------------------------------
module load Julia
module load bright shared mpich/ge/gcc/64

# Pure-MPI run: silence Extrae's "OMP_NUM_THREADS is set but OpenMP not
# supported" notice (run_extrae_example.sh also unsets it per rank).
unset OMP_NUM_THREADS || true

# --- make sure the project + Extrae are precompiled ONCE -------------------
# Doing this here (serially) avoids 32 ranks racing to precompile at launch,
# which is slow and memory-hungry.
echo "--- Instantiating + precompiling project (incl. Extrae) ---"
julia --project=. -e '
    import Pkg
    Pkg.instantiate()
    # add Extrae if it is not already a dependency (Linux-only binary)
    haskey(Pkg.project().dependencies, "Extrae") || Pkg.add("Extrae")
    Pkg.precompile()'
echo "--- Precompilation complete ---"

# --- set up the Extrae preload from the Julia artifact ---------------------
# EXTRAE_LIB     : the C MPI tracing library (libmpitrace.so, NOT libmpitracef.so)
# EXTRAE_LIBPATH : the jll's full dependency search path (libunwind, PAPI, ...)
ART="$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.artifact_dir)')"
export EXTRAE_LIB="${ART}/lib/libmpitrace.so"
export EXTRAE_LIBPATH="$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')"
export EXTRAE_CONFIG_FILE="${PROJECT_ROOT}/tools/Extrae/extrae.xml"

if [[ ! -f "${EXTRAE_LIB}" ]]; then
    echo "ERROR: ${EXTRAE_LIB} not found." >&2
    echo "       List what the artifact provides with:" >&2
    echo "         ls \"${ART}/lib\"/libmpitrace*.so" >&2
    exit 1
fi

echo "EXTRAE_LIB=${EXTRAE_LIB}"
echo "EXTRAE_CONFIG_FILE=${EXTRAE_CONFIG_FILE}"
echo "Launching on ${SLURM_NTASKS} ranks ..."

# --- run -------------------------------------------------------------------
# run_extrae_example.sh attaches LD_PRELOAD / LD_LIBRARY_PATH per rank and
# launches via MPI.jl's bundled mpiexec.
./tools/Extrae/run_extrae_example.sh "${SLURM_NTASKS}"

echo "--- Done. Paraver trace files in ${PROJECT_ROOT}: ---"
ls -la jexpresso-extrae.prv jexpresso-extrae.pcf jexpresso-extrae.row 2>/dev/null || \
    echo "  (no .prv yet — check the .err log; see tools/Extrae/README.md)"

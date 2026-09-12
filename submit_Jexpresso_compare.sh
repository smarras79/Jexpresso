#!/bin/bash -l
#SBATCH --job-name=smoothVortexScan
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks=128                 # = SV_NP × SV_JOBS = 16 × 4
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # avoid thread oversubscription
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
julia --project=. -e 'using MPI; using Jexpresso' 2>/dev/null || \
    julia --project=. -e 'include("src/Jexpresso.jl")' --warmup-only

echo "--- Setup complete, launching scan: 4 concurrent jobs x 16 ranks ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

export MPIEXEC="srun --exclusive --mpi=pmi2"

SV_NOPS="4 5 6 7" SV_NELX="16 32 64" SV_VISC="dsgs none" SV_NP=16 SV_JOBS=6 \
SV_PLOT_NOPS="4 6" tools/smooth_vortex_mpi_scan.sh

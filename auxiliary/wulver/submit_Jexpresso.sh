#!/bin/bash -l
#SBATCH --job-name=LESsmago
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=71:59:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # avoid 64x64 thread oversubscription
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
julia --project=. -e 'using MPI; using Jexpresso' 2>/dev/null || \
  julia --project=. -e 'include("src/Jexpresso.jl")' --warmup-only

echo "--- Setup complete, launching 64 ranks ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

mpirun -np 64 julia --project=. src/Jexpresso.jl CompEuler ffs_step

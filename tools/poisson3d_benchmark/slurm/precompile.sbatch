#!/bin/bash -l
#SBATCH --job-name=p3d-precompile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Instantiate and precompile once, before the benchmark jobs, so that they do
# not all precompile at the same time; then run the correctness checks
# (2D equivalence with Jexpresso, 3D convergence and solver agreement).
source "$P3D_DIR/common.sh"      # P3D_DIR: this directory, exported by submit_wulver3d.sh
export JULIA_PKG_OFFLINE=false
$JULIA --project="$REPO" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using Jexpresso; println("precompiled")'
$JULIA --project="$REPO" -t 4 tools/poisson3d_benchmark/verify_2d.jl
$JULIA --project="$REPO" -t 4 tools/poisson3d_benchmark/verify_3d.jl

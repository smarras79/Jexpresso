#!/bin/bash -l
#---------------------------------------------------------------------------
# The BDF-k study of the DynSGS residual sensor on Wulver: both smooth
# vortices, every case of tools/smooth_vortex_bdf_study.sh, on one node.
#
#   sbatch tools/smooth_vortex_bdf_study_wulver.sh
#   sbatch --export=ALL,STUDY_STEPS=dtcheck tools/smooth_vortex_bdf_study_wulver.sh   # the dt check alone
#   sbatch --export=ALL,STUDY_CASES=MHD/smoothVortex tools/smooth_vortex_bdf_study_wulver.sh
#
# Steps 1-3 are the setup of the fluxEmergence runs on this cluster,
# unchanged. The cases are SERIAL (SV_NP = 1), SV_JOBS = 16 at a time: each
# call of the sweep driver in the study runs at most 16 cases (4 meshes x 4
# orders, one variant), the largest about 15 min serial (MHD, 32^2, N = 7,
# DynSGS) and every one starting with a 2-3 min compile, so the whole study
# is about 3 h. SV_NP x SV_JOBS must equal --ntasks-per-node.
#
# MPI ranks per case instead (the 32^2 cases ~8x faster, the compile not) —
# the command line overrides the #SBATCH lines below, so nothing to edit:
#   sbatch --ntasks-per-node=128 --export=ALL,SV_NP=8,SV_JOBS=16 tools/smooth_vortex_bdf_study_wulver.sh
# Each case is then `mpirun -n 8`, sixteen at once; the driver runs a 2-rank
# hello through that launcher first, so one that does not work fails in a
# line, not in every log.
#
# Output (in the job's .out file): every error against Galerkin's, per case,
# and the dt check (L1 at dt and dt/2). Figures in
# output/<case>/output/bdf_study/, errors in problems/<case>/errors (dt/2 in
# errors.dt2). Stores already there are moved aside, never deleted.
#---------------------------------------------------------------------------
#SBATCH --job-name=sv-bdf-study
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:59:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # one core per serial case (or per rank)
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
julia --project=. -e 'using MPI; using Jexpresso' 2>/dev/null || \
  julia --project=. -e 'include("src/Jexpresso.jl")' --warmup-only

echo "--- Setup complete, launching the BDF-k study ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # the cases must never attempt to precompile

export SV_NP=${SV_NP:-1}
export SV_JOBS=${SV_JOBS:-16}
export MPIEXEC=${MPIEXEC:-mpirun}   # used only when SV_NP > 1

tools/smooth_vortex_bdf_study.sh

echo "--- done ---"

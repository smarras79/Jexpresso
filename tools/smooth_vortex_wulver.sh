#!/bin/bash -l
#---------------------------------------------------------------------------
# Smooth MHD vortex convergence sweep on Wulver.
#
#   sbatch tools/smooth_vortex_wulver.sh
#
# The setup steps (1-3) are the ones from the working orszagTang script on
# this cluster, unchanged; only the launch at the end differs, because this is
# a SWEEP of independent cases rather than one run.
#
# As written: orders 4 and 6 on 32² and 64² elements, DynSGS and plain
# Galerkin — the two panels of Dao & Nazarov (2022) Fig. 1 — with
# SV_NP = 16 ranks per case and SV_JOBS = 4 cases at a time, i.e. 64 cores
# and 8 cases in all.
#
#   SV_NP x SV_JOBS must equal --ntasks-per-node below. To use a whole
#   128-core node, either SV_JOBS=8 (same case size, twice as many at once)
#   or SV_NP=32 (twice the ranks per case).
#---------------------------------------------------------------------------
#SBATCH --job-name=smoothVortex
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # the parallelism here is MPI ranks
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: one real case, so the concurrent ones do not all compile at once ---"
mkdir -p logs
JEXPRESSO_SV_NOP=4 JEXPRESSO_SV_NELX=4 JEXPRESSO_SV_VISC=none \
JEXPRESSO_SV_SOLVER=vern9 JEXPRESSO_SV_TEND=0.01 \
    julia --project=. src/Jexpresso.jl MHD smoothVortex > logs/warmup.log 2>&1 \
    || { echo "warm-up FAILED — see logs/warmup.log"; tail -40 logs/warmup.log; exit 1; }
rm -f problems/MHD/smoothVortex/errors/*.dat   # it was a t = 0.01 case, not the sweep

echo "--- Setup complete, launching the sweep ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

#---------------------------------------------------------------------------
# HOW ONE CASE IS LAUNCHED. Four of them run at once inside this allocation,
# so each must be pinned to its own cores; that is what srun's step-level
# --exclusive does. MPICH here means PMI2 — `srun --mpi=list` confirms what
# this SLURM was built with.
#
# If srun cannot launch MPICH on this cluster, use the mpirun fallback: one
# case at a time on all 64 ranks. It is the same 8 cases and the same
# figures, just without the 4-way overlap.
#---------------------------------------------------------------------------
# MPICH on this cluster means PMI2. `srun --mpi=list` says what this SLURM
# was built with; the scan runs a 2-rank hello through this launcher before
# the sweep, so a wrong value fails in one line instead of eight log files.
export MPIEXEC=${MPIEXEC:-"srun --exclusive --mpi=pmi2"}
SV_NP=${SV_NP:-16}
SV_JOBS=${SV_JOBS:-4}

# --- fallback: uncomment these two lines and comment the three above ---
# export MPIEXEC="mpirun"
# SV_NP=64; SV_JOBS=1

SV_NOPS="4 6" \
SV_NELX="32 64" \
SV_VISC="dsgs none" \
SV_NP="$SV_NP" \
SV_JOBS="$SV_JOBS" \
SV_SOLVER=vern9 \
SV_PLOT_NOPS="4 6" \
SV_TEND=${SV_TEND:-1.0} \
    tools/smooth_vortex_mpi_scan.sh

echo "--- done ---"
echo "Nothing appears in THIS file while a case runs: each case writes to"
echo "logs/sv_<visc>_nop<N>_nelx<M>.log, and its error and the figures are"
echo "written when it reaches its final time. SV_TEND=0.05 first is a"
echo "ten-minute check of the whole pipeline before an eight-hour sweep."
echo "figures: output/MHD/smoothVortex/output/convergence_{dsgs,galerkin}*-it*.png"
echo "errors:  problems/MHD/smoothVortex/errors/"
echo "replot any subset without rerunning:"
echo "    julia --project=. tools/smooth_vortex_plot.jl 4,6"

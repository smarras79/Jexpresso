#!/bin/bash
#---------------------------------------------------------------------------
# Smooth MHD vortex convergence sweep on a SLURM cluster.
#
#   sbatch tools/smooth_vortex_slurm.sh
#
# As written: orders 4 and 6 on 32² and 64² elements, both panels of
# Dao & Nazarov (2022) Fig. 1 (DynSGS and plain Galerkin), 16 MPI ranks per
# case and 4 cases at a time — 64 cores, 8 cases in all.
#
#   SV_NP x SV_JOBS = 16 x 4 = 64 = --ntasks below. KEEP THE THREE IN STEP:
#   ask SLURM for exactly the product, or the steps will queue behind each
#   other (too few tasks) or leave cores idle (too many).
#
# Every case is one MPI job step (srun --exclusive -n 16), and four of them
# run at once inside the allocation. The sweep writes one error file per
# (order, mesh, viscosity) into problems/MHD/smoothVortex/errors and every
# case redraws the figures from the whole store as it finishes, atomically,
# so the last one to land leaves the complete comparison behind.
#---------------------------------------------------------------------------
#SBATCH --job-name=sv-scan
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=sv-scan-%j.out
#SBATCH --error=sv-scan-%j.err
##SBATCH --nodes=2                  # let SLURM place the 64 tasks, or pin it
##SBATCH --partition=<your-partition>
##SBATCH --exclusive                # if the node is shared and RAM is tight

set -u
cd "${SLURM_SUBMIT_DIR:-$PWD}"

#---------------------------------------------------------------------------
# Environment. Replace with your cluster's module names; Jexpresso's MPI.jl
# must be built against the SAME MPI these modules provide —
# tools/check_mpi_setup.sh checks exactly that and is worth running once,
# interactively, before the first sbatch.
#---------------------------------------------------------------------------
# module purge
# module load julia/1.11 openmpi/4.1

export JULIA=${JULIA:-julia}
export JULIA_NUM_THREADS=1          # the parallelism here is MPI ranks
export OMP_NUM_THREADS=1

mkdir -p logs

#---------------------------------------------------------------------------
# One warm-up case, serial, before the sweep. The first run of a session
# compiles the whole RHS and integrator; without this, four cases would do it
# at the same time and contend for the depot's precompile locks. It costs a
# couple of minutes and saves more.
#---------------------------------------------------------------------------
echo "=== warm-up $(date +%T)"
JEXPRESSO_SV_NOP=4 JEXPRESSO_SV_NELX=4 JEXPRESSO_SV_VISC=none \
JEXPRESSO_SV_SOLVER=vern9 JEXPRESSO_SV_TEND=0.01 \
    "$JULIA" --project=. src/Jexpresso.jl MHD smoothVortex > logs/warmup.log 2>&1 \
    || { echo "warm-up failed — see logs/warmup.log"; tail -30 logs/warmup.log; exit 1; }
# it wrote an error at t = 0.01 that does not belong to the sweep
rm -f problems/MHD/smoothVortex/errors/*.dat

#---------------------------------------------------------------------------
# The sweep. MPIEXEC is how one case is launched: inside an allocation that
# is what srun does, and --exclusive at STEP level is what stops the four
# concurrent steps from landing on the same cores. --mpi=pmix or pmi2
# depending on what your SLURM was built with (srun --mpi=list shows it).
#---------------------------------------------------------------------------
export MPIEXEC=${MPIEXEC:-"srun --exclusive --mpi=pmix"}

echo "=== sweep $(date +%T)"
SV_NOPS="4 6" \
SV_NELX="32 64" \
SV_VISC="dsgs none" \
SV_NP=16 \
SV_JOBS=4 \
SV_SOLVER=vern9 \
SV_PLOT_NOPS="4 6" \
    tools/smooth_vortex_mpi_scan.sh

echo "=== done $(date +%T)"
echo "    figures: output/MHD/smoothVortex/output/convergence_{dsgs,galerkin}*-it*.png"
echo "    errors:  problems/MHD/smoothVortex/errors/"
echo "    replot any subset without rerunning:"
echo "        julia --project=. tools/smooth_vortex_plot.jl 4,6"

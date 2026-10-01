#!/bin/bash -l
#SBATCH --job-name=ppb3d
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=71:59:00
#SBATCH --mem-per-cpu=4000M
# =============================================================================
#  3D periodic Poisson solver comparison on one Wulver node (128 cores, 500 GB):
#      sbatch tools/poisson3d_benchmark/slurm/run_wulver3d.sbatch
#  The seven solvers run side by side, one Julia process each with THREADS
#  threads, each sweeping the sizes below; then the figures are drawn.
#  Resubmitting resumes: finished configurations are skipped.
# =============================================================================
module load Julia/1.11.9
module load GCC MPICH
cd /project/smarras/smarras/Jexpresso/        # checkout of sm/elementLearning

export JULIA_NUM_THREADS=1          # each benchmark process sets its own (-t)
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only
# Wulver's MPICH uses UCX (ch4:ucx), whose SIGSEGV handler steals the signals
# Julia's GC uses for safepoints: threaded runs then die with "Segmentation
# fault: invalid permissions for mapped object" in ijl_gc_safepoint. Leave
# those signals to Julia (MPI.jl's documented setting for UCX):
export UCX_ERROR_SIGNALS=""

OUTDIR=ppb3d_wulver                 # results, figures, per-solver logs
SOLVERS="sem sem_amg sem_jacobi sc_direct sc_amg ps fft"
THREADS=16                          # per solver process: 7 x 16 = 112 of 128 cores
# Sweep as in the 2D benchmark: mesh levels of LEVELS^3 elements, each with
# the SEM orders NOPS (the Fourier solvers on the (ne*N)^3 grid, the same
# unknowns n = (ne*N)^3). Configurations above MAXN unknowns are not run:
# the direct solvers (sem, sc_direct) stop at 2.2e6 (the next sizes need
# 400 GB and more for the factor), the others at 1.7e7.
LEVELS="8 16 32 64"
NOPS="2 3 4 5 6 7 8"
MAXN=17000000;  MAXN_DIRECT=2200000
nops_for() {   # nops_for ne maxn -> the orders N of NOPS with (ne*N)^3 <= maxn, comma-separated
    local out=""
    for N in $NOPS; do [ $(( ($1*N)**3 )) -le $2 ] && out="$out${out:+,}$N"; done
    echo "$out"
}

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1
echo "--- 2. Serial precompile (one process, many cores internally) ---"
# resolve first: an existing Manifest.toml may predate packages added to Project.toml
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()' || exit 1
echo "--- 3. Serial warm-up and checks (2D against Jexpresso, 3D convergence) ---"
julia --project=. -e 'using MPI; using Jexpresso' || exit 1
julia --project=. -t 4 tools/poisson3d_benchmark/verify_2d.jl || exit 1
julia --project=. -t 4 tools/poisson3d_benchmark/verify_3d.jl || exit 1

echo "--- Setup complete, launching the solvers ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # the runs must never attempt to precompile
mkdir -p "$OUTDIR/logs"
for S in $SOLVERS; do
    case $S in sem|sc_direct) maxn=$MAXN_DIRECT ;; *) maxn=$MAXN ;; esac
    (   # one process per mesh level, orders in increasing order
        for ne in $LEVELS; do
            nops=$(nops_for $ne $maxn); [ -n "$nops" ] || continue
            julia --project=. -t $THREADS tools/poisson3d_benchmark/bench3d.jl \
                  --solver $S --nes $ne --nops $nops --blas-threads $THREADS \
                  --outdir "$OUTDIR/parts/$S" --resume
        done
    ) > "$OUTDIR/logs/$S.log" 2>&1 &
done
wait

echo "--- Merging results and drawing the figures ---"
python3 tools/poisson3d_benchmark/plot3d.py "$OUTDIR"
echo "results: $OUTDIR/results.md, figures: $OUTDIR/assets"

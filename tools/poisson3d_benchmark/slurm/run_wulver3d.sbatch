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

OUTDIR=ppb3d_wulver                 # results, figures, per-solver logs
SOLVERS="sem sem_amg sem_jacobi sc_direct sc_amg ps fft"
THREADS=16                          # per solver process: 7 x 16 = 112 of 128 cores
# Elements per direction; unknowns n = (ne*N)^3, from 4.1e3 to 1.7e7.
# The direct solvers (sem, sc_direct) stop at NES_DIRECT (n = 2.1e6): the next
# size needs about 400 GB for the factor alone, and 1.7e7 about 2.5 TB.
NES_N2="8,12,16,24,32,48,64,96,128";  NES_N2_DIRECT="8,12,16,24,32,48,64"
NES_N4="4,6,8,12,16,24,32,48,64";     NES_N4_DIRECT="4,6,8,12,16,24,32"

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
    case $S in sem|sc_direct) n2=$NES_N2_DIRECT; n4=$NES_N4_DIRECT ;; *) n2=$NES_N2; n4=$NES_N4 ;; esac
    (   # one process per N so that each sweeps its sizes in increasing order
        for sweep in "2 $n2" "4 $n4"; do
            set -- $sweep
            julia --project=. -t $THREADS tools/poisson3d_benchmark/bench3d.jl \
                  --solver $S --nop $1 --nes $2 --blas-threads $THREADS \
                  --outdir "$OUTDIR/parts/$S" --resume
        done
    ) > "$OUTDIR/logs/$S.log" 2>&1 &
done
wait

echo "--- Merging results and drawing the figures ---"
python3 tools/poisson3d_benchmark/plot3d.py "$OUTDIR"
echo "results: $OUTDIR/results.md, figures: $OUTDIR/assets"

#!/bin/bash -l
#SBATCH --job-name=ppb2d
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=4000M
# =============================================================================
#  The 2D periodic Poisson study of last week, rerun with the current code:
#      sbatch tools/periodic_poisson_benchmark/slurm/run_wulver2d.sbatch
#  Same sweep as the laptop run (README "high resolution"): the deck
#  problems/Elliptic/poisson_periodic_sem through Jexpresso.run_case, 16x16
#  elements refined L = 0..3 times (16², 32², 64², 128² elements), N = 2..8,
#  the six solvers, --warmup small, the same unknown caps. The six solvers
#  run side by side (one Julia process each, 5 threads), then the results
#  are merged, compared with earlier runs (compare2d.py) and plotted.
#  32 cores / 125 GB: the 2D sweep is serial per solver and needs < 10 GB.
#  Resubmitting resumes: finished configurations are skipped.
# =============================================================================
module load Julia/1.11.9
module load GCC MPICH
cd /project/smarras/smarras/Jexpresso/        # checkout of sm/elementLearning

export JULIA_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only
export UCX_ERROR_SIGNALS=""         # leave SIGSEGV to Julia's GC (Wulver's MPICH uses UCX)

OUTDIR=ppb2d_wulver                 # results, comparison, figures, per-solver logs
SOLVERS="sem sem_amg sc_direct sc_amg ps fft"
THREADS=5                           # per solver process: 6 x 5 = 30 of 32 cores
LEVELS=0:3
NOPS=2:8
CAPS=all=810000,sem_amg=300000      # last week's --max-unknowns
# Earlier results to compare with: the committed 16x16 run, plus last week's
# 4-level run if you copy its results.csv here (scp from the laptop:
# ppb_highres/results.csv -> $OUTDIR/reference_highres.csv)
REFERENCES="tools/periodic_poisson_benchmark/results.csv $OUTDIR/reference_highres.csv"

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1
echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()' || exit 1
echo "--- 3. Serial warm-up and the 2D check against the Kronecker cross-check ---"
julia --project=. -e 'using MPI; using Jexpresso' || exit 1
julia --project=. -t 4 tools/poisson3d_benchmark/verify_2d.jl || exit 1

echo "--- Setup complete, launching the solvers ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # the runs must never attempt to precompile
mkdir -p "$OUTDIR/logs"
for S in $SOLVERS; do
    julia --project=. -t $THREADS tools/periodic_poisson_benchmark/pipeline.jl \
          --solvers $S --levels $LEVELS --nops $NOPS --warmup small \
          --max-unknowns $CAPS --outdir "$OUTDIR/parts/$S" --resume --noplot \
          > "$OUTDIR/logs/$S.log" 2>&1 &
done
wait
for S in $SOLVERS; do grep "^PPB  [a-zA-Z]" "$OUTDIR/logs/$S.log" | grep -v " run \| skip "; done

echo "--- Merging, comparing with the earlier runs, drawing the figures ---"
python3 tools/periodic_poisson_benchmark/compare2d.py "$OUTDIR" $REFERENCES
status=$?
python3 tools/periodic_poisson_benchmark/plot.py "$OUTDIR/results.csv" "$OUTDIR/assets"
echo "results: $OUTDIR/results.csv, comparison: $OUTDIR/compare.md, figures: $OUTDIR/assets"
exit $status

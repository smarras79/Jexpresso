#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread
# One configuration (solver, ne, N) of the 3D benchmark. Submitted by
# submit_wulver3d.sh, which sets cores, memory and time per configuration and
# exports P3D_ENV P3D_SOLVER P3D_NE P3D_NOP P3D_THREADS. Protocol: a discarded
# warm-up run, then the recorded second run in the same Julia session; writes
#   $OUTDIR/parts/<solver>_ne<ne>_N<N>/results.csv
# (peak memory in it is this job's own: one configuration per process).
source "$P3D_DIR/common.sh"      # P3D_DIR: this directory, exported by submit_wulver3d.sh
PART="$OUTDIR/parts/${P3D_SOLVER}_ne${P3D_NE}_N${P3D_NOP}"
mkdir -p "$PART"
if [ -s "$PART/results.csv" ] && grep -q ",ok$" "$PART/results.csv"; then
    echo "already done: $PART"; exit 0     # requeued after preemption, or resubmitted
fi
export OPENBLAS_NUM_THREADS="$P3D_THREADS" OMP_NUM_THREADS="$P3D_THREADS"
echo "host $(hostname)  cpus ${SLURM_CPUS_PER_TASK:-?}  threads $P3D_THREADS  mem ${SLURM_MEM_PER_NODE:-?} MB"
srun --cpu-bind=cores $JULIA --project="$REPO" -t "$P3D_THREADS" tools/poisson3d_benchmark/bench3d.jl \
     --solver "$P3D_SOLVER" --ne "$P3D_NE" --nop "$P3D_NOP" --d "$D" --r "$R" \
     --ordering "$ORDERING" --rtol "$RTOL" --warmup "$WARMUP" --blas-threads "$P3D_THREADS" \
     --outdir "$PART"

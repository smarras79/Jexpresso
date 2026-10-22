#!/bin/bash -l
#---------------------------------------------------------------------------
# SWsphere trajectory ensemble on Wulver: one Slurm job per trajectory, all
# submitted at once as a JOB ARRAY, each on 12 MPI ranks.
#
#   bash problems/ShallowWater/SWsphere/submit_trajectories_wulver.sh
#
# Run it with BASH from the repository root on a login node, NOT with sbatch:
# it submits two jobs itself.
#
#   1. SWtraj-setup  (12 cores, one job)
#        MPIPreferences → Pkg.instantiate/precompile → a short 12-rank warm-up
#        run. The warm-up fills the package caches and the per-rank mesh/SEM
#        caches (.jexpresso_cache/ next to the case, keyed by rank-of-12), so
#        the array tasks only ever READ them. Without it, 25 jobs would compile
#        into the same depot and write the same cache files at the same time.
#
#   2. SWtraj        (job array, one task per trajectory, 12 cores each)
#        starts only if the setup job succeeded (--dependency=afterok). Task X
#        runs `run_trajectories.jl --only X` and writes <OUTDIR>/trajectory_X/.
#
# The trajectories are the ones run_trajectories.jl defines — by default the
# 10 trajectories picked from the 5 × 5 grid of amplitudes (40…200 m) and
# latitudes (30…60°): trajectory_1, 3, 7, 9, 13, 15, 16, 19, 22, 25. All 25
# with SWEEP_ARGS="--design grid" (the 10 already done are skipped). Each runs
# to tend = 6 days (SWEEP_ARGS="--tend <s>"
# changes it). Change the grid with SWEEP_ARGS (below); the
# array is sized from it automatically.
#
# KNOBS (environment variables, all optional):
#   SWEEP_ARGS     extra arguments for run_trajectories.jl, e.g.
#                  SWEEP_ARGS="--namp 7 --amp-halfwidth 90"   (default: none)
#                  SWEEP_ARGS="--amplitudes 100,120,140 --latitudes 40,50 --set μ=1.0e5"
#                  (no spaces inside a single value: SWEEP_ARGS is split on spaces)
#   OUTDIR         parent directory of trajectory_X/        (default: ./trajectories)
#   NP             MPI ranks per trajectory                  (default: 12)
#   TIME_LIMIT     wall time per trajectory                  (default: 02:00:00)
#   ARRAY          which trajectories to submit, Slurm syntax, e.g. "3,7,12-14"
#                  to rerun a few                            (default: all of them)
#   MAX_RUNNING    cap on simultaneously running trajectories (default: no cap)
#   SKIP_SETUP=1   submit the array without the setup job (caches already warm,
#                  e.g. when rerunning a few failed trajectories)
#   LAUNCHER       MPI launcher                     (default: "srun --mpi=pmi2")
#                  use LAUNCHER="mpirun" if srun cannot start MPICH here
#
# Examples:
#   bash problems/ShallowWater/SWsphere/submit_trajectories_wulver.sh
#   TIME_LIMIT=04:00:00 bash problems/ShallowWater/SWsphere/submit_trajectories_wulver.sh
#   SKIP_SETUP=1 ARRAY=4,17 bash problems/ShallowWater/SWsphere/submit_trajectories_wulver.sh
#
# A finished trajectory leaves trajectory_X/COMPLETED and is skipped if it is
# submitted again; a failed one leaves trajectory_X/FAILED and its array task
# ends FAILED in `sacct`. Logs: logs/SWtraj.<jobid>_<X>.out / .err
#---------------------------------------------------------------------------
#SBATCH --job-name=SWtraj
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000M

set -eo pipefail

SWEEP_SCRIPT=problems/ShallowWater/SWsphere/run_trajectories.jl

load_modules() {
    module load Julia/1.11.9
    module load GCC MPICH
    export JULIA_NUM_THREADS=1          # the parallelism is MPI ranks
}

#===========================================================================
# SUBMITTER: run on the login node with bash.
#===========================================================================
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    export JEXPRESSO_ROOT=${JEXPRESSO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
    cd "$JEXPRESSO_ROOT"
    load_modules

    export SWEEP_ARGS=${SWEEP_ARGS:-}
    export OUTDIR=$(realpath -m "${OUTDIR:-./trajectories}")
    export NP=${NP:-12}
    export LAUNCHER=${LAUNCHER:-"srun --mpi=pmi2"}
    TIME_LIMIT=${TIME_LIMIT:-02:00:00}
    mkdir -p logs "$OUTDIR"

    # Size the array from the sweep definition itself (--dry-run only parses
    # the arguments and prints the table; it does not load Jexpresso).
    if [[ -z "${ARRAY:-}" ]]; then
        # shellcheck disable=SC2086
        IDS=$(julia --project=. "$SWEEP_SCRIPT" $SWEEP_ARGS --dry-run | tail -n +2 | cut -d, -f1)
        [[ -n "$IDS" ]] || { echo "run_trajectories.jl --dry-run listed no trajectories" >&2; exit 1; }
        ARRAY=$(echo "$IDS" | paste -sd, -)      # ids need not be contiguous (--design subset)
    fi
    FIRST_ID=${ARRAY%%[-,]*}
    [[ -n "${MAX_RUNNING:-}" ]] && ARRAY="${ARRAY}%${MAX_RUNNING}"

    COMMON=(--account=smarras --partition=general --qos=standard
            --nodes=1 --ntasks="$NP" --cpus-per-task=1 --mem-per-cpu=4000M)

    DEPEND=()
    if [[ "${SKIP_SETUP:-0}" != 1 ]]; then
        SETUP_ID=$(sbatch --parsable "${COMMON[@]}" \
                          --job-name=SWtraj-setup --time=02:00:00 \
                          --output=logs/%x.%j.out --error=logs/%x.%j.err \
                          --export=ALL,STAGE=setup,WARMUP_ID="$FIRST_ID" \
                          "${BASH_SOURCE[0]}")
        echo "setup job:      $SETUP_ID"
        DEPEND=(--dependency=afterok:"$SETUP_ID")
    fi

    ARRAY_ID=$(sbatch --parsable "${COMMON[@]}" "${DEPEND[@]}" \
                      --job-name=SWtraj --time="$TIME_LIMIT" --array="$ARRAY" \
                      --output=logs/%x.%A_%a.out --error=logs/%x.%A_%a.err \
                      --export=ALL,STAGE=run \
                      "${BASH_SOURCE[0]}")
    echo "trajectory job: $ARRAY_ID  (array $ARRAY, $NP ranks each, time limit $TIME_LIMIT)"
    echo "output:         $OUTDIR/trajectory_<X>/"
    echo "monitor:        squeue -u \$USER ;  sacct -j $ARRAY_ID --format=JobID,State,Elapsed"
    exit 0
fi

#===========================================================================
# INSIDE A JOB.
#===========================================================================
: "${STAGE:?submit with bash (see the header of this script), not sbatch}"
cd "$JEXPRESSO_ROOT"
load_modules

# Load the caches built by the setup job, never write new ones from 12 ranks
# (or 25 × 12 ranks) at once. Older Julia without these flags just skips them.
JULIA_FLAGS=(--project=.)
if [[ "$STAGE" == run ]]; then
    export JULIA_PKG_PRECOMPILE_AUTO=0
    julia --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--compiled-modules=existing)
    julia --pkgimages=existing        -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--pkgimages=existing)
fi

# shellcheck disable=SC2086
mpi_launch() { $LAUNCHER -n "$NP" julia "${JULIA_FLAGS[@]}" "$@"; }

case "$STAGE" in
setup)
    export JULIA_PKG_PRECOMPILE_AUTO=1
    echo "--- 1. MPI preferences (before any compilation) ---"
    julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

    echo "--- 2. Serial precompile ---"
    julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

    echo "--- 3. ${NP}-rank warm-up: trajectory_${WARMUP_ID}, 1 simulated hour ---"
    WARMUP_DIR="$OUTDIR/.warmup"
    # shellcheck disable=SC2086
    mpi_launch "$SWEEP_SCRIPT" $SWEEP_ARGS --outdir "$WARMUP_DIR" \
               --only "$WARMUP_ID" --tend 3600 --force --stop-on-error
    rm -rf "$WARMUP_DIR"
    echo "--- setup complete; the trajectory array can start ---"
    ;;
run)
    X=$SLURM_ARRAY_TASK_ID
    echo "--- trajectory_${X} on ${NP} ranks ($(hostname), job ${SLURM_ARRAY_JOB_ID}_${X}) ---"
    # shellcheck disable=SC2086
    mpi_launch "$SWEEP_SCRIPT" $SWEEP_ARGS --outdir "$OUTDIR" --only "$X" --stop-on-error
    echo "--- trajectory_${X} done: $OUTDIR/trajectory_${X} ---"
    ;;
*)
    echo "unknown STAGE=$STAGE" >&2; exit 1
    ;;
esac

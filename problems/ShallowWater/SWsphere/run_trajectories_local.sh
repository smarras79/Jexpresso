#!/usr/bin/env bash
#---------------------------------------------------------------------------
# SWsphere trajectory ensemble on a workstation (e.g. a 10-core iMac): the
# local analogue of submit_trajectories_wulver.sh, without Slurm.
#
#   bash problems/ShallowWater/SWsphere/run_trajectories_local.sh
#
# Run it from anywhere; it works from the repository root. By default it runs
# the trajectories ONE AFTER ANOTHER, each on all NP = 10 MPI ranks, inside a
# single mpiexec launch — so Julia compiles once and every later trajectory
# starts immediately. Output: <OUTDIR>/trajectory_X/, as on Wulver.
#
# The trajectories are the ones run_trajectories.jl defines — by default
# trajectory_1 … trajectory_10 (see the header of that script). Same knobs
# as the Wulver script where they make sense:
#
#   SWEEP_ARGS   extra arguments for run_trajectories.jl, e.g.
#                SWEEP_ARGS="--design grid" for all 25, or "--tend 86400"
#                (no spaces inside a single value)         (default: none)
#   OUTDIR       parent directory of trajectory_X/         (default: ./trajectories)
#   NP           MPI ranks per trajectory                  (default: 10)
#   JOBS         trajectories run at the same time         (default: 1)
#                NP × JOBS should not exceed the core count; e.g. NP=5 JOBS=2
#   ONLY         run only these trajectory numbers, e.g. "3,7" (default: all)
#   JULIA        julia executable                          (default: julia)
#
# Examples:
#   bash problems/ShallowWater/SWsphere/run_trajectories_local.sh
#   NP=5 JOBS=2 bash problems/ShallowWater/SWsphere/run_trajectories_local.sh
#   ONLY=5 SWEEP_ARGS="--tend 21600" OUTDIR=./traj_test \
#       bash problems/ShallowWater/SWsphere/run_trajectories_local.sh
#
# MPI: the run is launched with the mpiexec that MPI.jl is configured for
# (MPI.mpiexec()), so it matches whatever route INSTALL.md §5 set up
# (system MPICH/OpenMPI or MPICH_jll). If MPI_Init fails on macOS, see
# INSTALL.md §5.5 (hostname) and §5.7 (FI_PROVIDER).
#
# Finished trajectories (trajectory_X/COMPLETED) are skipped on a rerun, so an
# interrupted sweep can simply be restarted. On macOS the machine is kept awake
# for the duration with `caffeinate`. Logs: logs/local_traj*.log
#---------------------------------------------------------------------------
set -eo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
cd "$ROOT"

SWEEP_SCRIPT=problems/ShallowWater/SWsphere/run_trajectories.jl
JULIA=${JULIA:-julia}
NP=${NP:-10}
JOBS=${JOBS:-1}
OUTDIR=${OUTDIR:-./trajectories}
mkdir -p logs "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

export JULIA_NUM_THREADS=1          # the parallelism is MPI ranks

# Keep a Mac from sleeping mid-run; a no-op elsewhere.
KEEPAWAKE=()
command -v caffeinate >/dev/null 2>&1 && KEEPAWAKE=(caffeinate -i)

NCORES=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 0)
if (( NCORES > 0 && NP * JOBS > NCORES )); then
    echo "WARNING: NP × JOBS = $((NP * JOBS)) ranks on $NCORES cores — the runs will oversubscribe the CPU." >&2
fi

# The trajectory numbers, from the sweep definition itself.
if [[ -n "${ONLY:-}" ]]; then
    IDS=$(echo "$ONLY" | tr ',' ' ')
else
    # shellcheck disable=SC2086
    IDS=$("$JULIA" --project=. "$SWEEP_SCRIPT" $SWEEP_ARGS --dry-run | tail -n +2 | cut -d, -f1 | tr '\n' ' ')
fi
[[ -n "${IDS// /}" ]] || { echo "no trajectories to run" >&2; exit 1; }
echo "--- trajectories: $IDS→ $OUTDIR  ($NP ranks each, $JOBS at a time) ---"

# mpiexec -n <np> julia --project=. <args…>, through MPI.jl's own mpiexec.
mpi_launch() {
    local np=$1; shift
    "${KEEPAWAKE[@]}" "$JULIA" --project=. -e '
        using MPI
        np = parse(Int, ARGS[1])
        run(`$(mpiexec()) -n $np $(Base.julia_cmd()) --project=. $(ARGS[2:end])`)' \
        "$np" "$@"
}

echo "--- 1. Precompile (once; fast when already done) ---"
"$JULIA" --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

if (( JOBS <= 1 )); then
    #-----------------------------------------------------------------------
    # One launch, every trajectory in turn, compiled once. A failed trajectory
    # is logged (trajectory_X/FAILED) and the sweep moves on to the next.
    #-----------------------------------------------------------------------
    echo "--- 2. Running the sweep (log: logs/local_traj.log) ---"
    # shellcheck disable=SC2086
    mpi_launch "$NP" "$SWEEP_SCRIPT" $SWEEP_ARGS --outdir "$OUTDIR" \
               --only "$(echo $IDS | tr ' ' ',')" 2>&1 | tee logs/local_traj.log || true
else
    #-----------------------------------------------------------------------
    # JOBS concurrent launches, the trajectories dealt round-robin. A short
    # warm-up first writes the per-rank mesh/SEM caches, so the concurrent
    # runs only read them instead of all writing the same files at once.
    #-----------------------------------------------------------------------
    FIRST=${IDS%% *}
    echo "--- 2. ${NP}-rank warm-up: trajectory_${FIRST}, 1 simulated hour (log: logs/local_warmup.log) ---"
    # shellcheck disable=SC2086
    mpi_launch "$NP" "$SWEEP_SCRIPT" $SWEEP_ARGS --outdir "$OUTDIR/.warmup" \
               --only "$FIRST" --tend 3600 --force --stop-on-error > logs/local_warmup.log 2>&1 \
        || { echo "warm-up FAILED — see logs/local_warmup.log"; tail -40 logs/local_warmup.log; exit 1; }
    rm -rf "$OUTDIR/.warmup"

    echo "--- 3. Running the sweep, $JOBS at a time (logs: logs/local_traj_<n>.log) ---"
    PIDS=()
    for (( g = 0; g < JOBS; g++ )); do
        GROUP=""; i=0
        for id in $IDS; do
            (( i % JOBS == g )) && GROUP="${GROUP:+$GROUP,}$id"
            i=$((i + 1))
        done
        [[ -n "$GROUP" ]] || continue
        echo "    launch $((g + 1)): trajectories $GROUP"
        # shellcheck disable=SC2086
        mpi_launch "$NP" "$SWEEP_SCRIPT" $SWEEP_ARGS --outdir "$OUTDIR" \
                   --only "$GROUP" > "logs/local_traj_$((g + 1)).log" 2>&1 &
        PIDS+=($!)
    done
    for pid in "${PIDS[@]}"; do wait "$pid" || true; done
fi

#---------------------------------------------------------------------------
# Summary, from the markers run_trajectories.jl leaves in each directory.
#---------------------------------------------------------------------------
echo "--- summary ($OUTDIR) ---"
NFAIL=0
for id in $IDS; do
    d="$OUTDIR/trajectory_$id"
    if   [[ -f "$d/COMPLETED" ]]; then echo "    trajectory_$id  COMPLETED"
    elif [[ -f "$d/FAILED"    ]]; then echo "    trajectory_$id  FAILED: $(head -c 200 "$d/FAILED")"; NFAIL=$((NFAIL + 1))
    else                               echo "    trajectory_$id  did not finish"; NFAIL=$((NFAIL + 1))
    fi
done
echo "compare them with:"
echo "    julia --project=. problems/ShallowWater/SWsphere/compare_trajectories.jl $OUTDIR"
exit $(( NFAIL > 0 ))

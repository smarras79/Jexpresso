#!/bin/bash -l
# ============================================================================
# submit_Jexpresso_precompile_pass.sh
#
# Large LES run on many ranks, with compilation paid for ONCE and up front
# instead of inside the production time loop.
#
# There are three distinct "compile" stages in a Julia MPI job, and they are
# routinely confused. This script deals with all three, in order:
#
#   (A) PACKAGE PRECOMPILATION  -- Pkg.precompile(), writes .ji / .so caches
#       into the depot. Must happen ONCE, SERIALLY, before any rank starts:
#       1536 ranks discovering the same missing cache file all try to build
#       it, into the same shared depot, at the same time. That is minutes to
#       hours of wall time, lock contention and half-written cache files.
#       -> phases 1-3 below, one process, no mpirun.
#
#   (B) LOAD-TIME WORK          -- `using Jexpresso` on every rank. Cheap only
#       if (A) already produced valid caches, which is why (A) runs first and
#       why the ranks are then forbidden from writing to the depot at all.
#
#   (C) JIT OF THE SOLVER       -- rhs!, the SciML integrator specialised on
#       the real CallbackSet, the MPI halo exchange, plus the first-touch
#       allocation of every per-rank working array. This one CANNOT be cached
#       across processes: it happens inside the run, on the first timestep.
#       Launched all at once on a big problem, it happens interleaved with
#       collective communication -- fast ranks block in MPI_Wait while slow
#       ranks compile, the first few hundred steps report a meaningless step
#       rate, and the long run inherits a heap full of compiler garbage.
#       -> JEXPRESSO_PRECOMPILE_PASS=1 below splits the run in two:
#             phase 1: ONE timestep, real problem / real callbacks / real
#                      kwargs -> everything compiles and allocates here
#             GC.gc() + MPI.Barrier
#             phase 2: the production solve, RESUMING from phase 1's state,
#                      reusing phase 1's specialisation, hot from step two.
#          Phase 1's step is kept -- it IS the simulation's first step.
#          See ENVIRONMENT_VARIABLES.md -> JEXPRESSO_PRECOMPILE_PASS.
#
# Submit with:
#   sbatch submit_Jexpresso_precompile_pass.sh
#   sbatch --export=ALL,CASE=LESICP2-coarse submit_Jexpresso_precompile_pass.sh
# ============================================================================

#SBATCH --job-name=LESsmago
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=71:59:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

JEXPRESSO_ROOT="${JEXPRESSO_ROOT:-/project/smarras/smarras/Jexpresso}"
EQS="${EQS:-CompEuler}"
CASE="${CASE:-LESICP2-coarse}"

cd "$JEXPRESSO_ROOT" || exit 1

# ---------------------------------------------------------------------------
# Threads. Jexpresso is essentially pure MPI: one core per rank. Without this,
# 128 ranks per node each spawn 128 BLAS/OpenMP threads and fight over the same
# 128 cores -- the classic "more nodes, slower run".
# ---------------------------------------------------------------------------
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# ===========================================================================
# (A) SERIAL SETUP -- exactly one process, no mpirun anywhere in this section.
# ===========================================================================
export JULIA_PKG_PRECOMPILE_AUTO=1   # allowed during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
# MPIPreferences rewrites LocalPreferences.toml, which INVALIDATES the
# precompile caches of MPI.jl and everything downstream of it. Doing this after
# step 2 would throw away everything step 2 just built.
julia --project=. --startup-file=no \
    -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. --startup-file=no \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
# This is a real check, not a no-op: it loads Jexpresso exactly the way a rank
# will and fails HERE, on one process with a readable backtrace, rather than on
# 1536 ranks at once. `using Jexpresso` hits the cache written in step 2; it
# does not run a case (src/Jexpresso.jl only auto-includes run.jl when it is
# the PROGRAM_FILE).
julia --project=. --startup-file=no -e '
    t = @elapsed (using MPI; using Jexpresso)
    println("    Jexpresso loaded in ", round(t, digits=2), " s")' || {
    echo "ERROR: Jexpresso failed to load serially. Fix that before launching ranks." >&2
    exit 1
}

# ===========================================================================
# (B) PARALLEL LAUNCH
# ===========================================================================

# Guard rails: nothing a rank does may write to the shared depot.
#   JULIA_PKG_PRECOMPILE_AUTO=0   Pkg never kicks off a precompile pass
#   --compiled-modules=existing   load existing .ji caches, never create them
#   --pkgimages=existing          same, for the native-code images
# A rank that somehow still misses a cache now compiles privately in its own
# memory instead of racing 1535 others for one file on a shared filesystem.
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=1

# (C) The two-phase run: 1 compile step, then the simulation resumes from it.
# Set to 0 to fall back to the historical single-phase launch.
export JEXPRESSO_PRECOMPILE_PASS="${JEXPRESSO_PRECOMPILE_PASS:-1}"

# Progress lines during the solve. Diagnostics for this case are thousands of
# seconds apart; without this, a healthy run and a hung one look identical.
export JEXPRESSO_STEP_HEARTBEAT="${JEXPRESSO_STEP_HEARTBEAT:-1}"

JULIA_FLAGS=(--project=. --startup-file=no)
# Probe rather than assume: `existing` was added to these flags in Julia 1.11,
# and an unrecognised flag would kill the launch on every rank.
if julia --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--compiled-modules=existing)
fi
if julia --pkgimages=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--pkgimages=existing)
fi

NTASKS="${SLURM_NTASKS:-64}"

echo "--- Setup complete, launching $NTASKS ranks ---"
echo "    case                      : $EQS / $CASE"
echo "    JEXPRESSO_PRECOMPILE_PASS : $JEXPRESSO_PRECOMPILE_PASS"
echo "    julia flags               : ${JULIA_FLAGS[*]}"
echo "    started                   : $(date)"

# MPICH/Hydra propagates the environment to every rank by default (-genvall).
# Under OpenMPI the exports above do NOT propagate -- pass them explicitly:
#   mpirun -x JULIA_PKG_PRECOMPILE_AUTO -x JEXPRESSO_PRECOMPILE_PASS ...
mpirun -np "$NTASKS" julia "${JULIA_FLAGS[@]}" src/Jexpresso.jl "$EQS" "$CASE"
rc=$?

echo "--- Finished: $(date) ---"

# Report the launcher's exit status. Without this the script printed
# "Finished" whether the run completed or every rank was killed, so an
# out-of-memory kill during the mesh read -- which produces no Julia
# backtrace at all, just a dead job -- was indistinguishable from success.
if [ "$rc" -ne 0 ]; then
    # Report on BOTH streams: stderr alone lands in the .err file while
    # "--- Finished ---" goes to .out, so reading .out on its own showed a run
    # that appeared to simply stop.
    #
    # Do NOT guess at the cause. An earlier version printed the empty-rank
    # explanation unconditionally on every failure, which reads as a diagnosis
    # and sent at least one debugging session down the wrong path. Instead,
    # look for the one symptom we can actually detect from here, and otherwise
    # say plainly that the cause is in the .err file.
    log="${SLURM_JOB_NAME:-job}.${SLURM_JOB_ID:-0}.out"
    for stream in /dev/stdout /dev/stderr; do
        {
            echo "--- FAILED: mpirun exited $rc ---"
            echo "    The run did NOT complete."
            echo "    The actual error is in the .err file for job ${SLURM_JOB_ID:-<id>}."
            if [ -f "$log" ] && grep -q "Min elements: 0" "$log"; then
                echo
                echo "    CONFIRMED in this run's output: 'Min elements: 0'."
                echo "    Some ranks own no elements, which is fatal downstream."
                echo "    ntasks exceeds what the mesh can be split into. Note that"
                echo "    the split is aspect-ratio aware: a non-square domain does"
                echo "    NOT give you nelemx*nelemy usable ranks. Check the nx x ny"
                echo "    the run chose against the element counts in each direction."
            fi
        } > "$stream"
    done
fi
exit $rc

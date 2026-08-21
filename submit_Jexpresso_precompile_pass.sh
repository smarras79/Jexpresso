#!/bin/bash -l

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
    # CRITICAL: write with >&2 / plain echo, NEVER "> /dev/stdout" or
    # "> /dev/stderr".
    #
    # Under SLURM, fd 1 and fd 2 are regular files (%x.%j.out / .err). On Linux
    # /dev/stdout is a symlink to /proc/self/fd/1, so "> /dev/stdout" RE-OPENS
    # that file with O_TRUNC -- it discards everything the job had written and
    # starts again at offset zero. A previous version of this handler did
    # exactly that to both streams, so every failed run destroyed its own log
    # and left only this message behind. The Julia backtrace, the mesh
    # diagnostics, the whole run: gone, and the failure looked like a job that
    # had produced no output at all. ("tail: file truncated" was the tell.)
    #
    # Duplicating the existing descriptor appends instead of truncating, which
    # is the whole point: this handler must never be able to lose evidence.
    _report() {
        echo "--- FAILED: mpirun exited $rc ---"
        echo "    The run did NOT complete."
        echo "    Julia's error and backtrace are ABOVE this line in the .err file"
        echo "    for job ${SLURM_JOB_ID:-<id>} (stderr from every rank lands there)."
        if [ -f "$log" ] && grep -q "Min elements: 0" "$log"; then
            echo
            echo "    CONFIRMED in this run's output: 'Min elements: 0'."
            echo "    Some ranks own no elements, which is fatal downstream."
            echo "    ntasks exceeds what this mesh can be split into. The split is"
            echo "    aspect-ratio aware, so nelemx*nelemy is an upper bound a"
            echo "    non-square domain does not reach. Run tools/pick_nranks.jl."
        fi
    }
    log="${SLURM_JOB_NAME:-job}.${SLURM_JOB_ID:-0}.out"
    _report          # -> .out, appended
    _report >&2      # -> .err, appended
fi
exit $rc

#!/bin/bash -l
#
# submit_Jexpresso_profile.sh -- submit_Jexpresso_precompile_pass.sh, unchanged
# except for a 200-step PROBE that measures where the time goes. Everything
# else in this file is the production script verbatim, comments included.
#
# THE DIFFERENCE, in full:
#   JEXPRESSO_HEVI_PROFILE=1   print the per-step cost breakdown
#   DBG_TEND=100.0             200 steps at the deck's dt = 0.5 s, then stop
#   CASE                       the production deck
#   SBATCH nodes / time / mem  production rank count, short walltime
#
# WHAT IT PRINTS
# --------------
# Every 50 steps (the first 10 discarded as JIT), rank 0 prints
#
#    | rhs!            4.00 /step   0.0xxxx s each   ...   xx.x %
#    | f_imp           4.00 /step   0.0xxxx s each   ...   xx.x %
#    | column solve    3.00 /step   0.0xxxx s each   ...   xx.x %
#    | refactorise     0.20 /step   0.0xxxx s each   ...   xx.x %
#    | accounted / measured step                     0.xxxx s/step
#
# For an IMEX3D run "column solve" IS the whole GMRES stage solve and "f_imp"
# is the 3D acoustic operator. The number the production plan turns on is
#
#     rho = T_rhs / A_vert = (rhs! s each) / ( (f_imp s each) / 3 )
#
# -- the 3 because the 3D operator does three derivative sweeps where the
# vertical operator does one. rho picks the column of the projection table in
# the deck header:
#
#     rho      10     15     20     30     50
#     hours   7.7    5.4    4.2    3.1    2.1
#
# Also worth reading off the same log:
#   * the setup report's "measured: N iterations ... from a COLD start", and
#     the :imex_monitor running average. Far above ~31 means raise
#     :imex_restart (DBG_RESTART) or lower DBG_DT.
#   * "measured step" minus "accounted" -- time outside rhs!/f_imp/solve, i.e.
#     filtering, diagnostics or MPI wait.
#
# RUN IT TWICE, AT TWO RANK COUNTS. Two five-minute jobs settle the scaling
# question that no estimate can: submit with --nodes=8 and again with
# --nodes=16 and compare "measured step". IMEX3D does ~93 halo exchanges and
# ~280 MPI reductions per step against the explicit scheme's handful, so it is
# far more sensitive to thin ranks than the explicit run this case is being
# migrated from. Do not assume the rank count that suited that run suits this
# one.

#SBATCH --job-name=LESprofile
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras
#
# 1024 ranks = a 32 x 32 rank grid over the 64 x 64 element columns, so each
# rank owns a 2 x 2 block of columns and 15 545 points -- the closest valid
# column decomposition to the per-rank load the explicit 1240-core run already
# carries. 1240 itself is NOT valid here: :lxy_partition, which HEVI and
# IMEX3D both require, cuts the mesh into vertical columns and the admissible
# rank counts are the ones this lists:
#
#     julia tools/pick_nranks.jl 64 64 60 4 2048 10240 10240
#
# --nodes=16 for the 2048-rank arm of the scaling comparison above.
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2000M

module load Julia/1.11.9
module load GCC MPICH

JEXPRESSO_ROOT="${JEXPRESSO_ROOT:-/project/smarras/smarras/Jexpresso}"
EQS="${EQS:-CompEuler}"
CASE="${CASE:-LESICP2-64x64x60-imex}"

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

# ---------------------------------------------------------------------------
# THE PROBE -- the only functional difference from the production script.
# ---------------------------------------------------------------------------
# Rank 0 only. _EVERY and _SKIP are the defaults, written out so the cadence is
# visible in the log: report every 50 steps, discard the first 10 because they
# are compiling rather than computing. JEXPRESSO_HEVI_PROFILE_ALLRANKS=1 prints
# every rank instead, which is how to find ONE slow rank -- worth it on a rank
# count pick_nranks flags as imbalanced, unnecessary at 1024 where it is exact.
export JEXPRESSO_HEVI_PROFILE=1
export JEXPRESSO_HEVI_PROFILE_EVERY=50
export JEXPRESSO_HEVI_PROFILE_SKIP=10

# 200 steps at the deck's dt = 0.5 s. The deck reads DBG_TEND, and when it is
# set it also turns off the ~640 MB initial VTK write -- minutes of I/O for a
# run whose output is thrown away. It does not distort s/step (it happens
# before the time loop), it just wastes the allocation.
export DBG_TEND="${DBG_TEND:-100.0}"

JULIA_FLAGS=(--project=. --startup-file=no)
# Probe rather than assume: `existing` was added to these flags in Julia 1.11,
# and an unrecognised flag would kill the launch on every rank.
if julia --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--compiled-modules=existing)
fi
if julia --pkgimages=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--pkgimages=existing)
fi

NTASKS="${SLURM_NTASKS:-1024}"

echo "--- Setup complete, launching $NTASKS ranks (PROFILE PROBE) ---"
echo "    case                      : $EQS / $CASE"
echo "    DBG_TEND                  : $DBG_TEND s  (200 steps at dt = 0.5 s)"
echo "    JEXPRESSO_HEVI_PROFILE    : $JEXPRESSO_HEVI_PROFILE  (every $JEXPRESSO_HEVI_PROFILE_EVERY steps, first $JEXPRESSO_HEVI_PROFILE_SKIP skipped)"
echo "    JEXPRESSO_PRECOMPILE_PASS : $JEXPRESSO_PRECOMPILE_PASS"
echo "    julia flags               : ${JULIA_FLAGS[*]}"
echo "    started                   : $(date)"

# MPICH/Hydra propagates the environment to every rank by default (-genvall).
# Under OpenMPI the exports above do NOT propagate -- pass them explicitly:
#   mpirun -x JULIA_PKG_PRECOMPILE_AUTO -x JEXPRESSO_PRECOMPILE_PASS \
#          -x JEXPRESSO_HEVI_PROFILE -x DBG_TEND ...
# Under OpenMPI without those, the probe silently prints nothing and the run is
# a full-length production run: DBG_TEND would not reach the ranks either.
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
            echo "    non-square domain does not reach. For this grid:"
            echo "      julia tools/pick_nranks.jl 64 64 60 4 2048 10240 10240"
        fi
    }
    log="${SLURM_JOB_NAME:-job}.${SLURM_JOB_ID:-0}.out"
    _report          # -> .out, appended
    _report >&2      # -> .err, appended
fi
exit $rc

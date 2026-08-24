#!/bin/bash -l
#
# submit_Jexpresso_profile.sh -- submit_Jexpresso_precompile_pass.sh plus a
# 200-step PROBE that measures where the time goes, and the launcher work the
# jump to multi-node turned out to need.
#
# THE DIFFERENCE, in full:
#   JEXPRESSO_HEVI_PROFILE=1   print the per-step cost breakdown
#   DBG_TEND=100.0             200 steps at the deck's dt = 0.5 s, then stop
#   CASE                       the production deck
#   SBATCH nodes / time / mem  production rank count, short walltime
#   launcher + preflight       srun on multi-node, and a 10-second MPI_Init
#                              check at the full rank count before the run
#
# START ON THE LADDER, NOT AT THE TOP. The script this was copied from had
# only ever run at --nodes=1 --ntasks-per-node=32, which never touches the
# network: MPICH keeps all 32 ranks in shared memory. Going straight to
# 8 x 128 changed both variables at once and died inside MPI_Init_thread,
# which cost a queue turnaround to attribute. Submit in this order and keep
# the first that works:
#
#     --nodes=1 --ntasks-per-node=32      32 ranks, one node, no network
#     --nodes=8 --ntasks-per-node=32      256 ranks -- the DEFAULT below, and
#         --cpus-per-task=2               what tools/pick_nranks.jl recommends
#     --nodes=8 --ntasks-per-node=128     1024 ranks, the scaling arm
#
# Every rung is a valid :lxy_partition rank count for this 64 x 64 mesh; the
# ladder is about the MPI stack, not the decomposition. The preflight below
# fails in seconds rather than after the mesh read, and its error message says
# which rung failed and what to change.
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
# 256 ranks = a 16 x 16 rank grid over the 64 x 64 element columns, so each
# rank owns a 4 x 4 block of columns and 62 179 points. This is what the repo's
# own tool picks for this mesh:
#
#     julia tools/pick_nranks.jl 64 64 60 4 2048 10240 10240
#     ...
#     256    16 x 16   4 x 4   16   62179   halo/elem 1.0   <== RECOMMENDED
#     1024   32 x 32   2 x 2    4   15545   halo/elem 2.0   (thin: comms-bound)
#
# 1240 -- the explicit run's core count -- is NOT valid here at all:
# :lxy_partition, which HEVI and IMEX3D both require, cuts the mesh into
# vertical columns, and only the counts that tool lists divide cleanly. Anything
# else leaves ranks owning zero elements.
#
# MEMORY PER RANK WHEN --mem-per-cpu IS CAPPED.
#
# Wulver enforces MaxMemPerCPU = 4000M, so --mem-per-cpu=8000M is rejected. That
# caps the memory per CPU, NOT the memory per RANK: --cpus-per-task=2 allocates
# two CPUs to each MPI rank, and the rank's cgroup budget is the sum, 8000M.
# Asking for fewer ranks per node does NOT help -- fewer tasks means fewer CPUs
# means proportionally less memory, and 4 GB per rank either way. cpus-per-task
# is the only lever the cap leaves.
#
#   nodes x ntasks-per-node x cpus-per-task = 8 x 32 x 2
#     -> 256 ranks, 64 CPUs and 256 GB per node, 8 GB per rank
#
# The second CPU of each pair sits idle: JULIA_NUM_THREADS/OMP_NUM_THREADS are 1
# below and Jexpresso is pure MPI. It is bought for its memory, and it doubles
# the core-hours charged. If the run comes up comfortably inside 8 GB (watch the
# banner and `seff` afterwards), drop back to --cpus-per-task=1 and halve that.
#
# 4 GB per rank died of OutOfMemoryError during setup, but that was WITHOUT
# --heap-size-hint, which is the actual bug (see the flag's block below). It is
# worth retrying at --cpus-per-task=1 once the hint is in place; the mesh
# broadcast is the peak and 15.9M global gridpoints may well fit in 4 GB when
# the GC is collecting against the right limit.
#
# The scaling arm of the comparison in the header is 1024 ranks
# (--nodes=8 --ntasks-per-node=128 --cpus-per-task=1). Run it SECOND, once 256
# has worked: at 128 tasks x 4000M it is 512 GB per node, so it needs nodes that
# large and it gets only 4 GB per rank.
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4000M

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

echo "--- 1b. Syntax check (1 second, before paying 133 for a precompile) ---"
# A one-character typo in src/ costs a full precompile (133 s) and then fails
# on EVERY rank -- that is how a stray `$` in krylov.jl burned a 256-rank job.
# Note that `Meta.parse` alone would not have caught it: `$` outside a string
# parses cleanly and only fails when the expression is LOWERED. This checks
# both, on every .jl file, without loading the package or needing MPI.
julia --startup-file=no tools/syntax_check.jl src test problems tools || {
    echo "ERROR: syntax error in the source tree -- fix it before submitting." >&2
    echo "       Nothing was compiled and no ranks were launched." >&2
    exit 1
}

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

# ===========================================================================
# --heap-size-hint -- THE ONE FLAG WITHOUT WHICH THIS RUN DIES OF OOM.
#
# Julia sizes its GC heap target from /proc/meminfo, i.e. from the NODE's RAM.
# It does not read the SLURM cgroup. So on a 512 GB node every rank believes
# it has ~400 GB to play with, while the cgroup gives it --mem-per-cpu. The GC
# therefore has no reason to collect until long past the real cap, and the
# first allocation that crosses it fails.
#
# That is what killed the 256-rank run, and the traceback is deceptive:
#
#     OutOfMemoryError()
#      [8] GMRESWorkspace(npoin=69649, nimp=5; m=30, ...)  krylov.jl:174
#
# GMRESWorkspace at that size is 90 MB (31 Arnoldi vectors + 3 scratch, each
# 69649 x 5 x 8 B). Nothing asks for 90 MB and fails on a rank holding 4 GB --
# unless the heap is already at the cap with garbage the GC never collected.
# The workspace is the last straw, not the load. Do not "fix" this by lowering
# :imex_restart; m=30 -> m=20 saves 26 MB.
#
# The hint is set to 85% of the cgroup budget, leaving room for the parts that
# are not Julia heap: the MPI library, UCX buffers, and the LAPACK banded
# factors.
# ===========================================================================
if [ -n "$SLURM_MEM_PER_CPU" ]; then
    RANK_MB=$(( SLURM_MEM_PER_CPU * ${SLURM_CPUS_PER_TASK:-1} ))
elif [ -n "$SLURM_MEM_PER_NODE" ] && [ -n "$SLURM_NTASKS_PER_NODE" ]; then
    RANK_MB=$(( SLURM_MEM_PER_NODE / SLURM_NTASKS_PER_NODE ))
else
    RANK_MB=0
fi
if [ "$RANK_MB" -gt 0 ] && julia --heap-size-hint=1G -e 'exit(0)' >/dev/null 2>&1; then
    HEAP_MB=$(( RANK_MB * 85 / 100 ))
    JULIA_FLAGS+=(--heap-size-hint=${HEAP_MB}M)
fi

NTASKS="${SLURM_NTASKS:-256}"
NNODES="${SLURM_NNODES:-1}"

# ===========================================================================
# LAUNCHER -- srun for multi-node, mpirun for single-node.
#
# WHY THIS EXISTS. The first multi-node run of this script died with
#
#     Abort(15) on node 814: Fatal error in internal_Init_thread: Other MPI error
#     [proxy:3@n0031] HYD_pmcd_pmip_control_cmd_cb: assert (!closed) failed
#
# -- a failure inside MPI_Init_thread, i.e. before a single line of Jexpresso
# runs. The precursor script this file was copied from had only ever been run
# at --nodes=1 --ntasks-per-node=32, where MPICH never touches the network at
# all: every rank is on one node, so ch4:ucx uses shared memory and the whole
# inter-node path is untested. This script was the first launch to use it.
#
# `mpirun` here is MPICH's Hydra, which under SLURM spawns one proxy per node
# and then runs its OWN PMI wire-up inside the allocation. The `assert
# (!closed)` lines are those proxies losing their control socket to a rank that
# has already aborted. Hydra-inside-SLURM is the fragile combination; srun
# talks to SLURM's PMI directly and is what the site supports.
#
# Override either way with JEXPRESSO_LAUNCHER=srun|mpirun. If srun rejects the
# PMI plugin, the list of ones this cluster actually has is printed below --
# set SRUN_MPI to one of them (pmix and pmi2 are the usual names).
# ===========================================================================
SRUN_MPI="${SRUN_MPI:-pmi2}"
if [ -n "$JEXPRESSO_LAUNCHER" ]; then
    LAUNCHER="$JEXPRESSO_LAUNCHER"
elif [ "$NNODES" -gt 1 ] && command -v srun >/dev/null 2>&1; then
    LAUNCHER="srun"
else
    LAUNCHER="mpirun"
fi

# -c is passed EXPLICITLY. Since Slurm 22.05 srun no longer inherits the batch
# job's --cpus-per-task on every site configuration, and without it each rank
# gets one CPU and one CPU's worth of memory -- silently undoing the
# cpus-per-task trick above and reinstating the OOM.
launch() {   # launch <julia args...>
    if [ "$LAUNCHER" = "srun" ]; then
        srun --mpi="$SRUN_MPI" -n "$NTASKS" -c "${SLURM_CPUS_PER_TASK:-1}" julia "$@"
    else
        mpirun -np "$NTASKS" julia "$@"
    fi
}

echo "--- MPI environment ---"
if [ "$LAUNCHER" = "srun" ]; then
    echo "    launcher                  : srun --mpi=$SRUN_MPI"
else
    echo "    launcher                  : mpirun (MPICH Hydra)"
fi
echo "    nodes x tasks x cpus      : $NNODES x $((NTASKS / NNODES)) x ${SLURM_CPUS_PER_TASK:-1} = $NTASKS ranks"
if [ "${RANK_MB:-0}" -gt 0 ]; then
    echo "    cgroup budget per rank    : ${RANK_MB} MB"
    echo "    --heap-size-hint          : ${HEAP_MB} MB  (Julia reads /proc/meminfo, NOT the cgroup)"
else
    echo "    cgroup budget per rank    : unknown (no SLURM_MEM_PER_CPU / _PER_NODE)"
    echo "    --heap-size-hint          : NOT SET -- Julia will size its heap from the"
    echo "                                whole node and can OOM against a smaller cgroup"
fi
echo "    node RAM (this node)      : $(awk '/MemTotal/{printf "%.0f MB", $2/1024}' /proc/meminfo)"
srun --mpi=list 2>&1 | sed 's/^/    srun --mpi=list: /' || true

# ---------------------------------------------------------------------------
# PREFLIGHT -- get MPI_Init_thread to succeed at the full rank count BEFORE
# paying for the mesh read and the operator build.
#
# This is seconds of walltime and it separates the two failures that look
# identical in a SLURM log: "the MPI stack cannot start 1024 ranks on this
# many nodes" and "Jexpresso failed". The original failure above was entirely
# the former, and it cost a full queue turnaround to learn that.
# ---------------------------------------------------------------------------
echo "--- Preflight: MPI_Init_thread + one Allreduce on $NTASKS ranks ---"
launch "${JULIA_FLAGS[@]}" -e '
    using MPI
    MPI.Init()
    c = MPI.COMM_WORLD
    n = MPI.Comm_size(c)
    s = MPI.Allreduce(1, MPI.SUM, c)
    if MPI.Comm_rank(c) == 0
        println("    MPI OK: ", n, " ranks, Allreduce = ", s,
                s == n ? "  (consistent)" : "  *** INCONSISTENT ***")
    end
    MPI.Barrier(c)
    MPI.Finalize()' || {
    rc=$?
    echo "ERROR: MPI itself cannot start $NTASKS ranks across $NNODES nodes (exit $rc)." >&2
    echo "       Nothing in Jexpresso has run. Bisect the two variables that" >&2
    echo "       changed from the known-good 1 x 32 launch, one at a time:" >&2
    echo "         1. --nodes=1 --ntasks-per-node=32    (one node, no network)" >&2
    echo "         2. --nodes=8 --ntasks-per-node=32    (adds the network)" >&2
    echo "         3. --nodes=8 --ntasks-per-node=128   (more ranks per node)" >&2
    echo "       If 1 passes and 2 fails, it is the inter-node transport:" >&2
    echo "       try JEXPRESSO_LAUNCHER=srun with SRUN_MPI from the list above," >&2
    echo "       or export UCX_TLS=self,sm,ud (ud scales further than rc at" >&2
    echo "       high rank counts; rc exhausts queue-pair memory)." >&2
    echo "       If 1 also fails, it is not the network at all: check /dev/shm" >&2
    echo "       room and that --mem-per-cpu x ntasks-per-node fits the node." >&2
    exit $rc
}

echo "--- Setup complete, launching $NTASKS ranks (PROFILE PROBE) ---"
echo "    case                      : $EQS / $CASE"
echo "    DBG_TEND                  : $DBG_TEND s  (200 steps at dt = 0.5 s)"
echo "    JEXPRESSO_HEVI_PROFILE    : $JEXPRESSO_HEVI_PROFILE  (every $JEXPRESSO_HEVI_PROFILE_EVERY steps, first $JEXPRESSO_HEVI_PROFILE_SKIP skipped)"
echo "    JEXPRESSO_PRECOMPILE_PASS : $JEXPRESSO_PRECOMPILE_PASS"
echo "    julia flags               : ${JULIA_FLAGS[*]}"
echo "    started                   : $(date)"

# ENVIRONMENT PROPAGATION. srun exports the submitting environment to every
# task by default, and MPICH/Hydra does the same (-genvall), so the exports
# above reach the ranks under either launcher here. Under OpenMPI's mpirun they
# would NOT -- pass them explicitly:
#   mpirun -x JULIA_PKG_PRECOMPILE_AUTO -x JEXPRESSO_PRECOMPILE_PASS \
#          -x JEXPRESSO_HEVI_PROFILE -x DBG_TEND ...
# Without them the probe silently prints nothing and the run is a full-length
# production run: DBG_TEND would not reach the ranks either.
launch "${JULIA_FLAGS[@]}" src/Jexpresso.jl "$EQS" "$CASE"
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

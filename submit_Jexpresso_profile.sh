#!/bin/bash -l
#
# submit_Jexpresso_profile.sh -- submit_Jexpresso_precompile_pass.sh plus a
# 200-step PROBE that measures where the time goes, and the launcher work the
# jump to multi-node turned out to need.
#
# THE DIFFERENCE, in full:
#   JEXPRESSO_HEVI_PROFILE=1   print the per-step cost breakdown
#   DBG_SCHUR                  which stage solve to run -- THE A/B, see below
#   DBG_TEND=45.0              75 steps at the deck's dt = 0.6 s, then stop
#   DBG_RTOL / DBG_RESTART     pinned so both arms solve to the same tolerance
#   CASE                       CompEuler/rtb3d_schur
#   SBATCH nodes / time / mem  25 ranks on one node, short walltime
#   launcher + preflight       srun on multi-node, and a 10-second MPI_Init
#                              check at the full rank count before the run
#
# THE A/B THIS EXISTS FOR. Submit it TWICE, changing ONE thing:
#
#     DBG_SCHUR=0 sbatch submit_Jexpresso_profile.sh    # five fields, 5*Np
#     DBG_SCHUR=1 sbatch submit_Jexpresso_profile.sh    # scalar Schur, Np
#
# Same rank count, same rtol, same restart, same step count -- all pinned below
# so neither arm can drift. The deck writes to output_imex_full/ and
# output_imex_schur/, so the second does not overwrite the first.
#
# WHAT TO COMPARE, in order of how much it tells you:
#   1. "measured step" from the profile block -- the number the whole exercise
#      is about.
#   2. The profile SPLIT. On the 64x64x60 baseline the stage solve was 87.4% of
#      the step, of which matvec 46.9%, banded 16.4%, gather/scatter 2.2%,
#      orthogonalise 23.5% and MPI reduce 6.8%. Everything but that last one
#      scales with the implicit field count, so the split is what shows WHERE
#      the saving comes from rather than that there is one.
#   3. Iterations/solve, from the setup report and :imex_monitor.
#
# AND EXPECT THE ITERATION COUNT TO GO THE WRONG WAY. Measured on the small
# variant of this deck (10x10x40, one rank): the scalar system took 61 cold
# iterations against the five-field system's 20 -- 3x MORE -- and was still
# 1.21x faster per step, because one implicit field instead of five makes each
# iteration much cheaper. A mock sweep had predicted 2.67x FEWER iterations at
# this anisotropy; it inverted on a real mesh. So a rise in iterations/solve
# under DBG_SCHUR=1 is the expected result here, not a fault -- read the step
# time, not the count.
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
# KEEP BOTH ARMS AT THE SAME RANK COUNT. IMEX3D does ~93 halo exchanges and
# ~280 MPI reductions per step, and the MPI reduce is precisely the part of the
# stage solve the Schur reduction CANNOT shrink. Comparing arms across rank
# counts therefore changes the one term that does not improve, and biases the
# answer -- thin ranks against Schur, fat ranks for it.
#
# If you want the scaling curve as well, run the PAIR again at another count
# (16 is the other efficient one for this mesh) rather than splitting a pair
# across two.

#SBATCH --job-name=rtb3dschur
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras
#
# 25 ranks = a 5 x 5 rank grid over rtb3d_schur's 20 x 20 element columns, so
# each rank owns a 4 x 4 block of columns and 84 243 of the mesh's 2 106 081
# points. This is what the repo's own tool picks:
#
#     julia tools/pick_nranks.jl 20 20 80 4 256 10000 10000
#     ...
#      16     4 x 4    5 x 5   25  131630   halo/elem 0.8
#      25     5 x 5    4 x 4   16   84243   halo/elem 1.0   <== RECOMMENDED
#      50     5 x 10   4 x 2    8   42122   halo/elem 1.5   (thin: comms-bound)
#
# 400 element columns cap the useful parallelism, and going past 25 is worse
# than merely wasteful HERE: the MPI reduce is the one part of the stage solve
# the Schur reduction cannot shrink, so thin ranks inflate the term that cannot
# improve and bias the A/B against it.
#
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
#   nodes x ntasks-per-node x cpus-per-task = 1 x 25 x 2
#     -> 25 ranks, 50 CPUs and 200 GB on the node, 8 GB per rank
#
# 8 GB per rank is generous for this mesh and deliberately so on the first run.
# Measured here: ONE rank on the full 20 x 20 x 80 peaked at 13.4 GB and was
# OOM-killed on a 16 GB machine, having got as far as the 3D operator. Spread
# over 25 that is roughly 0.5 GB of distributed data per rank on top of a
# ~1.7 GB Julia-plus-packages baseline, so ~2.2 GB -- comfortably inside 4 GB.
# But the MESH BROADCAST is the peak of the whole run, not the solve, and 4 GB
# is exactly what failed on LESICP2 before --heap-size-hint went in. Take 8 GB
# once, read `seff`, then drop to --cpus-per-task=1 and halve the core-hours.
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
# 25 ranks fit on one node, so this configuration never touches the network --
# MPICH keeps them all in shared memory. That makes the MPI reduce cheaper than
# it would be across nodes, which FLATTERS the Schur arm slightly. It is the
# right rank count for this mesh regardless (see the table above); just do not
# read the ratio as a multi-node result.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --cpus-per-task=2
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

JEXPRESSO_ROOT="${JEXPRESSO_ROOT:-/project/smarras/smarras/Jexpresso}"
EQS="${EQS:-CompEuler}"
CASE="${CASE:-rtb3d_schur}"

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

# PIN THE SOLVER KNOBS. The deck defaults :imex_rtol to 1e-6, but the 16.4
# s/step baseline this probe is meant to explain was measured at 1e-4 (the
# monitor showed residuals of ~9e-05 and 26.2 iterations/solve; at 1e-6 it is
# 54.4 and the step is half as fast again). Leaving it to the deck default
# would produce a profile of a DIFFERENT run and the split would not be
# comparable to the number it is supposed to break down. Both are overridable.
# ---------------------------------------------------------------------------
# THE A/B SWITCH. 0 = the five-field stage solve, 1 = the scalar Schur one.
#
# It is exported with an explicit default rather than left unset so that BOTH
# arms are stated: a run whose log says DBG_SCHUR=0 cannot later be mistaken
# for one that simply predated the flag.
#
# Everything else in this block is pinned so the two arms differ in ONE thing.
# That matters more than the values: an A/B where the tolerance also moved
# measures nothing, and :imex_rtol is the single largest lever on this scheme's
# cost.
#
# NOTE THAT DBG_SCHUR=1 IS A DIFFERENT SPLITTING, not just a different solver.
# It forces the advective Theta row, without which the elimination does not
# close on one scalar. The two rows differ by 0.06% of the flux form
# (test/hevi/test_theta_advective.jl). That is the right comparison for wall
# clock and the wrong one for reading a state difference between the two output
# trees as an error.
# ---------------------------------------------------------------------------
export DBG_SCHUR="${DBG_SCHUR:-0}"

# The rtb3d_schur deck's own defaults, restated here so both arms are pinned to
# them and neither can pick up a different tolerance from an edited deck. 1e-8
# is tight enough that the stage solve is not the leading error term against a
# third-order tableau, and loose enough not to buy digits the step discards.
export DBG_RTOL="${DBG_RTOL:-1.0e-8}"
export DBG_RESTART="${DBG_RESTART:-20}"

# 75 steps at rtb3d_schur's dt = 0.6 s. The profile reports once n_step reaches
# JEXPRESSO_HEVI_PROFILE_EVERY, and n_step only counts steps AFTER the _SKIP
# JIT steps -- so it needs 10 + 50 = 60 before it prints anything, which is
# 36 s of model time, and 45 leaves margin for one clean block. There is no
# reason to run to the deck's 1000 s: the breakdown is a per-step average and a
# second block only confirms it.
#
# RAISE THIS, NOT LOWER IT, if the profile block never appears -- a run that
# stops at 59 steps prints no breakdown at all and looks like a silent failure.
export DBG_TEND="${DBG_TEND:-45.0}"

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
# ---------------------------------------------------------------------------
# SLURM_MEM_PER_CPU vs SLURM_MEM_PER_NODE -- srun refuses when both are set:
#
#     srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE
#           are mutually exclusive.
#
# and it is fatal to the STEP, so it kills the run after the whole serial setup
# and precompile have already been paid for -- the most expensive place to
# discover it.
#
# Both are set because this job asks for --mem-per-cpu (which sets
# SLURM_MEM_PER_CPU) on a site whose defaults also put SLURM_MEM_PER_NODE in the
# batch environment. srun inherits both and refuses.
#
# SLURM_MEM_PER_NODE is the one to drop. The batch job's --mem-per-cpu=4000M is
# what actually sized the cgroup, it is what the cpus-per-task trick above
# depends on, and it is what RANK_MB and --heap-size-hint are derived from.
# Unsetting the per-node variable changes NOTHING about the allocation -- the
# cgroup is already established -- it only stops srun seeing a contradiction it
# cannot resolve.
#
# This runs before `srun --mpi=list` below as well as before any launch, since
# that probe hits the same check.
# ---------------------------------------------------------------------------
if [ -n "$SLURM_MEM_PER_CPU" ] && [ -n "$SLURM_MEM_PER_NODE" ]; then
    echo "    NOTE: both SLURM_MEM_PER_CPU ($SLURM_MEM_PER_CPU) and" \
         "SLURM_MEM_PER_NODE ($SLURM_MEM_PER_NODE) are set;"
    echo "          unsetting SLURM_MEM_PER_NODE -- srun treats them as mutually" \
         "exclusive and this job is sized per CPU."
    unset SLURM_MEM_PER_NODE
fi

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
# WHY AN ALLREDUCE AND NOT JUST MPI_Init.
#
# MPI_Init_thread returning 0 on every rank proves only that each process
# initialised ITS OWN library. It says nothing about whether those processes
# found each other. The failure that hides in that gap is the one worth a
# preflight: with the wrong PMI plugin, MPICH does not error -- it starts N
# INDEPENDENT ONE-RANK JOBS, each with its own COMM_WORLD of size 1. Every
# rank runs, every rank succeeds, and the job then behaves as 256 serial
# Jexpressos writing over each other's output. Nothing in MPI_Init detects it.
#
# A collective is what forces the ranks to actually reach each other, and one
# whose answer is known in advance is what turns "it did not hang" into a
# verified result: summing 1 over the communicator must give Comm_size, and
# Comm_size must be the rank count SLURM was asked for. Both are checked below,
# and either being wrong exits non-zero. The Barrier that follows costs nothing
# and confirms a second collective completes after the first.
export JEXPRESSO_EXPECT_RANKS="$NTASKS"
echo "--- Preflight: MPI_Init_thread + one Allreduce on $NTASKS ranks ---"
launch "${JULIA_FLAGS[@]}" -e '
    using MPI
    MPI.Init()
    c = MPI.COMM_WORLD
    n = MPI.Comm_size(c)
    r = MPI.Comm_rank(c)
    s = MPI.Allreduce(1, MPI.SUM, c)
    want = parse(Int, ENV["JEXPRESSO_EXPECT_RANKS"])
    ok = (n == want) && (s == n)
    if r == 0
        println("    MPI ", ok ? "OK" : "**BROKEN**", ": Comm_size = ", n,
                " (asked for ", want, "), Allreduce(1) = ", s)
        if n == 1 && want > 1
            println("    Every process got a COMM_WORLD of size 1: the ranks never found")
            println("    each other, so the job would have run as ", want, " independent")
            println("    serial Jexpressos writing over each other. That is the PMI plugin --")
            println("    try another SRUN_MPI from the --mpi=list above.")
        elseif n != want
            println("    The launcher started ", n, " ranks, not ", want,
                    ". The allocation and the launch disagree;")
            println("    check -n against SLURM_NTASKS and whether --cpus-per-task changed")
            println("    how many tasks fit.")
        elseif s != n
            println("    Allreduce disagrees with Comm_size: the collective is returning")
            println("    the wrong answer, which is a broken transport, not a broken launch.")
        end
    end
    MPI.Barrier(c)
    MPI.Finalize()
    exit(ok ? 0 : 1)' || {
    rc=$?
    echo "ERROR: MPI cannot run $NTASKS ranks across $NNODES nodes (exit $rc)." >&2
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
echo "    stage solve               : DBG_SCHUR=$DBG_SCHUR  ($([ "$DBG_SCHUR" = "1" ] && echo "SCALAR Schur, Np unknowns" || echo "five-field, 5*Np unknowns"))"
echo "    DBG_TEND / rtol / restart  : $DBG_TEND s (75 steps at dt=0.6) / $DBG_RTOL / $DBG_RESTART"
echo "    JEXPRESSO_HEVI_PROFILE    : $JEXPRESSO_HEVI_PROFILE  (every $JEXPRESSO_HEVI_PROFILE_EVERY steps, first $JEXPRESSO_HEVI_PROFILE_SKIP skipped)"
echo "    JEXPRESSO_PRECOMPILE_PASS : $JEXPRESSO_PRECOMPILE_PASS"
echo "    julia flags               : ${JULIA_FLAGS[*]}"
echo "    started                   : $(date)"

# ENVIRONMENT PROPAGATION. srun exports the submitting environment to every
# task by default, and MPICH/Hydra does the same (-genvall), so the exports
# above reach the ranks under either launcher here. Under OpenMPI's mpirun they
# would NOT -- pass them explicitly:
#   mpirun -x JULIA_PKG_PRECOMPILE_AUTO -x JEXPRESSO_PRECOMPILE_PASS \
#          -x JEXPRESSO_HEVI_PROFILE -x DBG_TEND -x DBG_SCHUR \
#          -x DBG_RTOL -x DBG_RESTART ...
# Without them the probe silently prints nothing and the run is a full-length
# production run: DBG_TEND would not reach the ranks either.
#
# DBG_SCHUR IS THE WORST ONE TO LOSE, because losing it fails QUIETLY and
# PLAUSIBLY. The banner above is printed by this script, from its own
# environment, so it would still say DBG_SCHUR=1 while every rank ran the
# five-field solve -- and the A/B would come back a dead heat with nothing
# anywhere saying why. Check the SETUP REPORT in the run's own output instead:
# it prints "Stage solve: preconditioned GMRES on the SCALAR SCHUR system"
# against "on all 5 fields", from inside the ranks. That line is the one that
# proves which arm actually ran.
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

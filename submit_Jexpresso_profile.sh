#!/bin/bash
#=============================================================================
#  Jexpresso profile run.
#
#      sbatch submit_Jexpresso_profile.sh
#
#  NO ARGUMENTS. Which stage solve runs is the DECK's business, not this
#  file's: it is `limex_schur` in the case's user_inputs.jl, on by default.
#  This script sets resources and runs the case.
#
#  EDIT THE TWO BLOCKS MARKED "EDIT ME". Nothing else in this file needs
#  changing to run a different problem, on a different rank count, on a
#  different cluster.
#
#  For a one-off A/B without editing the deck, set the variable in the
#  environment and it is passed through:
#
#      DBG_SCHUR=0 sbatch submit_Jexpresso_profile.sh     # five-field arm
#
#  WHAT STILL GUARDS THE A/B. This used to demand `schur` or `full` as an
#  argument so that a forgotten word could not silently run the wrong arm. The
#  guard is not lost by dropping it: the ranks print which stage solve they
#  built, from inside the run --
#
#      Stage solve: preconditioned GMRES on the SCALAR SCHUR system
#      Stage solve: preconditioned GMRES on all 5 fields
#
#  -- and that line, not anything this script echoes, is the authoritative
#  record of what was measured.
#=============================================================================

#--------------------------------------------------------------- EDIT ME (1) -
#  RESOURCES.  nodes x ntasks-per-node = your MPI RANK COUNT.
#
#  DO NOT pick a rank count freely. With :lxy_partition the mesh is cut into
#  vertical columns, so the usable parallelism is bounded by nelemx*nelemy and
#  anything else leaves ranks owning zero elements (fatal). Ask first:
#
#      julia tools/pick_nranks.jl <nelemx> <nelemy> <nelemz> <nop> <max_cores>
#
#  Answers for the cases here:
#      LESICP2-64x64x60-imex   256  (16 x 16 over 64 x 64 columns)  <- set below
#      rtb3d_schur              25  (5 x 5 over 20 x 20)
#      rtb2d_schur               4  (a slab has only 20 columns)
#
#  cpus-per-task BUYS MEMORY, NOT SPEED. This is pure MPI -- the job exports
#  JULIA_NUM_THREADS=1 -- so the second core per task sits idle and its only
#  effect is that each rank gets mem-per-cpu x cpus-per-task. 64 x 2 = 128
#  fills a Wulver node; 3500M x 2 = 7000M/rank is what setup_les_run.sh found
#  this domain needs, and 64 x 7000M = 448 GB/node leaves the node headroom.
#-----------------------------------------------------------------------------
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3500M
#
#  Site settings -- set once for your cluster, then forget.
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras
#SBATCH --job-name=jexpresso
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

#--------------------------------------------------------------- EDIT ME (2) -
#  WHAT TO RUN.
#-----------------------------------------------------------------------------
EQS="CompEuler"
CASE="LESICP2-64x64x60-imex"
TEND="10800.0"
MESH=""
ROOT="/project/smarras/smarras/Jexpresso"
MODULES=(Julia/1.11.9 GCC MPICH)

#=============================================================================
#  Nothing below here needs editing.
#=============================================================================
set -u

# An argument is no longer accepted. Refuse one rather than ignore it: a
# leftover `sbatch ... schur` from muscle memory should not look like it did
# something.
if [ $# -gt 0 ]; then
    echo "ERROR: this script takes no arguments (got '$1')." >&2
    echo "       Which stage solve runs is set in the deck --" >&2
    echo "       limex_schur in problems/$EQS/$CASE/user_inputs.jl, on by" >&2
    echo "       default. For a one-off five-field run:" >&2
    echo "           DBG_SCHUR=0 sbatch $(basename "$0")" >&2
    exit 2
fi

for m in "${MODULES[@]}"; do module load "$m" || exit 1; done
cd "$ROOT" || exit 1

#-- 2. deck settings, passed through the environment ----------------------
# DBG_SCHUR is NOT set here -- the deck owns it. Passed through only when the
# caller put it in the environment, on the same terms as rtol/restart/maxiter
# below.
[ -n "${DBG_SCHUR:-}" ] && export DBG_SCHUR
[ -n "$TEND" ] && export DBG_TEND="$TEND"
[ -n "$MESH" ] && export DBG_MESH="$MESH"
# THE LAUNCHER DOES NOT SET rtol OR restart. It used to default them to 1.0e-8
# and 20 -- the rtb values -- which silently OVERRODE any deck with its own
# tuning. LESICP2 ships 1.0e-6 and 30 for a reason (it runs at CFL_h 2.77), and
# a launcher quietly imposing a 100x tighter tolerance at a shorter restart is
# how a case that converges becomes a case that does not. Pass them only when
# the user asked for them; otherwise the deck decides.
[ -n "${DBG_RTOL:-}" ]    && export DBG_RTOL
[ -n "${DBG_RESTART:-}" ] && export DBG_RESTART
[ -n "${DBG_MAXITER:-}" ] && export DBG_MAXITER

export JEXPRESSO_HEVI_PROFILE=1
export JEXPRESSO_HEVI_PROFILE_EVERY=50
export JEXPRESSO_HEVI_PROFILE_SKIP=10
export JEXPRESSO_PRECOMPILE_PASS="${JEXPRESSO_PRECOMPILE_PASS:-1}"
export JEXPRESSO_STEP_HEARTBEAT=1

# One BLAS/Julia thread per rank: the ranks already fill the node, and nested
# threading oversubscribes it.
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

#-- 3. serial setup, before any rank is launched --------------------------
# Each step here fails cheaply and says why. Compiling inside the parallel
# launch instead means every rank compiles the same code at once.
export JULIA_PKG_PRECOMPILE_AUTO=1

echo "--- MPI preferences (must precede any compilation) ---"
julia --project=. --startup-file=no \
    -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1

echo "--- Syntax check (seconds, before paying minutes for a precompile) ---"
julia --startup-file=no tools/syntax_check.jl src test problems tools || {
    echo "ERROR: syntax error in the source tree. Nothing compiled, no ranks launched." >&2
    exit 1; }

echo "--- Precompile ---"
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' || exit 1

echo "--- Serial load (exercises the real load path) ---"
julia --project=. --startup-file=no -e '
    t = @elapsed (using MPI; using Jexpresso)
    println("    Jexpresso loaded in ", round(t, digits=2), " s")' || {
    echo "ERROR: Jexpresso failed to load serially. Fix that before launching ranks." >&2
    exit 1; }

export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=1

#-- 4. julia flags --------------------------------------------------------
JULIA_FLAGS=(--project=. --startup-file=no)
julia --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--compiled-modules=existing)
julia --pkgimages=existing       -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--pkgimages=existing)

# Julia sizes its heap from /proc/meminfo, i.e. the whole NODE, not the cgroup
# this job actually gets -- so without a hint it happily grows past the limit
# and is OOM-killed with no Julia-level error.
if [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    RANK_MB=$(( SLURM_MEM_PER_CPU * ${SLURM_CPUS_PER_TASK:-1} ))
elif [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ -n "${SLURM_NTASKS_PER_NODE:-}" ]; then
    RANK_MB=$(( SLURM_MEM_PER_NODE / SLURM_NTASKS_PER_NODE ))
else
    RANK_MB=0
fi
HEAP_MB=0
if [ "$RANK_MB" -gt 0 ] && julia --heap-size-hint=1G -e 'exit(0)' >/dev/null 2>&1; then
    HEAP_MB=$(( RANK_MB * 85 / 100 ))
    JULIA_FLAGS+=(--heap-size-hint=${HEAP_MB}M)
fi

#-- 5. launcher -----------------------------------------------------------
NTASKS="${SLURM_NTASKS:-1}"
NNODES="${SLURM_NNODES:-1}"

# srun rejects a job that carries both, and this one is sized per CPU.
if [ -n "${SLURM_MEM_PER_CPU:-}" ] && [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    echo "    (unsetting SLURM_MEM_PER_NODE: srun treats it as exclusive with _PER_CPU)"
    unset SLURM_MEM_PER_NODE
fi

SRUN_MPI="${SRUN_MPI:-pmi2}"
if [ -n "${JEXPRESSO_LAUNCHER:-}" ]; then LAUNCHER="$JEXPRESSO_LAUNCHER"
elif [ "$NNODES" -gt 1 ] && command -v srun >/dev/null 2>&1; then LAUNCHER="srun"
else LAUNCHER="mpirun"; fi

launch() {
    if [ "$LAUNCHER" = "srun" ]; then
        srun --mpi="$SRUN_MPI" -n "$NTASKS" -c "${SLURM_CPUS_PER_TASK:-1}" julia "$@"
    else
        mpirun -np "$NTASKS" julia "$@"
    fi
}

#-- 6. preflight: can MPI actually run these ranks? -----------------------
# A launcher that gives every process a COMM_WORLD of size 1 produces N
# independent serial runs writing over each other, which looks like a slow
# success. Catch it before the mesh read.
export JEXPRESSO_EXPECT_RANKS="$NTASKS"
echo "--- Preflight: MPI_Init_thread + one Allreduce on $NTASKS ranks ---"
launch "${JULIA_FLAGS[@]}" -e '
    using MPI; MPI.Init()
    c = MPI.COMM_WORLD; n = MPI.Comm_size(c); r = MPI.Comm_rank(c)
    s = MPI.Allreduce(1, MPI.SUM, c)
    want = parse(Int, ENV["JEXPRESSO_EXPECT_RANKS"])
    ok = (n == want) && (s == n)
    r == 0 && println("    MPI ", ok ? "OK" : "**BROKEN**", ": Comm_size = ", n,
                      " (asked for ", want, "), Allreduce(1) = ", s)
    MPI.Barrier(c); MPI.Finalize(); exit(ok ? 0 : 1)' || {
    rc=$?
    echo "ERROR: MPI cannot run $NTASKS ranks across $NNODES nodes (exit $rc)." >&2
    echo "       Nothing in Jexpresso has run. If Comm_size came back 1, it is the" >&2
    echo "       PMI plugin: try another SRUN_MPI (srun --mpi=list). If it came" >&2
    echo "       back wrong-but-not-1, the allocation and the launch disagree." >&2
    echo "       Single node failing too means it is not the network: check /dev/shm" >&2
    echo "       and that --mem-per-cpu x --ntasks-per-node fits the node." >&2
    exit $rc; }

#-- 7. run ----------------------------------------------------------------
# This shell cannot know what the deck chose, so it does not pretend to: it
# reports an override if one was given and defers otherwise. The ranks' own
# setup report is the authoritative line -- see the header.
case "${DBG_SCHUR:-}" in
    1) arm_desc="SCALAR Schur, Np unknowns   (forced by DBG_SCHUR=1)" ;;
    0) arm_desc="five-field, 5*Np unknowns   (forced by DBG_SCHUR=0)" ;;
    *) arm_desc="<deck decides -- see 'Stage solve:' in the run's own output>" ;;
esac

echo "--- Launching $NTASKS ranks ---"
echo "    case          : $EQS / $CASE${MESH:+  (mesh $MESH)}"
echo "    stage solve   : $arm_desc"
echo "    resources     : $NNODES node(s) x $((NTASKS / NNODES)) tasks x ${SLURM_CPUS_PER_TASK:-1} cpus"
[ "$LAUNCHER" = srun ] && LAUNCHER_DESC="srun --mpi=$SRUN_MPI" || LAUNCHER_DESC="mpirun"
echo "    launcher      : $LAUNCHER_DESC"
if [ "$HEAP_MB" -gt 0 ]; then
    echo "    memory/rank   : ${RANK_MB} MB budget, --heap-size-hint=${HEAP_MB}M"
else
    echo "    memory/rank   : no heap hint set -- Julia sizes its heap from the whole"
    echo "                    node and can be OOM-killed against a smaller cgroup"
fi
echo "    tend / rtol   : ${DBG_TEND:-<deck>} s / ${DBG_RTOL:-<deck>}, restart ${DBG_RESTART:-<deck>}, maxiter ${DBG_MAXITER:-<deck>}"
echo "    started       : $(date)"
echo
# The banner above is printed from THIS shell. The authoritative statement of
# which arm the ranks ran is in their own setup report:
#     "Stage solve: preconditioned GMRES on the SCALAR SCHUR system"
#     "    matvec: bespoke scalar sweeps ..."
# versus "... on all 5 fields". Read that, not this.

launch "${JULIA_FLAGS[@]}" src/Jexpresso.jl "$EQS" "$CASE"
rc=$?
echo "--- Finished: $(date) ---"

#-- 8. report the exit status ---------------------------------------------
# Write with plain echo / >&2, NEVER "> /dev/stdout": under SLURM those are
# regular files and re-opening them truncates the whole job log.
if [ "$rc" -ne 0 ]; then
    log="${SLURM_JOB_NAME:-job}.${SLURM_JOB_ID:-0}.out"
    _report() {
        echo "--- FAILED: launcher exited $rc ---"
        echo "    Julia's error and backtrace are ABOVE this line, in the .err for"
        echo "    job ${SLURM_JOB_ID:-<id>}."
        if [ -f "$log" ] && grep -q "Min elements: 0" "$log"; then
            echo "    'Min elements: 0' in this run: some ranks own no elements."
            echo "    The rank count exceeds what this mesh can be split into --"
            echo "    ask tools/pick_nranks.jl and use a count it lists."
        fi
    }
    _report; _report >&2
fi
exit $rc

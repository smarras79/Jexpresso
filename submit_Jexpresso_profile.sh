#!/bin/bash
#=============================================================================
#  Jexpresso profile run.
#
#      sbatch submit_Jexpresso_profile.sh schur      scalar Schur stage solve
#      sbatch submit_Jexpresso_profile.sh full       five-field stage solve
#
#  EDIT THE TWO BLOCKS MARKED "EDIT ME". Nothing else in this file needs
#  changing to run a different problem, on a different rank count, on a
#  different cluster.
#
#  The argument is REQUIRED and a typo exits 2: half an A/B that turns out to
#  be two copies of the same arm looks exactly like a result.
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
#  25 is that tool's answer for rtb3d_schur's 20x20x80. For rtb2d_schur use 4.
#-----------------------------------------------------------------------------
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --cpus-per-task=2
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=4000M
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
CASE="rtb3d_schur"          # rtb3d_schur | rtb2d_schur | LESICP2-64x64x60-imex
TEND="45.0"                 # model seconds. 45 = 75 steps at dt 0.6 = one
                            # profile block. Use the deck's own tend (blank)
                            # for physics rather than timing.
MESH=""                     # override the deck's mesh, e.g. 10x10x40. Blank =
                            # whatever the deck picks.
ROOT="/project/smarras/smarras/Jexpresso"
MODULES=(Julia/1.11.9 GCC MPICH)

#=============================================================================
#  Nothing below here needs editing.
#=============================================================================
set -u

#-- 1. which stage solve --------------------------------------------------
case "${1:-}" in
    schur) SCHUR=1 ;;
    full)  SCHUR=0 ;;
    *)
        echo "ERROR: give exactly one of: schur | full" >&2
        echo "       (got '${1:-<nothing>}'). Refusing to guess: running the" >&2
        echo "       wrong arm looks exactly like a result." >&2
        exit 2 ;;
esac

for m in "${MODULES[@]}"; do module load "$m" || exit 1; done
cd "$ROOT" || exit 1

#-- 2. deck settings, passed through the environment ----------------------
export DBG_SCHUR="$SCHUR"
[ -n "$TEND" ] && export DBG_TEND="$TEND"
[ -n "$MESH" ] && export DBG_MESH="$MESH"
export DBG_RTOL="${DBG_RTOL:-1.0e-8}"
export DBG_RESTART="${DBG_RESTART:-20}"

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
[ "$SCHUR" = 1 ] && arm_desc="SCALAR Schur, Np unknowns" \
                 || arm_desc="five-field, 5*Np unknowns"

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
echo "    tend / rtol   : ${DBG_TEND:-<deck>} s / $DBG_RTOL, restart $DBG_RESTART"
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

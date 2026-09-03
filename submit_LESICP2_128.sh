#!/bin/bash -l
#=============================================================================
#  Jexpresso profile run.
#
#      sbatch submit_Jexpresso_profile.sh
#
#  NO ARGUMENTS. Which stage solve runs is the DECK's business, not this
#  file's: it is `use_imex` / `use_schur` in the case's user_inputs.jl, both
#  on by default. On the 64x64x60 decks `use_schur = !_vdiff`, so DBG_VDIFF=0
#  is what selects the Schur arm (with explicit vertical diffusion).
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
#      LESICP2-128x128x60-imex 1024 (32 x 32 over 128 x 128 columns) <- set below
#      LESICP2-64x64x60-imex    256 (16 x 16 over 64 x 64 columns)
#      1280 ranks DOES NOT divide 128x128 = 16384 columns (12.8 each).
#      rtb3d_schur              25  (5 x 5 over 20 x 20)
#      rtb2d_schur               4  (a slab has only 20 columns)
#
#  cpus-per-task BUYS MEMORY, NOT SPEED. This is pure MPI -- the job exports
#  JULIA_NUM_THREADS=1 -- so the extra cores per task sit idle and their only
#  effect is that each rank gets mem-per-cpu x cpus-per-task.
#
#  4000M IS THE NODE'S OWN RATIO: a Wulver node is 128 cores x 4000M = 512 GB,
#  and hw59's submit_high_LESICP2.sh ran 128 x 4000M, so the scheduler accepts
#  it. It replaces the 3500M this script was copied with -- that figure came
#  from setup_les_run.sh sizing the 64x64x60 DOMAIN, not from any property of
#  the hardware, and at 64 ranks/node it left 64 GB/node unclaimed while we
#  were busy suspecting a memory wall.
#
#  MEMORY PER RANK IS SET BY RANKS-PER-NODE, NOT BY mem-per-cpu. Keep
#  ntasks-per-node x cpus-per-task = 128 and the node total is 512 GB either
#  way; what moves is how it is divided:
#
#    nodes  tasks/node  cpus/task   mem/rank   ranks
#      16       64          2         8 GB      1024   <- here
#      32       32          4        16 GB      1024
#      64       16          8        32 GB      1024
#
#  All three are the same 1024 ranks and the same 32 x 32 rank grid; they
#  differ only in how thinly each node is packed. THE GLOBAL MESH IS THE
#  REASON THIS MATTERS: :lxy_partition uses the rank-0-read + MPI.bcast path,
#  so every rank transiently holds the WHOLE 63.4 M-point mesh -- ~1.5 GB of
#  coordinates plus ~1.0 GB of connectivity plus Gridap's structures -- before
#  it is partitioned. That transient is 4x what the 64x64x60 deck needed at
#  the same 7000M/rank, which is the leading suspect for the SIGBUS on job
#  1219979. If 8 GB/rank still dies in the mesh read, go to the 32-node row:
#  it is a three-line change and nothing else in the deck moves.
#-----------------------------------------------------------------------------
#SBATCH --exclusive
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
# 71:59:00, not 2 hours. At dt = 0.5 the deck's tend of 10800 s is 21,600 steps
# and the measured rate is ~3.7 s/step, i.e. around 22 h -- a 2 h limit killed
# it about 9% in. Drop this back to something short only for a TEND-limited
# probe run.
#SBATCH --time=71:59:00
#SBATCH --mem-per-cpu=4000M
#
#  Site settings -- set once for your cluster, then forget.
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --job-name=ICP2_128
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

#--------------------------------------------------------------- EDIT ME (2) -
#  WHAT TO RUN.
#-----------------------------------------------------------------------------
EQS="CompEuler"
CASE="LESICP2-128x128x60-imex"
TEND="10800.0"
MESH=""
# ROOT: the repo checkout. Defaults to the directory you ran `sbatch` FROM,
# which SLURM records in SLURM_SUBMIT_DIR -- so a second clone at a different
# path needs no edit here. Override with JEXPRESSO_ROOT=... if you sbatch from
# somewhere else.
#
# NOT `dirname $0` / BASH_SOURCE: SLURM runs a COPY of this script out of its
# spool directory, so inside the job those point at the spool, not the repo.
ROOT="${JEXPRESSO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
# Where the VTK trees go. INFERRED from ROOT by analogy -- hw59's old script
# did not set it, the deck defaults to another user's scratch, and ~543 GB
# lands here. CHECK IT EXISTS AND HAS THE QUOTA before a long run.
export JEXPRESSO_OUTDIR="${JEXPRESSO_OUTDIR:-/scratch/smarras/hw59/output_new}"
# Taken verbatim from hw59's working submit_high_LESICP2.sh. NOTE: it loads
# no Julia module -- julia comes from the login environment, which is why
# the shebang above is `bash -l`. If `module load bright` fails on its own
# (this script loads them one at a time), collapse them back into a single
# `module load bright shared mpich/ge/gcc/64` here.
MODULES=(bright shared mpich/ge/gcc/64)

# JULIA SELECTOR. `+1.11.9` is juliaup syntax: it picks that channel for the
# invocation. INSTALL.md recommends 1.11.9 and Project.toml's floor is 1.11.5.
# Every julia call below goes through this array, so it is the ONE place to
# change -- and, more to the point, the serial precompile and the 1024 ranks
# are then guaranteed to be the same Julia. They must be: a mismatch means
# every rank finds an unusable cache and falls back to compiling privately,
# which is the failure the two-phase setup exists to prevent.
#
# If julia here is NOT juliaup, `+1.11.9` is passed through as a script
# argument and Julia will complain -- set JULIA=(julia) in that case, or point
# JEXPRESSO_JULIA at an absolute path.
JULIA=("${JEXPRESSO_JULIA:-julia}" +1.11.9)

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
    echo "       use_imex / use_schur in problems/$EQS/$CASE/user_inputs.jl," >&2
    echo "       default. For a one-off five-field run:" >&2
    echo "           DBG_SCHUR=0 sbatch $(basename "$0")" >&2
    exit 2
fi

# LOCKED MEMORY, BEFORE ANYTHING TOUCHES THE FABRIC.
#
# InfiniBand registers ("pins") the buffers it sends from, and pinned pages
# count against RLIMIT_MEMLOCK -- `ulimit -l` -- not against the node's RAM.
# When that limit is small, MPI_Init fails at scale with
#
#     Unable to alloc send buffer MR on mlx5_0: Cannot allocate memory
#     Unable to allocate UD send buffer pool
#
# which is what killed job 1220165 in the PREFLIGHT (a bare MPI_Init plus one
# Allreduce) on 1024 ranks, while the node had ~395 GiB of RAM free. It is a
# limit, not a shortage, and it scales with ranks-per-node: 64 ranks each
# pinning their own send pools is 64x one rank''s demand.
#
# Raising it is allowed only if the hard limit permits, so this is best-effort
# and never fatal -- and the value actually in force is printed either way,
# because "we asked" and "we got" are different things and the difference is
# exactly what this failure looks like.
ulimit -l unlimited 2>/dev/null || true
echo "--- Limits: locked memory (ulimit -l) = $(ulimit -l), hard = $(ulimit -H -l) ---"

for m in "${MODULES[@]}"; do module load "$m" || exit 1; done
cd "$ROOT" || exit 1
# Fail here, with the path in hand, rather than 200 lines later on a missing
# deck -- a wrong ROOT is the single most likely thing to be wrong in a fresh
# clone or a moved checkout.
[ -f "$ROOT/src/Jexpresso.jl" ] || {
    echo "ERROR: ROOT=$ROOT is not a Jexpresso checkout (no src/Jexpresso.jl)." >&2
    echo "       sbatch this script FROM the repo root, or set JEXPRESSO_ROOT." >&2
    exit 1; }
[ -d "$ROOT/problems/$EQS/$CASE" ] || {
    echo "ERROR: no deck at $ROOT/problems/$EQS/$CASE" >&2
    echo "       Are you on the branch that carries it? git log --oneline -1" >&2
    exit 1; }

#-- 2. deck settings, passed through the environment ----------------------
# DBG_SCHUR is NOT set here -- the deck owns it. Passed through only when the
# caller put it in the environment, on the same terms as rtol/restart/maxiter
# below.
[ -n "${DBG_SCHUR:-}" ] && export DBG_SCHUR
# DBG_VDIFF picks the IMEX ARM on the LESICP2 64x64x60 decks: the deck sets
# `use_schur = !_vdiff`, so DBG_VDIFF=0 is what selects the Schur stage solve
# with EXPLICIT vertical diffusion. Without this line a caller's DBG_VDIFF is a
# plain shell variable that never reaches Julia, the deck keeps its default,
# and the run silently takes the other arm.
[ -n "${DBG_VDIFF:-}" ]   && export DBG_VDIFF
[ -n "${DBG_FILTER:-}" ]  && export DBG_FILTER
[ -n "${DBG_SCHEME:-}" ]  && export DBG_SCHEME
# `${DBG_TEND:-$TEND}`, not `$TEND`: a bare assignment here OVERRODE the
# caller, so `DBG_TEND=5 sbatch ...` silently ran the deck's full tend. The
# EDIT ME value is the default; the environment wins, as it does for DBG_SCHUR
# and the Krylov knobs below.
export DBG_TEND="${DBG_TEND:-$TEND}"
[ -n "${DBG_MESH:-$MESH}" ] && export DBG_MESH="${DBG_MESH:-$MESH}"
# THE LAUNCHER DOES NOT SET rtol OR restart. It used to default them to 1.0e-8
# and 20 -- the rtb values -- which silently OVERRODE any deck with its own
# tuning. LESICP2 ships 1.0e-6 and 30 for a reason (it runs at CFL_h 2.77), and
# a launcher quietly imposing a 100x tighter tolerance at a shorter restart is
# how a case that converges becomes a case that does not. Pass them only when
# the user asked for them; otherwise the deck decides.
[ -n "${DBG_RTOL:-}" ]    && export DBG_RTOL
[ -n "${DBG_RESTART:-}" ] && export DBG_RESTART
[ -n "${DBG_MAXITER:-}" ] && export DBG_MAXITER

# `:-` on all three: bare assignments clobbered the caller, so
# JEXPRESSO_HEVI_PROFILE_EVERY=1 on the sbatch line was ignored and the
# breakdown still came every 50 steps. Same rule as DBG_TEND and DBG_VDIFF --
# the value here is the default, the environment wins.
export JEXPRESSO_HEVI_PROFILE="${JEXPRESSO_HEVI_PROFILE:-1}"
export JEXPRESSO_HEVI_PROFILE_EVERY="${JEXPRESSO_HEVI_PROFILE_EVERY:-50}"
export JEXPRESSO_HEVI_PROFILE_SKIP="${JEXPRESSO_HEVI_PROFILE_SKIP:-10}"
export JEXPRESSO_PRECOMPILE_PASS="${JEXPRESSO_PRECOMPILE_PASS:-1}"
export JEXPRESSO_STEP_HEARTBEAT=1
# THE IMPLICIT-VDIFF ARM IS THIS SCRIPT'S DEFAULT, AND THAT IS THE CHANGE FROM
# submit_jexpresso_profile.sh, WHICH DEFAULTS TO 0. This is a production
# launcher, not a profiling one, and the two arms are not interchangeable:
#
#   DBG_VDIFF=1  vertical SGS diffusion IMPLICIT -> five-field stage solve.
#                Slower per step. RUNS.
#   DBG_VDIFF=0  vertical SGS diffusion EXPLICIT -> scalar Schur stage solve,
#                ~3.56x faster per step. DIES AROUND t = 500 s.
#
# One variable sets both because the deck derives `use_schur = !_vdiff`: a
# Schur reduction cannot see the diffusion operator, so the cheap solve and
# the stable diffusion are mutually exclusive. The t = 500 s failure is
# recorded at the top of the deck -- with diffusion explicit, nu_eff/dz^2 is
# back in the explicit budget and ARS343 loses it once nu_t reaches 20-40
# m^2/s, which this boundary layer reaches in the first few hundred seconds.
#
# A bare 0 here (as in the profile script) OVERRIDES the deck's own default of
# true, so an sbatch would silently take the arm that blows up. `:-1` keeps
# `DBG_VDIFF=0 sbatch ...` working for a deliberate A/B or a timing run.
export DBG_VDIFF="${DBG_VDIFF:-1}"

# One BLAS/Julia thread per rank: the ranks already fill the node, and nested
# threading oversubscribes it.
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

#-- 3. serial setup, before any rank is launched --------------------------
# Each step here fails cheaply and says why. Compiling inside the parallel
# launch instead means every rank compiles the same code at once.
export JULIA_PKG_PRECOMPILE_AUTO=1

# Resolve the selector before anything depends on it. juliaup will try to
# INSTALL a missing channel, which needs the network and is not something to
# discover 16 nodes into an allocation.
echo "--- Julia ---"
"${JULIA[@]}" --version || {
    echo "ERROR: ${JULIA[*]} does not run. If julia here is not juliaup, drop" >&2
    echo "       the +1.11.9 selector (JULIA=(julia)); if it is, install the" >&2
    echo "       channel on a login node first: juliaup add 1.11.9" >&2
    exit 1; }

echo "--- MPI preferences (must precede any compilation) ---"
"${JULIA[@]}" --project=. --startup-file=no \
    -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1

echo "--- Syntax check (seconds, before paying minutes for a precompile) ---"
"${JULIA[@]}" --startup-file=no tools/syntax_check.jl src test problems tools || {
    echo "ERROR: syntax error in the source tree. Nothing compiled, no ranks launched." >&2
    exit 1; }

echo "--- Precompile ---"
"${JULIA[@]}" --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' || exit 1

echo "--- Serial load (exercises the real load path) ---"
"${JULIA[@]}" --project=. --startup-file=no -e '
    t = @elapsed (using MPI; using Jexpresso)
    println("    Jexpresso loaded in ", round(t, digits=2), " s")' || {
    echo "ERROR: Jexpresso failed to load serially. Fix that before launching ranks." >&2
    exit 1; }

export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=1

#-- 4. julia flags --------------------------------------------------------
JULIA_FLAGS=(--project=. --startup-file=no)
"${JULIA[@]}" --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--compiled-modules=existing)
"${JULIA[@]}" --pkgimages=existing       -e 'exit(0)' >/dev/null 2>&1 && JULIA_FLAGS+=(--pkgimages=existing)

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
if [ "$RANK_MB" -gt 0 ] && "${JULIA[@]}" --heap-size-hint=1G -e 'exit(0)' >/dev/null 2>&1; then
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

# LAUNCHER: mpirun with an explicit hostfile, which is what hw59's
# submit_high_LESICP2.sh used and is therefore the form KNOWN to work with
# this site's mpich/ge module. srun --mpi=pmi2 is the generic path and is
# available with JEXPRESSO_LAUNCHER=srun if the site prefers it.
SRUN_MPI="${SRUN_MPI:-pmi2}"
if [ -n "${JEXPRESSO_LAUNCHER:-}" ]; then LAUNCHER="$JEXPRESSO_LAUNCHER"
else LAUNCHER="mpirun"; fi

# Built once, here, rather than inline in the run command -- the preflight
# below has to use the SAME launch line as the run, or it proves nothing.
NODEFILE="$ROOT/nodefile.$SLURM_JOB_ID"
if [ "$LAUNCHER" = "mpirun" ] && [ -n "${SLURM_NODELIST:-}" ]; then
    scontrol show hostname "$SLURM_NODELIST" > "$NODEFILE"
    echo "    hostfile: $NODEFILE ($(wc -l < "$NODEFILE") nodes)"
fi
trap '[ -n "${NODEFILE:-}" ] && rm -f "$NODEFILE"' EXIT

# RANKS PER NODE, STATED RATHER THAN ASSUMED. --ntasks-per-node shapes the
# ALLOCATION; what actually places ranks is mpirun, and its default mapping is
# its own business. Our hostfile lists each node ONCE, so `--map-by node`
# round-robins and happens to land 1024/16 = 64 per node -- correct by
# arithmetic, not by instruction, and silently wrong the moment NTASKS is not
# an exact multiple of NNODES. The flag that says it outright differs by
# family, so pick it from the launcher's own banner.
PPN=$(( NTASKS / NNODES ))
MPIRUN_BANNER="$(mpirun --version 2>&1 | head -3)"
case "$MPIRUN_BANNER" in
    *Hydra*|*HYDRA*|*MPICH*)      MAP_FLAGS=(-ppn "$PPN") ;;         # MPICH/Hydra
    *"Open MPI"*|*OpenRTE*|*prte*) MAP_FLAGS=(--map-by "ppr:${PPN}:node") ;;
    *)                            MAP_FLAGS=(--map-by node) ;;        # unknown: as before
esac
echo "    ranks/node    : $PPN  (mpirun flags: ${MAP_FLAGS[*]})"

launch() {
    if [ "$LAUNCHER" = "srun" ]; then
        # srun honours --ntasks-per-node from the allocation natively.
        srun --mpi="$SRUN_MPI" -n "$NTASKS" -c "${SLURM_CPUS_PER_TASK:-1}" "${JULIA[@]}" "$@"
    elif [ -s "${NODEFILE:-/nonexistent}" ]; then
        mpirun -np "$NTASKS" -hostfile "$NODEFILE" "${MAP_FLAGS[@]}" "${JULIA[@]}" "$@"
    else
        mpirun -np "$NTASKS" "${MAP_FLAGS[@]}" "${JULIA[@]}" "$@"
    fi
}

#-- 6. preflight: can MPI actually run these ranks? -----------------------
# A launcher that gives every process a COMM_WORLD of size 1 produces N
# independent serial runs writing over each other, which looks like a slow
# success. Catch it before the mesh read.
export JEXPRESSO_EXPECT_RANKS="$NTASKS"
echo "--- Preflight: MPI_Init_thread + one Allreduce on $NTASKS ranks ---"
export JEXPRESSO_EXPECT_PPN="$PPN"
launch "${JULIA_FLAGS[@]}" -e '
    using MPI; MPI.Init()
    c = MPI.COMM_WORLD; n = MPI.Comm_size(c); r = MPI.Comm_rank(c)
    s = MPI.Allreduce(1, MPI.SUM, c)
    want = parse(Int, ENV["JEXPRESSO_EXPECT_RANKS"])
    ok = (n == want) && (s == n)

    # WHERE THE RANKS ACTUALLY LANDED. Asking for 64/node and getting 128 on
    # half the nodes is how a job OOMs or SIGBUSes with every resource request
    # looking correct in the SLURM accounting -- so count them rather than
    # trust the flag. Gathered as fixed-width bytes: hostnames vary in length
    # and Gather wants a uniform count.
    hn = gethostname(); W = 64
    buf = zeros(UInt8, W); b = codeunits(hn); copyto!(buf, 1, b, 1, min(length(b), W))
    all_hn = MPI.Gather(buf, 0, c)
    if r == 0
        hosts = [replace(String(all_hn[(i-1)*W+1 : i*W]), "\0" => "") for i in 1:n]
        tally = Dict{String,Int}(); for h in hosts; tally[h] = get(tally, h, 0) + 1; end
        wantppn = parse(Int, ENV["JEXPRESSO_EXPECT_PPN"])
        counts  = sort(collect(values(tally)))
        even    = all(==(wantppn), counts)
        println("    placement     : ", length(tally), " node(s), ",
                even ? "$wantppn ranks each -- as asked" :
                       "UNEVEN $(minimum(counts))..$(maximum(counts)) ranks/node (wanted $wantppn)")
        even || for (h, k) in sort(collect(tally)); println("        ", h, "  ", k); end
        ok &= even
    end
    ok_arr = [ok]; MPI.Bcast!(ok_arr, 0, c); ok = ok_arr[1]
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
    echo "       'Unable to alloc send buffer MR on mlx5_x: Cannot allocate memory'" >&2
    echo "       or 'Unable to allocate UD send buffer pool' is NOT a shortage of" >&2
    echo "       RAM -- it is RLIMIT_MEMLOCK. See the 'Limits:' line above; if it" >&2
    echo "       is not unlimited the site caps it, and the fix is fewer ranks per" >&2
    echo "       node (--nodes=32 --ntasks-per-node=32 --cpus-per-task=4 keeps the" >&2
    echo "       same 1024 ranks) or a raised hard limit from the admins." >&2
    exit $rc; }

#-- 7. run ----------------------------------------------------------------
# This shell cannot know what the deck chose, so it does not pretend to: it
# reports an override if one was given and defers otherwise. The ranks' own
# setup report is the authoritative line -- see the header.
case "${DBG_SCHUR:-}" in
    1) arm_desc="SCALAR Schur, Np unknowns   (forced by DBG_SCHUR=1)" ;;
    0) arm_desc="five-field, 5*Np unknowns   (forced by DBG_SCHUR=0)" ;;
    *) case "${DBG_VDIFF:-}" in
           # On the 64x64x60 decks use_schur = !_vdiff, so DBG_VDIFF settles it
           # when DBG_SCHUR is not given. Other decks ignore DBG_VDIFF, hence
           # "expected" rather than a claim.
           0) arm_desc="SCALAR Schur, Np unknowns   (expected: DBG_VDIFF=0, explicit diffusion)" ;;
           1) arm_desc="five-field, 5*Np unknowns   (expected: DBG_VDIFF=1, implicit diffusion)" ;;
           *) arm_desc="<deck decides -- see 'Stage solve:' in the run's own output>" ;;
       esac ;;
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
echo "    vert. diff    : ${DBG_VDIFF:+DBG_VDIFF=}${DBG_VDIFF:-<deck>}"
# Print the diagnostic switches the ranks will actually see. A monitor that
# does not fire is otherwise indistinguishable from a variable that never
# arrived, and that ambiguity has cost a debugging cycle already.
echo "    diagnostics   : DSGS_MONITOR=${JEXPRESSO_DSGS_MONITOR:-0}  HEVI_PROFILE=${JEXPRESSO_HEVI_PROFILE:-0} every ${JEXPRESSO_HEVI_PROFILE_EVERY}/skip ${JEXPRESSO_HEVI_PROFILE_SKIP}  imex_monitor_every=${DBG_IMEXMONEVERY:-<deck>}"
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

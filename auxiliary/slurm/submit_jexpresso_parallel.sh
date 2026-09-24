#!/bin/bash -l
#
# ============================================================================
# submit_jexpresso_parallel.sh — run Jexpresso on a large Slurm cluster
# ============================================================================
#
# Submit with:
#
#   sbatch --export=ALL,EQS=CompEuler,CASE=theta auxiliary/slurm/submit_jexpresso_parallel.sh
#
# or edit the EQS/CASE defaults below and just `sbatch` it.
#
# ---------------------------------------------------------------------------
# THE POINT OF THIS SCRIPT
# ---------------------------------------------------------------------------
# A Julia MPI job on a big machine has one failure mode that dwarfs all the
# others: N ranks starting at once against a cold precompile cache. Every rank
# discovers the same missing `.ji`/`.so` files, every rank tries to build them,
# and they all write into the same depot. On 512 ranks that is 512 compilers
# hammering one shared filesystem — minutes to hours of wall time, lock
# contention, half-written cache files, and occasionally a job that never
# starts at all.
#
# So the run is split into two phases:
#
#   PHASE 1 (serial, ONE process on ONE node)
#       Pkg.instantiate() + Pkg.precompile(), and optionally a sysimage build.
#       Nothing else in the allocation is running. All the compilation happens
#       exactly once.
#
#   PHASE 2 (parallel, all ranks)
#       Launched with precompilation switched OFF at the Julia level:
#           JULIA_PKG_PRECOMPILE_AUTO=0    — Pkg never auto-precompiles
#           --compiled-modules=existing    — load existing caches, write none
#           --pkgimages=existing           — same, for the native-code images
#       A rank that somehow still misses a cache falls back to compiling in
#       its own memory; it can no longer write to the shared depot and it can
#       no longer race the other 511 ranks for the same file.
#
# The idle-node cost of phase 1 is real (you pay N-nodes x precompile-time),
# and it is still far cheaper than the alternative. If your precompile is very
# slow, run phase 1 once from a login node or a small interactive allocation
# and submit with PRECOMPILE=0.
#
# ---------------------------------------------------------------------------
# KNOBS (all overridable via `sbatch --export=ALL,NAME=value,...`)
# ---------------------------------------------------------------------------
#   EQS, CASE          problem directory and case, i.e. problems/$EQS/$CASE
#   JEXPRESSO_ROOT     repo root (default: the directory you submitted from)
#   PRECOMPILE         1 = run phase 1 (default), 0 = skip it
#   PRECOMPILE_TASKS   Julia's internal precompile worker count (default 8)
#   BUILD_SYSIMAGE     1 = build jexpresso.so in phase 1 (default 0)
#   SYSIMAGE           sysimage to use if present (default jexpresso.so)
#   LAUNCHER           auto | srun | mpirun    (default auto)
#   SRUN_MPI           value for srun --mpi=   (e.g. pmix, pmi2; default: site)
#   TAG_OUTPUT         1 = prefix each line with its rank (default 1)
#   JULIA              julia to use (default: whatever is on PATH)
#
# ============================================================================

# ---------------------------------------------------------------------------
# Slurm resource request — EDIT THIS BLOCK for your site and your case.
# ---------------------------------------------------------------------------
#SBATCH --job-name=jexpresso
#SBATCH --output=%x.%j.out          # %x = job name, %j = job id
#SBATCH --error=%x.%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=64        # = MPI ranks per node
#SBATCH --cpus-per-task=1           # Jexpresso is pure MPI: 1 core per rank
#SBATCH --time=24:00:00             # HH:MM:SS (or D-HH:MM:SS)
#SBATCH --exclusive
#SBATCH --hint=nomultithread        # ranks on physical cores, not SMT siblings

## Site-specific — uncomment and fill in:
##SBATCH --account=CHANGE_ME
##SBATCH --partition=CHANGE_ME
##SBATCH --qos=CHANGE_ME
##SBATCH --mem-per-cpu=4000M

# Fail loudly and early: an unset variable or a failed precompile must not
# silently fall through into a 512-rank launch.
set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Site environment — EDIT: load the same modules you built Jexpresso with.
#
#    The MPI module here MUST be the one MPIPreferences is bound to. If they
#    differ, the ranks never form a common world and the job hangs inside
#    MPI_Init with no diagnostic (INSTALL.md section 5.2).
# ---------------------------------------------------------------------------
module load Julia            2>/dev/null || true
module load openmpi          2>/dev/null || true

EQS="${EQS:-CompEuler}"
CASE="${CASE:-theta}"
JEXPRESSO_ROOT="${JEXPRESSO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
PRECOMPILE="${PRECOMPILE:-1}"
PRECOMPILE_TASKS="${PRECOMPILE_TASKS:-8}"
BUILD_SYSIMAGE="${BUILD_SYSIMAGE:-0}"
SYSIMAGE="${SYSIMAGE:-jexpresso.so}"
LAUNCHER="${LAUNCHER:-auto}"
TAG_OUTPUT="${TAG_OUTPUT:-1}"
JULIA="${JULIA:-julia}"

cd "$JEXPRESSO_ROOT"

if [ ! -f src/Jexpresso.jl ]; then
    echo "ERROR: $JEXPRESSO_ROOT is not a Jexpresso checkout (no src/Jexpresso.jl)." >&2
    echo "       Submit from the repo root, or pass JEXPRESSO_ROOT=/path/to/Jexpresso." >&2
    exit 1
fi

if [ ! -d "problems/$EQS/$CASE" ]; then
    echo "ERROR: no such case: problems/$EQS/$CASE" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 1. Threads. Jexpresso is essentially pure MPI (one @threads loop, in matrix
#    assembly), so each rank gets exactly one core. Without this, N ranks per
#    node each spawn N BLAS threads and fight over the same cores — the
#    classic "more nodes, slower run" oversubscription.
# ---------------------------------------------------------------------------
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# ---------------------------------------------------------------------------
# 2. Depot. Every rank reads the precompile caches from here, so it must be on
#    a filesystem that all compute nodes can see, and phase 1 must be able to
#    write to it.
#
#    Do NOT try to stage the depot to node-local scratch: Julia's cache files
#    record the absolute paths of the sources they were built from, so a depot
#    at a different path is a depot with no valid caches — every rank would
#    recompile, which is precisely what this script exists to prevent.
#
#    JULIA_CPU_TARGET matters on heterogeneous partitions: caches built on one
#    microarchitecture are rejected on another, and the rejection happens
#    per-rank at load time. Pin a multiversioned target if your nodes differ.
# ---------------------------------------------------------------------------
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"
# export JULIA_CPU_TARGET="generic;skylake-avx512,clone_all;znver3,clone_all"

# Resolve the real julia binary. A juliaup shim is a shell wrapper, and some
# launchers cannot exec it on the compute nodes.
JULIA_BIN="$("$JULIA" -e 'print(joinpath(Sys.BINDIR, "julia"))')"

echo "============================================================"
echo " Jexpresso parallel run"
echo "   job            : ${SLURM_JOB_NAME:-?} (${SLURM_JOB_ID:-?})"
echo "   case           : $EQS / $CASE"
echo "   root           : $JEXPRESSO_ROOT"
echo "   nodes x ranks  : ${SLURM_JOB_NUM_NODES:-?} x ${SLURM_NTASKS:-?}"
echo "   julia          : $JULIA_BIN ($("$JULIA_BIN" --version))"
echo "   depot          : $JULIA_DEPOT_PATH"
echo "   started        : $(date)"
echo "============================================================"

# ===========================================================================
# PHASE 1 — SERIAL PRECOMPILATION (one process, one node)
#
# This runs in the batch step itself, i.e. as a single process on the first
# node of the allocation. No srun, no ranks: there is exactly one compiler.
# PRECOMPILE_TASKS bounds Julia's own worker pool for the pass, so it stays a
# handful of cores on one node rather than the whole machine.
# ===========================================================================
if [ "$PRECOMPILE" = "1" ]; then
    echo
    echo "--- Phase 1: serial instantiate + precompile (single process) ---"
    t0=$SECONDS

    # Never ask for more precompile workers than the node has cores.
    if [ "${SLURM_CPUS_ON_NODE:-0}" -gt 0 ] && \
       [ "$PRECOMPILE_TASKS" -gt "$SLURM_CPUS_ON_NODE" ]; then
        PRECOMPILE_TASKS="$SLURM_CPUS_ON_NODE"
    fi
    export JULIA_NUM_PRECOMPILE_TASKS="$PRECOMPILE_TASKS"
    echo "    (single process, JULIA_NUM_PRECOMPILE_TASKS=$PRECOMPILE_TASKS)"

    # Auto-precompile off during instantiate so the pass below is the only
    # one, and so a failure is attributable (INSTALL.md step 3a).
    JULIA_PKG_PRECOMPILE_AUTO=0 "$JULIA_BIN" --project=. --startup-file=no \
        -e 'using Pkg; Pkg.instantiate()'

    "$JULIA_BIN" --project=. --startup-file=no \
        -e 'using Pkg; Pkg.precompile()'

    echo "--- Phase 1 complete in $((SECONDS - t0)) s ---"

    # Optional: bake everything into a sysimage. This removes the remaining
    # per-rank JIT cost (tens of seconds x every rank x every launch), at the
    # price of several minutes here. See INSTALL.md section 6.1 — the build
    # requires `using Revise` to be commented out in src/Jexpresso.jl.
    if [ "$BUILD_SYSIMAGE" = "1" ]; then
        echo "--- Phase 1b: building sysimage $SYSIMAGE (several minutes) ---"
        t0=$SECONDS
        "$JULIA_BIN" --project=. --startup-file=no create_Jexpresso_sysimage.jl
        echo "--- Sysimage built in $((SECONDS - t0)) s ---"
    fi
else
    echo
    echo "--- Phase 1 skipped (PRECOMPILE=0): assuming the depot is already warm ---"
fi

# ===========================================================================
# PHASE 2 — PARALLEL RUN
# ===========================================================================

# The guard rail. Anything that would make a rank compile-and-write is off:
#   JULIA_PKG_PRECOMPILE_AUTO=0  Pkg never kicks off a precompile pass
#   --compiled-modules=existing  load existing .ji caches, never create them
#   --pkgimages=existing         same for the native-code images
# A rank with a missing cache now compiles privately in RAM instead of racing
# the others for the same file on a shared filesystem.
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=1

JULIA_FLAGS=(--project=. --startup-file=no)

# Probe rather than assume: these two flags gained the `existing` value in
# Julia 1.11, and an unrecognised flag would kill the launch on every rank.
if "$JULIA_BIN" --compiled-modules=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--compiled-modules=existing)
else
    echo "WARNING: this Julia does not support --compiled-modules=existing;" >&2
    echo "         relying on JULIA_PKG_PRECOMPILE_AUTO=0 alone." >&2
fi
if "$JULIA_BIN" --pkgimages=existing -e 'exit(0)' >/dev/null 2>&1; then
    JULIA_FLAGS+=(--pkgimages=existing)
fi

if [ -f "$SYSIMAGE" ]; then
    echo "==> Using sysimage $SYSIMAGE"
    JULIA_FLAGS+=(--sysimage "$SYSIMAGE")
else
    echo "==> No sysimage ($SYSIMAGE not found). Every rank JIT-compiles the"
    echo "    solver on startup. Set BUILD_SYSIMAGE=1 to build one in phase 1."
fi

# ---------------------------------------------------------------------------
# Launcher choice.
#
# srun is the right answer when MPI.jl is bound to the site's MPI ("system"),
# which is what you want on a cluster: it inherits the allocation, binds ranks
# and cleans up after a failure.
#
# If MPI.jl is bound to a bundled JLL binary, srun's PMI cannot talk to it —
# that build must be launched by its own mpiexec, with an explicit hostfile.
# It also means MPI.jl is NOT using the cluster's high-speed fabric, so treat
# this branch as a warning sign rather than a supported configuration.
# ---------------------------------------------------------------------------
MPI_BINARY="$("$JULIA_BIN" --project=. --startup-file=no \
    -e 'using MPIPreferences; print(MPIPreferences.binary)' 2>/dev/null || echo unknown)"
echo "==> MPIPreferences.binary = $MPI_BINARY"

if [ "$LAUNCHER" = "auto" ]; then
    if [ "$MPI_BINARY" = "system" ]; then
        LAUNCHER=srun
    else
        LAUNCHER=mpirun
        echo "WARNING: MPI.jl is not bound to the system MPI ($MPI_BINARY)." >&2
        echo "         Falling back to mpiexec + hostfile. On a cluster you almost" >&2
        echo "         certainly want MPIPreferences.use_system_binary() instead —" >&2
        echo "         see INSTALL.md section 5.2." >&2
    fi
fi

NTASKS="${SLURM_NTASKS:?SLURM_NTASKS unset — are you running under sbatch?}"

echo
echo "--- Phase 2: parallel run on $NTASKS ranks via $LAUNCHER ---"
echo "    $(date)"
t0=$SECONDS

case "$LAUNCHER" in
    srun)
        SRUN_FLAGS=(--ntasks="$NTASKS" --cpu-bind=cores)
        [ -n "${SRUN_MPI:-}" ] && SRUN_FLAGS+=(--mpi="$SRUN_MPI")
        # --label prefixes every line with its rank. Under a launcher stdout is
        # a block-buffered pipe, so an untagged parallel run looks like a hang
        # and then dumps everything at once.
        [ "$TAG_OUTPUT" = "1" ] && SRUN_FLAGS+=(--label)

        srun "${SRUN_FLAGS[@]}" \
            "$JULIA_BIN" "${JULIA_FLAGS[@]}" src/Jexpresso.jl "$EQS" "$CASE"
        ;;

    mpirun)
        NODEFILE="nodefile.${SLURM_JOB_ID:-$$}"
        scontrol show hostnames "${SLURM_JOB_NODELIST:?SLURM_JOB_NODELIST unset}" > "$NODEFILE"
        trap 'rm -f "$NODEFILE"' EXIT

        MPIEXEC="${MPIEXEC:-$("$JULIA_BIN" --project=. --startup-file=no -e '
            using MPI
            try
                print(first(MPI.mpiexec().exec))
            catch
                print("")
            end' 2>/dev/null)}"
        [ -z "$MPIEXEC" ] && MPIEXEC="$(command -v mpiexec || command -v mpirun || true)"
        if [ -z "$MPIEXEC" ]; then
            echo "ERROR: no mpiexec/mpirun found. Set MPIEXEC=/path/to/mpiexec." >&2
            exit 1
        fi
        echo "    mpiexec: $MPIEXEC"

        MPI_FLAGS=(-np "$NTASKS" -hostfile "$NODEFILE")
        if [ "$TAG_OUTPUT" = "1" ]; then
            case "$("$MPIEXEC" --version 2>&1 | head -2 | tr '[:upper:]' '[:lower:]')" in
                *"open mpi"*|*openmpi*|*prterun*) MPI_FLAGS+=(--output tag) ;;
                *mpich*|*hydra*)                  MPI_FLAGS+=(-prepend-rank) ;;
            esac
        fi

        "$MPIEXEC" "${MPI_FLAGS[@]}" \
            "$JULIA_BIN" "${JULIA_FLAGS[@]}" src/Jexpresso.jl "$EQS" "$CASE"
        ;;

    *)
        echo "ERROR: unknown LAUNCHER='$LAUNCHER' (expected auto, srun or mpirun)." >&2
        exit 1
        ;;
esac

echo
echo "============================================================"
echo " Phase 2 finished in $((SECONDS - t0)) s"
echo " Ended: $(date)"
echo "============================================================"

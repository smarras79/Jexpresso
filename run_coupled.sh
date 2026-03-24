#!/bin/bash
# run_coupled.sh — Launch Alya + Jexpresso in MPMD mode locally.
#
# USAGE
#   ./run_coupled.sh [alya_procs] [julia_procs]
#
# EXAMPLES
#   ./run_coupled.sh               # 1 Alya rank + 2 Julia ranks (default)
#   ./run_coupled.sh 1 4           # 1 Alya rank + 4 Julia ranks
#   REBUILD_SYSIMAGE=1 ./run_coupled.sh   # force sysimage rebuild first
#
# REQUIRES
#   ./AlyaProxy/Alya.x  — compiled Alya proxy  (cd AlyaProxy && bash compilef90.sh)
#   ./jexpresso.so      — sysimage (auto-built on first run)
#
set -e

EQUATIONS="CompEuler"
CASE="thetaAlya"
ALYA_PROCS="${1:-1}"
JULIA_PROCS="${2:-2}"
SYSIMAGE="./jexpresso.so"
ALYA_EXE="./AlyaProxy/Alya.x"

if [ ! -f "$ALYA_EXE" ]; then
    echo "ERROR: Alya proxy not compiled. Run:  cd AlyaProxy && bash compilef90.sh" >&2
    exit 1
fi

# ── Step 1: Ensure all packages are installed and precompiled ─────────────────
# Without this, PackageCompiler spends hours downloading/compiling packages
# that are missing, and stdlib deps (SparseArrays etc.) may not resolve.
echo "==> Instantiating Julia project (installs any missing packages) ..."
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
echo "==> Project ready."

# ── Step 2: Build sysimage once (or on demand) ────────────────────────────────
# The warmup runs the standalone theta case for 1 timestep.  This precompiles
# sem_setup, initialize, params_setup, and the ODE solver — the four heaviest
# JIT targets.  Coupling-specific functions (je_receive_alya_data, etc.) add
# only a small fraction of compilation time and are omitted from the warmup.
#
# make_sysimage.jl exits 0 in two cases:
#   a) sysimage built successfully  → jexpresso.so now exists
#   b) not enough RAM to build      → jexpresso.so NOT created, run without it
# It exits 1 only on a hard error (build crashed, bad args, etc.).
if [ ! -f "$SYSIMAGE" ] || [ "${REBUILD_SYSIMAGE:-0}" = "1" ]; then
    echo "==> Building sysimage for $EQUATIONS/theta (one-time, ~5-15 min) ..."
    # Use standalone 'theta' for the warmup — it exercises the same code paths
    # (same mesh, same equation type, same backends) without needing Alya to
    # be running.  The sysimage is valid for both standalone and coupled runs.
    #
    # Use 'if' so that a graceful skip (exit 0, no file) doesn't abort the
    # script under set -e.
    if julia --project=. scripts/make_sysimage.jl "$EQUATIONS" theta; then
        if [ -f "$SYSIMAGE" ]; then
            echo "==> Sysimage ready: $SYSIMAGE"
        fi
        # If script exited 0 but no file was created, it means RAM was too low
        # and the build was deliberately skipped.  We continue without sysimage.
    else
        echo "==> WARNING: Sysimage build failed — continuing without --sysimage (slower startup)"
    fi
fi

# ── Step 3: Launch MPMD ───────────────────────────────────────────────────────
# The ':' separator creates one MPI_COMM_WORLD spanning both executables.
# Alya splits with color=1, Jexpresso splits with color=2 (APPID default),
# so each code sees only its own ranks via its local communicator.

# Prevent Julia MPI ranks from re-precompiling on startup.
# Pkg.precompile() above already built a valid cache; the sysimage covers
# everything else.  Without this, all Julia ranks race to recompile
# Jexpresso simultaneously, adding ~25 s of latency and causing Alya to
# hang at the first MPI collective waiting for slow Julia ranks.
export JULIA_PKG_PRECOMPILE_AUTO=0

# Use sysimage only if it was successfully built.
if [ -f "$SYSIMAGE" ]; then
    SYSIMAGE_FLAG="--sysimage $SYSIMAGE"
    echo "==> Launching: $ALYA_PROCS Alya rank(s) + $JULIA_PROCS Julia rank(s) (with sysimage)"
else
    SYSIMAGE_FLAG=""
    echo "==> Launching: $ALYA_PROCS Alya rank(s) + $JULIA_PROCS Julia rank(s) (no sysimage — first run will be slower)"
fi

# shellcheck disable=SC2086  # SYSIMAGE_FLAG is intentionally word-split
mpirun \
    -n "$ALYA_PROCS" "$ALYA_EXE" \
    : \
    -n "$JULIA_PROCS" \
        julia $SYSIMAGE_FLAG --project=. \
        run_jexpresso.jl "$EQUATIONS" "$CASE"

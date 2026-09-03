#!/bin/bash
#=============================================================================
#  LESICP2-4col-imex on 4 MPI ranks, on a laptop. Run from the REPO ROOT:
#
#      ./problems/CompEuler/LESICP2-4col-imex/run_4ranks.sh
#      DBG_TEND=10 ./problems/CompEuler/LESICP2-4col-imex/run_4ranks.sh
#      DBG_VDIFF=0 ./problems/CompEuler/LESICP2-4col-imex/run_4ranks.sh   # Schur arm
#
#  WHY THIS SCRIPT EXISTS RATHER THAN A BARE `mpiexec -n 4 julia ...`.
#  MPI.jl here is bound to MPICH_jll, not to a system MPI, and the JLL's
#  mpiexec cannot be run straight off its artifact path -- it fails with
#
#      dyld: Library not loaded: @rpath/libhwloc.15.dylib
#
#  because the launcher needs the JLL's own LIBPATH. MPI.mpiexec() sets that
#  environment up and hands you the executable inside a do-block, which is what
#  the julia -e below is doing. A homebrew mpiexec may happen to work (both are
#  MPICH) but it is a different install from the library the ranks load.
#=============================================================================
set -eu

cd "$(dirname "$0")/../../.." || exit 1

# One thread per rank: 4 ranks already fill a laptop, and nested threading
# oversubscribes it.
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=0

NR="${NRANKS:-4}"

echo "--- Precompile (serial, so 4 ranks do not compile the same code at once) ---"
JULIA_PKG_PRECOMPILE_AUTO=1 julia --project=. --startup-file=no \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- Launching $NR ranks: CompEuler / LESICP2-4col-imex ---"
julia --project=. --startup-file=no -e '
    using MPI
    nr = ARGS[1]
    MPI.mpiexec() do exe
        run(`$exe -n $nr julia --project=. --startup-file=no src/Jexpresso.jl CompEuler LESICP2-4col-imex`)
    end' "$NR"

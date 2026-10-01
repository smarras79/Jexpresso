#!/bin/bash -l
#SBATCH --job-name=ppb3d_mpi
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=71:59:00
#SBATCH --mem-per-cpu=4000M
# =============================================================================
#  MPI version of the 3D periodic Poisson solver comparison (distributed
#  memory): MUMPS (sparse direct), hypre BoomerAMG + CG and Jacobi-CG, on all
#  MPI ranks, one solver after the other:
#      sbatch tools/poisson3d_benchmark/slurm/run_wulver3d_mpi.sbatch
#  More nodes: raise --nodes (and NP); the MUMPS sizes can then grow.
#  Resubmitting resumes: finished configurations are skipped.
# =============================================================================
module load Julia/1.11.9
module load GCC MPICH
cd /project/smarras/smarras/Jexpresso/        # checkout of sm/elementLearning

export JULIA_NUM_THREADS=1          # one thread per MPI rank
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

MPIENV=tools/poisson3d_benchmark/mpi          # Julia environment: MPI, HYPRE, MUMPS (no Jexpresso)
OUTDIR=ppb3d_wulver_mpi             # results, figures, per-solver logs
NP=${SLURM_NTASKS:-128}             # MPI ranks
SOLVERS="boomeramg jacobi mumps"
# Elements per direction; unknowns n = (ne*N)^3, from 4.1e3 to 5.7e7.
# MUMPS stops at n = 2.1e6 on one node (the next size needs about 400 GB for
# the factor alone); the iterative solvers go on to 5.7e7.
NES_N2="8,12,16,24,32,48,64,96,128,192";  NES_N2_DIRECT="8,12,16,24,32,48,64"
NES_N4="4,6,8,12,16,24,32,48,64,96";      NES_N4_DIRECT="4,6,8,12,16,24,32"

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1
# Install the packages of $MPIENV first, without compiling, so that the
# preference is still set before any compilation. Its manifest is resolved
# afresh every time (it is generated, not in git): MPICH_jll must stay at 4.3
# (Project.toml compat) to match Wulver's MPICH 4.3; a stale manifest with
# MPICH_jll 5 fails with "libmpifort.so: undefined symbol: MPIR_fortran_false".
rm -f $MPIENV/Manifest.toml
JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=$MPIENV -e 'using Pkg; Pkg.instantiate()' || exit 1
julia --project=$MPIENV -e 'using MPIPreferences; MPIPreferences.use_system_binary()' || exit 1
echo "--- 2. Serial precompile (one process, many cores internally) ---"
# resolve first: an existing Manifest.toml may predate packages added to Project.toml
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()' || exit 1
julia --project=$MPIENV -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()' || exit 1
echo "--- 3. Serial warm-up, then checks against the serial (Jexpresso) benchmark ---"
julia --project=$MPIENV -e 'using MPI, HYPRE, MUMPS' || exit 1
julia --project=. tools/poisson3d_benchmark/bench3d.jl --solvers sem --nes 4,6 --nops 2,3,4 \
      --outdir "$OUTDIR/reference" || exit 1
mpirun -np 8 julia --project=$MPIENV tools/poisson3d_benchmark/mpi/verify_mpi.jl \
       --reference "$OUTDIR/reference/results.csv" || exit 1

echo "--- Setup complete, launching $NP ranks ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile
mkdir -p "$OUTDIR/logs"
for S in $SOLVERS; do
    case $S in mumps) n2=$NES_N2_DIRECT; n4=$NES_N4_DIRECT ;; *) n2=$NES_N2; n4=$NES_N4 ;; esac
    for sweep in "2 $n2" "4 $n4"; do
        set -- $sweep
        mpirun -np $NP julia --project=$MPIENV tools/poisson3d_benchmark/mpi/bench3d_mpi.jl \
               --solver $S --nop $1 --nes $2 --outdir "$OUTDIR/parts/$S" --resume \
               2>&1 | tee -a "$OUTDIR/logs/$S.log" | grep --line-buffered "^P3DMPI  $S"
        # (full output in the log; the result line of each configuration also here)
    done
done

echo "--- Merging results and drawing the figures ---"
python3 tools/poisson3d_benchmark/plot3d.py "$OUTDIR"
echo "--- Summary ---"
cat "$OUTDIR/results.md"
missing=0                           # configurations without a row: the run was killed (memory, time)
for S in $SOLVERS; do
    case $S in mumps) n2=$NES_N2_DIRECT; n4=$NES_N4_DIRECT ;; *) n2=$NES_N2; n4=$NES_N4 ;; esac
    for sweep in "2 $n2" "4 $n4"; do
        set -- $sweep
        for ne in ${2//,/ }; do
            awk -F, -v s=$S -v ne=$ne -v N=$1 'NR>1 && $1==s && $4==ne && $5==N {f=1} END {exit !f}' "$OUTDIR/results.csv" ||
                { echo "MISSING: $S ne=$ne N=$1 (no result: see $OUTDIR/logs/$S.log)"; missing=$((missing+1)); }
        done
    done
done
echo "missing configurations: $missing"
echo "results: $OUTDIR/results.md, figures: $OUTDIR/assets"

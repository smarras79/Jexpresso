#!/bin/bash -l
#SBATCH --job-name=smoothVortexScan
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks=128                 # = SV_NP × SV_JOBS = 16 × 4
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4000M

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # avoid thread oversubscription
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
julia --project=. -e 'using MPI; using Jexpresso' 2>/dev/null || \
    julia --project=. -e 'include("src/Jexpresso.jl")' --warmup-only

echo "--- Setup complete, launching scan: 4 concurrent jobs x 16 ranks ---"
export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

export MPIEXEC="srun --exclusive --mpi=pmi2"

# ---------------------------------------------------------------------------
# THE BOX. SV_L is its WIDTH, and it decides how far the study can go.
#
# The vortex is a Gaussian, not compactly supported, so on [-L/2, L/2]^2 the
# exact solution is not periodic: its velocity perturbation at the middle of
# an edge is (L/2)exp((1-(L/2)^2)/2)/2pi, with the OPPOSITE SIGN on the
# opposite edge, so the initial data jumps across the periodic seam by twice
# that. The jump is a discontinuity in the DATA, and no scheme of any order
# converges below it:
#
#     SV_L=10   [-5,5]^2      seam jump 4.9e-06   every order flattens at ~1e-5
#     SV_L=20   [-10,10]^2    seam jump 5.1e-22   machine zero        <-- default
#
# SV_L=20 is the box Dao & Nazarov run this accuracy test on. Set SV_L=10
# only to reproduce that floor on purpose.
# ---------------------------------------------------------------------------
export SV_L=${SV_L:-20}

# The meshes of that box. The repository ships 4,8,16,32,64 for both boxes,
# so this is a no-op unless a resolution is missing (and it needs gmsh).
SV_NELX="16 32 64" tools/smooth_vortex_mesh.sh || true

SWEEP='SV_NOPS=4 5 6 7 | SV_NELX=16 32 64 | dsgs and Galerkin | 16 ranks x 6 jobs'

echo "=== 1/2  ideal GLM-MHD smooth vortex   [$SWEEP]"
SV_CASE=MHD/smoothVortex \
SV_NOPS="4 5 6 7" SV_NELX="16 32 64" SV_VISC="dsgs none" SV_NP=16 SV_JOBS=6 \
SV_PLOT_NOPS="4 6" tools/smooth_vortex_mpi_scan.sh

# The same vortex without the magnetic field, at the amplitude that matches
# the MHD one (JEXPRESSO_EV_BETA=1): what the MHD error would be if div(B)
# cost nothing. JEXPRESSO_EV_BETA=5 is the classical Shu vortex instead; the
# two are stored and plotted separately, they never share a curve.
echo "=== 2/2  Euler control, same vortex without B   [$SWEEP]"
SV_CASE=CompEuler/smoothVortex JEXPRESSO_EV_BETA=1 \
SV_NOPS="4 5 6 7" SV_NELX="16 32 64" SV_VISC="dsgs none" SV_NP=16 SV_JOBS=6 \
SV_PLOT_NOPS="4 6" tools/smooth_vortex_mpi_scan.sh

echo "=== both sweeps finished. The figures:"
echo "    output/MHD/smoothVortex/output/convergence[_L1|_L2|_Linf]-it<n>.png"
echo "    output/CompEuler/smoothVortex/output/convergence[_L1|_L2|_Linf]-it<n>.png"

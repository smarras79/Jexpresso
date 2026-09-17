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

# ---------------------------------------------------------------------------
# THE SWEEP. These are the settings the reference figures were measured with,
# extended from P4/P6 to P4..P7. Override any of them from the command line:
#
#   SV_NOPS="4 6" sbatch submit_Jexpresso_compare.sh
#
# Δt is NOT pinned: the driver computes one step for the whole sweep from its
# finest (nelx, nop), so every order and mesh is integrated with the same
# step and the comparison is of the SPACE discretization. At nop 7 / 64
# elements that rule gives 2.86e-04.
#
# t = 0.5 is half a crossing of the box. It is not the paper's t = 0.05: this
# is a longer integration, so the errors are larger than theirs and the
# orders are what is being compared, not the absolute numbers.
# ---------------------------------------------------------------------------
export SV_NOPS=${SV_NOPS:-"3 4 5 6"}
export SV_NELX=${SV_NELX:-"16 32 64"}
export SV_VISC=${SV_VISC:-"dsgs none"}
export SV_TEND=${SV_TEND:-0.5}
# Δt PER CASE, from the deck's CFL rule, refined WITH the mesh. One step for
# the whole sweep is right when only the space discretization is compared —
# a plain Galerkin study — but wrong with the residual viscosity on: nu is
# C_R h^2 R, and at a FIXED step R stops falling once it reaches its time
# floor, so nu stalls and every RV curve bends to slope h^2 as soon as its
# spatial error drops under it. That is what caps P6 and P7 at p = 2.6 and
# 2.4 while their Galerkin curves hold 5.7 and 9.6.
#   SV_DT=2.857e-4 sbatch ...   pins one step again (the old behaviour)
export SV_DT=${SV_DT:-auto}
export SV_SOLVER=${SV_SOLVER:-ck54}
export SV_NP=${SV_NP:-16}
export SV_JOBS=${SV_JOBS:-6}
export SV_PLOT_NOPS=${SV_PLOT_NOPS:-"3 6"}   # the extremes, on their own axes

# ---------------------------------------------------------------------------
# THE DynSGS KNOBS, and how to add a curve without losing the one you have.
#
#   SV_HOLD    the startup hold: the number of steps at the beginning of a run
#              over which the viscosity is held off (default 2, the minimum the
#              BDF2 residual needs). The DynSGS excess at P6/P7 is mostly the
#              startup transient, so a longer hold is the knob that recovers
#              the rate; measured at P6/32x32, hold 20 removed 96% of it.
#   SV_CUTOFF  the smoothness cutoff on nu (default 0, off).
#
# Either one makes the run a DIFFERENT experiment, so its records are tagged
# dsgs_hold / dsgs_cut and drawn as their own curve beside the default DynSGS
# one rather than overwriting it. SV_KEEP=1 ADDS to the store instead of
# starting a fresh one, which is what puts both on the same figure:
#
#   SV_HOLD=20 SV_KEEP=1 SV_VISC=dsgs SV_NOPS="5 6" sbatch submit_Jexpresso_compare.sh
#
# Without SV_KEEP=1 the store is cleared first and the figure shows this
# sweep alone.
# ---------------------------------------------------------------------------
export SV_HOLD=${SV_HOLD:-2}
export SV_CUTOFF=${SV_CUTOFF:-0}
export SV_KEEP=${SV_KEEP:-0}

# The meshes of that box. The repository ships 4, 8, 16, 32 and 64 elements
# per side for both boxes, so this is a no-op unless a resolution is missing
# (and only then does it need gmsh).
tools/smooth_vortex_mesh.sh || true

echo "=== 1/2  ideal GLM-MHD smooth vortex"
echo "===      nops [$SV_NOPS] x nelx [$SV_NELX] x [$SV_VISC], t = $SV_TEND, box width $SV_L, dt $SV_DT"
SV_CASE=MHD/smoothVortex tools/smooth_vortex_mpi_scan.sh

# The same vortex without the magnetic field, at the amplitude that matches
# the MHD one (JEXPRESSO_EV_BETA=1): what the MHD error would be if div(B)
# cost nothing. JEXPRESSO_EV_BETA=5 is the classical Shu vortex instead; the
# two are stored and plotted separately, they never share a curve.
echo "=== 2/2  Euler control, the same vortex without B"
SV_CASE=CompEuler/smoothVortex JEXPRESSO_EV_BETA=1 tools/smooth_vortex_mpi_scan.sh

echo "=== both sweeps finished. The figures:"
echo "    output/MHD/smoothVortex/output/convergence[_L1|_L2|_Linf]-it<n>.png"
echo "    output/CompEuler/smoothVortex/output/convergence[_L1|_L2|_Linf]-it<n>.png"

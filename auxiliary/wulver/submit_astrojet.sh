#!/bin/bash -l
#SBATCH --job-name=astrojet
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:59:00
#SBATCH --mem-per-cpu=4000M
#
# problems/MHD/astroJetWuShu2018 — the classical magnetized astrophysical jet.
#
#   sbatch auxiliary/wulver/submit_astrojet.sh                      # the deck as committed
#   sbatch --export=ALL,AJ_RUNG=B auxiliary/wulver/submit_astrojet.sh
#   sbatch --export=ALL,AJ_CMAX=1.0 auxiliary/wulver/submit_astrojet.sh
#
# WHY A RUNG SELECTOR AND NOT LOOSE VARIABLES. The low-Mach rungs of
# README.md §6 need JEXPRESSO_AJ_TEND and JEXPRESSO_AJ_DT moved TOGETHER with the
# injection speed: the head speed scales with u_jet, so a rung run at the default
# tend simply integrates less of the flow, and one run at the default Δt has a
# different Courant number. Setting three variables by hand and getting one of them
# wrong gives a run that looks like a result. The case statement below moves them
# as a unit, so a rung is one word.
#
#   RUNG | u_jet | B_a²  | β_a   | ambient thermal margin | what it tests
#   -----+-------+-------+-------+------------------------+-----------------------
#     B  |  800  |     2 | 1     | 71 %                   | is the FIELD the problem?
#     1  |   20  |   200 | 1e-2  | 2.4 %                  | geometry, BCs, mesh
#     2  |   80  |   200 | 1e-2  | 2.4 %                  | a real jet
#     3  |  250  |   200 | 1e-2  | 2.4 %                  | p/ρE = 5.7e-5
#     4  |  800  |   200 | 1e-2  | 2.4 %                  | the paper's case (i)  <- default
#     5  |  800  |  2000 | 1e-3  | 0.249 %                | the paper's case (ii)
#     6  |  800  | 20000 | 1e-4  | 0.025 %                | the paper's case (iii)
#
# Rung 4 sets NOTHING: it is the deck exactly as committed, so this script with no
# arguments reproduces problems/MHD/astroJetWuShu2018/user_inputs.jl and nothing
# else. Rung B is the one to run when the repair report says the damage is in the
# ambient gas rather than in the beam — see README.md §11(c) and §12.
#
# Other knobs, all optional, all defaulting to the deck:
#   AJ_MESH   40x60 (default) | 100x150
#   AJ_CMAX   DynSGS cap coefficient; 0.5 is the deck, ~3.5 is the viscous limit
#   AJ_SMOOTH nozzle-lip transition half-width; 0 = the paper's top hat (aborts)
#   AJ_TRAMP  inflow turn-on time; 0 = the paper's impulsive start (aborts)
#   NRANKS    MPI ranks (default 64; may go to 128 without touching #SBATCH above)
#
# NOTE on --ntasks-per-node=128 with NRANKS=64: the allocation is deliberately
# larger than the launch, which doubles the memory available per rank and leaves
# room to raise NRANKS from the command line — #SBATCH directives are static
# comments and cannot read a shell variable.

module load Julia/1.11.9
module load GCC MPICH

cd /project/smarras/smarras/Jexpresso/

export JULIA_NUM_THREADS=1          # avoid 64x64 thread oversubscription
export JULIA_PKG_PRECOMPILE_AUTO=1  # allow precompile during setup only

echo "--- 1. MPI preferences FIRST (before any compilation) ---"
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

echo "--- 2. Serial precompile (one process, many cores internally) ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "--- 3. Serial warm-up: exercise the real load path ---"
julia --project=. -e 'using MPI; using Jexpresso' 2>/dev/null || \
  julia --project=. -e 'include("src/Jexpresso.jl")' --warmup-only

#------------------------------------------------------------------------------
# The case configuration.
#
# These are read AT INCLUDE TIME inside each rank's Julia process
# (problems/MHD/astroJetWuShu2018/user_flux.jl and user_inputs.jl), so they have
# to be in the RANKS' environment: exported here, before mpirun. MPICH's Hydra
# propagates the whole environment by default, so `export` is enough. OpenMPI does
# NOT — there it would need `mpirun -x JEXPRESSO_AJ_BA2 -x ...`, and the
# MPICH-explicit form is `-genv JEXPRESSO_AJ_BA2 2`.
#------------------------------------------------------------------------------
AJ_RUNG="${AJ_RUNG:-4}"

case "$AJ_RUNG" in
  B)  export JEXPRESSO_AJ_BA2=2                                             ;;
  1)  export JEXPRESSO_AJ_UJET=20   JEXPRESSO_AJ_TEND=8e-2  JEXPRESSO_AJ_DT=1e-5   ;;
  2)  export JEXPRESSO_AJ_UJET=80   JEXPRESSO_AJ_TEND=2e-2  JEXPRESSO_AJ_DT=4e-6   ;;
  3)  export JEXPRESSO_AJ_UJET=250  JEXPRESSO_AJ_TEND=6.3e-3 JEXPRESSO_AJ_DT=1.5e-6 ;;
  4)  : ;;                          # the deck as committed — set nothing
  5)  export JEXPRESSO_AJ_BA2=2000                                          ;;
  6)  export JEXPRESSO_AJ_BA2=20000                                         ;;
  *)  echo "ERROR: AJ_RUNG='$AJ_RUNG' is not one of B 1 2 3 4 5 6" >&2; exit 1 ;;
esac

# Optional overrides. Unset by default so the deck decides; each is exported only
# when the submitter asked for it, which keeps the log honest about what changed.
[ -n "${AJ_MESH:-}"   ] && export JEXPRESSO_AJ_MESH="$AJ_MESH"
[ -n "${AJ_CMAX:-}"   ] && export JEXPRESSO_AJ_CMAX="$AJ_CMAX"
[ -n "${AJ_SMOOTH:-}" ] && export JEXPRESSO_AJ_SMOOTH="$AJ_SMOOTH"
[ -n "${AJ_TRAMP:-}"  ] && export JEXPRESSO_AJ_TRAMP="$AJ_TRAMP"

# DELIBERATELY NOT SET: JEXPRESSO_DSGS_DEBUG.
#
# Its @printf in kernel/physics/SGS.jl has no rank guard, so on 64 ranks it is 64
# interleaved streams of RANK-LOCAL numbers that read like global ones. And the one
# thing it is usually wanted for — whether ν has saturated its cap — is already in
# every CFL block as "max ν", which soundSpeed.jl MPI.Allreduce(..., MAX)es and
# prints from rank 0. Compare it with C_max·Δ·(|v|+c_f): equal means the cap binds
# and C_R is irrelevant. Use DSGS_DEBUG on a short 1-rank run when the per-equation
# residual breakdown and the argmax node are actually needed.

#------------------------------------------------------------------------------
# What ran, in the log. This matters more than it looks: these runs differ only by
# environment, and a .out file that does not name its own configuration cannot be
# compared with the next one.
#------------------------------------------------------------------------------
NRANKS="${NRANKS:-64}"
echo "--- astroJetWuShu2018: rung $AJ_RUNG on $NRANKS ranks ---"
env | grep '^JEXPRESSO_' | sort || echo "  (no JEXPRESSO_* set: the deck's own defaults)"
echo "-------------------------------------------------------"

export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

mpirun -np "$NRANKS" julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018

#
# READING THE RESULT (problems/MHD/astroJetWuShu2018/README.md §10-§12):
#
#   POSITIVITY REPAIR ENGAGED ... GLOBAL first repair at (x, y) on RHS call N
#     N <= 5            the FIRST time step: a boundary datum, not propagation
#     y = 0             the inlet itself
#     y = h, x = ±0.05  the lip's shear layer
#     energy-RAISED     dominating means ρE is dipping below ½|B|², i.e. the
#                       ambient thermal margin — run rung B
#     injected energy   compare with 153.75, the domain's initial total energy.
#                       Much larger means the repair wrote the answer.
#   Viscous CFL / max ν   equal to C_max·Δ·(|v|+c_f) means the cap binds
#

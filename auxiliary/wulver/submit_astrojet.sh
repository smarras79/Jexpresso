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
# problems/MHD/astroJetWuShu2018 as committed:  sbatch [--export=ALL,NRANKS=32] auxiliary/wulver/submit_astrojet.sh
# Variants are edited in the case: B_a (aj_Ba) in user_flux.jl, mesh and Δt in user_inputs.jl.
#------------------------------------------------------------------------------
NRANKS="${NRANKS:-64}"
echo "--- astroJetWuShu2018 on $NRANKS ranks ---"

export JULIA_PKG_PRECOMPILE_AUTO=0  # ranks must never attempt to precompile

mpirun -np "$NRANKS" julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018

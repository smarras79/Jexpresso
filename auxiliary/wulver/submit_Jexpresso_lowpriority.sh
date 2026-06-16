#!/bin/bash -l
#SBATCH --job-name=LESsmago
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=low
#SBATCH --account=smarras # Replace PI_ucid which the NJIT UCID of PI
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=128
#SBATCH --time=71:59:00  # D-HH:MM:SS

#SBATCH --mem-per-cpu=4000M
ml Julia
ml bright shared mpich/ge/gcc/64
cd /project/smarras/smarras/Jexpresso/
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > nodefile

echo "--- Instantiating Julia Project ---"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
echo "--- Instantiation Complete ---"

# Same launch approach as the laptop (jexp_mpich.sh / INSTALL.md sec 5.6,
# Route A/B system MPI): the script form `src/Jexpresso.jl <EQ> <CASE>`,
# NOT `-e 'push!(ARGS,...); include(...)'`. Runtime env defaults match
# jexp_mpich.sh; override in the job environment if needed.
: "${JEXPRESSO_STEP_HEARTBEAT:=0}"
: "${JEXPRESSO_ALLOC_SUMMARY:=0}"
: "${JEXPRESSO_PRECOMPILE_WARMUP:=1}"
export JEXPRESSO_STEP_HEARTBEAT JEXPRESSO_ALLOC_SUMMARY JEXPRESSO_PRECOMPILE_WARMUP

mpirun -np 300 -hostfile nodefile julia --project=. src/Jexpresso.jl CompEuler LESsmago

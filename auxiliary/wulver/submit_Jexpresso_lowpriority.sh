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

mpirun -np 300 -hostfile nodefile julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "LESsmago"); include("src/Jexpresso.jl")'

#!/bin/bash -l
#SBATCH --job-name=LESsmago100
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=high_smarras
#SBATCH --account=smarras # Replace PI_ucid which the NJIT UCID of PI
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=71:59:00  # D-HH:MM:SS

#SBATCH --mem-per-cpu=4000M

ml Julia
module load bright shared mpich/ge/gcc/64
cd /project/smarras/smarras/Jexpresso/
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > nodefile
mpirun -np 200 -hostfile nodefile julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "LESsmago"); include("src/Jexpresso.jl")'

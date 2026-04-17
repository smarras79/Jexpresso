#!/bin/bash
#SBATCH --job-name=Cuyo
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err
#SBATCH --partition=32crs
#SBATCH --ntasks=32
#SBATCH --time=24:00:00  # D-HH:MM:SS
##SBATCH --mem=100G

module load mpich-3.1
cd /home/esalcedo/Jexpresso

JULIA="/home/esalcedo/julia-1.12.5/bin/julia"
/app/mpich-3.1/bin/mpirun -np 32 $JULIA --project=. -e 'push!(empty!(ARGS), "CompEuler", "Cuyo"); include("src/Jexpresso.jl")'

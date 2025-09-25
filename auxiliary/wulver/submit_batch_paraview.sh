#!/bin/bash -l
#SBATCH --job-name=LESsmago100
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=high_smarras
#SBATCH --account=smarras # Replace PI_ucid which the NJIT UCID of PI
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:59:00  # D-HH:MM:SS

#SBATCH --mem-per-cpu=4000M


module load bright shared mpich/ge/gcc/64 foss/2024a ParaView
cd /scratch/smarras/smarras/output/64x64x24/CompEuler/LESsmago/output/
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > nodefile
python3 batch_paraview_analysis.py


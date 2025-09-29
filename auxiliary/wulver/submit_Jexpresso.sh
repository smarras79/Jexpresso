#!/bin/bash -l
#SBATCH --job-name=LESsmago
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=high_smarras
#SBATCH --account=smarras # Replace PI_ucid which the NJIT UCID of PI
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=71:59:00  # D-HH:MM:SS

#SBATCH --mem-per-cpu=4000M

module load Julia
module load bright shared mpich/ge/gcc/64
cd /project/smarras/smarras/Jexpresso/
nodelist=$(scontrol show hostname $SLURM_NODELIST)

# Monitor in background
(
  while true; do
    echo "=== $(date) ==="
    ps aux --sort=-%mem | grep julia | head -10
    sleep 60
  done
) > memory_monitor.log 2>&1 &


printf "%s\n" "${nodelist[@]}" > nodefile
mpirun -np 256 -hostfile nodefile julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "LESsmago"); include("src/Jexpresso.jl")'

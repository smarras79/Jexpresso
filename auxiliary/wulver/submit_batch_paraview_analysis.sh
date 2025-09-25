#!/bin/bash -l
#SBATCH --job-name=ParaView_Simple
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10  # Adjust based on your needs
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=4000M

module load bright shared mpich/ge/gcc/64 foss/2024a ParaView Julia
cd /scratch/smarras/smarras/output/64x64x24/CompEuler/LESsmago/output/
#cd /scratch/smarras/smarras/output/64x64x36_5kmX5kmX3km/CompEuler/LESsmago/output/

echo "Starting simple parallel ParaView processing..."
echo "Job ID: $SLURM_JOB_ID, CPUs: $SLURM_NTASKS_PER_NODE"

#=============================================
# CUSTOMIZE THESE RANGES BASED ON YOUR FILES:
# First, run this to see what files you have:
# python3 batch_paraview_analysis.py --dry-run
#=============================================
# Then split the work across processes:
#python3 batch_paraview_analysis.py --range 100 199 1 --process-id 1 &
#python3 batch_paraview_analysis.py --range 200 299 1 --process-id 2 &
#python3 batch_paraview_analysis.py --range 300 399 1 --process-id 3 &
#python3 batch_paraview_analysis.py --range 400 499 1 --process-id 4 &
#python3 batch_paraview_analysis.py --range 500 582 1 --process-id 5 &


echo "Ensuring Julia packages are installed..."
julia -e 'using Pkg; Pkg.add(["ArgParse", "Dates"])' #<-- EXAMPLE: Add ArgParse

#julia --project batch_paraview_analysis.jl --range 1000 199 1 --process-id 1 &
#julia --project batch_paraview_analysis.jl --range 200 299 1 --process-id 2 &
#julia --project batch_paraview_analysis.jl --range 300 399 1 --process-id 3 &
#julia --project batch_paraview_analysis.jl --range 400 499 1 --process-id 4 &
#julia --project batch_paraview_analysis.jl --range 500 582 1 --process-id 5 &
julia --project batch_paraview_analysis.jl --range 500 582 1 --process-id 5 &


# Wait for all processes to complete
wait

echo "All processes completed! Check batch_output/ for results."

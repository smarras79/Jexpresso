#!/bin/bash -l
#SBATCH --job-name=ParaView_Parallel
#SBATCH --output=%x.%j.out # %x.%j expands to slurm JobName.JobID
#SBATCH --error=%x.%j.err
#SBATCH --partition=general
#SBATCH --qos=high_smarras
#SBATCH --account=smarras
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16  # Increased for parallel processing
#SBATCH --time=23:59:00  # D-HH:MM:SS
#SBATCH --mem-per-cpu=4000M

# Load required modules
module load bright shared mpich/ge/gcc/64 foss/2024a ParaView

# Change to your data directory
cd /scratch/smarras/smarras/output/64x64x24/CompEuler/LESsmago/output/

# Create nodefile for reference
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > nodefile

echo "Starting ParaView Parallel Batch Processing"
echo "============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs available: $SLURM_NTASKS_PER_NODE"
echo "Start time: $(date)"
echo "============================================="

# First, get suggestions for parallel processing
echo "Analyzing files for optimal parallel processing..."
python3 batch_paraview_analysis.py --suggest-parallel --num-processes $SLURM_NTASKS_PER_NODE

echo ""
echo "Starting parallel processing with $SLURM_NTASKS_PER_NODE processes..."

# Method 1: Auto-suggested ranges (recommended)
# Uncomment this section if you want the script to auto-determine ranges
# python3 batch_paraview_analysis.py --suggest-parallel --num-processes $SLURM_NTASKS_PER_NODE > parallel_suggestions.txt
# # Parse suggestions and launch processes (would need additional scripting)

# Method 2: Manual range specification (more predictable)
# Adjust these ranges based on your actual file numbers
# You can run with --dry-run first to see what files exist

# Example: Process files in parallel chunks
# Adjust START, END, STEP based on your actual file range
START_FILE=82
END_FILE=582
STEP=5

# Calculate range per process
TOTAL_RANGE=$((END_FILE - START_FILE))
RANGE_PER_PROCESS=$((TOTAL_RANGE / SLURM_NTASKS_PER_NODE))

echo "File range: $START_FILE to $END_FILE (step $STEP)"
echo "Range per process: ~$RANGE_PER_PROCESS"
echo ""

# Launch parallel processes
pids=()
for ((i=0; i<SLURM_NTASKS_PER_NODE; i++)); do
    process_start=$((START_FILE + i * RANGE_PER_PROCESS))
    if [ $i -eq $((SLURM_NTASKS_PER_NODE - 1)) ]; then
        # Last process gets remaining files
        process_end=$END_FILE
    else
        process_end=$((process_start + RANGE_PER_PROCESS - 1))
    fi
    
    echo "Launching Process $((i+1)): files $process_start to $process_end"
    
    # Launch background process
    python3 batch_paraview_analysis.py \
        --range $process_start $process_end $STEP \
        --process-id $((i+1)) \
        --log-level INFO &
    
    # Store process ID
    pids+=($!)
    
    # Small delay to avoid resource conflicts
    sleep 2
done

echo ""
echo "All $SLURM_NTASKS_PER_NODE processes launched. Process IDs: ${pids[*]}"
echo "Waiting for all processes to complete..."

# Wait for all background processes
for pid in "${pids[@]}"; do
    wait $pid
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "Process $pid completed successfully"
    else
        echo "Process $pid failed with exit code $exit_code"
    fi
done

echo ""
echo "============================================="
echo "All parallel processing completed!"
echo "End time: $(date)"

# Collect summary information
echo ""
echo "Processing Summary:"
echo "==================="

# Count output files
if [ -d "batch_output" ]; then
    total_outputs=$(find batch_output -name "*.png" | wc -l)
    echo "Total PNG files created: $total_outputs"
    
    # Show some examples
    echo "Sample output files:"
    ls batch_output/*.png | head -5
    
    if [ $total_outputs -gt 5 ]; then
        echo "... and $((total_outputs - 5)) more files"
    fi
else
    echo "No batch_output directory found"
fi

# Show log files
echo ""
echo "Log files created:"
ls -la batch_processing_proc_*.log 2>/dev/null || echo "No parallel log files found"

# Check for any failures
failed_logs=$(grep -l "Failed:" batch_processing_proc_*.log 2>/dev/null | wc -l)
if [ $failed_logs -gt 0 ]; then
    echo ""
    echo "WARNING: $failed_logs processes reported failures"
    echo "Check individual log files for details"
fi

echo ""
echo "Job completed. Check output files in batch_output/ directory."
echo "============================================="

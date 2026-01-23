#!/bin/bash
# test_continuous_coupling.sh
# Test continuous execution with both codes running in parallel

echo "=========================================="
echo "Testing Continuous Coupled Execution"
echo "=========================================="
echo ""
echo "This test runs both codes for 10 seconds"
echo "Both codes will:"
echo "  1. Participate in MPI_Comm_split"
echo "  2. Participate in MPI_Gather"
echo "  3. Continue running simulation steps"
echo "  4. Synchronize at barriers"
echo "  5. Finalize together"
echo ""

# Check if Fortran continuous test exists
if [ ! -f "./alya_gather_continuous" ]; then
    echo "Compiling Fortran continuous test..."
    mpif90 -cpp -DUSEMPIF08 alya_gather_continuous.f90 -o alya_gather_continuous
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to compile alya_gather_continuous.f90"
        exit 1
    fi
    echo "✓ Compiled successfully"
    echo ""
fi

echo "Running coupled simulation:"
echo "  - 2 ranks: Alya (ranks 0-1)"
echo "  - 2 ranks: Jexpresso (ranks 2-3)"
echo "  - Duration: 10 seconds"
echo ""

# Run the coupled simulation
# Both codes will run for 10 seconds then exit cleanly
mpirun --tag-output \
    -np 2 ./alya_gather_continuous : \
    -np 2 julia --project=. jexpresso_gather_continuous.jl --code-name "Jexpresso" --max-time 10.0

echo ""
echo "=========================================="
echo "Test completed"
echo "=========================================="

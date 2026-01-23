#!/bin/bash
# test_jexpresso_intercomm.sh
# Test Jexpresso with intercommunicator coupling mode

echo "=========================================="
echo "Testing Jexpresso with Intercommunicator Coupling"
echo "=========================================="
echo ""
echo "This test demonstrates truly independent execution:"
echo "- Alya and Jexpresso run independently"
echo "- Each has its own communicator for internal MPI operations"
echo "- Communication only when needed via intercommunicator"
echo ""

# Check if alya_mini_coupler exists
if [ ! -f "./alya_mini_coupler" ]; then
    echo "Compiling alya_mini_coupler.f90..."
    mpif90 alya_mini_coupler.f90 -o alya_mini_coupler
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to compile alya_mini_coupler.f90"
        exit 1
    fi
    echo "✓ Compiled successfully"
    echo ""
fi

echo "Test 1: Minimal intercommunicator test (Jexpresso-mini)"
echo "--------------------------------------------------------"
echo "Running:"
echo "  - 2 ranks: Alya (APPID=0)"
echo "  - 2 ranks: Jexpresso-mini (APPID=1)"
echo ""

mpirun --tag-output \
    -np 2 -x APPID=0 ./alya_mini_coupler : \
    -np 2 -x APPID=1 julia --project=. Jexpresso-mini-coupled.jl

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Minimal test failed"
    exit 1
fi

echo ""
echo "✓ Minimal test passed"
echo ""
echo "=========================================="
echo ""

echo "Test 2: Full Jexpresso with intercommunicator coupling"
echo "-------------------------------------------------------"
echo "Running:"
echo "  - 2 ranks: Alya (APPID=0)"
echo "  - 2 ranks: Jexpresso full simulation (APPID=1)"
echo ""

mpirun --tag-output \
    -np 2 -x APPID=0 ./alya_mini_coupler : \
    -np 2 -x APPID=1 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --intercomm-coupling --code-name "Jexpresso"

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Full Jexpresso test failed"
    exit 1
fi

echo ""
echo "✓ Full Jexpresso test passed"
echo ""
echo "=========================================="
echo "All tests completed successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Modify your Alya code to use intercommunicator pattern"
echo "2. Run with real Alya: mpirun -np N -x APPID=0 ./alya/Alya.x : \\"
echo "                               -np M -x APPID=1 julia src/Jexpresso.jl ... --intercomm-coupling"
echo ""

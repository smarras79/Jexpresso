#!/bin/bash
# test_jexpresso_alya_gather.sh
# Test Jexpresso with gather-based coupling mode
# Usage: ./test_jexpresso_alya_gather.sh

echo "=========================================="
echo "Testing Jexpresso with gather-based coupling"
echo "=========================================="

# Check if Alya.x exists
if [ ! -f "./alya/Alya.x" ]; then
    echo "ERROR: Alya.x not found in ./alya/"
    echo "Please compile your Alya code first or update the path"
    exit 1
fi

# Check if Jexpresso is accessible
if [ ! -f "src/Jexpresso.jl" ]; then
    echo "ERROR: src/Jexpresso.jl not found"
    echo "Please run this script from the Jexpresso root directory"
    exit 1
fi

echo ""
echo "Running coupled simulation:"
echo "  - 2 ranks: Alya (ranks 0-1)"
echo "  - 2 ranks: Jexpresso (ranks 2-3)"
echo ""

# Run the coupled simulation
# Alya gets ranks 0-1, Jexpresso gets ranks 2-3
# Using --coupling-test-only to exit after coupling initialization
mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false --gather-coupling --coupling-test-only --code-name "Jexpresso"

echo ""
echo "=========================================="
echo "Test completed"
echo "=========================================="

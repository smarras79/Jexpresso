#!/usr/bin/env julia
#
# Test launcher: Couple two Jexpresso instances together
#
# This script demonstrates and tests MPI coupling by running two separate
# Jexpresso instances that communicate with each other.
#
# Usage:
#   mpirun -np 4 julia --project=. test_jexpresso_coupling.jl
#
# This will run:
#   - Ranks 0-1: Jexpresso-A (code_id=1)
#   - Ranks 2-3: Jexpresso-B (code_id=2)
#

using MPI

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

# Verify we have at least 2 ranks
if world_size < 2
    if world_rank == 0
        println("ERROR: Need at least 2 MPI ranks to test coupling")
        println("Usage: mpirun -np 4 julia --project=. test_jexpresso_coupling.jl")
    end
    MPI.Finalize()
    exit(1)
end

# Split ranks: first half = Jexpresso-A, second half = Jexpresso-B
jexpresso_a_ranks = world_size ÷ 2
code_id = (world_rank < jexpresso_a_ranks) ? 1 : 2
code_name = (code_id == 1) ? "Jexpresso-A" : "Jexpresso-B"

# Print startup info
println("═"^60)
println("World Rank $world_rank: Starting $code_name (code_id=$code_id)")
println("  Total world ranks: $world_size")
println("  Jexpresso-A ranks: 0-$(jexpresso_a_ranks-1)")
println("  Jexpresso-B ranks: $(jexpresso_a_ranks)-$(world_size-1)")
println("═"^60)
flush(stdout)

# Small delay to ensure output is ordered
sleep(0.1 * world_rank)

# Set command-line arguments for Jexpresso
# Both instances run the same case but with coupling enabled
push!(empty!(ARGS),
      "--coupling",
      "--code-id", string(code_id),
      "--n-codes", "2",
      "--code-name", code_name,
      "CompEuler",  # Change to your preferred equation type
      "theta")      # Change to your preferred test case

# Load and run Jexpresso
println("[$code_name, World Rank $world_rank] Loading Jexpresso module...")
flush(stdout)

using Jexpresso

# Note: Jexpresso will handle MPI.Finalize() in its cleanup

#!/usr/bin/env julia
#
# Example launcher for coupled Jexpresso with external code
#
# This script demonstrates how to launch Jexpresso coupled with an external
# code from a single mpirun command. Ranks are split between the two codes.
#
# Usage:
#   mpirun -np 8 julia --project=. launch_coupled.jl
#
# This will run:
#   - Ranks 0-3: Jexpresso
#   - Ranks 4-7: External code (you need to implement this part)
#

using MPI

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

# Configure rank split (adjust these values for your needs)
jexpresso_ranks = world_size ÷ 2  # First half for Jexpresso
code_id = (world_rank < jexpresso_ranks) ? 1 : 2

if code_id == 1
    # ===== JEXPRESSO CODE =====
    if world_rank == 0
        println("Launching Jexpresso on ranks 0-$(jexpresso_ranks-1)")
    end

    # Set command-line arguments for Jexpresso
    push!(empty!(ARGS),
          "--coupling",
          "--code-id", "1",
          "--n-codes", "2",
          "--code-name", "Jexpresso",
          "CompEuler",  # Change to your equation type
          "theta")      # Change to your case name

    # Load and run Jexpresso
    using Jexpresso

else
    # ===== EXTERNAL CODE =====
    if world_rank == jexpresso_ranks
        println("Launching ExternalCode on ranks $(jexpresso_ranks)-$(world_size-1)")
    end

    # Set arguments for external code
    push!(empty!(ARGS),
          "--code-id", "2",
          "--n-codes", "2")

    # TODO: Replace this with your external code
    # Example: include("path/to/external_code.jl")

    # For demonstration, we'll just use the example coupling code
    include("src/mpi/example_coupling.jl")
end

# Note: MPI.Finalize() is called by Jexpresso and should be called by external code too

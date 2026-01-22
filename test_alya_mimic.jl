#!/usr/bin/env julia
#
# Simple test program that mimics Alya's coupling behavior
# Use this to test MPMD launch before trying with real Alya
#
# Usage with Jexpresso:
#   mpirun -np 4 julia --project=. run_jexpresso.jl \
#       --coupling --code-id 1 --n-codes 2 CompEuler theta : \
#       -np 4 julia --project=. test_alya_mimic.jl
#

using MPI

println("Starting Alya mimic...")

MPI.Init()

world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

println("[Alya mimic, rank $world_rank] Started with world_size=$world_size")

# Load JexpressoCoupling
include("src/mpi/JexpressoCoupling.jl")
using .JexpressoCoupling

# Initialize coupling as code 2
println("[Alya mimic, rank $world_rank] Initializing coupling...")

ctx = initialize_coupling(
    MPI.COMM_WORLD,
    2,  # code_id for Alya
    2;  # n_codes
    code_name="AlyaMimic"
)

println("[Alya mimic, rank $world_rank] Coupling initialized:")
println("  Local rank: $(ctx.local_rank) / $(ctx.local_size)")
println("  World rank: $(ctx.world_rank)")
println("  Is root: $(ctx.is_root)")

# Simulate Alya running
for step in 1:5
    if ctx.is_root
        println("\n[Alya mimic] Simulation step $step/5")
    end

    # Synchronize with Jexpresso
    synchronize_coupling(ctx)

    # Simulate work
    sleep(0.5)

    # Optional: Test data exchange
    if ctx.is_root && step == 3
        println("[Alya mimic] Testing data exchange...")

        # Send dummy pressure to Jexpresso
        pressure_data = rand(10) .* 1000.0
        send_field_array!(ctx, pressure_data, 1; tag=201)
        println("[Alya mimic] Sent pressure data to Jexpresso")

        # Receive dummy velocity from Jexpresso
        velocity_data = zeros(10)
        recv_field_array!(ctx, velocity_data, 1; tag=200)
        println("[Alya mimic] Received velocity data from Jexpresso")
        println("  Velocity range: $(minimum(velocity_data)) to $(maximum(velocity_data))")
    end
end

if ctx.is_root
    println("\n[Alya mimic] Simulation finished successfully")
end

# Cleanup
finalize_coupling(ctx)
MPI.Finalize()

println("[Alya mimic, rank $world_rank] Done")

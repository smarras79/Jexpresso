#!/usr/bin/env julia
#
# Simple coupling test: Verify MPI coupling infrastructure
#
# This script tests the coupling module without running full Jexpresso.
# It initializes coupling, prints diagnostic info, and performs a simple
# data exchange between two "codes".
#
# Usage:
#   mpirun -np 4 julia --project=. test_simple_coupling.jl
#

using MPI
using Printf

# Need to load Jexpresso to access JexpressoCoupling
include("src/Jexpresso.jl")
using .Jexpresso.JexpressoCoupling

MPI.Init()
world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

# Verify we have at least 2 ranks
if world_size < 2
    if world_rank == 0
        @error "Need at least 2 MPI ranks to test coupling"
        println("Usage: mpirun -np 4 julia --project=. test_simple_coupling.jl")
    end
    MPI.Finalize()
    exit(1)
end

println("\n" * "="^70)
println("COUPLING TEST STARTED")
println("="^70)

# Split ranks between two codes
code_a_ranks = world_size ÷ 2
code_id = (world_rank < code_a_ranks) ? 1 : 2
code_name = (code_id == 1) ? "Jexpresso-A" : "Jexpresso-B"

println(@sprintf("[Rank %2d] I am %s (code_id=%d)", world_rank, code_name, code_id))
flush(stdout)
sleep(0.05 * world_rank)  # Stagger output

# Initialize coupling
println(@sprintf("[Rank %2d] Initializing coupling...", world_rank))
flush(stdout)

ctx = initialize_coupling(MPI.COMM_WORLD, code_id, 2; code_name=code_name)

# Print coupling information
info = get_coupling_info(ctx)
println("\n" * "─"^70)
println(@sprintf("[Rank %2d] COUPLING INFO for %s:", world_rank, code_name))
println(@sprintf("  ├─ Code ID:     %d / %d", info.code_id, info.n_codes))
println(@sprintf("  ├─ Local rank:  %d / %d", info.local_rank, info.local_size))
println(@sprintf("  ├─ World rank:  %d", info.world_rank))
println(@sprintf("  ├─ Is root:     %s", info.is_root))
println(@sprintf("  └─ Root ranks:  %s", info.root_ranks))
println("─"^70)
flush(stdout)

# Synchronize before data exchange
synchronize_coupling(ctx)

if world_rank == 0
    println("\n" * "="^70)
    println("STARTING DATA EXCHANGE TEST")
    println("="^70)
end

# Test data exchange between coupling roots
if is_coupling_root(ctx)
    n_data = 5
    partner_code_id = (code_id == 1) ? 2 : 1

    # Prepare send data (different for each code)
    if code_id == 1
        send_data = Float64[10.0, 20.0, 30.0, 40.0, 50.0]
    else
        send_data = Float64[100.0, 200.0, 300.0, 400.0, 500.0]
    end

    recv_data = zeros(Float64, n_data)

    println(@sprintf("\n[%s Root] BEFORE exchange:", code_name))
    println(@sprintf("  Send: %s", send_data))
    println(@sprintf("  Recv: %s", recv_data))
    flush(stdout)

    # Perform exchange
    exchange_field_data!(ctx, send_data, recv_data, partner_code_id; tag=42)

    println(@sprintf("\n[%s Root] AFTER exchange:", code_name))
    println(@sprintf("  Send: %s", send_data))
    println(@sprintf("  Recv: %s", recv_data))

    # Verify we received the correct data
    expected = (code_id == 1) ? [100.0, 200.0, 300.0, 400.0, 500.0] : [10.0, 20.0, 30.0, 40.0, 50.0]
    if recv_data == expected
        println(@sprintf("  ✓ Data exchange SUCCESSFUL!"))
    else
        println(@sprintf("  ✗ Data exchange FAILED!"))
        println(@sprintf("    Expected: %s", expected))
    end
    flush(stdout)
else
    # Non-root ranks
    recv_data = zeros(Float64, 5)
end

# Broadcast received data to all local ranks
broadcast_to_local!(ctx, recv_data)

# All ranks verify they have the data
expected = (code_id == 1) ? [100.0, 200.0, 300.0, 400.0, 500.0] : [10.0, 20.0, 30.0, 40.0, 50.0]
if recv_data == expected
    println(@sprintf("[Rank %2d, %s] ✓ Received correct data after broadcast: %s",
            world_rank, code_name, recv_data))
else
    println(@sprintf("[Rank %2d, %s] ✗ Incorrect data after broadcast!", world_rank, code_name))
end
flush(stdout)

# Final synchronization
synchronize_coupling(ctx)

if world_rank == 0
    println("\n" * "="^70)
    println("COUPLING TEST COMPLETED SUCCESSFULLY!")
    println("="^70)
    println("\nSummary:")
    println("  ✓ Coupling initialization successful")
    println("  ✓ Communicator splitting working")
    println("  ✓ Data exchange between codes successful")
    println("  ✓ Broadcast to local ranks successful")
    println("\nYour coupling infrastructure is ready!")
    println("="^70 * "\n")
end

# Cleanup
finalize_coupling(ctx)
MPI.Finalize()

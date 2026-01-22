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
# Run long enough to cover Jexpresso's simulation
# Jexpresso: tinit=0.0, tend=3.0, dt=0.001 means ~3000-6000 steps depending on solver
# We'll run for a generous duration with frequent synchronization

max_duration = 300.0  # seconds - generous time for Jexpresso to complete
start_time = time()

if ctx.is_root
    println("[Alya mimic] Starting simulation loop for up to $(max_duration)s")
end

step = 0
while time() - start_time < max_duration
    global step
    step += 1

    if ctx.is_root && (step % 500 == 0)
        elapsed = time() - start_time
        println("[Alya mimic] Step $step (elapsed: $(round(elapsed, digits=1))s)")
    end

    # Synchronize with Jexpresso on coupling communicator
    synchronize_coupling(ctx)

    # Also sync on COMM_WORLD to catch any stray collective operations
    MPI.Barrier(MPI.COMM_WORLD)

    # Small sleep to avoid busy-waiting
    sleep(0.01)
end

if ctx.is_root
    println("\n[Alya mimic] Simulation loop ended")
end

# Cleanup
finalize_coupling(ctx)

# Don't call MPI.Finalize() - let Jexpresso finish first
# In real coupling, both codes would coordinate their finalization
if ctx.is_root
    println("[Alya mimic] Waiting for final synchronization...")
end
MPI.Barrier(MPI.COMM_WORLD)

MPI.Finalize()

println("[Alya mimic, rank $world_rank] Done")

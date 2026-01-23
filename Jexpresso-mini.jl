#!/usr/bin/env julia
#
# Jexpresso-mini.jl - Minimal coupling test (Jexpresso side)
#
# Usage with Alya:
#   mpirun -np 2 julia --project=. Jexpresso-mini.jl : -np 2 ./alya/Alya.x
#

using MPI

println("Jexpresso-mini starting...")

MPI.Init()

world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
world_size = MPI.Comm_size(MPI.COMM_WORLD)

println("[Jexpresso-mini rank $world_rank] World size: $world_size")

# Load coupling module
include("src/mpi/JexpressoCoupling.jl")
using .JexpressoCoupling

# Initialize coupling
println("[Jexpresso-mini rank $world_rank] Initializing coupling...")

ctx = initialize_coupling(
    MPI.COMM_WORLD,
    1,  # code_id for Jexpresso
    2;  # n_codes = 2 (Jexpresso + Alya)
    code_name="Jexpresso"
)

println("[Jexpresso-mini rank $world_rank] Coupling initialized!")
println("  Local rank: $(ctx.local_rank) / $(ctx.local_size)")
println("  World rank: $(ctx.world_rank)")
println("  Is root: $(ctx.is_root)")

# Test: Exchange application names
if ctx.is_root
    println("\n[Jexpresso-mini] Testing name exchange...")

    # Send "Jexpresso" to Alya
    my_name = "Jexpresso"
    partner_name = exchange_app_names!(ctx, my_name, 2; max_length=128)

    println("[Jexpresso-mini] Exchanged names successfully!")
    println("  My name: $my_name")
    println("  Partner name: $partner_name")
end

# Simulate some work with synchronization
for i in 1:5
    if ctx.is_root
        println("\n[Jexpresso-mini] Step $i/5")
    end

    synchronize_coupling(ctx)
    sleep(0.5)
end

println("\n[Jexpresso-mini rank $world_rank] Test complete!")

# Cleanup
finalize_coupling(ctx)
MPI.Finalize()

println("[Jexpresso-mini rank $world_rank] Done")

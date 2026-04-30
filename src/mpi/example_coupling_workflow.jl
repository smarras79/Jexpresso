"""
    example_coupling_workflow.jl

Practical example showing how to use coupling helpers in Jexpresso.

This demonstrates the typical workflow for coupling Jexpresso with an
external code that expects:
1. Application name exchange (MPI_Gather compatible)
2. Field data exchange at each timestep

You can adapt this pattern for your user_bc.jl or initialization files.
"""

using MPI
using Printf

# These would normally be included automatically via Jexpresso module
include("JexpressoCoupling.jl")
include("coupling_helpers.jl")
using .JexpressoCoupling

"""
Example 1: Send application name to external code during initialization
"""
function example_send_app_name(ctx::CouplingContext)
    println("\n" * "="^70)
    println("EXAMPLE 1: Sending Application Name")
    println("="^70)

    # This would be called during Jexpresso initialization
    partner_code_id = 2  # Assuming external code has code_id=2

    if is_coupling_root(ctx)
        # Send Jexpresso's name to partner
        # max_length=128 matches the external code's MPI_CHARACTER buffer size
        send_app_name_to_partner!(ctx, "Jexpresso", partner_code_id; max_length=128)

        println("[$(ctx.code_name)] Sent application name to partner code")
    end

    # Synchronize before continuing
    synchronize_coupling(ctx)
end

"""
Example 2: Two-way name exchange (both codes send and receive)
"""
function example_exchange_names(ctx::CouplingContext)
    println("\n" * "="^70)
    println("EXAMPLE 2: Exchanging Application Names")
    println("="^70)

    partner_code_id = 2

    partner_name = ""
    if is_coupling_root(ctx)
        # Exchange names with partner
        partner_name = exchange_app_names!(ctx, ctx.code_name, partner_code_id)
        println("[$(ctx.code_name)] Coupled with: '$partner_name'")
    end

    # Broadcast partner name to all local ranks
    if is_coupling_root(ctx)
        name_bytes = Vector{UInt8}(rpad(partner_name, 128))
    else
        name_bytes = Vector{UInt8}(undef, 128)
    end

    broadcast_from_root_to_all_local!(ctx, name_bytes)
    partner_name = String(name_bytes) |> rstrip

    println(@sprintf("[Rank %d, %s] Partner code name: %s",
                    ctx.local_rank, ctx.code_name, partner_name))
end

"""
Example 3: Send boundary field data to external code
This would typically be called at each timestep or coupling iteration
"""
function example_send_boundary_data(ctx::CouplingContext, timestep::Int)
    println("\n" * "="^70)
    println("EXAMPLE 3: Sending Boundary Field Data (timestep $timestep)")
    println("="^70)

    partner_code_id = 2
    n_boundary_points = 100

    # Step 1: Each rank computes its local boundary data
    # (In real code, this would extract actual boundary values)
    local_boundary_size = n_boundary_points ÷ ctx.local_size
    local_velocity = zeros(local_boundary_size)

    for i in 1:local_boundary_size
        # Simulate some boundary velocity data
        global_idx = ctx.local_rank * local_boundary_size + i
        local_velocity[i] = sin(2π * global_idx / n_boundary_points + timestep * 0.1)
    end

    # Step 2: Gather to coupling root
    if is_coupling_root(ctx)
        global_velocity = zeros(n_boundary_points)
    else
        global_velocity = zeros(0)  # Non-roots don't need this
    end

    gather_to_coupling_root!(ctx, local_velocity, global_velocity)

    # Step 3: Coupling root sends to partner
    if is_coupling_root(ctx)
        send_field_array!(ctx, global_velocity, partner_code_id; tag=200)
        println(@sprintf("[%s Root] Sent velocity field (size=%d) to code %d",
                        ctx.code_name, length(global_velocity), partner_code_id))
    end

    synchronize_coupling(ctx)
end

"""
Example 4: Receive boundary data from external code
"""
function example_receive_boundary_data(ctx::CouplingContext, timestep::Int)
    println("\n" * "="^70)
    println("EXAMPLE 4: Receiving Boundary Field Data (timestep $timestep)")
    println("="^70)

    partner_code_id = 2
    n_boundary_points = 100

    # Step 1: Coupling root receives from partner
    if is_coupling_root(ctx)
        global_pressure = zeros(n_boundary_points)
        recv_field_array!(ctx, global_pressure, partner_code_id; tag=201)
        println(@sprintf("[%s Root] Received pressure field from code %d",
                        ctx.code_name, partner_code_id))
        println(@sprintf("  Min: %.4f, Max: %.4f, Mean: %.4f",
                        minimum(global_pressure), maximum(global_pressure),
                        sum(global_pressure)/length(global_pressure)))
    else
        global_pressure = zeros(0)
    end

    # Step 2: Scatter to all local ranks
    local_pressure_size = n_boundary_points ÷ ctx.local_size
    local_pressure = zeros(local_pressure_size)

    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)

    # Step 3: Each rank applies the boundary condition
    println(@sprintf("[Rank %d, %s] Got local pressure data (size=%d)",
                    ctx.local_rank, ctx.code_name, length(local_pressure)))

    # In real code, you would now apply this as a boundary condition
    # apply_pressure_boundary_condition!(mesh, state, local_pressure)
end

"""
Example 5: Complete bidirectional exchange workflow
This is the most common pattern: send one field, receive another
"""
function example_bidirectional_exchange(ctx::CouplingContext, timestep::Int)
    println("\n" * "="^70)
    println("EXAMPLE 5: Bidirectional Exchange (timestep $timestep)")
    println("="^70)

    partner_code_id = 2
    n_interface = 100

    # Step 1: Prepare data to send
    local_size = n_interface ÷ ctx.local_size
    local_velocity = rand(local_size) .* 10.0  # Simulate velocity data

    # Gather to root
    if is_coupling_root(ctx)
        global_velocity = zeros(n_interface)
        global_pressure = zeros(n_interface)
    else
        global_velocity = zeros(0)
        global_pressure = zeros(0)
    end

    gather_to_coupling_root!(ctx, local_velocity, global_velocity)

    # Step 2: Root exchanges with partner
    if is_coupling_root(ctx)
        println(@sprintf("[%s Root] Starting exchange...", ctx.code_name))

        # Send velocity, receive pressure (different tags!)
        coupling_exchange_workflow!(ctx, global_velocity, global_pressure,
                                    partner_code_id;
                                    tag_send=300, tag_recv=301)

        println(@sprintf("[%s Root] Exchange complete", ctx.code_name))
        println(@sprintf("  Sent velocity: min=%.2f, max=%.2f",
                        minimum(global_velocity), maximum(global_velocity)))
        println(@sprintf("  Recv pressure: min=%.2f, max=%.2f",
                        minimum(global_pressure), maximum(global_pressure)))
    end

    # Step 3: Scatter received data to all ranks
    local_pressure = zeros(local_size)
    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)

    println(@sprintf("[Rank %d, %s] Applied pressure BC (mean=%.2f)",
                    ctx.local_rank, ctx.code_name, sum(local_pressure)/length(local_pressure)))
end

"""
Example 6: Alternative - Broadcast instead of Scatter
When all ranks need the full global array (not just their portion)
"""
function example_broadcast_workflow(ctx::CouplingContext)
    println("\n" * "="^70)
    println("EXAMPLE 6: Broadcast Workflow")
    println("="^70)

    partner_code_id = 2
    n_data = 50

    # Receive on root
    if is_coupling_root(ctx)
        global_data = zeros(n_data)
        recv_field_array!(ctx, global_data, partner_code_id; tag=400)
        println(@sprintf("[%s Root] Received data from partner", ctx.code_name))
    else
        global_data = zeros(n_data)
    end

    # Broadcast to ALL local ranks (everyone gets full array)
    broadcast_from_root_to_all_local!(ctx, global_data)

    # Now every rank has the full array
    println(@sprintf("[Rank %d, %s] Have full global data (size=%d)",
                    ctx.local_rank, ctx.code_name, length(global_data)))
end

"""
Main demonstration function
"""
function run_coupling_workflow_demo()
    MPI.Init()

    world_rank = MPI.Comm_rank(MPI.COMM_WORLD)
    world_size = MPI.Comm_size(MPI.COMM_WORLD)

    # Split ranks
    code_a_ranks = world_size ÷ 2
    code_id = (world_rank < code_a_ranks) ? 1 : 2
    code_name = (code_id == 1) ? "Jexpresso" : "ExternalCode"

    # Initialize coupling
    ctx = initialize_coupling(MPI.COMM_WORLD, code_id, 2; code_name=code_name)

    println("\n" * "━"^70)
    println("COUPLING WORKFLOW EXAMPLES")
    println("Code: $(ctx.code_name), Local Rank: $(ctx.local_rank)/$(ctx.local_size)")
    println("━"^70)

    # Run examples (comment out the ones you don't want to test)
    example_exchange_names(ctx)

    # Simulate a few timesteps with data exchange
    for timestep in 1:2
        example_bidirectional_exchange(ctx, timestep)
        synchronize_coupling(ctx)
    end

    # Cleanup
    finalize_coupling(ctx)
    MPI.Finalize()

    if world_rank == 0
        println("\n" * "━"^70)
        println("All examples completed successfully!")
        println("━"^70 * "\n")
    end
end

# Uncomment to run as standalone script
# run_coupling_workflow_demo()

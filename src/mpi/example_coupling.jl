# example_coupling.jl
# Example: Coupling Jexpresso with an atmospheric boundary condition code

using MPI
include("JexpressoCoupling.jl")
using .JexpressoCoupling

function run_jexpresso_coupled()
    MPI.Init()
    
    comm_world = MPI.COMM_WORLD
    world_rank = MPI.Comm_rank(comm_world)
    world_size = MPI.Comm_size(comm_world)
    
    # Configuration: first half = Jexpresso (code 1), second half = AtmosBCs (code 2)
    # In practice, this would come from a configuration file or command line
    jexpresso_ranks = world_size ÷ 2
    code_id = (world_rank < jexpresso_ranks) ? 1 : 2
    code_name = (code_id == 1) ? "Jexpresso" : "AtmosBCs"
    
    # Initialize coupling
    ctx = initialize_coupling(comm_world, code_id, 2; code_name=code_name)
    
    info = get_coupling_info(ctx)
    println("Rank $world_rank: $(info.code_name), local_rank=$(info.local_rank)/$(info.local_size)")
    
    # Simulation parameters
    n_interface_points = 100  # Points on the coupling interface
    n_timesteps = 10
    dt = 0.01
    
    # Allocate interface data arrays (on coupling roots only, but allocate for all)
    # In Jexpresso: we send velocity, receive pressure
    # In AtmosBCs: we send pressure, receive velocity  
    send_buffer = zeros(n_interface_points)
    recv_buffer = zeros(n_interface_points)
    
    # Register buffers for efficient repeated exchanges
    if is_coupling_root(ctx)
        partner_id = (code_id == 1) ? 2 : 1
        register_exchange_buffer!(ctx, :interface, n_interface_points, partner_id; tag=100)
    end
    
    synchronize_coupling(ctx)
    
    # Time stepping loop
    for step in 1:n_timesteps
        t = step * dt
        
        if code_id == 1  # Jexpresso
            # 1. Compute local solution (simplified)
            local_solution = sin.(range(0, π, n_interface_points) .+ t) .* (1.0 + 0.1 * info.local_rank)
            
            # 2. Gather interface data to root (if distributed)
            # For this example, assume interface data is already on root
            
            # 3. Exchange with partner code
            if is_coupling_root(ctx)
                # Prepare send data (e.g., velocity at interface)
                send_buffer .= local_solution
                
                # Async exchange - overlap with computation
                start_buffered_exchange!(ctx, :interface, send_buffer)
                
                # ... could do other work here ...
                
                # Complete exchange
                recv_data = finish_buffered_exchange!(ctx, :interface)
                recv_buffer .= recv_data
                
                println("[$code_name] Step $step: sent max=$(maximum(send_buffer)), " *
                        "recv max=$(maximum(recv_buffer))")
            end
            
            # 4. Broadcast received BCs to all local ranks
            broadcast_to_local!(ctx, recv_buffer)
            
            # 5. Apply boundary conditions and continue...
            
        else  # AtmosBCs code
            # 1. Receive velocity from Jexpresso, compute pressure response
            if is_coupling_root(ctx)
                # Prepare atmospheric pressure data
                send_buffer .= 101325.0 .+ 100.0 .* cos.(range(0, 2π, n_interface_points) .+ t)
                
                start_buffered_exchange!(ctx, :interface, send_buffer)
                recv_data = finish_buffered_exchange!(ctx, :interface)
                recv_buffer .= recv_data
                
                println("[$code_name] Step $step: sent max=$(maximum(send_buffer)), " *
                        "recv max=$(maximum(recv_buffer))")
            end
            
            broadcast_to_local!(ctx, recv_buffer)
        end
        
        # Synchronize codes before next timestep (optional, for tightly coupled schemes)
        synchronize_coupling(ctx)
    end
    
    # Cleanup
    finalize_coupling(ctx)
    
    MPI.Finalize()
end

# Run
run_jexpresso_coupled()

"""
Mutable struct to store timing statistics for a function that's called multiple times
"""
mutable struct MPIFunctionTimer
    total_time::Float64
    call_count::Int64
    min_time::Float64
    max_time::Float64
    comm::MPI.Comm
    skip_first_n::Int64
    skipped_count::Int64
    
    function MPIFunctionTimer(comm::MPI.Comm=MPI.COMM_WORLD; skip_first_n::Int=1)
        new(0.0, 0, Inf, 0.0, comm, skip_first_n, 0)
    end
end

"""
    reset_timer!(timer::MPIFunctionTimer)

Reset the timer to initial state
"""
function reset_timer!(timer::MPIFunctionTimer)
    timer.total_time = 0.0
    timer.call_count = 0
    timer.min_time = Inf
    timer.max_time = 0.0
    timer.skipped_count = 0
end

"""
    time_function!(timer::MPIFunctionTimer, f::Function, args...)

Time a single function call and accumulate statistics.
Skips the first N calls (configurable) to avoid JIT compilation overhead.
"""
function time_function!(timer::MPIFunctionTimer, f::Function, args...)
    start = time()
    result = f(args...)
    elapsed = time() - start
    
    # Skip the first N calls to avoid compilation overhead
    if timer.skipped_count < timer.skip_first_n
        timer.skipped_count += 1
        return result
    end
    
    timer.total_time += elapsed
    timer.call_count += 1
    timer.min_time = min(timer.min_time, elapsed)
    timer.max_time = max(timer.max_time, elapsed)
    
    return result
end

"""
    report_timer(timer::MPIFunctionTimer; name::String="Function")

Gather timing statistics from all MPI ranks and print report on rank 0
"""
function report_timer(timer::MPIFunctionTimer; name::String="Function")
    rank = MPI.Comm_rank(timer.comm)
    nprocs = MPI.Comm_size(timer.comm)
    
    # Gather statistics from all ranks
    all_total_times = MPI.Gather(timer.total_time, 0, timer.comm)
    all_call_counts = MPI.Gather(timer.call_count, 0, timer.comm)
    all_min_times = MPI.Gather(timer.min_time, 0, timer.comm)
    all_max_times = MPI.Gather(timer.max_time, 0, timer.comm)
    
    if rank == 0
        println("\n" * "="^70)
        println("MPI Timing Report: $name")
        println("="^70)
        println("Number of MPI processes: $nprocs")
        println("\nPer-process statistics:")
        println("-"^70)
        println("Rank | Calls      | Total Time (s) | Avg (ms)  | Min (ms)  | Max (ms)")
        println("-"^70)
        
        for i in 1:nprocs
            rank_id = i - 1
            calls = all_call_counts[i]
            total = all_total_times[i]
            avg = calls > 0 ? (total / calls) * 1000 : 0.0
            min_t = all_min_times[i] == Inf ? 0.0 : all_min_times[i] * 1000
            max_t = all_max_times[i] * 1000
            
            @printf("%4d | %10d | %14.6f | %9.6f | %9.6f | %9.6f\n", 
                    rank_id, calls, total, avg, min_t, max_t)
        end
        
        println("-"^70)
        total_calls = sum(all_call_counts)
        total_time_all = sum(all_total_times)
        
        println("\nGlobal statistics:")
        println("  Total calls across all ranks: $total_calls")
        println("  Total time across all ranks:  $(round(total_time_all, digits=6))s")
        println("  Average time per call:         $(round(total_time_all/total_calls*1000, digits=6))ms")
        println("  Min time (single call):        $(round(minimum(all_min_times)*1000, digits=6))ms")
        println("  Max time (single call):        $(round(maximum(all_max_times)*1000, digits=6))ms")
        println("="^70 * "\n")
        
        return (
            total_calls = total_calls,
            total_time = total_time_all,
            all_call_counts = all_call_counts,
            all_total_times = all_total_times,
            all_min_times = all_min_times,
            all_max_times = all_max_times
        )
    else
        return nothing
    end
end

# ============================================================================
# Multiple function timing
# ============================================================================

"""
    create_timer_dict(function_names::Vector{String}, comm::MPI.Comm=MPI.COMM_WORLD; skip_first_n::Int=1)

Create a dictionary of timers for multiple functions
"""
function create_timer_dict(function_names::Vector{String}, comm::MPI.Comm=MPI.COMM_WORLD; skip_first_n::Int=1)
    return Dict(name => MPIFunctionTimer(comm, skip_first_n=skip_first_n) for name in function_names)
end

"""
    report_all_timers(timers::Dict{String, MPIFunctionTimer})

Report timing for all functions in the dictionary
"""
function report_all_timers(timers::Dict{String, MPIFunctionTimer})
    for name in sort(collect(keys(timers)))
        report_timer(timers[name], name=name)
    end
end

# ============================================================================
# Example usage with ODEProblem
# ============================================================================

"""
Example: Lorenz system with timed form_rhs function
"""
function example_lorenz_with_timing()
    
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    # Create timer for form_rhs
    rhs_timer = MPIFunctionTimer(comm)
    
    # Lorenz parameters
    σ = 10.0
    ρ = 28.0
    β = 8.0/3.0
    
    # Define the actual physics computation
    function form_rhs(u, p, t)
        σ, ρ, β = p
        du1 = σ * (u[2] - u[1])
        du2 = u[1] * (ρ - u[3]) - u[2]
        du3 = u[1] * u[2] - β * u[3]
        return [du1, du2, du3]
    end
    
    # Wrapper function that times form_rhs
    function lorenz_timed!(du, u, p, t)
        result = time_function!(rhs_timer, form_rhs, u, p, t)
        du .= result
        return nothing
    end
    
    # Initial conditions (different for each rank for demonstration)
    u0 = [1.0 + rank, 0.0, 0.0]
    tspan = (0.0, 10.0)
    p = (σ, ρ, β)
    
    # Create and solve ODE problem
    prob = ODEProblem(lorenz_timed!, u0, tspan, p)
    
    if rank == 0
        println("Solving Lorenz system with timing...")
    end
    
    MPI.Barrier(comm)
    solve_start = time()
    
    sol = solve(prob, Tsit5(), saveat=0.1)
    
    MPI.Barrier(comm)
    solve_time = time() - solve_start
    
    # Report timing statistics
    report_timer(rhs_timer, name="form_rhs")
    
    # Also report total solve time
    all_solve_times = MPI.Gather(solve_time, 0, comm)
    if rank == 0
        println("Total ODE solve time per rank:")
        for (i, t) in enumerate(all_solve_times)
            println("  Rank $(i-1): $(round(t, digits=6))s")
        end
        println()
    end
    
    MPI.Finalize()
end




# ============================================================================
# Example with your specific structure
# ============================================================================

"""
Example showing how to integrate with your existing code structure with multiple timed functions
"""
function your_solver_with_timing(inputs, prob_params, comm=MPI.COMM_WORLD)
    # Create multiple timers - one for each function
    timers = create_timer_dict(["form_rhs_1", "form_rhs_2", "form_rhs_3"], comm)
    
    # OR create them individually:
    # rhs1_timer = MPIFunctionTimer(comm)
    # rhs2_timer = MPIFunctionTimer(comm)
    # rhs3_timer = MPIFunctionTimer(comm)
    
    # Your original functions
    function form_rhs_1(u, p, t)
        # First RHS computation
        return du1
    end
    
    function form_rhs_2(u, p, t)
        # Second RHS computation
        return du2
    end
    
    function form_rhs_3(u, p, t)
        # Third RHS computation
        return du3
    end
    
    # Create timed wrapper that calls all functions
    function your_ode_function!(du, u, p, t)
        # Time each function call
        result1 = time_function!(timers["form_rhs_1"], form_rhs_1, u, p, t)
        result2 = time_function!(timers["form_rhs_2"], form_rhs_2, u, p, t)
        result3 = time_function!(timers["form_rhs_3"], form_rhs_3, u, p, t)
        
        # Combine results (adjust based on your needs)
        du .= result1 .+ result2 .+ result3
        return nothing
    end
    
    # Create callbacks (your existing callbacks)
    cb = DiscreteCallback(condition, affect!)
    cb_restart = DiscreteCallback(condition_restart, affect_restart!)
    dosetimes = inputs[:dosetimes]
    
    # Create and solve problem
    u0 = prob_params[:u0]
    tspan = (inputs[:tinit], inputs[:tend])
    prob = ODEProblem(your_ode_function!, u0, tspan, prob_params[:p])
    
    solution = solve(prob,
                     inputs[:ode_solver], 
                     dt=Float32(inputs[:Δt]),
                     callback = CallbackSet(cb, cb_restart), 
                     tstops = dosetimes,
                     save_everystep = false,
                     adaptive=inputs[:ode_adaptive_solver],
                     saveat = range(inputs[:tinit],
                                   inputs[:tend],
                                   length=inputs[:ndiagnostics_outputs]))
    
    # Report timing for all functions
    report_all_timers(timers)
    
    # OR report individually:
    # report_timer(timers["form_rhs_1"], name="form_rhs_1")
    # report_timer(timers["form_rhs_2"], name="form_rhs_2")
    # report_timer(timers["form_rhs_3"], name="form_rhs_3")
    
    return solution
end

# Uncomment to run example:
# example_lorenz_with_timing()
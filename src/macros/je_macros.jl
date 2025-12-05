
"""
    @missinginput "Error message"
Macro used to raise an error, when something is not implemented.
"""
macro missinginput(message=" # Missing input variable in user_input")
    quote
        error($(esc(message)))
    end
end


macro mystop(message=" MY STOP HERE")
    quote
        error($(esc(message)))
    end
end

"""
    @timers func(args...)

Macro to time a function call. Expects a variable named `timers` in scope.

# Example
```julia
timers = create_timer_dict([...])

@timers DSS_nc_gather_rhs!(params.RHS, SD, QT, ...)
```
"""
macro timers(expr)
    if expr.head == :call
        func_name = string(expr.args[1])
    else
        error("@time_function expects a function call")
    end
    
    return quote
        local timers = $(esc(:timers))
        local func_name = $func_name
        
        # Create timer if it doesn't exist
        if !haskey(timers, func_name)
            timers[func_name] = MPIFunctionTimer(MPI.COMM_WORLD; skip_first_n=10)
        end
        
        local timer_obj = timers[func_name]
        local start = time()
        local result = $(esc(expr))
        local elapsed = time() - start
        
        if timer_obj.skipped_count < timer_obj.skip_first_n
            timer_obj.skipped_count += 1
        else
            timer_obj.total_time += elapsed
            timer_obj.call_count += 1
            timer_obj.min_time = min(timer_obj.min_time, elapsed)
            timer_obj.max_time = max(timer_obj.max_time, elapsed)
        end
        
        result
    end
end

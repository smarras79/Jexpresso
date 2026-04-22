
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
    @outputrootonly expr

Macro that runs `expr` normally on MPI rank 0, and silences stdout on all other ranks.
Requires a variable named `rank` (integer MPI rank) to be in scope at the call site.

# Example
```julia
model = @outputrootonly GmshDiscreteModel(parts, filename, renumber=true)
```
"""
macro outputrootonly(expr)
    quote
        if $(esc(:rank)) != 0
            redirect_stdout(open("/dev/null", "w")) do
                $(esc(expr))
            end
        else
            $(esc(expr))
        end
    end
end

"""
    @mpi_time expr

Drop-in replacement for `@time` in MPI code. Prints max wall time, allocations,
and bytes only from rank 0. Requires `comm` and `rank` in scope.

# Example
```julia
@mpi_time prob, partitioned_model = amr_strategy!(...)
```
"""
macro mpi_time(expr)
    label = string(expr)
    quote
        local _stats = @timed $(esc(expr))
        local _t_max = MPI.Allreduce(_stats.time, MPI.MAX, $(esc(:comm)))
        if $(esc(:rank)) == 0
            local _allocs = _stats.gcstats.poolalloc + _stats.gcstats.malloc
            @printf("  %s\n  %.6f seconds (max across ranks, %d allocs, %.3f MiB)\n",
                    $label, _t_max, _allocs, _stats.bytes / 1024^2)
        end
        _stats.value
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

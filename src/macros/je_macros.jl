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
            # Use a fresh per-call sink instead of the global Base.devnull:
            # 88438d9 swapped open("/dev/null","w") -> devnull for
            # portability, but on macOS arm64 with Open MPI 5 + Gridap +
            # GridapGmsh the shared devnull stream interacts badly with
            # Gmsh's C-level stdout buffering inside the parallel
            # GmshDiscreteModel collective, triggering a Bus error 10
            # in _platform_memmove right after "Done reading *.msh". Going
            # back to a fresh IOStream sink avoids the shared-handle issue.
            local _sink = open("/dev/null", "w")
            try
                redirect_stdout(_sink) do
                    $(esc(expr))
                end
            finally
                close(_sink)
            end
        else
            $(esc(expr))
        end
    end
end


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
            # get_mpi_comm() returns Julia's local sub-comm in coupled
            # mode (where rank 0 of COMM_WORLD is Alya, not Julia) and
            # COMM_WORLD otherwise. Using it here keeps the timer's
            # internal MPI.Gather (timing.jl:101-104) collectives
            # restricted to Julia ranks.
            timers[func_name] = MPIFunctionTimer(get_mpi_comm(); skip_first_n=10)
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

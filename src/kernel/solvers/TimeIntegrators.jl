using BenchmarkTools

function time_loop!(inputs, params, u)

    println(" # Solving ODE  ................................ ")
    
    prob = ODEProblem(rhs!,
                      u,
                      params.tspan,
                      params);
    
    #------------------------------------------------------------------------
    # Callback to plot on the run
    #------------------------------------------------------------------------
    dosetimes = inputs[:diagnostics_at_times]
    idx_ref   = Ref{Int}(0)
    c         = Float64(0.0)
    function condition(u, t, integrator)
        idx  = findfirst(x -> x == t, dosetimes)
        
        if idx !== nothing
            idx_ref[] = idx
            return true
        else
            return false
        end
    end
    function affect!(integrator)
        idx  = idx_ref[]

        println(" #  t=", integrator.t)
        computeCFL(params.mesh.npoin, inputs[:Δt], params.mesh.Δeffective_s, integrator, params.SD)
        
        write_output(params.SD, integrator.u, integrator.t, idx,
                     params.mesh,
                     inputs[:output_dir], inputs,
                     params.qp.qvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe, case="rtb")

    
    end
    cb = DiscreteCallback(condition, affect!)
    
    #------------------------------------------------------------------------
    # END Callback to plot on the run
    #------------------------------------------------------------------------
    
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                           callback = cb, tstops = dosetimes,
                           #saveat = dosetimes,
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    
    println(" # Solving ODE  ................................ DONE")
    
    return solution
end

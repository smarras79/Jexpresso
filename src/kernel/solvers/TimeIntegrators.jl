function time_loop!(inputs, params, u)

    println(" # Solving ODE  ................................ ")
    
    prob = ODEProblem(rhs!,
                      u,
                      params.tspan,
                      params);
    
    #------------------------------------------------------------------------
    # Runtime callbacks
    #------------------------------------------------------------------------
    dosetimes = inputs[:diagnostics_at_times]
    idx_ref   = Ref{Int}(0)
    c         = Float64(0.0)
    rad_time  = inputs[:radiation_time_step]
    lnew_mesh = true   
    function two_stream_condition(u, t, integrator)
        if (rem(t,rad_time) < 1e-3)
            return true
        else
            return false
        end
    end

    function do_radiation!(integrator)
        println(" doing two stream radiation heat flux calculations at t=", integrator.t)
        @info "doing rad test"
        compute_radiative_fluxes!(lnew_mesh, params.mesh, params.uaux, params.qp.qe, params.mp, params.phys_grid, params.inputs[:backend], params.SOL_VARS_TYPE)
    end
    
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
        
        #CFL
        #computeCFL(params.mesh.npoin, params.neqs, inputs[:Δt], params.mesh.Δeffective_s, integrator, params.SD; visc=inputs[:μ])
        
        #Write results to file
        write_output(params.SD, integrator.u, integrator.t, idx,
                     params.mesh, params.mp,
                     inputs[:output_dir], inputs,
                     params.qp.qvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe, case=inputs[:case])

    
    end
    cb = DiscreteCallback(condition, affect!)
    cb_rad = DiscreteCallback(two_stream_condition, do_radiation!)
    #------------------------------------------------------------------------
    # END runtime callbacks
    #------------------------------------------------------------------------

    
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                           callback = cb, tstops = dosetimes,
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    
    println(" # Solving ODE  ................................ DONE")
    
    return solution
end

function time_loop!(inputs, params, u)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
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
    ret_dosetime_ref  = Ref{Bool}(false)

    
    # #------------------------------------------------------------------------
    # # AMR config
    # #------------------------------------------------------------------------
    # amr_freq        = inputs[:amr_freq]
    # Δt_amr          = amr_freq * inputs[:Δt]
    # ret_amrtime_ref = Ref{Bool}(false)

    function condition(u, t, integrator)
        idx  = findfirst(x -> x == t, dosetimes)
        if idx !== nothing
            idx_ref[] = idx
            ret_dosetime_ref[] = true
        else
            ret_dosetime_ref[] = false
        end

        tol             = 1e-6
        # ret_amrtime_ref = abs(mod(t, Δt_amr)) < tol
        # return (ret_dosetime_ref[] || ret_amrtime_ref[])
        return ret_dosetime_ref[]
    end
    function affect!(integrator)
        idx          = idx_ref[]
        ret_dosetime = ret_dosetime_ref[]
        # ret_amrtime  = ret_amrtime_ref[]
        if ret_dosetime == true
            println_rank(" #  t=", integrator.t; msg_rank = rank)

            #CFL
        #   @time  computeCFL(params.mesh.npoin, inputs[:Δt], params.mesh.Δeffective_s, integrator, params.SD; visc=inputs[:μ])

            #Write results to file
            write_output(integrator.p.SD, integrator.u, params.uaux, integrator.t, idx,
                        integrator.p.mesh, integrator.p.mp,
                        inputs[:output_dir], inputs,
                        integrator.p.qp.qvars,
                        inputs[:outformat];
                        nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe, case="rtb")

        end

        # if ret_amrtime == true
        #     println_rank(" # mesh adapt at t=", integrator.t; msg_rank = rank)
        #     # nothing
        #     @info "u1",  size(integrator.u,1), integrator.p[38].npoin
        #     amr_strategy!(inputs,integrator, integrator.p, integrator.u)
        #     @info "u2",  size(integrator.u,1)
        # end
    end
    cb_rad = DiscreteCallback(two_stream_condition, do_radiation!)
    cb = DiscreteCallback(condition, affect!)    
    cb_amr = DiscreteCallback(condition, affect!)    
    #------------------------------------------------------------------------
    # END runtime callbacks
    #------------------------------------------------------------------------

    
    solution = solve(prob,
                           inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                           callback = cb, tstops = dosetimes,
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    # @info fieldnames(typeof(prob))
    if inputs[:ladapt] == true
        while solution.t[end] < inputs[:tend]
            prob = amr_strategy!(inputs, prob.p, solution.u[end][:], solution.t[end])
            
            solution = solve(prob,
                                inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                                callback = cb_amr, tstops = dosetimes,
                                save_everystep = false,
                                adaptive=inputs[:ode_adaptive_solver],
                                saveat = []);
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return solution
end

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
        compute_radiative_fluxes!(lnew_mesh,
                                  params.mesh,
                                  params.uaux,
                                  params.qp.qe,
                                  params.mp,
                                  params.phys_grid,
                                  params.inputs[:backend],
                                  params.SOL_VARS_TYPE)
    end
    ret_dosetime_ref  = Ref{Bool}(false)


    
    #------------------------------------------------------------------------
    # AMR config
    #------------------------------------------------------------------------
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
        if ret_dosetime == true
            println_rank(" #  t=", integrator.t; msg_rank = rank)

            #CFL
            #@time  computeCFL(params.mesh.npoin, inputs[:Δt], params.mesh.Δeffective_s, integrator, params.SD; visc=inputs[:μ])
            
            write_output(integrator.p.SD, integrator.u, params.uaux, integrator.t, idx,
                         integrator.p.mesh, integrator.p.mp, params.μdsgs,
                         integrator.p.connijk_original, integrator.p.poin_in_bdy_face_original,
                         integrator.p.x_original, integrator.p.y_original, integrator.p.z_original,
                         inputs[:output_dir], inputs,
                         integrator.p.qp.qvars,
                         integrator.p.qp.qoutvars,                         
                         inputs[:outformat];
                         nvar=integrator.p.qp.neqs,
                         qexact=integrator.p.qp.qe)
            
        end
    end
    cb_rad = DiscreteCallback(two_stream_condition, do_radiation!)
    cb     = DiscreteCallback(condition, affect!)    
    cb_amr = DiscreteCallback(condition, affect!)
    CallbackSet(cb)#,cb_rad)
    #------------------------------------------------------------------------
    # END runtime callbacks
    #------------------------------------------------------------------------

    #
    # Write initial conditions:
    #
    if rank == 0 println(" # Write initial condition to ",  typeof(inputs[:outformat]), " .........") end
    write_output(params.SD, u, params.uaux, inputs[:tinit], 0,
                 params.mesh, params.mp, params.μdsgs,
                 params.connijk_original, params.poin_in_bdy_face_original,
                 params.x_original, params.y_original, params.z_original,
                 inputs[:output_dir], inputs,
                 params.qp.qvars, params.qp.qoutvars,
                 inputs[:outformat];
                 nvar=params.qp.neqs, qexact=params.qp.qe)
    if rank == 0  println(" # Write initial condition to ",  typeof(inputs[:outformat]), " ......... END") end
    
    #
    # Simulation
    #
    sol = solve(prob,
                inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                callback = CallbackSet(cb), tstops = dosetimes,
                save_everystep = false,
                adaptive=inputs[:ode_adaptive_solver],
                saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    
    
    #------------------------------------------------------------------------
    # Initialize a buffer to store the previous two sol values
    #------------------------------------------------------------------------
   #= previous_two_u = Vector{typeof(sol.u[1])}(undef, 2)
    fill!(previous_two_u, sol.u[1]) # Initialize with the first sol

    println("Time | Current u | Previous u (t-dt) | Previous u (t-2dt)")
    println("-----|-----------|-------------------|---------------------")

    for (i, t) in enumerate(sol.t)
        current_u = sol.u[i]
        previous_u_t_minus_dt = previous_two_u[2]
        previous_u_t_minus_2dt = previous_two_u[1]

        println("$(round(t, digits=2))  | $(round(current_u[1], digits=8))   | $(round(previous_u_t_minus_dt[1], digits=8))        | $(round(previous_u_t_minus_2dt[1], digits=8))")
        
        
        # Update the buffer for the next iteration
        previous_two_u[1] = previous_two_u[2]
        previous_two_u[2] = current_u
    end
    =#
    
    if inputs[:ladapt] == true
        while sol.t[end] < inputs[:tend]
            prob = amr_strategy!(inputs, prob.p, sol.u[end][:], sol.t[end])
            
            sol = solve(prob,
                        inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                        callback = cb_amr, tstops = dosetimes,
                        save_everystep = false,
                        adaptive=inputs[:ode_adaptive_solver],
                        saveat = []);
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return sol
end

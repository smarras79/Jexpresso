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
    dosetimes    = inputs[:diagnostics_at_times]
    idx_ref      = Ref{Int}(0)
    c            = Float64(0.0)
    restart_time = inputs[:restart_time]
    rad_time     = inputs[:radiation_time_step]
    lnew_mesh    = true   
    function two_stream_condition(u, t, integrator)
        if (rem(t,rad_time) < 1e-3)
            return true
        else
            return false
        end
    end

    function do_radiation!(integrator)
        println(" doing two stream radiation heat flux calculations at t=", integrator.t)
        #@info "doing rad test"
        compute_radiative_fluxes!(lnew_mesh, params.mesh, params.uaux, params.qp.qe, params.mp, params.phys_grid, params.inputs[:backend], params.SOL_VARS_TYPE)
    end

    function restart_condition(u, t, integrator)
        if restart_time ≠ 0.0 && (rem(t,restart_time) < 1e-3)
            return true
        else
            return false
        end
    end
    function do_restart!(integrator)
        idx         = idx_ref[]
        res_fortmat = HDF5()
        println_rank(" #  writing restart ........................", round(integrator.t,digits=2); msg_rank = rank)
        tmp_restart_path = joinpath(inputs[:output_dir],"tmp_restart")
        if (rank == 0)
            if !isdir(tmp_restart_path)
                mkpath(tmp_restart_path)
            end
        end
        MPI.Barrier(comm)
        write_output(integrator.p.SD, integrator.u, params.uaux, integrator.t, idx,
                        integrator.p.mesh, integrator.p.mp,
                        integrator.p.connijk_original, integrator.p.poin_in_bdy_face_original,
                        integrator.p.x_original, integrator.p.y_original, integrator.p.z_original,
                        tmp_restart_path, inputs,
                        integrator.p.qp.qvars,
                        integrator.p.qp.qoutvars,
                        res_fortmat;
                        nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe)
        MPI.Barrier(comm)
        if rank == 0
            cp(tmp_restart_path, inputs[:restart_output_file_path]; force=true)
            rm(tmp_restart_path; recursive=true, force=true)
        end

        println_rank(" #  writing restart ........................ DONE"; msg_rank = rank)
    end
    # #------------------------------------------------------------------------
    # #  config
    # #------------------------------------------------------------------------
    ret_dosetime_ref  = Ref{Bool}(false)
    function condition(u, t, integrator)
        idx  = findfirst(x -> x == t, dosetimes)
        if idx !== nothing
            idx_ref[] = idx
            ret_dosetime_ref[] = true
        else
            ret_dosetime_ref[] = false
        end

        tol = 1e-6
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
            computeCFL(params.mesh.npoin, integrator.p.qp.neqs,
                       integrator.p.mp, integrator.p.uaux[:,end], inputs[:Δt],
                       params.mesh.Δeffective_s,
                       integrator,
                       params.SD; visc=inputs[:μ])
            
            write_output(integrator.p.SD, integrator.u, params.uaux, integrator.t, idx,
                         integrator.p.mesh, integrator.p.mp,
                         integrator.p.connijk_original, integrator.p.poin_in_bdy_face_original,
                         integrator.p.x_original, integrator.p.y_original, integrator.p.z_original,
                         inputs[:output_dir], inputs,
                         integrator.p.qp.qvars,
                         integrator.p.qp.qoutvars,
                         inputs[:outformat];
                         nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe)
        end
    end
    cb_rad     = DiscreteCallback(two_stream_condition, do_radiation!)
    cb         = DiscreteCallback(condition, affect!)    
    cb_amr     = DiscreteCallback(condition, affect!)
    cb_restart = DiscreteCallback(restart_condition, do_restart!)
    CallbackSet(cb)#,cb_rad)
    #------------------------------------------------------------------------
    # END runtime callbacks
    #------------------------------------------------------------------------

    #
    # Write initial conditions:
    #
    idx  = (inputs[:tinit] == 0.0) ? 0 : findfirst(x -> x == inputs[:tinit], dosetimes)
    if idx ≠ nothing
        if rank == 0 println(" # Write initial condition to ",  typeof(inputs[:outformat]), " .........") end
        write_output(params.SD, u, params.uaux, inputs[:tinit], idx,
                    params.mesh, params.mp,
                    params.connijk_original, params.poin_in_bdy_face_original,
                    params.x_original, params.y_original, params.z_original,
                    inputs[:output_dir], inputs,
                    params.qp.qvars, params.qp.qoutvars,
                    inputs[:outformat];
                    nvar=params.qp.neqs, qexact=params.qp.qe)
        if rank == 0  println(" # Write initial condition to ",  typeof(inputs[:outformat]), " ......... END") end
    end
    
    #
    # Simulation
    #
    limex = false
    if limex
        ntime_steps = floor(Int32, inputs[:tend]/inputs[:Δt])
        
        # Basic usage
        u_final = imex_integration_simple_2d!(u, params, params.mesh.connijk, params.qp.qe, params.mesh.coords, 
                                           inputs[:Δt], ntime_steps, inputs[:lsource])
        
        # Or step-by-step
        for n = 1:ntime_steps
            imex_time_step_simple_2d!(u, params, params.mesh.connijk,  params.qp.qe,  params.mesh.coords, inputs[:Δt], inputs[:lsource])
        end
        println(" IMEX RAN IT SEEMS. IS IT CORRECT? WHO KNOWS?")
        @mystop("TimeIntegrators.jl WIP on IMEX by Simone (obsolete)")
    else
        solution = solve(prob,
                         inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                         #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                         callback = CallbackSet(cb, cb_restart), tstops = dosetimes,
                         save_everystep = false,
                         adaptive=inputs[:ode_adaptive_solver],
                         saveat = range(inputs[:tinit],
                                        inputs[:tend],
                                        length=inputs[:ndiagnostics_outputs]));
    end
    
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

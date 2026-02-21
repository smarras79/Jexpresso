function time_loop!(inputs, params, u, args...)
    
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    partitioned_model = length(args) >= 1 ? args[1] : nothing
    is_coupled        = length(args) >= 2 ? args[2] : nothing
    coupling          = length(args) >= 3 ? args[3] : nothing
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    prob = ODEProblem(rhs!, u, params.tspan, params)

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
        compute_radiative_fluxes!(lnew_mesh, params.mesh, params.uaux, params.qp.qe, params.mp,
                                  params.phys_grid, params.inputs[:backend], params.SOL_VARS_TYPE)
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

    #------------------------------------------------------------------------
    # Diagnostics callback
    #------------------------------------------------------------------------
    ret_dosetime_ref  = Ref{Bool}(false)
    function condition(u, t, integrator)
        idx  = findfirst(x -> x == t, dosetimes)
        if idx !== nothing
            idx_ref[] = idx
            ret_dosetime_ref[] = true
        else
            ret_dosetime_ref[] = false
        end
        return ret_dosetime_ref[]
    end

    function affect!(integrator)
        idx          = idx_ref[]
        ret_dosetime = ret_dosetime_ref[]
        if ret_dosetime == true
            println_rank(" #  t=", integrator.t; msg_rank = rank)

            if inputs[:ladapt] == false
                computeCFL(integrator.p.mesh.npoin, integrator.p.qp.neqs,
                           integrator.p.mp, integrator.p.uaux[:,end], inputs[:Δt],
                           integrator.p.mesh.Δeffective_s,
                           integrator,
                           integrator.p.SD; visc=inputs[:μ])
            end
            
            write_output(integrator.p.SD, integrator.u, integrator.p.uaux, integrator.t, idx,
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

    #------------------------------------------------------------------------
    # Coupling callback: RECEIVE METADATA + EXCHANGE SOLUTION at every step
    #------------------------------------------------------------------------
    #------------------------------------------------------------------------
    # Coupling callback: EXCHANGE at every step (if enabled)
    #------------------------------------------------------------------------
    coupling_enabled = (is_coupled !== false)

    if coupling_enabled

        # Pull coupling object prepared in setup_coupling_and_mesh
        cpg = params.coupling
        @assert cpg !== nothing "params.coupling must be set during setup."
        
        # Condition: couple at every step after initial time
        t0   = params.tspan[1]
        tol0 = get(inputs, :couple_time_tol, 1e-12)
        
        @inline function coupling_condition(u_state, t, integrator) t > t0 + tol0 end
        
        basis = hasproperty(params, :basis) ? params.basis : params.SD.basis
        ξ     = params.ξ
        neqs  = params.neqs
        mesh  = params.mesh
        
        function do_coupling_exchange!(integrator)
            je_perform_coupling_exchange(integrator.u, integrator.p.uaux, integrator.t,
                                         cpg, mesh, basis, inputs, ξ, neqs)
        end
        
        cb_coupling = DiscreteCallback(coupling_condition, do_coupling_exchange!)
    end
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
    # Build callbacks
    #
    callbacks = coupling_enabled ? CallbackSet(cb, cb_restart, cb_coupling) : CallbackSet(cb, cb_restart)
    #callbacks = coupling_enabled ? CallbackSet(cb, cb_restart) : CallbackSet(cb, cb_restart)
    tstops_all = dosetimes
    
    #------------------------------------------------------------------------
    # TIME INTEGRATION
    #------------------------------------------------------------------------
    solution = solve(prob,
                     inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                     callback = callbacks, tstops = tstops_all,
                     save_everystep = false,
                     adaptive=false,
                     saveat = range(inputs[:tinit],
                                    inputs[:tend],
                                    length=inputs[:ndiagnostics_outputs]))

    MPI.Barrier(comm)
    report_all_timers(params.timers)
    MPI.Barrier(comm)
        
    if inputs[:lamr] == true
        while solution.t[end] < inputs[:tend]
            @time prob, partitioned_model = amr_strategy!(inputs, prob.p, solution.u[end][:], solution.t[end], partitioned_model)

            @time solution = solve(prob,
                                   inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                                   #callback = coupling_enabled ? CallbackSet(cb_amr, cb_restart)
                                   callback = coupling_enabled ? CallbackSet(cb_amr, cb_restart, cb_coupling)
                                   : CallbackSet(cb_amr, cb_restart),
                                   tstops = dosetimes,
                                   save_everystep = false,
                                   adaptive=false,
                                   saveat = [])
            MPI.Barrier(comm)
            report_all_timers(prob.p.timers)
            MPI.Barrier(comm)
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)

    MPI.Finalize
    
    return solution
end

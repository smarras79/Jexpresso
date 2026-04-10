function time_loop!(inputs, params, u, args...)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    partitioned_model = args[1]
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    prob = ODEProblem(rhs!,
                      u,
                      params.tspan,
                      params);
    
    #------------------------------------------------------------------------
    # Runtime callbacks
    #------------------------------------------------------------------------
    dosetimes    = inputs[:diagnostics_at_times]
    les_stat_t         = inputs[:statistics_time]
    tstops_all   = sort(unique(vcat(collect(Float64, dosetimes), collect(Float64, les_stat_t))))
    idx_ref      = Ref{Int}(0)
    c            = Float64(0.0)
    restart_time = inputs[:restart_time]
    stats_online_start = get(inputs, :statistics_online_start, Inf)
    rad_time           = inputs[:radiation_time_step]
    lnew_mesh    = true   
    lwrite_time  = (inputs[:outformat] == VTK()) && (rank == 0)

    if (lwrite_time == true)
        pvd_path = joinpath(inputs[:output_dir], "simulation.pvd")
        if get(inputs, :lrestart_vtk, false) && isfile(pvd_path)
            # VTK restart: preserve existing simulation.pvd; continue appending
        else
            init_pvd_file(pvd_path)
        end
    end

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


    # LES statistics callback:
    function les_stat_condition(u, t, integrator)
        # return les_stat_t ≠ 0.0 && rem(t, les_stat_t) < 1e-3
        idx  = findfirst(x -> x == t, les_stat_t)
        if idx !== nothing
            return true
        else
            return false
        end
    end
    function do_les_statistics!(integrator)
        println_rank(" # LES Statistics at t=", integrator.t; msg_rank = rank)
        les_statistics(integrator.u, integrator.p, integrator.t)
    end

    # Online statistics accumulation callback (Approach 2): fires every interval, no MPI
    stats_online_interval = Float64(get(inputs, :statistics_online_interval, inputs[:Δt]))
    online_last_t         = Ref{Float64}(-Inf)
    function les_online_condition(u, t, integrator)
        return !isnothing(integrator.p.les_stat_cache) &&
               t >= stats_online_start &&
               t - online_last_t[] >= stats_online_interval - eps(t)
    end
    function do_les_online!(integrator)
        online_last_t[] = integrator.t
        println_rank(" # LES online accumulation at t=", integrator.t; msg_rank = rank)
        les_accumulate_online!(integrator.u, integrator.p)
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
            # if inputs[:ladapt] == false

                ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl), MPI.MAX, comm)
                dt         = Float32(inputs[:Δt]/(2.0^(ad_lvl_max)))
                computeCFL(integrator.p.mesh.npoin, integrator.p.qp.neqs,
                        integrator.p.mp, integrator.p.uaux[:,end], Float32(dt),
                        integrator.p.mesh.Δeffective_s,
                        integrator,
                        integrator.p.SD; visc=inputs[:μ])
            # end
            write_output(integrator.p.SD, integrator.u, integrator.p.uaux, integrator.t, idx,
                         integrator.p.mesh, integrator.p.mp,
                         integrator.p.connijk_original, integrator.p.poin_in_bdy_face_original,
                         integrator.p.x_original, integrator.p.y_original, integrator.p.z_original,
                         inputs[:output_dir], inputs,
                         integrator.p.qp.qvars,
                         integrator.p.qp.qoutvars,
                         inputs[:outformat];
                         nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe)
            if (lwrite_time == true) 
                append_pvd_entry(pvd_path, integrator.t, "iter_$(idx).pvtu")
            end
        end
    end
    cb_les_stat    = DiscreteCallback(les_stat_condition, do_les_statistics!)
    cb_les_online  = DiscreteCallback(les_online_condition, do_les_online!)

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
    # Skipped on VTK restart — the snapshot already exists in the output dir
    # and its entry is already in simulation.pvd from the previous run.
    #
    idx  = (inputs[:tinit] == 0.0) ? 0 : findfirst(x -> x == inputs[:tinit], dosetimes)
    if idx ≠ nothing && !get(inputs, :lrestart_vtk, false)
        if rank == 0 println(" # Write initial condition to ",  typeof(inputs[:outformat]), " .........") end
        write_output(params.SD, u, params.uaux, inputs[:tinit], idx,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe)
        if (lwrite_time == true)
            append_pvd_entry(pvd_path, inputs[:tinit], "iter_$(idx).pvtu")
        end
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
        @mystop()
    else
        ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl), MPI.MAX, comm)
        dt         = Float32(inputs[:Δt]/(2.0^(ad_lvl_max)))
        solution = solve(prob,
                         inputs[:ode_solver], dt=dt,
                         #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                         callback = CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online), tstops = tstops_all,
                         save_everystep = false,
                         adaptive=inputs[:ode_adaptive_solver],
                         saveat = range(inputs[:tinit],
                                        inputs[:tend],
                                        length=inputs[:ndiagnostics_outputs]));
    end
    
    MPI.Barrier(comm)
    report_all_timers(params.timers)
    MPI.Barrier(comm)

    if inputs[:lamr] == true
        while solution.t[end] < inputs[:tend]
            @time prob, partitioned_model = amr_strategy!(inputs, prob.p, solution.u[end][:], solution.t[end], partitioned_model)
            ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl), MPI.MAX, comm)
            dt         = Float32(inputs[:Δt]/(2.0^(ad_lvl_max)))
            println_rank(" #  dt=", dt; msg_rank = rank)
            @time solution = solve(prob,
                                inputs[:ode_solver], dt=Float32(dt),
                                callback = CallbackSet(cb_amr, cb_restart), tstops = dosetimes,
                                save_everystep = false,
                                adaptive=inputs[:ode_adaptive_solver],
                                saveat = []);
            MPI.Barrier(comm)
            report_all_timers(prob.p.timers)
            MPI.Barrier(comm)
        end
    end

    # Finalize online statistics: single Allreduce + write output
    if stats_online_start < Inf
        les_finalize_online!(params, solution.t[end])
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)

    return solution
end

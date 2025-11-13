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

    # Time averaging variables
    tavg_timestep_counter = Ref{Int}(0)
    tavg_start_time       = inputs[:tavg_start_time]
    tavg_end_time         = inputs[:tavg_end_time]
    tavg_every_n          = inputs[:tavg_every_n_timesteps]
    ltavg                 = inputs[:ltavg]
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

    # Time averaging callback
    # Accumulate whenever condition is true (within averaging window)
    function tavg_condition(u, t, integrator)
        if !ltavg
            return false
        end
        # Check if we're within the averaging window
        in_window = (t >= tavg_start_time && t <= tavg_end_time)
        if !in_window
            return false
        end

        # Apply sampling frequency: only trigger every N-th call
        tavg_timestep_counter[] += 1
        should_sample = (mod(tavg_timestep_counter[], tavg_every_n) == 0)

        return should_sample
    end

    function do_tavg!(integrator)
        # Initialize averaging window if first sample
        if integrator.p.sample_count[] == 0
            integrator.p.t_start[] = integrator.t
            println_rank(" # Starting time averaging at t=", integrator.t; msg_rank = rank)
        end

        # Accumulate time-averaged quantities
        # Read directly from ODE state vector u with proper indexing
        npoin = integrator.p.mesh.npoin
        neqs = integrator.p.qp.neqs

        for ieq = 1:neqs
            # u is stored as: [all points for eq1, all points for eq2, ...]
            # So for equation ieq, the data starts at index (ieq-1)*npoin + 1
            idx_start = (ieq-1)*npoin
            for ip = 1:npoin
                val = integrator.u[idx_start + ip]
                integrator.p.q_tavg[ip, ieq] += val
                integrator.p.q2_tavg[ip, ieq] += val * val
            end
        end

        integrator.p.sample_count[] += 1
        integrator.p.t_end[] = integrator.t

        # Print progress every 100 samples
        if mod(integrator.p.sample_count[], 100) == 0
            println_rank(" # Time-averaging: ", integrator.p.sample_count[], " samples at t=", integrator.t; msg_rank = rank)
        end
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
    cb_tavg    = DiscreteCallback(tavg_condition, do_tavg!)
    CallbackSet(cb, cb_tavg)#,cb_rad)
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
        @mystop()
    else
        solution = solve(prob,
                         inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                         #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                         callback = CallbackSet(cb, cb_restart, cb_tavg), tstops = dosetimes,
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

    # Write time-averaged output if enabled
    if ltavg && params.sample_count[] > 0
        println_rank(" # Writing time-averaged output .................... "; msg_rank = rank)
        println_rank(" #   Samples collected: ", params.sample_count[]; msg_rank = rank)
        println_rank(" #   Averaging window: t=", params.t_start[], " to t=", params.t_end[]; msg_rank = rank)

        # Compute the mean by dividing accumulated values by sample count
        n_samples = params.sample_count[]
        println_rank(" #   Dividing by n_samples = ", n_samples; msg_rank = rank)

        q_tavg_mean_2d  = params.q_tavg  ./ n_samples
        q2_tavg_mean_2d = params.q2_tavg ./ n_samples

        # Debug: print sample values to verify averaging
        if rank == 0 && params.mesh.npoin > 0
            println_rank(" #   Sample check - accumulated[1,1] = ", params.q_tavg[1,1]; msg_rank = rank)
            println_rank(" #   Sample check - mean[1,1] = ", q_tavg_mean_2d[1,1]; msg_rank = rank)
            println_rank(" #   Sample check - ratio = ", params.q_tavg[1,1] / q_tavg_mean_2d[1,1]; msg_rank = rank)
        end

        # Compute variance: Var(X) = E[X²] - E[X]²
        q_tavg_var_2d = q2_tavg_mean_2d .- (q_tavg_mean_2d .^ 2)

        # Convert from 2D (npoin, neqs) to flat (npoin*neqs) format for write_output
        npoin = params.mesh.npoin
        neqs = params.qp.neqs
        q_tavg_mean_flat = zeros(npoin * neqs)
        q_tavg_var_flat  = zeros(npoin * neqs)

        uaux2u!(q_tavg_mean_flat, q_tavg_mean_2d, neqs, npoin)
        uaux2u!(q_tavg_var_flat, q_tavg_var_2d, neqs, npoin)

        # Create output directory for time-averaged data
        tavg_output_dir = joinpath(inputs[:output_dir], "time_averaged")
        if rank == 0 && !isdir(tavg_output_dir)
            mkpath(tavg_output_dir)
        end
        MPI.Barrier(comm)

        # Write mean (pass flat array as sol, 2D array as uaux)
        write_output(params.SD, q_tavg_mean_flat, q_tavg_mean_2d, params.t_end[], 1,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     tavg_output_dir, inputs,
                     params.qp.qvars,
                     params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe, case="mean")

        # Write variance (pass flat array as sol, 2D array as uaux)
        write_output(params.SD, q_tavg_var_flat, q_tavg_var_2d, params.t_end[], 2,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     tavg_output_dir, inputs,
                     params.qp.qvars,
                     params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe, case="variance")

        println_rank(" # Writing time-averaged output .................... DONE"; msg_rank = rank)
    end
    # MPI.Barrier(comm)
    # report_all_timers(params.timers)
    # MPI.Barrier(comm)

    return solution
end

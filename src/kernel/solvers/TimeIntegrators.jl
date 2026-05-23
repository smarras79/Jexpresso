using TrixiBase
using TimerOutputs

"""
    alloc_summary_enabled(inputs) -> Bool

Decide whether the end-of-run per-function timing & allocation summary
table is printed. Precedence (highest first):

 1. ENV variable `JEXPRESSO_ALLOC_SUMMARY` — set it from the REPL before
    running, e.g. `ENV["JEXPRESSO_ALLOC_SUMMARY"] = "false"`.
 2. Command-line flag — pass `--no-alloc-summary` (or `no-alloc-summary`)
    in Julia's `ARGS`, e.g. `julia run.jl ... --no-alloc-summary`.
 3. `:lalloc_summary` in the case's `user_inputs.jl`.
 4. Default: `true` (on).

Accepted truthy ENV values: 1/true/yes/on; falsy: 0/false/no/off
(case-insensitive).
"""
function alloc_summary_enabled(inputs)
    e = get(ENV, "JEXPRESSO_ALLOC_SUMMARY", nothing)
    if e !== nothing
        v = lowercase(strip(e))
        v in ("0", "false", "no", "off")  && return false
        v in ("1", "true", "yes", "on")   && return true
    end
    if any(a -> a in ("--no-alloc-summary", "no-alloc-summary"), ARGS)
        return false
    end
    return get(inputs, :lalloc_summary, true) == false
end

function time_loop!(inputs, params, u, args...)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    partitioned_model = args[1]
    # Optional coupled-mode positional args: args[2] = is_coupled::Bool,
    # args[3] = coupling::CouplingData. Both default to "off" when callers
    # use the historical 1-arg form (drivers.jl standalone path).
    is_coupled = length(args) >= 2 ? args[2] : false
    coupling   = length(args) >= 3 ? args[3] : nothing
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)

    #------------------------------------------------------------------------
    # Diagnostic: identify any AbstractDict-typed field inside params before
    # SciMLBase introspects it.  Any such field will trigger the
    # "Using arrays or dicts to store parameters of different types" warning
    # from SciMLBase/performance_warnings.jl — print it so the offender is
    # named, then leave the run untouched.  Set :diag_params_dicts => false
    # in user_inputs to silence.  Rank 0 only.
    #------------------------------------------------------------------------
    if rank == 0 && get(inputs, :diag_params_dicts, true) == true
        offenders = Tuple{Symbol, DataType}[]
        for (k, v) in pairs(params)
            v isa AbstractDict && push!(offenders, (k, typeof(v)))
        end
        if !isempty(offenders)
            println(" # [params-dict-scan] AbstractDict fields in params (these trigger the SciMLBase perf warning):")
            for (k, t) in offenders
                println("     • params.$(k) :: $(t)")
            end
            flush(stdout)
        end
    end

    # FullSpecialize: no type-erasing FunctionWrapper around rhs!. The
    # default AutoSpecialize wrapper type-erases p=params, so every
    # params.* access inside rhs!/inviscid_rhs_el!/viscous_rhs_el! boxes
    # (the ~2 KiB/call seen in the allocation summary). FullSpecialize
    # makes the integrator specialize on the concrete rhs!/params/du
    # types. Requires RHStoDU! to scalar-assign du so low-storage RK's
    # ArrayFuse du works without the wrapper.
    prob = ODEProblem{true, FullSpecialize}(rhs!,
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
    lwrite_time  = (inputs[:outformat] == VTK()) && (rank == 0)

    if (lwrite_time == true) 
        pvd_path = joinpath(inputs[:output_dir], "simulation.pvd")
        init_pvd_file(pvd_path)
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
            if (lwrite_time == true) 
                append_pvd_entry(pvd_path, integrator.t, "iter_$(idx).pvtu")
            end
        end
    end
    cb_rad     = DiscreteCallback(two_stream_condition, do_radiation!)
    cb         = DiscreteCallback(condition, affect!)    
    cb_amr     = DiscreteCallback(condition, affect!)
    cb_restart = DiscreteCallback(restart_condition, do_restart!)
    cb_coupling = is_coupled ? setup_coupling_callback(is_coupled, params, inputs) : nothing
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
        if (lwrite_time == true) 
            append_pvd_entry(pvd_path, inputs[:tinit], "iter_$(idx).pvtu")
        end
        if rank == 0  println(" # Write initial condition to ",  typeof(inputs[:outformat]), " ......... END") end
    end

    TimerOutputs.reset_timer!(TrixiBase.timer())
    #
    # Simulation
    #
    callbacks_main = (is_coupled && cb_coupling !== nothing) ?
        CallbackSet(cb, cb_restart, cb_coupling) :
        CallbackSet(cb, cb_restart)

    # Warm-up: run the full rhs! path once, untimed, so JIT compilation
    # of rhs!/_build_rhs!/the inviscid·viscous barriers/expansions/BCs/DSS
    # is excluded from the summary. Without this the @timeit table counts
    # the first (compiling) call, inflating per-call allocations by the
    # one-time compile cost / n_calls (huge on short runs). Skipped when
    # the summary is disabled so we don't add an eval needlessly.
    if alloc_summary_enabled(inputs)
        let du_warmup = similar(u)
            rhs!(du_warmup, u, prob.p, params.tspan[1])
        end
        MPI.Barrier(comm)
    end

    TimerOutputs.reset_timer!(JEXPRESSO_TIMER)
    rank == 0 && println(" # Simulation timing and allocations (steady state; compile warm-up excluded):")
    solution = @time solve(prob,
                           inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                           #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                           callback = callbacks_main, tstops = dosetimes,
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit],
                                          inputs[:tend],
                                          length=inputs[:ndiagnostics_outputs]));

    # End-of-simulation per-function timing & allocation summary table.
    # On by default. Disable via ENV["JEXPRESSO_ALLOC_SUMMARY"]="false"
    # (REPL), the --no-alloc-summary CLI flag, or :lalloc_summary => false
    # in user_inputs.jl. See alloc_summary_enabled above for precedence.
    if rank == 0 && alloc_summary_enabled(inputs)
        println("\n # Per-function timing & allocation summary (disable via ENV[\"JEXPRESSO_ALLOC_SUMMARY\"]=\"false\", --no-alloc-summary, or :lalloc_summary => false):")
        show(stdout, JEXPRESSO_TIMER; allocations=true, sortby=:firstexec)
        println()
    end
    MPI.Barrier(comm)
    report_all_timers(params.timers)
    MPI.Barrier(comm)
    
    if inputs[:lamr] == true
        while solution.t[end] < inputs[:tend]
            @time prob, partitioned_model = amr_strategy!(inputs, prob.p, solution.u[end][:], solution.t[end], partitioned_model)
            
            @time solution = solve(prob,
                                   inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                                   callback = CallbackSet(cb_amr, cb_restart), tstops = dosetimes,
                                   save_everystep = false,
                                   adaptive=inputs[:ode_adaptive_solver],
                                   saveat = []);
            MPI.Barrier(comm)
            report_all_timers(prob.p.timers)
            MPI.Barrier(comm)
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)

    return solution
end

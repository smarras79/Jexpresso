using Logging: with_logger, NullLogger, current_logger

"""
    alloc_summary_enabled(inputs) -> Bool

Decide whether the end-of-run per-function timing & allocation summary
table is printed at end of `time_loop!`. Precedence (highest first):

 1. ENV variable `JEXPRESSO_ALLOC_SUMMARY` — set from shell or REPL
    before running, e.g.
      `ENV["JEXPRESSO_ALLOC_SUMMARY"] = "1"`     (REPL, enable)
      `JEXPRESSO_ALLOC_SUMMARY=1 julia ...`      (shell, enable)
      `JEXPRESSO_ALLOC_SUMMARY=false julia ...`  (shell, disable)
    For coupled-mode mpirun, pass via mpirun's `-x`:
      `mpirun -np 2 ./AlyaProxy/Alya.x : \\`
      `       -x JEXPRESSO_COUPLED=1 -x JEXPRESSO_ALLOC_SUMMARY=1 \\`
      `       -np 2 julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya`
 2. Command-line flag — pass `--no-alloc-summary` (or `no-alloc-summary`)
    in Julia's `ARGS`, e.g. `julia run.jl ... --no-alloc-summary`.
 3. `:lalloc_summary` in the case's `user_inputs.jl` (Bool).
 4. Default: `false` (off) — the summary adds an extra full-RHS warm-up
    call before the real solve, so it's opt-in for performance runs.

Accepted truthy ENV values: 1/true/yes/on; falsy: 0/false/no/off
(case-insensitive).

Ported from origin/sm/alyacouple-merge:src/kernel/solvers/TimeIntegrators.jl
after the giga_les wholesale swap dropped it.
"""
function alloc_summary_enabled(inputs)
    e = get(ENV, "JEXPRESSO_ALLOC_SUMMARY", nothing)
    if e !== nothing
        v = lowercase(strip(e))
        v in ("0", "false", "no", "off") && return false
        v in ("1", "true", "yes", "on")  && return true
    end
    if any(a -> a in ("--no-alloc-summary", "no-alloc-summary"), ARGS)
        return false
    end
    return get(inputs, :lalloc_summary, false) == true
end

"""
    precompile_warmup_enabled(inputs) -> Bool

Decide whether `time_loop!` runs a one-shot RHS warm-up call before
the real `solve(...)`. The warm-up triggers JIT compilation of
`rhs!`, `_build_rhs!`, and the rest of the per-step kernel chain on
an arbitrary (initial-condition) state, then the real run starts
already compiled.

Most useful when launching from the command line, where every run
incurs JIT-compile cost on the first timestep that REPL users avoid
by re-running inside the same session. Default: ON.

Precedence (highest first):

 1. ENV variable `JEXPRESSO_PRECOMPILE_WARMUP`. Examples:
      `JEXPRESSO_PRECOMPILE_WARMUP=0 julia ...`           (disable)
      `JEXPRESSO_PRECOMPILE_WARMUP=false julia ...`       (disable)
    For coupled-mode mpirun, pass via mpirun's `-x`:
      `mpirun -np 2 ./AlyaProxy/Alya.x : \\`
      `       -x JEXPRESSO_COUPLED=1 -x JEXPRESSO_PRECOMPILE_WARMUP=0 \\`
      `       -np 2 julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya`
 2. Command-line flag — pass `--no-precompile-warmup` in Julia's `ARGS`.
 3. `:lprecompile_warmup` in the case's `user_inputs.jl` (Bool).
 4. Default: `true` (on).

Note: when `alloc_summary_enabled(inputs)` is true, the warm-up runs
unconditionally (the alloc summary needs the post-JIT measurement
window to be meaningful). The flag here only controls the
warm-up-without-alloc-summary case.

Accepted truthy ENV values: 1/true/yes/on; falsy: 0/false/no/off
(case-insensitive).
"""
function precompile_warmup_enabled(inputs)
    e = get(ENV, "JEXPRESSO_PRECOMPILE_WARMUP", nothing)
    if e !== nothing
        v = lowercase(strip(e))
        v in ("0", "false", "no", "off") && return false
        v in ("1", "true", "yes", "on")  && return true
    end
    if any(a -> a in ("--no-precompile-warmup", "no-precompile-warmup"), ARGS)
        return false
    end
    return get(inputs, :lprecompile_warmup, true) == true
end

function time_loop!(inputs, params, u, args...)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    partitioned_model = args[1]
    # Optional coupled-mode positional args: args[2] = is_coupled::Bool,
    # args[3] = coupling::CouplingData. Both default to "off" when callers
    # use the historical 1-arg form (standalone non-coupled runs).
    is_coupled = length(args) >= 2 ? args[2] : false
    coupling   = length(args) >= 3 ? args[3] : nothing
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    # FullSpecialize: SciMLBase's default AutoSpecialize wraps rhs! in a
    # FunctionWrapper that type-erases `params` to ::Any. That defeats
    # type inference inside the entire RHS chain - every `params.field`
    # access in rhs!/_build_rhs!/inviscid_rhs_el!/viscous_rhs_el! boxes,
    # adding ~2 KiB / RK stage on this case. FullSpecialize keeps the
    # concrete params type all the way down. Paired with the typed
    # function barriers below in inviscid_rhs_el!/viscous_rhs_el!.
    prob = ODEProblem{true, FullSpecialize}(rhs!,
                      u,
                      params.tspan,
                      params);
    
    #------------------------------------------------------------------------
    # Runtime callbacks
    #------------------------------------------------------------------------
    dosetimes    = inputs[:diagnostics_at_times]
    les_stat_t   = inputs[:statistics_time]
    tstops_all   = sort(unique(vcat(collect(Float64, dosetimes), collect(Float64, les_stat_t))))
    idx_ref      = Ref{Int}(0)
    c            = Float64(0.0)
    restart_time = inputs[:restart_time]
    stats_online_start = get(inputs, :statistics_online_start, Inf)
    rad_time           = inputs[:radiation_time_step]
    lnew_mesh    = true   
    lwrite_time  = (inputs[:outformat] == VTK()) && (rank == 0)
    lwrite_init  = !(inputs[:lrestart] || inputs[:lrestart_vtk] || inputs[:lrestart_amr]) 

    if (lwrite_time == true)
        pvd_path = joinpath(inputs[:output_dir], "simulation.pvd")
        if !lwrite_init && isfile(pvd_path)
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
        #println(" # doing rad test")
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
        # println_rank(" # LES online accumulation at t=", integrator.t; msg_rank = rank)
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

                ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl; init=0), MPI.MAX, comm)
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
            # Save p4est forest checkpoint alongside VTK for AMR restart support.
            # p8est_save is MPI-collective: all ranks call together.
            # Julia closures capture `partitioned_model` by binding — it always
            # reflects the current forest after each AMR iteration.
            if get(inputs, :lamr, false)
                write_p4est_checkpoint(inputs[:output_dir], idx, partitioned_model)
            end
        end
    end
    cb_les_stat    = DiscreteCallback(les_stat_condition, do_les_statistics!)
    cb_les_online  = DiscreteCallback(les_online_condition, do_les_online!)

    cb_rad     = DiscreteCallback(two_stream_condition, do_radiation!)
    cb         = DiscreteCallback(condition, affect!)
    cb_amr     = DiscreteCallback(condition, affect!)
    cb_restart = DiscreteCallback(restart_condition, do_restart!)
    # Coupled-mode exchange callback: fires once per accepted timestep,
    # sends Julia's interpolated solution to Alya so its MPI.Waitall in
    # the time loop can advance. Without this Alya hangs and never
    # writes its VTS output.
    cb_coupling = is_coupled ? setup_coupling_callback(is_coupled, params, inputs) : nothing
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
    if idx ≠ nothing && lwrite_init
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
        ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl; init=0), MPI.MAX, comm)
        dt         = Float32(inputs[:Δt]/(2.0^(ad_lvl_max)))
        # Include cb_coupling in coupled mode so Julia's per-timestep
        # send to Alya actually fires; without it Alya's MPI.Waitall
        # blocks and its VTS output never gets written.
        callbacks_main = (is_coupled && cb_coupling !== nothing) ?
            CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online, cb_coupling) :
            CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online)

        # Precompile warm-up: run the full rhs! path once UNTIMED so JIT
        # compilation of rhs!, _build_rhs!, inviscid/viscous kernels,
        # BCs, DSS, etc. is excluded from the real run's wall time. From
        # the command line, the first timestep of every fresh process
        # otherwise pays the full JIT bill (~5–40s depending on the
        # case), which dominates short runs and inflates per-call alloc
        # measurements. Skipped only when explicitly disabled via
        # JEXPRESSO_PRECOMPILE_WARMUP=0 / --no-precompile-warmup /
        # :lprecompile_warmup => false.
        #
        # Coupled-mode safety: this is a LOCAL Julia call (rhs! + an
        # MPI.Barrier on the Julia sub-comm). It does NOT exchange data
        # with Alya - the coupling callback `cb_coupling` only fires
        # inside `solve(...)`, not during the warm-up. Safe to run
        # before Alya's first send/recv.
        if precompile_warmup_enabled(inputs) || alloc_summary_enabled(inputs)
            rank == 0 && (print(" # Precompile warm-up (1 rhs! call) ......... "); flush(stdout))
            t0_warmup = time_ns()
            let du_warmup = similar(u)
                rhs!(du_warmup, u, prob.p, params.tspan[1])
            end
            MPI.Barrier(comm)
            rank == 0 && @printf("%.2f s\n", (time_ns() - t0_warmup) / 1e9)
            # If the alloc summary is also on, reset JEXPRESSO_TIMER so
            # the warm-up's allocations aren't counted in the summary.
            if alloc_summary_enabled(inputs)
                TimerOutputs.reset_timer!(JEXPRESSO_TIMER)
                rank == 0 && println(" # Simulation timing and allocations (steady state; compile warm-up excluded):")
            end
        end

        # Silence SciMLBase's per-call "Using arrays or dicts to store
        # parameters of different types can hurt performance" warning on
        # non-root ranks - the warning is identical from every rank so
        # printing it nparts times is pure noise.  Root rank still sees
        # the warning once, which is the right amount.
        solve_logger = rank == 0 ? current_logger() : NullLogger()
        solution = with_logger(solve_logger) do
            solve(prob,
                  inputs[:ode_solver], dt=dt,
                  #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                  callback = callbacks_main, tstops = tstops_all,
                  save_everystep = false,
                  adaptive=inputs[:ode_adaptive_solver],
                  saveat = range(inputs[:tinit],
                                 inputs[:tend],
                                 length=inputs[:ndiagnostics_outputs]))
        end

        # End-of-simulation per-function timing & allocation summary.
        # Set ENV["JEXPRESSO_ALLOC_SUMMARY"]="1" or :lalloc_summary => true
        # to enable; see alloc_summary_enabled above for precedence.
        if rank == 0 && alloc_summary_enabled(inputs)
            println("\n # Per-function timing & allocation summary (set ENV[\"JEXPRESSO_ALLOC_SUMMARY\"]=\"0\" to disable):")
            show(stdout, JEXPRESSO_TIMER; allocations=true, sortby=:firstexec)
            println()
        end
    end
    
    MPI.Barrier(comm)
    report_all_timers(params.timers)
    MPI.Barrier(comm)

    if inputs[:lamr] == true
        while solution.t[end] < inputs[:tend]
            @mpi_time prob, partitioned_model = amr_strategy!(inputs, prob.p, solution.u[end][:], solution.t[end], partitioned_model)
            ad_lvl_max = MPI.Allreduce(maximum(prob.p.mesh.ad_lvl; init=0), MPI.MAX, comm)
            dt         = Float32(inputs[:Δt]/(2.0^(ad_lvl_max)))
            println_rank(" #  dt=", dt; msg_rank = rank)
            @mpi_time solution = solve(prob,
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
    # Finalize Approach 1 statistics: write final time-and-space averages
    les_finalize!(params, solution.t[end])
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)

    return solution
end

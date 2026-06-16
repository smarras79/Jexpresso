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

"""
    precompile_warmup_run!(inputs, params, u, partitioned_model, is_coupled, coupling)

Run one timestep of the real solve as a JIT warm-up, snapshotting and
restoring mutable state so the actual time loop starts from a clean
initial condition. Returns `nothing`.

Called from drivers.jl BEFORE the `@time time_loop!(...)` wrapper so the
warm-up's JIT-compile allocations are NOT counted in the displayed
total. A 1-step `solve(...)` (not just an `rhs!` call) compiles the
SciML integrator's `step!`, the callback dispatch machinery, and
`save_*` paths in addition to the entire RHS chain - matches the
"REPL second run" warmth.

No-op when both `precompile_warmup_enabled(inputs) == false` AND
`alloc_summary_enabled(inputs) == false`.

Coupled-mode safety: the warm-up `solve` uses `CallbackSet()` (empty),
so the coupling exchange callback (`cb_coupling`) and the diagnostic
output callbacks are NOT registered for this 1-step run. Julia does
not exchange data with Alya during warm-up; only the local Julia
sub-comm is touched (via an `MPI.Barrier` at the end).

State preserved across warm-up: `u`, `params.qp.qnm1`, `params.qp.qnm2`.
Other params working arrays (RHS, uaux, fluxes, ...) get overwritten on
the real solve's first `rhs!` call.
"""
function precompile_warmup_run!(inputs, params, u,
                                partitioned_model, is_coupled, coupling)
    (precompile_warmup_enabled(inputs) || alloc_summary_enabled(inputs)) || return nothing

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    rank == 0 && (print(" # Precompile warm-up (1 step solve) ......... "); flush(stdout))
    t0 = time_ns()

    # Snapshot mutable state that the warm-up step would advance.
    u_snapshot    = copy(u)
    qnm1_snapshot = copy(params.qp.qnm1)
    qnm2_snapshot = copy(params.qp.qnm2)

    # 1-step problem; same params and same FullSpecialize as the real
    # solve, so the compiled code is reused.
    Δt_warmup    = Float32(inputs[:Δt])
    t0_warmup    = params.tspan[1]
    warmup_tspan = (t0_warmup, t0_warmup + Δt_warmup)
    warmup_prob  = ODEProblem{true, FullSpecialize}(rhs!, u, warmup_tspan, params)

    # Silence all log output during the warm-up step. Also silences the
    # SciMLBase "Using arrays or dicts..." warning so it doesn't appear
    # twice (warm-up + real solve).
    with_logger(NullLogger()) do
        try
            solve(warmup_prob,
                  inputs[:ode_solver], dt=Δt_warmup,
                  callback = CallbackSet(),
                  save_everystep = false,
                  save_start = false,
                  save_end = false,
                  adaptive = false)
        catch e
            # Warm-up failure must not prevent the real run. Surface a
            # one-line warning on rank 0 and continue uncompiled - the
            # real solve will JIT on its first step as before.
            rank == 0 && @warn "precompile warm-up failed; continuing without it" exception=e
        end
    end

    # Restore u and the qnm1/qnm2 history slots so the real solve sees
    # the same initial condition the caller passed in.
    u .= u_snapshot
    params.qp.qnm1 .= qnm1_snapshot
    params.qp.qnm2 .= qnm2_snapshot

    # The VTK write path is JIT-compiled on first use by the IC write in
    # time_loop! (when :lwrite_initial is true) or by the first diagnostic
    # callback (when false). No need for a separate throw-away write here.

    MPI.Barrier(comm)
    rank == 0 && @printf("%.2f s\n", (time_ns() - t0) / 1e9)

    # Reset JEXPRESSO_TIMER so the alloc summary table only reflects
    # steady-state allocations from the real solve. The
    # @timeit_debug compile-time gate is set at module load time in
    # src/auxiliary/timing.jl based on JEXPRESSO_ALLOC_SUMMARY - it
    # MUST happen before any rhs!-containing function is JIT-compiled
    # (calling enable_debug_timings here would be too late).
    if alloc_summary_enabled(inputs)
        TimerOutputs.reset_timer!(JEXPRESSO_TIMER)
    end
    return nothing
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

    #------------------------------------------------------------------------
    # Write the initial condition.
    #
    # This runs FIRST, before the integrator warm-up below: the warm-up
    # temporarily swaps inputs[:output_dir] to a throw-away tempdir, and the
    # IC write must never interact with that. It also doubles as the JIT
    # warm-up of the write_output/write_vtk path.
    #
    # File number: the slot of tinit in dosetimes (1 for a t=0 start). The
    # diagnostic callback never fires at the initial time, so this slot is
    # otherwise unused and the sequence is contiguous: iter_1 (IC), iter_2
    # (first diagnostic hit), ...
    #
    # When skipped, say so and why — never silently.
    #------------------------------------------------------------------------
    idx = something(findfirst(x -> x == inputs[:tinit], dosetimes), 1)
    if !lwrite_init
        println_rank(" # Skipping initial-condition write (restart run)"; msg_rank = rank)
    elseif get(inputs, :lwrite_initial, true) != true
        println_rank(" # Skipping initial-condition write (:lwrite_initial => false)"; msg_rank = rank)
    else
        println_rank(" # Write initial condition to ", typeof(inputs[:outformat]),
                     " in ", inputs[:output_dir], " ........."; msg_rank = rank)
        write_output(params.SD, u, params.uaux, inputs[:tinit], idx,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe,
                     μ_dsgs_pnode = (params.VT == DSGS()) ? params.μ_dsgs_pnode : nothing)
        if (lwrite_time == true)
            append_pvd_entry(pvd_path, inputs[:tinit], "iter_$(idx).pvtu")
        end
        println_rank(" # Write initial condition to ", typeof(inputs[:outformat]),
                     " ......... END"; msg_rank = rank)
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
                         nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe,
                         μ_dsgs_pnode = (integrator.p.VT == DSGS()) ? integrator.p.μ_dsgs_pnode : nothing)
            # The DSGS viscosity panel is rendered by the 1D PNG writer
            # itself (write_output -> plot_results, fed by μ_dsgs_pnode
            # above) so that the whole output time is a single GR render:
            # a separate plot_dsgs_1d call here would either flash its own
            # gksqt window or, via the silent export path, close the GKS
            # session and with it the live plot-matrix window.
            if (lwrite_time == true)
                append_pvd_entry(pvd_path, integrator.t, "iter_$(idx).pvtu")
            end
            # Save p4est forest checkpoint alongside VTK for AMR restart support.
            # p8est_save is MPI-collective: all ranks call together.
            # Julia closures capture `partitioned_model` by binding — it always
            # reflects the current forest after each AMR iteration.
            if get(inputs, :lamr, false)
                # PERF: GridapP4est is lazy-loaded; ensure it's in
                # scope before the p8est_save underlying this call.
                _ensure_amr_loaded!()
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
        # DEBUG: per-step heartbeat so the user can tell whether solve()
        # is making progress between diagnostic writes. Diagnostic
        # writes only fire at dosetimes (e.g. every 10 time units for
        # city2d), so for Δt=0.004 the integrator can be silently doing
        # 2500 steps between two user-visible prints — looks identical
        # to a hang.
        #
        # Throttled: every step for the first 5, then every 100. With
        # 150 000-step runs that's ~1 500 lines of output instead of
        # 150 000. DEFAULT IS OFF (debugging-only); opt in with
        # `:lstep_heartbeat => true` in user_inputs.jl or via env
        # `JEXPRESSO_STEP_HEARTBEAT=1`. The env var, if set, takes
        # precedence over the user_inputs.jl flag.
        _step_count = Ref{Int}(0)
        function step_heartbeat_condition(u, t, integrator)
            _step_count[] += 1
            n = _step_count[]
            return n <= 5 || (n % 100 == 0)
        end
        function step_heartbeat_affect!(integrator)
            rank == 0 && (@printf(" #   step %d   t = %.6f\n", _step_count[], integrator.t); flush(stdout))
        end
        # Default OFF. Env var, if set, wins over user_inputs.jl.
        _env_hb = lowercase(strip(get(ENV, "JEXPRESSO_STEP_HEARTBEAT", "")))
        _heartbeat_on = if _env_hb in ("1", "true", "yes", "on")
            true
        elseif _env_hb in ("0", "false", "no", "off")
            false
        else
            get(inputs, :lstep_heartbeat, false) == true
        end
        cb_heartbeat = _heartbeat_on ?
            DiscreteCallback(step_heartbeat_condition, step_heartbeat_affect!) :
            nothing

        callbacks_main = if is_coupled && cb_coupling !== nothing
            cb_heartbeat === nothing ?
                CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online, cb_coupling) :
                CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online, cb_coupling, cb_heartbeat)
        else
            cb_heartbeat === nothing ?
                CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online) :
                CallbackSet(cb, cb_restart, cb_les_stat, cb_les_online, cb_heartbeat)
        end

        # PERF: SciML integrator warmup with the REAL callback set.
        #
        # The outer precompile_warmup_run! (called from drivers.jl) uses
        # `CallbackSet()` because it does not have access to `cb`, `cb_restart`,
        # … which are constructed inside time_loop!. That warmup pre-JITs
        # `rhs!` but the integrator's `step!` is specialised on the
        # CallbackSet *type*, and CallbackSet{Tuple{…}} ≠ CallbackSet{}, so
        # the first real `solve(…)` below was recompiling the entire
        # integrator on first hit — the user's "long wait after initial-VTK
        # END" wall.
        #
        # Run one throw-away step here with the same prob/callbacks_main type
        # the real solve uses, then snapshot/restore u and qnm1/qnm2 so the
        # production run sees the original IC. The throw-away step's
        # diagnostic-VTK output, if any, goes to a per-rank mktempdir that's
        # removed right after.
        if precompile_warmup_enabled(inputs)
            rank == 0 && (print(YELLOW_FG(" # Integrator warm-up with real callbacks (PATIENCE: ONLY DONE ON 1st RUN!) ......... ")); flush(stdout))
            _t_wm = time_ns()
            u_snap    = copy(u)
            qnm1_snap = copy(params.qp.qnm1)
            qnm2_snap = copy(params.qp.qnm2)
            # If a callback ends up actually writing during the warmup
            # (only possible for cases whose first dosetime falls inside
            # [t0, t0+Δt]), redirect that output to a per-rank tempdir.
            # When `inputs` is a NamedTuple (rare; user_inputs flag
            # :use_named_tuples) we can't mutate it — that's fine, only
            # a few cases hit the affect!() inside this single step,
            # and the temp redirect is best-effort anyway.
            warm_outdir  = mktempdir(; prefix = "jexpresso_intwarmup_")
            saved_outdir = inputs isa Dict ? inputs[:output_dir] : nothing
            inputs isa Dict && (inputs[:output_dir] = warm_outdir)
            t0_w = params.tspan[1]
            Δt_w = Float32(inputs[:Δt])
            # PERF: build the warmup problem with `remake(prob, …)` so it
            # has the IDENTICAL Julia type as `prob` — same `tspan`
            # container type (Vector vs Tuple was the previous mismatch),
            # same `u0` field type, same params type. SciML specialises
            # the integrator on the prob type + kwarg types, so any
            # type-level difference (even Vector vs Tuple for tspan) makes
            # the real solve recompile from scratch despite the warmup.
            warmup_prob = remake(prob; tspan = [t0_w, t0_w + Δt_w])
            # PERF: kwargs also match the real `solve(...)` exactly. The
            # only deviation is `tstops = [t0_w + Δt_w]` (just one point)
            # to keep the warmup cheap.
            warm_saveat = range(t0_w, t0_w + Δt_w, length = inputs[:ndiagnostics_outputs])
            with_logger(NullLogger()) do
                try
                    solve(warmup_prob,
                          inputs[:ode_solver], dt = dt,
                          callback = callbacks_main,
                          tstops = [t0_w + Δt_w],
                          save_everystep = false,
                          adaptive = inputs[:ode_adaptive_solver],
                          saveat = warm_saveat)
                catch e
                    rank == 0 && @warn "integrator warm-up failed; continuing without it" exception=e
                end
            end
            u .= u_snap
            params.qp.qnm1 .= qnm1_snap
            params.qp.qnm2 .= qnm2_snap
            inputs isa Dict && saved_outdir !== nothing && (inputs[:output_dir] = saved_outdir)
            try; rm(warm_outdir; recursive = true, force = true); catch; end
            # Reset the heartbeat counter so the real solve gets its
            # own first-5-steps detail (the warmup just consumed one).
            _step_count[] = 0
            MPI.Barrier(comm)
            #rank == 0 && (print(YELLOW_FG(@sprintf("DONE (%.2f s)\n", (time_ns() - _t_wm) / 1e9))); flush(stdout))
        end

        if alloc_summary_enabled(inputs)
            rank == 0 && println(" # Simulation timing and allocations (steady state; compile warm-up excluded):")
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

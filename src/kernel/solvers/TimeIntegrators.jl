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
    precompile_pass_enabled(inputs) -> Bool

Decide whether `time_loop!` runs the simulation as a **pre-compilation pass**
(exactly ONE timestep) followed by the production loop, instead of launching
the whole thing as a single `solve(...)` call.

Why this exists
---------------
On a large problem launched at full size, the first timestep is not a
timestep: it is Julia's JIT compiling the entire RHS chain, the SciML
integrator specialised on the real `CallbackSet`, the MPI halo-exchange
paths and the diagnostic/output path, plus the first-touch allocation of
every per-rank working array. On 1 500+ ranks that work happens *inside* the
production time loop, interleaved with collective communication: the ranks
that compiled quickly sit in `MPI_Wait` while the slow ones finish, the
measured step rate for the first hundreds of steps is meaningless, and the
long run inherits a heap full of compilation garbage.

How the split is done
---------------------
`solve(prob, alg; kw...)` is, by definition, `solve!(init(prob, alg; kw...))`.
The pass takes that apart:

    integrator = init(prob, alg; kwargs...)
    step!(integrator)      # PHASE 1 — one timestep: all JIT, all first touches
    GC.gc(); MPI.Barrier()  # compacted heap, all ranks aligned
    solve!(integrator)     # PHASE 2 — the rest of the run, hot

Phase 2 continues the **same integrator object**. This is not a restart:
there is no second `init`, and the step-size controller, the callback caches,
the FSAL history and `u` all carry across untouched. The trajectory is
therefore bit-for-bit the one a single `solve(...)` would have produced —
adaptive stepping included — and not one method is compiled twice.

Phase 1's step IS the simulation's first step: nothing is snapshotted,
restored or thrown away.

Relationship to `precompile_warmup_enabled`
-------------------------------------------
The older warm-ups (`precompile_warmup_run!` in drivers.jl and the
in-`time_loop!` integrator warm-up) run a throw-away step and then *restore*
the initial condition, so a run pays for two extra full-mesh RHS evaluations
whose results are discarded. When the pre-compilation pass is on, both are
skipped: the pass compiles the same code and keeps the work.

Failure is not fatal
--------------------
If `init` or the first `step!` throws, the advanced history is rolled back
and `time_loop!` falls through to the historical single-phase `solve(...)`,
JIT-ing on its first step exactly as it always did. That decision is made
collectively (an `Allreduce` over the per-rank outcome): one rank failing
sends every rank down the same path, because half a job in `solve!` and half
in a fresh `solve` is a hang, not a fallback.

Precedence (highest first)
--------------------------
 1. ENV `JEXPRESSO_PRECOMPILE_PASS`. Examples:
      `JEXPRESSO_PRECOMPILE_PASS=1 mpirun -np 1536 julia ...`   (enable)
      `JEXPRESSO_PRECOMPILE_PASS=0 julia ...`                   (disable)
    MPICH/Hydra propagates the environment by default; under OpenMPI pass it
    explicitly with `mpirun -x JEXPRESSO_PRECOMPILE_PASS ...`.
 2. Command-line flag — `--precompile-pass` / `--no-precompile-pass` in ARGS.
 3. `:lprecompile_pass` in the case's `user_inputs.jl` (Bool).
 4. Default: `false` (historical single-phase behaviour).

Accepted truthy ENV values: 1/true/yes/on; falsy: 0/false/no/off
(case-insensitive).
"""
function precompile_pass_enabled(inputs)
    e = get(ENV, "JEXPRESSO_PRECOMPILE_PASS", nothing)
    if e !== nothing
        v = lowercase(strip(e))
        v in ("0", "false", "no", "off") && return false
        v in ("1", "true", "yes", "on")  && return true
    end
    if any(a -> a in ("--no-precompile-pass", "no-precompile-pass"), ARGS)
        return false
    elseif any(a -> a in ("--precompile-pass", "precompile-pass"), ARGS)
        return true
    end
    return get(inputs, :lprecompile_pass, false) == true
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

No-op when `precompile_pass_enabled(inputs) == true` (the two-phase
pre-compilation pass in `time_loop!` subsumes this warm-up and, unlike it,
keeps the step), or when both `precompile_warmup_enabled(inputs) == false`
AND `alloc_summary_enabled(inputs) == false`.

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
    # The two-phase pre-compilation pass compiles the SAME rhs! chain, with the
    # real callback set, and keeps the resulting step. Running this throw-away
    # solve first would just buy a second full-mesh RHS evaluation for nothing.
    precompile_pass_enabled(inputs) && return nothing
    (precompile_warmup_enabled(inputs) || alloc_summary_enabled(inputs)) || return nothing

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    rank == 0 && (print(" # Precompile warm-up (1 step solve) ......... "); flush(stdout))
    t0 = time_ns()

    # Snapshot mutable state that the warm-up step would advance.
    u_snapshot    = copy(u)
    qnm1_snapshot = copy(params.qp.qnm1)
    qnm2_snapshot = copy(params.qp.qnm2)
    # DynSGS-MHD carries its own step-cadenced history plus the time stamp
    # that gates it; the warm-up step would advance both.
    dsgs_qnm1_snapshot = copy(params.dsgs_qnm1)
    dsgs_qnm2_snapshot = copy(params.dsgs_qnm2)
    dsgs_thist_snapshot = params.dsgs_thist[]

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
    params.dsgs_qnm1 .= dsgs_qnm1_snapshot
    params.dsgs_qnm2 .= dsgs_qnm2_snapshot
    params.dsgs_thist[] = dsgs_thist_snapshot

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

"""
    flatten_times(x) -> Vector{Float64}

Every time list a deck can write, flattened to plain seconds: a range, a
vector, a tuple of numbers, or -- the case this exists for -- a tuple that
MIXES numbers and ranges.

WHY THIS IS NOT JUST `collect(Float64, x)`. Decks specify these as tuples of
splatted ranges:

    :diagnostics_at_times => (0.0:100.0:1000.0..., 1000.0:500.0:9000.0...,
                              9000.0:10.0:tend...)

and leaving the `...` off ONE of them is silent at parse time: it builds a
tuple whose last element is a range instead of a number. `collect(Float64, ...)`
then dies with

    MethodError: Cannot `convert` an object of type StepRangeLen{...} to Float64

and it dies HERE, in time_loop!, i.e. after the mesh read, after the operator
build and after the setup self-check -- on the full rank count, minutes into a
job. That has now cost two runs of a 256-rank case on two different keys.

Flattening costs nothing on the correct input and removes the failure mode, so
a missing `...` becomes a no-op rather than a dead job.
"""
flatten_times(x) = (out = Float64[]; _push_times!(out, x); out)
_push_times!(out::Vector{Float64}, x::Number) = (push!(out, Float64(x)); out)
_push_times!(out::Vector{Float64}, x) = (for y in x; _push_times!(out, y); end; out)

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
    tstops_all   = sort(unique(vcat(flatten_times(dosetimes), flatten_times(les_stat_t))))
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
                     μ_dsgs_pnode = (params.VT == DSGS() || params.VT == DSGS_MHD()) ? params.μ_dsgs_pnode : nothing,
                     schlieren = maybe_compute_schlieren(inputs, params, u))
        if (lwrite_time == true)
            append_pvd_entry(pvd_path, inputs[:tinit], "iter_$(idx).pvtu")
        end
        println_rank(" # Write initial condition to ", typeof(inputs[:outformat]),
                     " ......... END"; msg_rank = rank)
    end

    #------------------------------------------------------------------------
    # HEVI wiring check, and the direction-wise stability report.
    #
    # The report costs one RHS evaluation and is the only thing that tells
    # "Δt is limited by vertical sound", which HEVI removes, apart from
    # "Δt is limited by the SGS diffusion", which it does not. That RHS call
    # is also what fills the SGS viscosity cache: without it ν_t is zero
    # everywhere and the viscous limits would read as infinite, which is the
    # one answer guaranteed to be wrong.
    #
    # rhs! writes through u (Dirichlet BCs, the filter) and rolls the qnm /
    # DynSGS histories, so everything it touches is snapshotted and put back
    # -- same set the integrator warm-up below restores.
    #------------------------------------------------------------------------
    if inputs[:ode_solver] isa HEVI_ARK
        hasproperty(params, :hevi) ||
            error("The deck asks for HEVI_ARK but params carries no HEVI cache. ",
                  "params_setup builds it only when hevi_enabled(inputs) is true, so ",
                  "set :ode_solver => HEVI_ARK(:ARS343) in user_inputs.jl rather than ",
                  "swapping the integrator in later.")
        inputs[:ode_adaptive_solver] == true &&
            error("HEVI_ARK is fixed-step: its tableaux carry no embedded error ",
                  "estimate. Set :ode_adaptive_solver => false.")
    end

    if inputs[:ode_solver] isa IMEX_ARK
        hasproperty(params, :imex) ||
            error("The deck asks for IMEX_ARK but params carries no IMEX3D cache. ",
                  "params_setup builds it only when imex3d_enabled(inputs) is true, so ",
                  "set :ode_solver => IMEX_ARK(:ARS343) in user_inputs.jl rather than ",
                  "swapping the integrator in later.")
        inputs[:ode_adaptive_solver] == true &&
            error("IMEX_ARK is fixed-step: its tableaux carry no embedded error ",
                  "estimate, and an adaptive controller would disagree between ranks ",
                  "about the step size -- which, since rhs! contains MPI collectives, ",
                  "deadlocks rather than fails. Set :ode_adaptive_solver => false.")
    end

    if get(inputs, :lcfl_report, false) == true
        _u    = copy(u)
        _qnm1 = copy(params.qp.qnm1);   _qnm2 = copy(params.qp.qnm2)
        _dq1  = copy(params.dsgs_qnm1); _dq2  = copy(params.dsgs_qnm2)
        _dth  = params.dsgs_thist[]
        rhs!(similar(u), u, params, inputs[:tinit])
        u .= _u
        params.qp.qnm1 .= _qnm1;   params.qp.qnm2 .= _qnm2
        params.dsgs_qnm1 .= _dq1;  params.dsgs_qnm2 .= _dq2
        params.dsgs_thist[] = _dth
        cfl_report(params, u, inputs[:tinit]; dt = inputs[:Δt])
    end

    function rad_condition(u, t, integrator)
        if (rem(t,rad_time) < 1e-3)
            return true
        else
            return false
        end
    end

    function do_radiation!(integrator)
        if (params.inputs[:RT_atmos_coupling])
            println(" doing full 3D RT solve + heat flux calculations at t=", integrator.t)
            get_RT_heat_fluxes!(params.uaux, params.qp.qe, params.mesh, params.mp, params.metrics, params.atmos_data, params, params.basis.dψ, params.basis.ψ, params.ω, PhysicalConst{TFloat}(), params.inputs)
        else
            println(" doing two stream radiation heat flux calculations at t=", integrator.t)
            compute_radiative_fluxes!(lnew_mesh, params.mesh, params.uaux, params.qp.qe, params.mp, params.phys_grid, params.inputs[:backend], params.SOL_VARS_TYPE)
        end
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
                         μ_dsgs_pnode = (integrator.p.VT == DSGS() || integrator.p.VT == DSGS_MHD()) ? integrator.p.μ_dsgs_pnode : nothing,
                         schlieren = maybe_compute_schlieren(inputs, integrator.p, integrator.u))
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

    cb_rad     = DiscreteCallback(rad_condition, do_radiation!)
    cb         = DiscreteCallback(condition, affect!)
    cb_amr     = DiscreteCallback(condition, affect!)
    cb_restart = DiscreteCallback(restart_condition, do_restart!)
    # Coupled-mode exchange callback: fires once per accepted timestep,
    # sends Julia's interpolated solution to Alya so its MPI.Waitall in
    # the time loop can advance. Without this Alya hangs and never
    # writes its VTS output.
    cb_coupling = is_coupled ? setup_coupling_callback(is_coupled, params, inputs) : nothing
    lrad = inputs[:RT_atmos_coupling] || inputs[:lphysics_grid]
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
        #
        # Every heartbeat carries the wall-clock cost of the interval it
        # closes, which turns the throttled trace into a performance meter:
        #
        #   #   step 200   t = 4.000000   wall 00:03:22   1.004 s/step   ETA 3d 00:52:16
        #
        #   wall    time since the FIRST step of this loop. Under the
        #           pre-compilation pass that first step is the pass itself,
        #           so `wall` measures hot code only - JIT and first-touch
        #           allocation are already excluded, which is what you want
        #           when comparing rank counts, mesh sizes or machines.
        #   s/step  measured over the steps since the PREVIOUS heartbeat, not
        #           averaged from t0, so a slowdown shows up in the line where
        #           it happens instead of being diluted by earlier history.
        #   ETA     the interval's simulated-seconds-per-wall-second carried
        #           out to :tend. A straight-line extrapolation of the last
        #           100 steps, nothing more.
        #
        # Rank 0 does the timing and the printing. Nothing here is collective:
        # a heartbeat must never be able to hang a 1500-rank job.
        _step_count   = Ref{Int}(0)
        _hb_wall0     = Ref{UInt64}(0)      # time_ns() at the first step
        _hb_wall_last = Ref{UInt64}(0)      # time_ns() at the previous heartbeat
        _hb_step_last = Ref{Int}(0)         # step index of the previous heartbeat
        _hb_tsim_last = Ref{Float64}(NaN)   # integrator.t at the previous heartbeat
        _hb_tend      = Float64(inputs[:tend])
        # The wall-clock anchors live and die with the step counter: every
        # site that rewinds `_step_count` (integrator warm-up, failed
        # pre-compilation pass) must rewind these too, or the production
        # loop's first heartbeat would bill it for the warm-up's seconds.
        function reset_heartbeat!()
            _step_count[]   = 0
            _hb_wall0[]     = 0
            _hb_wall_last[] = 0
            _hb_step_last[] = 0
            _hb_tsim_last[] = NaN
        end
        # hh:mm:ss, with a leading day count once past 24 h - these runs are
        # measured in days and "72:41:09" is harder to read than "3d 00:41:09".
        function _hb_hms(secs)
            (isfinite(secs) && secs >= 0) || return "--"
            s    = round(Int, secs)
            d, s = divrem(s, 86400)
            h, s = divrem(s, 3600)
            m, s = divrem(s, 60)
            return d > 0 ? @sprintf("%dd %02d:%02d:%02d", d, h, m, s) :
                           @sprintf("%02d:%02d:%02d", h, m, s)
        end
        function step_heartbeat_condition(u, t, integrator)
            _step_count[] += 1
            n = _step_count[]
            if n == 1
                # The clock starts here, at the end of step 1. Anchoring the
                # interval on the same instant makes step 1's own rate
                # unmeasurable (it prints `--`) rather than fake.
                now             = time_ns()
                _hb_wall0[]     = now
                _hb_wall_last[] = now
                _hb_step_last[] = n
                _hb_tsim_last[] = Float64(integrator.t)
            end
            return n <= 5 || (n % 100 == 0)
        end
        function step_heartbeat_affect!(integrator)
            now    = time_ns()
            n      = _step_count[]
            tsim   = Float64(integrator.t)
            dwall  = (now - _hb_wall_last[]) / 1e9
            dsteps = n - _hb_step_last[]
            if rank == 0
                rate = dsteps > 0 && dwall > 0 ?
                       @sprintf("%.3f s/step", dwall / dsteps) : "-- s/step"
                # Extrapolate on simulated time per wall second, so the
                # estimate stays honest under adaptive stepping, where a
                # step is not a fixed amount of physics.
                dtsim = tsim - _hb_tsim_last[]
                eta   = (dwall > 0 && isfinite(dtsim) && dtsim > 0) ?
                        _hb_hms((_hb_tend - tsim) * dwall / dtsim) : "--"
                @printf(" #   step %d   t = %.6f   wall %s   %s   ETA %s\n",
                        n, tsim, _hb_hms((now - _hb_wall0[]) / 1e9), rate, eta)
                flush(stdout)
            end
            _hb_wall_last[] = now
            _hb_step_last[] = n
            _hb_tsim_last[] = tsim
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

        #---------------------------------------------------------------------
        # PERIODIC CFL REPORT -- `:lcfl_report_every => N` (steps). Default off.
        #
        # `:lcfl_report` prints the stability table ONCE, at startup, and at
        # startup its viscous row can be actively misleading. ν_t is whatever
        # the SGS closure returns on the initial sounding, which for a laminar
        # profile is ~0, so `dt_viscous` comes back enormous and the report
        # says diffusion never limits Δt. That stops being true as soon as the
        # flow is turbulent -- on a convective boundary layer, hundreds of
        # seconds in -- and by then nothing prints it again.
        #
        # That gap is exactly where a blow-up at a fixed MODEL time comes from,
        # and it is invisible in the worst way: the term that killed the run is
        # the one the startup report called harmless. It matters more the more
        # implicit the scheme is, because removing the acoustic limit EXPOSES
        # whatever is next, and the diffusive rate goes as ν/h_z² -- which on a
        # vertically refined LES mesh is the term waiting underneath.
        #
        # UNLIKE THE HEARTBEAT ABOVE, THIS IS COLLECTIVE. `cfl_limits` does
        # about fifteen Allreduces, so every rank has to reach it: the
        # condition counts steps, every rank increments that counter on the
        # same steps, and nothing here branches on the rank. (`cfl_report`
        # itself returns early on non-root only AFTER the reduction, which is
        # the correct order and why it is safe to call from all ranks.)
        #
        # Cost is one pass over the local mesh plus those reductions -- no RHS
        # evaluation, and `cfl_limits` refreshes `uaux` from `u` itself, so it
        # needs nothing set up first and mutates nothing the integrator owns.
        # At a cadence of a few hundred steps it is free.
        #
        # The CFL numbers are quoted against the deck's nominal `:Δt` rather
        # than `integrator.dt`, so that a report landing on a tstop-shortened
        # step does not read as a sudden improvement. The point is to watch a
        # RATE move against a fixed line.
        #---------------------------------------------------------------------
        _cfl_count = Ref{Int}(0)
        _env_cfl   = strip(get(ENV, "JEXPRESSO_CFL_REPORT_EVERY", ""))
        _cfl_every = isempty(_env_cfl) ? Int(get(inputs, :lcfl_report_every, 0)) :
                                         parse(Int, _env_cfl)
        function cfl_periodic_condition(u, t, integrator)
            _cfl_count[] += 1
            return _cfl_count[] % _cfl_every == 0
        end
        function cfl_periodic_affect!(integrator)
            if rank == 0
                println()
                @printf(" # CFL report at t = %.4f  (step %d)\n",
                        Float64(integrator.t), _cfl_count[])
                flush(stdout)
            end
            cfl_report(params, integrator.u, integrator.t; dt = inputs[:Δt])
            return nothing
        end
        cb_cfl = _cfl_every > 0 ?
            DiscreteCallback(cfl_periodic_condition, cfl_periodic_affect!) :
            nothing

        _cbs = Any[cb, cb_restart, cb_les_stat, cb_les_online]
        lrad                                  && push!(_cbs, cb_rad)
        is_coupled && cb_coupling !== nothing  && push!(_cbs, cb_coupling)
        cb_heartbeat !== nothing               && push!(_cbs, cb_heartbeat)
        cb_cfl !== nothing                     && push!(_cbs, cb_cfl)
        callbacks_main = CallbackSet(_cbs...)

        # The `saveat` grid is built ONCE and shared by the pre-compilation
        # pass and the production solve. SciML specialises the integrator on
        # the *types* of the solve kwargs, so handing the two phases two
        # different `saveat` objects of two different types would make phase 2
        # recompile the very thing phase 1 just compiled.
        saveat_main = range(inputs[:tinit], inputs[:tend],
                            length = inputs[:ndiagnostics_outputs])

        # Two-phase run: 1 compile step, then the production solve resuming
        # from it. See precompile_pass_enabled() for the full rationale.
        _precompile_pass = precompile_pass_enabled(inputs)

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
        #
        # Skipped entirely under the pre-compilation pass: that pass compiles
        # the identical prob/callbacks_main specialisation and, unlike this
        # warm-up, keeps the step instead of restoring the IC.
        if !_precompile_pass && precompile_warmup_enabled(inputs)
            rank == 0 && (print(YELLOW_FG(" # Integrator warm-up with real callbacks (PATIENCE: ONLY DONE ON 1st RUN!) ......... ")); flush(stdout))
            _t_wm = time_ns()
            u_snap    = copy(u)
            qnm1_snap = copy(params.qp.qnm1)
            qnm2_snap = copy(params.qp.qnm2)
            dsgs_qnm1_snap = copy(params.dsgs_qnm1)
            dsgs_qnm2_snap = copy(params.dsgs_qnm2)
            dsgs_thist_snap = params.dsgs_thist[]
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
            params.dsgs_qnm1 .= dsgs_qnm1_snap
            params.dsgs_qnm2 .= dsgs_qnm2_snap
            params.dsgs_thist[] = dsgs_thist_snap
            inputs isa Dict && saved_outdir !== nothing && (inputs[:output_dir] = saved_outdir)
            try; rm(warm_outdir; recursive = true, force = true); catch; end
            # Reset the heartbeat counter and its wall clock so the real
            # solve gets its own first-5-steps detail and its own timing
            # baseline (the warmup just consumed one step's worth of both).
            reset_heartbeat!()
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

        # The production solve, as ONE call. `solve(prob, alg; kw...)` is by
        # definition `solve!(init(prob, alg; kw...))`, which is what the
        # two-phase branch below takes apart.
        run_single_phase() = with_logger(solve_logger) do
            solve(prob,
                  inputs[:ode_solver], dt=dt,
                  #callback = CallbackSet(cb,cb_rad), tstops = dosetimes,
                  callback = callbacks_main, tstops = tstops_all,
                  save_everystep = false,
                  adaptive=inputs[:ode_adaptive_solver],
                  saveat = saveat_main)
        end

        if !_precompile_pass
            solution = run_single_phase()
        else
            #----------------------------------------------------------------
            # PRE-COMPILATION PASS (see precompile_pass_enabled).
            #
            # `solve` is split into its two halves along the SciML integrator
            # interface:
            #
            #     integrator = init(prob, alg; kwargs...)
            #     step!(integrator)      <-- PHASE 1: exactly one timestep
            #     GC.gc(); MPI.Barrier()
            #     solve!(integrator)     <-- PHASE 2: the rest of the run
            #
            # Phase 1 is where the whole cost of "first timestep on a large
            # problem" is paid: `rhs!` and the entire kernel chain, the
            # integrator specialised on the REAL CallbackSet, the MPI halo
            # exchange, the diagnostic/output path, and the first touch of
            # every per-rank working array. Doing it as its own call means
            # that cost is measured, reported, and finished BEFORE the
            # production loop starts — instead of happening inside it, where
            # every rank that compiled quickly sits in a collective waiting
            # for one that did not.
            #
            # Phase 2 then continues the SAME integrator object. This is not a
            # restart: there is no second `init`, the step-size controller,
            # the callback caches, the FSAL history and `u` are carried
            # across untouched. The trajectory is therefore bit-for-bit the
            # one a single `solve(...)` would have produced — including under
            # adaptive stepping — and not one method is compiled twice.
            #
            # Phase 1's step IS the simulation's first step. Nothing is
            # snapshotted, restored or thrown away (contrast the warm-ups
            # above, which pay for a full-mesh RHS evaluation and then discard
            # its result).
            #----------------------------------------------------------------
            t0_pass = params.tspan[1]
            if rank == 0
                println(YELLOW_FG(" # ┌─ Pre-compilation pass: 1 timestep on the production problem"))
                println(YELLOW_FG(" # │  JIT + first-touch allocations happen HERE, not inside the time loop."))
                @printf(" # │  t = %.6f → %.6f  (dt = %.6g).  PATIENCE: this is the slow step.\n",
                        t0_pass, t0_pass + dt, dt)
                flush(stdout)
            end

            # Restored only if the pass fails: a partially-advanced history
            # must never leak into the single-phase fallback. `u` itself needs
            # no snapshot — the integrator works on `recursivecopy(prob.u0)`,
            # so `u` is untouched until phase 2 completes.
            qnm1_snap       = copy(params.qp.qnm1)
            qnm2_snap       = copy(params.qp.qnm2)
            dsgs_qnm1_snap  = copy(params.dsgs_qnm1)
            dsgs_qnm2_snap  = copy(params.dsgs_qnm2)
            dsgs_thist_snap = params.dsgs_thist[]

            _t_pass = time_ns()
            integrator = with_logger(solve_logger) do
                try
                    integ = init(prob,
                                 inputs[:ode_solver], dt=dt,
                                 callback = callbacks_main, tstops = tstops_all,
                                 save_everystep = false,
                                 adaptive=inputs[:ode_adaptive_solver],
                                 saveat = saveat_main)
                    # One step. `step!` honours tstops, so if a diagnostic time
                    # falls inside [t0, t0+Δt] the step is shortened to land on
                    # it and its callback fires for real — this is production
                    # time stepping, not a rehearsal.
                    step!(integ)
                    integ
                catch e
                    rank == 0 && @warn "pre-compilation pass failed; falling back to a single-phase run" exception=e
                    nothing
                end
            end

            _t_pass_s   = (time_ns() - _t_pass) / 1e9
            _t_pass_max = MPI.Allreduce(_t_pass_s, MPI.MAX, comm)
            _t_pass_min = MPI.Allreduce(_t_pass_s, MPI.MIN, comm)

            # The fallback decision MUST be collective. A `catch` is per-rank,
            # so if one rank's init/step! threw and the others' did not, half
            # the job would enter `solve!` and half would enter a fresh
            # `solve` — two different sequences of collectives, i.e. a hang
            # with no error message. One rank failing therefore sends everyone
            # down the single-phase path.
            _pass_ok = MPI.Allreduce(integrator === nothing ? 0 : 1,
                                     MPI.MIN, comm) == 1

            # Hand phase 2 a compacted heap: the compiler's garbage is dead by
            # now, and a 1500-rank run should not drag it through hours of time
            # stepping. The barrier makes every rank enter the production loop
            # together — otherwise the ranks that compiled fastest immediately
            # block in the first halo exchange of the real loop.
            GC.gc(true)
            MPI.Barrier(comm)

            if rank == 0
                @printf(" # │  DONE in %.2f s   (slowest rank %.2f s, fastest %.2f s)\n",
                        _t_pass_s, _t_pass_max, _t_pass_min)
                if !_pass_ok
                    println(" # └─ Pass unusable on at least one rank; running single-phase from t0, JIT-ing as usual.")
                else
                    @printf(" # └─ Hot code. Simulation continues from t = %.6f to t = %.6f\n",
                            integrator.t, params.tspan[2])
                end
                flush(stdout)
            end

            # The pass absorbed the compile-time allocations, so the summary
            # table below reflects steady state only.
            if alloc_summary_enabled(inputs)
                TimerOutputs.reset_timer!(JEXPRESSO_TIMER)
            end

            if !_pass_ok
                # Roll back the history this rank's step advanced (if it got
                # that far) so the single-phase run below starts from the same
                # initial condition it would have had without the pass. `u`
                # needs nothing: the integrator only ever touched its own copy.
                integrator = nothing
                params.qp.qnm1     .= qnm1_snap
                params.qp.qnm2     .= qnm2_snap
                params.dsgs_qnm1   .= dsgs_qnm1_snap
                params.dsgs_qnm2   .= dsgs_qnm2_snap
                params.dsgs_thist[] = dsgs_thist_snap
                reset_heartbeat!()
                solution = run_single_phase()
            else
                solution = with_logger(solve_logger) do
                    solve!(integrator)
                end
            end
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

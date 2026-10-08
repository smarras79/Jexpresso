"""
Accumulates shallow-water clipping diagnostics between accepted timesteps.

The counters are deliberately stored in `params`, rather than as global state,
so separate simulations and MPI ranks cannot interfere with one another.  MPI
reductions are deferred to the accepted-step diagnostic callback.
"""
mutable struct ShallowWaterLimiterDiagnostics
    accepted_steps::Int
    stage_events::Int
    corrected_nodes::Int
    min_h_before::Float64
    max_abs_hu_before::Float64
    max_abs_hv_before::Float64
    mass_added::Float64
    first_time::Float64
    last_time::Float64
    first_owned_time::Float64
    first_node::Int
    first_x::Float64
    first_y::Float64
    first_h::Float64
    first_hu::Float64
    first_hv::Float64
end

ShallowWaterLimiterDiagnostics() = ShallowWaterLimiterDiagnostics(
    0, 0, 0, Inf, 0.0, 0.0, 0.0, Inf, -Inf,
    Inf, 0, NaN, NaN, NaN, NaN, NaN)

function reset_shallow_water_limiter_diagnostics!(diag;
                                                   reset_step_counter=false)
    reset_step_counter && (diag.accepted_steps = 0)
    diag.stage_events = 0
    diag.corrected_nodes = 0
    diag.min_h_before = Inf
    diag.max_abs_hu_before = 0.0
    diag.max_abs_hv_before = 0.0
    diag.mass_added = 0.0
    diag.first_time = Inf
    diag.last_time = -Inf
    diag.first_owned_time = Inf
    diag.first_node = 0
    diag.first_x = NaN
    diag.first_y = NaN
    diag.first_h = NaN
    diag.first_hu = NaN
    diag.first_hv = NaN
    return nothing
end

@inline function shallow_water_limiter_diagnostic_enabled(inputs)
    return get(inputs, :swe_limiter_diagnostic, false) == true
end

@inline function _shallow_water_mass_weight(M, ip)
    ndims(M) == 1 && return Float64(M[ip])
    # For a consistent mass matrix, integrating a nodal basis function gives
    # its row sum.  The common lumped-mass path above is allocation-free.
    return Float64(sum(@view M[ip, :]))
end

function _record_shallow_water_limiter_event!(diag, H, Hu, Hv, M, x, y,
                                               gip2owner, rank, hmin, time)
    event = false
    @inbounds for ip in eachindex(H)
        h = H[ip]
        if h < hmin
            event = true
            diag.min_h_before = min(diag.min_h_before, Float64(h))
            diag.max_abs_hu_before = max(diag.max_abs_hu_before,
                                         abs(Float64(Hu[ip])))
            diag.max_abs_hv_before = max(diag.max_abs_hv_before,
                                         abs(Float64(Hv[ip])))
            # Shared interface nodes occur on more than one rank. Count and
            # integrate only the owning copy to avoid double counting.
            if gip2owner[ip] == rank
                diag.corrected_nodes += 1
                diag.mass_added += _shallow_water_mass_weight(M, ip) *
                                   Float64(hmin - h)
                if diag.first_node == 0
                    diag.first_owned_time = Float64(time)
                    diag.first_node = Int(ip)
                    diag.first_x = Float64(x[ip])
                    diag.first_y = Float64(y[ip])
                    diag.first_h = Float64(h)
                    diag.first_hu = Float64(Hu[ip])
                    diag.first_hv = Float64(Hv[ip])
                end
            end
        end
    end
    if event
        diag.stage_events += 1
        t = Float64(time)
        diag.first_time = min(diag.first_time, t)
        diag.last_time = max(diag.last_time, t)
    end
    return nothing
end

"""
    apply_shallow_water_limiter!(u, params, time)

Apply the configured shallow-water limiter directly to the prognostic state.
The state uses Jexpresso's field-major layout `[H, Hu, Hv]`, with `npoin`
entries per field.

The initial implementation is deliberately simple: depths below `:swe_hmin`
are replaced by that floor and their two momentum components are set to zero.
This is not mass conservative.  Keeping the operation behind this single
interface lets a conservative positivity limiter replace it later without
changing any time-integrator wiring.
"""
function apply_shallow_water_limiter!(u, params, time)
    npoin = Int(params.mesh.npoin)
    params.neqs >= 3 ||
        error("The shallow-water limiter requires at least the fields [H, Hu, Hv].")
    length(u) >= 3npoin ||
        error("The shallow-water state is too short for [H, Hu, Hv].")

    hmin = convert(eltype(u), params.inputs[:swe_hmin])
    hmin >= zero(hmin) || error(":swe_hmin must be nonnegative.")

    # Change momentum first: after H is raised to hmin the dry-node mask is no
    # longer available. Broadcasted assignments work for both host and device
    # arrays and avoid allocating a Boolean mask.
    @views begin
        H  = u[1:npoin]
        Hu = u[npoin + 1:2npoin]
        Hv = u[2npoin + 1:3npoin]

        if shallow_water_limiter_diagnostic_enabled(params.inputs)
            rank = MPI.Comm_rank(get_mpi_comm())
            if params.backend isa CPU
                _record_shallow_water_limiter_event!(
                    params.swe_limiter_diag, H, Hu, Hv, params.M,
                    params.mesh.x, params.mesh.y,
                    params.mesh.gip2owner, rank, hmin, time)
            elseif minimum(H) < hmin
                # Diagnostics are intentionally opt-in.  Copy only when a
                # device-side reduction shows that clipping is actually
                # needed; the normal GPU limiter path below stays on-device.
                _record_shallow_water_limiter_event!(
                    params.swe_limiter_diag, Array(H), Array(Hu), Array(Hv),
                    Array(params.M), Array(params.mesh.x), Array(params.mesh.y),
                    Array(params.mesh.gip2owner), rank, hmin, time)
            end
        end

        Hu .= ifelse.(H .< hmin, zero(eltype(u)), Hu)
        Hv .= ifelse.(H .< hmin, zero(eltype(u)), Hv)
        H  .= max.(H, hmin)
    end

    return nothing
end


function shallow_water_limiter_diagnostic_condition(u, time, integrator)
    diag = integrator.p.swe_limiter_diag
    diag.accepted_steps += 1
    interval = max(1, Int(get(integrator.p.inputs,
                              :swe_limiter_diagnostic_interval, 1)))
    return diag.accepted_steps % interval == 0
end

function report_shallow_water_limiter_diagnostics!(params, time)
    diag = params.swe_limiter_diag
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    # MAX avoids multiplying the number of stage calls by the number of MPI
    # ranks. It is an event-count estimate if different ranks first encounter
    # clipping in different stages.
    global_events = MPI.Allreduce(diag.stage_events, MPI.MAX, comm)
    if global_events > 0
        global_corrections = MPI.Allreduce(diag.corrected_nodes, MPI.SUM, comm)
        global_min_h = MPI.Allreduce(diag.min_h_before, MPI.MIN, comm)
        global_max_hu = MPI.Allreduce(diag.max_abs_hu_before, MPI.MAX, comm)
        global_max_hv = MPI.Allreduce(diag.max_abs_hv_before, MPI.MAX, comm)
        global_mass_added = MPI.Allreduce(diag.mass_added, MPI.SUM, comm)
        global_first_time = MPI.Allreduce(diag.first_time, MPI.MIN, comm)
        global_last_time = MPI.Allreduce(diag.last_time, MPI.MAX, comm)

        # Select one deterministic, owned node from the earliest clipping
        # stage in this reporting interval. If several ranks clip at exactly
        # the same stage time, the lowest rank wins.
        global_first_owned_time = MPI.Allreduce(diag.first_owned_time,
                                                MPI.MIN, comm)
        rank_candidate = diag.first_owned_time == global_first_owned_time ?
                         rank : typemax(Int)
        first_rank = MPI.Allreduce(rank_candidate, MPI.MIN, comm)
        selected = rank == first_rank
        first_node = MPI.Allreduce(selected ? diag.first_node : 0,
                                   MPI.SUM, comm)
        first_x = MPI.Allreduce(selected ? diag.first_x : 0.0, MPI.SUM, comm)
        first_y = MPI.Allreduce(selected ? diag.first_y : 0.0, MPI.SUM, comm)
        first_h = MPI.Allreduce(selected ? diag.first_h : 0.0, MPI.SUM, comm)
        first_hu = MPI.Allreduce(selected ? diag.first_hu : 0.0, MPI.SUM, comm)
        first_hv = MPI.Allreduce(selected ? diag.first_hv : 0.0, MPI.SUM, comm)

        if rank == 0
            @printf(" # SWE limiter diagnostics: step=%d t=%.6e events=%d corrections=%d min(H)_before=%.16e max|Hu|=%.6e max|Hv|=%.6e mass_added=%.6e event_times=[%.6e, %.6e]\n",
                    diag.accepted_steps, time, global_events,
                    global_corrections, global_min_h, global_max_hu,
                    global_max_hv, global_mass_added, global_first_time,
                    global_last_time)
            @printf(" #   first owned clip: stage_t=%.6e rank=%d node=%d x=%.6e y=%.6e H=%.16e Hu=%.6e Hv=%.6e\n",
                    global_first_owned_time, first_rank, first_node,
                    first_x, first_y, first_h, first_hu, first_hv)
            flush(stdout)
        end
    end

    reset_shallow_water_limiter_diagnostics!(diag)
    return nothing
end

function shallow_water_limiter_diagnostic_affect!(integrator)
    report_shallow_water_limiter_diagnostics!(integrator.p, integrator.t)
    return nothing
end


# OrdinaryDiffEq limiter adapters.  The numerical operation above intentionally
# has no dependency on an integrator type.
function shallow_water_stage_limiter!(u, integrator, params, time)
    apply_shallow_water_limiter!(u, params, time)
    return nothing
end

function shallow_water_step_limiter_affect!(integrator)
    apply_shallow_water_limiter!(integrator.u, integrator.p, integrator.t)
    # A callback is outside the algorithm's normal stage/cache update. Tell
    # SciML that FSAL values and other cached derivatives must be refreshed.
    u_modified!(integrator, true)
    return nothing
end

shallow_water_step_limiter_condition(u, time, integrator) = true


# Compose rather than discard a limiter already supplied by the user.
struct ShallowWaterStageLimiter{F}
    inner::F
end

function (limiter::ShallowWaterStageLimiter)(u, integrator, params, time)
    limiter.inner(u, integrator, params, time)
    shallow_water_stage_limiter!(u, integrator, params, time)
    return nothing
end


@inline function shallow_water_limiter_strategy(inputs)
    return lowercase(string(get(inputs, :swe_limiter, "none")))
end

@inline function shallow_water_limiter_schedule(inputs)
    return lowercase(string(get(inputs, :swe_limiter_schedule, "none")))
end

@inline function shallow_water_limiter_enabled(inputs)
    return shallow_water_limiter_strategy(inputs) != "none" &&
           shallow_water_limiter_schedule(inputs) != "none"
end


"""
    configure_shallow_water_limiter(alg, inputs)

Return `(configured_alg, stage_active, step_callback_active)`.  Algorithms
with a `stage_limiter!` field receive the shallow-water stage adapter while
preserving their other settings and composing any existing stage limiter.
Every algorithm can use the accepted-step callback fallback.
"""
function configure_shallow_water_limiter(alg, inputs)
    shallow_water_limiter_enabled(inputs) || return alg, false, false

    strategy = shallow_water_limiter_strategy(inputs)
    strategy == "clip" ||
        error("Unsupported :swe_limiter = $(repr(strategy)); supported values are \"none\" and \"clip\".")

    schedule = shallow_water_limiter_schedule(inputs)
    schedule in ("stage", "step", "stage_and_step", "auto") ||
        error("Unsupported :swe_limiter_schedule = $(repr(schedule)).")

    wants_stage = schedule in ("stage", "stage_and_step", "auto")
    stage_active = false
    configured_alg = alg

    if wants_stage && hasfield(typeof(alg), :stage_limiter!)
        names = fieldnames(typeof(alg))
        values = ntuple(length(names)) do i
            name = names[i]
            value = getfield(alg, i)
            name == :stage_limiter! ? ShallowWaterStageLimiter(value) : value
        end
        kwargs = NamedTuple{names}(values)
        Alg = Base.typename(typeof(alg)).wrapper
        try
            # OrdinaryDiffEq algorithm families do not share a positional
            # constructor order (Rosenbrock methods differ from explicit RK,
            # for example), but their keyword constructors use the stored
            # field names. Rebuilding from those names preserves every option.
            configured_alg = Alg(; kwargs...)
            stage_active = true
        catch err
            schedule == "stage" && rethrow(err)
            @warn "Could not install the shallow-water stage limiter; using accepted-step limiting" algorithm=typeof(alg) exception=err
        end
    elseif wants_stage && schedule == "stage"
        error("$(nameof(typeof(alg))) has no stage_limiter! hook. Use :swe_limiter_schedule => \"step\" or \"auto\".")
    elseif wants_stage && schedule == "stage_and_step"
        @warn "The selected time integrator has no stage_limiter! hook; shallow-water limiting will run after accepted steps only" algorithm=typeof(alg)
    end

    step_callback_active = schedule in ("step", "stage_and_step") ||
                           (schedule == "auto" && !stage_active)
    return configured_alg, stage_active, step_callback_active
end

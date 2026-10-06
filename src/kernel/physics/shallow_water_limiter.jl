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

        Hu .= ifelse.(H .< hmin, zero(eltype(u)), Hu)
        Hv .= ifelse.(H .< hmin, zero(eltype(u)), Hv)
        H  .= max.(H, hmin)
    end

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

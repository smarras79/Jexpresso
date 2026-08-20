#---------------------------------------------------------------------------------
# test/test_precompile_pass.jl — the two-phase pre-compilation pass.
#
#   julia --project=. test/test_precompile_pass.jl
#
# `time_loop!` can run a large case as ONE compile timestep followed by the
# production loop, instead of a single `solve(...)` that JITs inside its own
# hot loop (ENVIRONMENT_VARIABLES.md -> JEXPRESSO_PRECOMPILE_PASS). It does
# that by splitting `solve` along the SciML integrator interface:
#
#     integrator = init(prob, alg; kwargs...)
#     step!(integrator)      # phase 1
#     solve!(integrator)     # phase 2
#
# The whole design rests on ONE claim: that this is not a restart, and produces
# exactly what a single `solve(...)` would have produced. Part 2 below pins that
# claim down — bitwise — for both fixed and adaptive stepping, on a toy in-place
# ODE with the same solve kwargs shape `time_loop!` uses (a CallbackSet, a
# tstops vector, `save_everystep = false` and a `saveat` range).
#
# A note on why this test exists at all: an earlier version of the pass ran
# phase 1 as its own `solve(...)` over [t0, t0+Δt] and read the advanced state
# out of the returned solution. That silently did nothing, because SciML's
# `save_end` defaults to `prob.tspan[2] in saveat` — false for a one-step
# tspan — so the returned solution held only the INITIAL state.
#
# Nothing here touches the filesystem, a mesh, or an MPI collective.
#---------------------------------------------------------------------------------
using Test
using MPI
using Jexpresso

using Jexpresso: precompile_pass_enabled
using Jexpresso: ODEProblem, FullSpecialize, CallbackSet, DiscreteCallback,
                 init, step!, solve!, solve
using OrdinaryDiffEq: CarpenterKennedy2N54, Tsit5

MPI.Initialized() || MPI.Init()

#-------------------------------------------------------------------------------
# 1. The switch: env var > CLI flag > user_inputs.jl > default(false)
#-------------------------------------------------------------------------------
@testset "precompile_pass_enabled precedence" begin
    saved_env  = get(ENV, "JEXPRESSO_PRECOMPILE_PASS", nothing)
    saved_args = copy(ARGS)
    try
        delete!(ENV, "JEXPRESSO_PRECOMPILE_PASS")
        empty!(ARGS)

        # 4. default
        @test precompile_pass_enabled(Dict{Symbol,Any}()) == false

        # 3. user_inputs.jl — Dict and NamedTuple decks alike
        @test precompile_pass_enabled(Dict(:lprecompile_pass => true))  == true
        @test precompile_pass_enabled(Dict(:lprecompile_pass => false)) == false
        @test precompile_pass_enabled((; lprecompile_pass = true))      == true

        # 2. CLI flag beats the deck
        push!(ARGS, "--precompile-pass")
        @test precompile_pass_enabled(Dict(:lprecompile_pass => false)) == true
        empty!(ARGS); push!(ARGS, "--no-precompile-pass")
        @test precompile_pass_enabled(Dict(:lprecompile_pass => true))  == false
        empty!(ARGS)

        # 1. env beats both
        for v in ("1", "true", "YES", "on")
            ENV["JEXPRESSO_PRECOMPILE_PASS"] = v
            @test precompile_pass_enabled(Dict(:lprecompile_pass => false)) == true
        end
        for v in ("0", "false", "No", "off")
            ENV["JEXPRESSO_PRECOMPILE_PASS"] = v
            @test precompile_pass_enabled(Dict(:lprecompile_pass => true))  == false
        end
        # Garbage falls through to the next level rather than being truthy.
        ENV["JEXPRESSO_PRECOMPILE_PASS"] = "maybe"
        @test precompile_pass_enabled(Dict(:lprecompile_pass => true))  == true
        @test precompile_pass_enabled(Dict{Symbol,Any}())               == false
    finally
        delete!(ENV, "JEXPRESSO_PRECOMPILE_PASS")
        saved_env === nothing || (ENV["JEXPRESSO_PRECOMPILE_PASS"] = saved_env)
        append!(empty!(ARGS), saved_args)
    end
end

#-------------------------------------------------------------------------------
# 2. init + step! + solve!  ==  solve, bit for bit.
#-------------------------------------------------------------------------------
const NEQ = 200

function _pass_rhs!(du, u, p, t)
    @inbounds for i in eachindex(u)
        du[i] = -0.3 * u[i] + 0.1 * sin(t + i) * u[mod1(i + 1, length(u))]
    end
    p.nrhs[] += 1
    return nothing
end

function _run_toy(alg, adaptive::Bool, twophase::Bool, hits::Vector{Float64})
    params = (nrhs = Ref(0),)
    u0     = [sin(i / 7) + 2.0 for i in 1:NEQ]
    u      = copy(u0)

    tinit, tend, dtf = 0.0, 2.0, Float32(0.02)
    dosetimes   = [0.5, 1.0, 2.0]
    tstops_all  = sort(unique(vcat(dosetimes, [1.5])))
    saveat_main = range(tinit, tend, length = 5)

    empty!(hits)
    cb  = DiscreteCallback((u, t, integ) -> any(x -> x == t, dosetimes),
                           integ -> push!(hits, integ.t))
    cb2 = DiscreteCallback((u, t, integ) -> false, integ -> nothing)
    cbs = CallbackSet(cb, cb2)

    # `tspan` as a Vector — that is what build_tspan hands params_setup, and
    # a Vector-vs-Tuple mismatch is enough to force a full recompile.
    prob = ODEProblem{true, FullSpecialize}(_pass_rhs!, u, [tinit, tend], params)
    kw   = (dt = dtf, callback = cbs, tstops = tstops_all,
            save_everystep = false, adaptive = adaptive, saveat = saveat_main)

    t_phase1 = NaN
    sol = if twophase
        integ = init(prob, alg; kw...)
        step!(integ)                      # PHASE 1 — the one compile step
        t_phase1 = integ.t
        solve!(integ)                     # PHASE 2 — resumes the same integrator
    else
        solve(prob, alg; kw...)
    end

    return (u_end = copy(sol.u[end]), t = copy(sol.t), hits = copy(hits),
            nrhs = params.nrhs[], retcode = sol.retcode,
            u0_untouched = (u == u0), t_phase1 = t_phase1, dt = Float64(dtf))
end

@testset "two-phase run == single-phase run" begin
    hits = Float64[]
    for (name, alg, adaptive) in (("fixed dt",    CarpenterKennedy2N54(), false),
                                  ("adaptive dt", Tsit5(),               true))
        @testset "$name" begin
            one = _run_toy(alg, adaptive, false, hits)
            two = _run_toy(alg, adaptive, true,  hits)

            # Phase 1 really did advance exactly one step.
            @test two.t_phase1 ≈ two.dt

            # …and phase 2 landed on the identical trajectory.
            @test two.u_end == one.u_end          # BITWISE, not ≈
            @test two.t     == one.t              # same saved output times
            @test two.hits  == one.hits           # same callback firings
            @test two.retcode == one.retcode

            # No extra work: the split costs zero additional RHS evaluations.
            @test two.nrhs == one.nrhs

            # The integrator works on a copy of u0, so a failed pass can always
            # fall back to a clean single-phase run (time_loop! relies on this).
            @test one.u0_untouched
            @test two.u0_untouched
        end
    end
end

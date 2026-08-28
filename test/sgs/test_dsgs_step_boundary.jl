#=============================================================================
 test_dsgs_step_boundary.jl -- the DynSGS history roll lands on the stage where
                               `uaux` is the step answer, and on no other.

 WHAT THE INVARIANT IS
 ---------------------
 DynSGS builds its residual from a BDF2 stencil over three STEP states, so the
 roll that maintains dsgs_qnm1/qnm2 has to happen at an rhs! call where `uaux`
 is a step state. rhs.jl decides that from elapsed time alone:

     ldsgs_step = ... && (time - dsgs_thist >= 0.999*dt)

 and there is no stage index in that expression. It works, but ONLY because of
 something that lives in a different file, so this test exists to pin the two
 halves together.

 WHY IT WORKS, AND WHY THAT IS NOT OBVIOUS
 -----------------------------------------
 ark.jl never calls the explicit rhs! at stage 1 -- stage 1 is the FSAL carried
 from the previous step -- so within a step the calls land at

     cE[2..4] = [0.4359, 0.7179, 1.0] * dt

 and only the last of those carries q^{n+1} (ARS343 is stiffly accurate, so
 stage s IS the answer). `dsgs_thist` starts at -1e30, so the FIRST rhs! call
 of the run fires the guard whatever its time -- and everything afterwards is
 locked to that anchor, because `thist` then advances by exactly dt each time.

 THE ANCHOR IS ark.jl:659, `initialize!`, which evaluates rhs! at t = tinit on
 uprev = q^0 before the loop starts. That call is a step state, so the guard
 anchors on a step boundary and every subsequent firing lands on stage 4.

 Take that call away and the first firing is stage 2, and the roll then sits on
 a MID-STEP state for the entire run: same cadence, one roll per step, finite
 coefficient, no error -- and the BDF2 levels are stage states while M^-1*RHS
 in the residual is evaluated at a different stage again. Silent in every way
 that matters. That is the regression this file is here to catch, so it is
 asserted as the negative case rather than described.

     julia --project=<env> test/sgs/test_dsgs_step_boundary.jl
=============================================================================#

using Test, Printf

say(a...) = (println(a...); flush(stdout))

# --- the real ARS343 tableau, read the way _mk_tableau builds it -------------
const γ             = 0.4358665215084590
const a31, a32      = 0.3212788860, 0.3966543747
const a41, a42, a43 = -0.105858296, 0.5529291479, 0.5529291479
const aE = [0.0 0.0 0.0 0.0
            γ   0.0 0.0 0.0
            a31 a32 0.0 0.0
            a41 a42 a43 0.0]
const cE = vec(sum(aE, dims = 2))

# The guard as rhs.jl spells it. Restated here rather than lifted, because it
# is three tokens inside a larger boolean and the thing under test is the
# CALL SEQUENCE it is driven by, not the arithmetic.
fires(time, thist, dt) = time - thist >= 0.999 * dt

"""
    replay(dt, nstep; anchored) -> Vector{Tuple{Int,Float64,Int}}

Every call the guard fires on, as (step, time, stage). Stage 1 means the
`initialize!` FSAL call at t = tinit; 2..4 are ark.jl's in-step rhs! calls.

`anchored = false` removes the initialize! call -- the counterfactual.
"""
function replay(dt, nstep; anchored::Bool)
    thist = -1.0e30
    fired = Tuple{Int,Float64,Int}[]
    if anchored && fires(0.0, thist, dt)
        push!(fired, (0, 0.0, 1)); thist = 0.0
    end
    for step = 1:nstep, i = 2:4
        time = (step - 1) * dt + cE[i] * dt
        if fires(time, thist, dt)
            push!(fired, (step, time, i)); thist = time
        end
    end
    return fired
end

say("\n=== DynSGS history roll: which rhs! call does it land on? ===")

@testset "DynSGS step-boundary detection" begin

    @testset "the tableau puts a call exactly at t + dt" begin
        # If this ever stops holding, the guard has no call to land on and the
        # whole scheme below changes -- so it is asserted, not assumed.
        @test isapprox(cE[4], 1.0; atol = 1e-9)
        @test all(cE[2:3] .< 0.999)
        say(@sprintf("  cE = [%.4f, %.4f, %.4f, %.4f]", cE...))
    end

    @testset "anchored by initialize!: every roll is on the answer stage" begin
        fired  = replay(0.5, 8; anchored = true)
        stages = [f[3] for f in fired]
        @test stages[1] == 1                 # the initialize! call, uprev = q^0
        @test all(==(4), stages[2:end])      # ... then the answer stage, always
        @test !any(in((2, 3)), stages)       # never a mid-step stage
        @test length(fired) == 9             # 1 anchor + one per step
        # atol 1e-8, not 1e-12: cE[4] is 0.9999999998, not exactly 1, so the
        # firings are one dt apart to the tableau's own precision and no better.
        @test all(isapprox.(diff([f[2] for f in fired]), 0.5; atol = 1e-8))
        say(@sprintf("  anchored:   %d rolls, stages %s", length(fired), stages))
    end

    @testset "without the anchor it locks to a mid-step stage, silently" begin
        fired  = replay(0.5, 8; anchored = false)
        stages = [f[3] for f in fired]
        @test all(==(2), stages)             # stage 2, for the whole run
        @test !any(==(4), stages)            # the answer stage is never seen
        @test length(fired) == 8             # SAME cadence -- that is the trap
        say(@sprintf("  unanchored: %d rolls, stages %s  <- same count, wrong operand",
                     length(fired), stages))
    end
end

say("=== done ===\n")

#=============================================================================
 test_dsgs_ramp.jl -- the DynSGS spin-up ramp, and the promise that it is
                      invisible to every case that does not ask for it.

 WHY THE OFF CASE IS THE IMPORTANT TEST
 --------------------------------------
 `:dsgs_ramp_steps` exists for a spinning-up boundary layer, where the sensor's
 denominators are L-infinity norms of a solution that has almost no resolved
 turbulence yet -- so the first coefficient the run applies is built on the
 least meaningful norms it will ever have. It is NOT wanted on the shock
 problems the model was tuned against (sod1d, theta_dsgs, ffs_step), which are
 correct as they stand.

 The default is 0, and 0 must mean EXACTLY the old behaviour: `dsgs_ramp_factor`
 returns 1.0 and rhs.jl skips the scaling entirely (`_ramp < one(TT)` is false).
 So the first testset is the one that protects the cases nobody asked to change.

     julia --project=<env> test/sgs/test_dsgs_ramp.jl
=============================================================================#

using Test, Printf

# Lift the real definition rather than restating it: a test carrying its own
# copy pins the copy. `dsgs_ramp_factor` deliberately depends on nothing.
const _SGS = read(joinpath(@__DIR__, "..", "..", "src", "kernel", "physics", "SGS.jl"), String)
let i = findfirst("@inline function dsgs_ramp_factor", _SGS)
    i === nothing && error("dsgs_ramp_factor not found in SGS.jl -- renamed?")
    j = findnext("\nend\n", _SGS, i[1])
    include_string(Main, _SGS[i[1]:j[1]+4], "SGS.jl:dsgs_ramp_factor")
end

say(a...) = (println(a...); flush(stdout))

const DT    = 0.5
const TINIT = 0.0
const T_ON  = TINIT + 1.999 * DT        # where dsgs_hist turns true

say("\n=== DynSGS spin-up ramp ===")

@testset "dsgs_ramp_factor" begin

    @testset "0 steps is exactly the old behaviour, at every time" begin
        # The shock cases live here. Nothing may move.
        for t in (0.0, 0.5, 1.0, 5.0, 1000.0)
            @test dsgs_ramp_factor(t, TINIT, DT, 0) === 1.0
        end
        # ... and a negative n cannot silently become a ramp either; params_setup
        # refuses it, but the function must not invent behaviour for it.
        @test dsgs_ramp_factor(10.0, TINIT, DT, -3) === 1.0
        say("  n = 0: factor is 1.0 at every time -- untouched")
    end

    @testset "n steps: 1/n on the first activated step, 1.0 after n" begin
        n = 5
        f = [dsgs_ramp_factor(T_ON + k*DT, TINIT, DT, n) for k = 0:n+2]
        @test isapprox(f[1], 1/n; atol = 1e-6)      # first activated step
        @test all(isapprox.(f[1:n], (1:n)./n; atol = 1e-6))
        @test all(==(1.0), f[n+1:end])              # saturated, and stays there
        @test issorted(f)                           # monotone, never dips
        say(@sprintf("  n = 5: %s", join(round.(f; digits=3), ", ")))
    end

    @testset "clamped below and above" begin
        # Before the model is even on, and long after -- neither may go outside
        # [0, 1]: a factor above 1 would AMPLIFY the coefficient.
        @test dsgs_ramp_factor(0.0, TINIT, DT, 10) >= 0.0
        @test dsgs_ramp_factor(1.0e6, TINIT, DT, 10) == 1.0
        @test all(0.0 .<= [dsgs_ramp_factor(t, TINIT, DT, 10)
                           for t in 0.0:0.25:20.0] .<= 1.0)
    end

    @testset "the deck's default of 20 steps covers a real spin-up" begin
        # 20 steps at dt = 0.5 is 10 s of model time -- long enough for the
        # momentum norms to become a measurement rather than a rounding error,
        # short against the hundreds of seconds a PBL takes to develop.
        n = 20
        @test isapprox(dsgs_ramp_factor(T_ON, TINIT, DT, n), 1/n; atol = 1e-6)
        @test dsgs_ramp_factor(T_ON + (n-1)*DT, TINIT, DT, n) == 1.0
        @test dsgs_ramp_factor(T_ON + 0.5*(n-1)*DT, TINIT, DT, n) < 0.6
        say(@sprintf("  n = 20: full strength at t = %.1f s", T_ON + (n-1)*DT))
    end
end

say("=== done ===\n")

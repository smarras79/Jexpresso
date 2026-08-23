#=============================================================================
 test_joint_stability.jl -- the IMEX figure of merit, pinned.

 The explicit half's imaginary radius is the number a Runge-Kutta tableau is
 normally judged by, and using it to pick a HEVI tableau is what put :ARS343
 in this package as the default. A HEVI acoustic mode with horizontal
 wavenumber kx and vertical kz loads BOTH halves at once -- zE = Δt·c·kx
 explicit, zI = Δt·c·kz implicit -- so the whole (zE, zI) rectangle the mesh
 supports has to sit inside the unit disc, not just its two edges. The two
 criteria rank the ARS family differently, and the joint one is the one that
 decides whether a run survives.

 These are regression tests in the strict sense: each number here was measured
 once and would have caught the bug before it cost a week.
=============================================================================#

using Test
using Jexpresso
const J = Jexpresso

@testset "IMEX joint stability" begin

    ars232 = J.ark_tableau(:ARS232)
    ars343 = J.ark_tableau(:ARS343)
    ars443 = J.ark_tableau(:ARS443)

    @testset "the explicit-only radius ranks ARS343 first" begin
        # This is the trap, stated as a test so it cannot be un-learned: on the
        # explicit axis alone ARS343 genuinely is the best of the three.
        @test J.ark_imaginary_radius(ars343) > J.ark_imaginary_radius(ars232)
        @test J.ark_imaginary_radius(ars343) > J.ark_imaginary_radius(ars443)
        @test J.ark_imaginary_radius(ars343) ≈ 2.828 atol = 0.01
    end

    @testset "the joint region ranks it last" begin
        # rtb_hevi: explicit half 26.6 1/s, implicit half 69.5 1/s (the latter
        # MEASURED from the assembled column operator, not estimated).
        rE, rI = 26.6, 69.5

        @test J.ark_joint_dt_max(ars232, rE, rI) > 0.045
        @test J.ark_joint_dt_max(ars443, rE, rI) > 0.045
        @test J.ark_joint_dt_max(ars343, rE, rI) < 0.005

        # At the step size the case actually ran at, ARS343 amplifies and the
        # other two are exactly neutral -- which is what a hyperbolic operator
        # must be.
        Δt = 0.03
        @test J.ark_joint_amplification(ars343, Δt*rE, Δt*rI) > 1.02
        @test J.ark_joint_amplification(ars232, Δt*rE, Δt*rI) ≤ 1 + 1e-9
        @test J.ark_joint_amplification(ars443, Δt*rE, Δt*rI) ≤ 1 + 1e-9
    end

    @testset "growth rate predicts the observed failure" begin
        # 2.4% per step at Δt = 0.03 is 0.79 1/s, so ~17 s of model time to
        # amplify by 1e6. The case blew up at t ≈ 15 s. That the rate is in 1/s
        # rather than per step is why two different Δt failed at the SAME
        # physical time, which is the clue that pointed here.
        rE, rI, Δt = 26.6, 69.5, 0.03
        amp = J.ark_joint_amplification(ars343, Δt*rE, Δt*rI)
        rate = log(amp) / Δt
        @test 0.5 < rate < 1.2
        @test 10.0 < log(1e6) / rate < 30.0
    end

    @testset "zI = 0 is the edge that hides it" begin
        # ARS343 is fine on the explicit axis and unstable just off it. Any
        # check that only samples zI = 0 passes.
        @test J.ark_joint_amplification(ars343, 2.5, 0.0) ≤ 1 + 1e-9
        @test J.ark_joint_amplification(ars343, 2.5, 1.0) > 1.0
    end

    @testset "the shipped default is a joint-stable tableau" begin
        @test J.HEVI_ARK().tab.name === :ARS232
        rE, rI = 26.6, 69.5
        @test J.ark_joint_amplification(J.HEVI_ARK().tab, 0.05*rE, 0.05*rI) ≤ 1 + 1e-9
    end

    @testset "the rectangle is HEVI's set, and not IMEX3D's" begin
        # This repository ships TWO defaults that disagree about ARS343:
        #
        #     HEVI_ARK()  -> :ARS232      (everything above says why)
        #     IMEX_ARK()  -> :ARS343
        #
        # and a reader who finds one without the other will reasonably conclude
        # the second is a mistake. It is not, and the difference is not a
        # matter of judgement -- it is which (zE, zI) pairs the mesh can
        # produce.
        #
        # HEVI splits the acoustics between the halves: zE is HORIZONTAL sound
        # and zI VERTICAL sound, set by independent wavenumbers, so the whole
        # rectangle is reachable and ARS343's unstable corner is in it. Move
        # ALL the acoustics implicit and zE becomes ADVECTION at the same
        # wavenumber that sets zI, so zE/zI = |v|/c and the reachable set is a
        # WEDGE that ARS343's unstable region lies outside of.
        #
        # See src/kernel/solvers/hevi/README_IMEX3D.md, and
        # test/imex3d/test_wedge_stability.jl, which recomputes the whole
        # comparison without needing this dependency stack.
        @test J.IMEX_ARK().tab.name === :ARS343

        rE, rI = 26.6, 69.5          # the same mesh as above
        Δt = 0.03
        # unusable on the rectangle...
        @test J.ark_joint_amplification(ars343, Δt*rE, Δt*rI) > 1.02
        # ...and neutral on the wedge, at any subsonic Mach number.
        for mach in (0.05, 0.1, 0.3)
            @test J.ark_wedge_amplification(ars343, Δt*rE, Δt*rI, mach)[1] ≤ 1 + 1e-6
        end

        # And on the wedge it is the BEST of the family, which is the reverse
        # of the ranking twenty lines above. Reported as the largest zE each
        # stays neutral to, since that is what converts into Δt.
        wedge_lim(tab; mach = 0.1) = J.ark_wedge_dt_max(tab, 1.0, 1.0/mach, mach)
        @test wedge_lim(ars343) > wedge_lim(ars443) > wedge_lim(ars232)
        @test wedge_lim(ars343) ≈ J.ark_imaginary_radius(ars343) rtol = 0.02
    end
end

# The bisection in ark_joint_dt_max used to run inside a fixed hi = 0.5 and
# return that cap whenever the true limit was above it -- indistinguishable, in
# the setup report, from a real limit, and erring towards calling a Δt safe.
# At rate_exp = rate_imp = 1 both ARS232 and ARS443 sit well above 0.5, so a
# capped implementation returns 0.4999 for both.
@testset "joint Δt_max is not capped by the search bracket" begin
    ars232 = J.ark_tableau(:ARS232)
    ars443 = J.ark_tableau(:ARS443)
    @test J.ark_joint_dt_max(ars232, 1.0, 1.0) > 1.0
    @test J.ark_joint_dt_max(ars443, 1.0, 1.0) > 1.0
    # ARS443 is uncoupled: its joint limit is its explicit imaginary radius and
    # does not move with the stiffness ratio. ARS232 is mildly coupled and
    # degrades. Both must be found without a caller-supplied bracket.
    @test J.ark_joint_dt_max(ars443, 1.0, 40.0) ≈ J.ark_imaginary_radius(ars443) rtol=0.01
    @test J.ark_joint_dt_max(ars232, 1.0, 40.0) < J.ark_joint_dt_max(ars232, 1.0, 1.0)
end

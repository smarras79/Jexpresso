#---------------------------------------------------------------------------------
# test/test_plastic_cylinder.jl — quadrature-point history on the quasi-static
# solid path (src/kernel/solvers/quasistatic_solid.jl), with J2 plasticity.
#
#   julia --project=. test/test_plastic_cylinder.jl
#
# Needs the package (`using Jexpresso`). About three minutes.
#
#   1. THE MATERIAL POINT (solid_constitutive_update!, model :j2_plasticity)
#      a. A cyclic uniaxial-strain path (ε_yy = ε_zz = 0) with isotropic AND
#         kinematic hardening, against an independent analytic integration:
#         along a fixed deviatoric direction with linear hardening the response
#         is piecewise linear, so it is integrated exactly, event by event.
#      b. The consistent tangent against central finite differences, at
#         plastic states, with and without hardening.
#   2. THE CYLINDER (problems/Elasticity/plastic_cylinder2d), bore displaced
#      radially, perfectly plastic:
#      a. before first yield the bore pressure is Lamé's and ēᵖ = 0;
#      b. at every increment the bore pressure p = R_inner,x / a matches the 1D
#         radial reference (same constitutive routine, same increments);
#      c. it levels off at the exact limit p_L = (2/√3) σy0 ln(b/a);
#      d. Newton converges in a few iterations per increment.
#
# The quadrature-point path with the LINEAR viscoelastic model is checked
# against the nodal path in test_viscoelastic_cylinder.jl.
#---------------------------------------------------------------------------------

using Test
using Printf
using DelimitedFiles
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using Jexpresso

const CSV_PATH = joinpath(ROOT, "output", "Elasticity", "plastic_cylinder2d", "output", "solid_history.csv")

# Functions the case defines are `include`d by run_case into a newer world.
case_fn(f::Symbol, args...; kw...) = Base.invokelatest(getfield(Jexpresso, f), args...; kw...)

#
# Uniaxial strain, analytically. e = dev ε = εxx (2/3, -1/3, -1/3, 0) = e_N N with
# N = (2, -1, -1, 0)/√6 and e_N = √(2/3) εxx. Everything stays along N:
# eᵖ = p N, X = (2G(e_N - p) - Hkin p) N, R = √(2/3)(σy0 + Hiso ēᵖ). A monotone
# strain segment is elastic until s·X reaches R and then plastic at the rate
# dp = s·2G de_N/(2G + Hkin + 2Hiso/3), dēᵖ = √(2/3)|dp|.
#
function uniaxial_exact(m, εpath)
    G, K = m.G, m.K
    Hm = 2G + m.Hkin + (2/3)*m.Hiso
    p = 0.0; ebar = 0.0
    σxx = zeros(length(εpath))
    for n = 2:length(εpath)
        ea = sqrt(2/3)*εpath[n-1]; eb = sqrt(2/3)*εpath[n]
        s  = sign(eb - ea)
        if s != 0
            R  = sqrt(2/3)*(m.σy0 + m.Hiso*ebar)
            es = p + (s*R + m.Hkin*p)/(2G)                 # strain at which s·X reaches R
            if s*(eb - es) > 0                             # the segment ends plastic
                Δp    = s*2G*abs(eb - es)/Hm
                p    += Δp
                ebar += sqrt(2/3)*abs(Δp)
            end
        end
        σxx[n] = K*εpath[n] + 2G*(eb - p)*(2/sqrt(6))
    end
    return σxx
end

@testset verbose = true "J2 plasticity on quadrature-point history" begin

    #---------------------------------------------------------------------------
    @testset "material point: cyclic uniaxial strain, isotropic + kinematic hardening" begin
        m = Jexpresso.SolidJ2Plasticity(1000.0, 0.3, 1000/(3*0.4), 1000/2.6, 10.0, 50.0, 100.0)
        # load to +5%, reverse to -5%, reload to +5%: the reversals sit on step ends
        εpath = vcat(range(0, 0.05, length = 41), range(0.05, -0.05, length = 81)[2:end],
                     range(-0.05, 0.05, length = 81)[2:end])
        hold = zeros(5); hnew = zeros(5); σ = zeros(4); D = zeros(3, 3)
        σnum = zeros(length(εpath))
        for n = 2:length(εpath)
            Jexpresso.solid_constitutive_update!(σ, D, hnew, m, εpath[n], 0.0, 0.0, hold, 0.0)
            σnum[n] = σ[1]
            hold .= hnew
        end
        σex = uniaxial_exact(m, εpath)
        err = maximum(abs, σnum .- σex)/maximum(abs, σex)
        @printf("   max |σxx − exact| / max|σxx| = %.2e over %d steps (Bauschinger loop)\n", err, length(εpath) - 1)
        @test err < 1.0e-12
        @test hold[5] > 0 && abs(hold[1] + hold[2] + hold[3]) < 1.0e-14   # ēᵖ grew; eᵖ traceless
    end

    #---------------------------------------------------------------------------
    @testset "material point: consistent tangent = finite differences" begin
        for (Hiso, Hkin) in ((0.0, 0.0), (50.0, 100.0))
            m = Jexpresso.SolidJ2Plasticity(1000.0, 0.3, 1000/(3*0.4), 1000/2.6, 10.0, Hiso, Hkin)
            hold = [0.004, -0.001, -0.003, 0.002, 0.006]       # a plastic history
            ε0   = (0.031, -0.012, 0.017)                      # εxx, εyy, εxy: a plastic step
            σ = zeros(4); D = zeros(3, 3); hn = zeros(5)
            Jexpresso.solid_constitutive_update!(σ, D, hn, m, ε0..., hold, 0.0)
            @test hn[5] > hold[5]                              # it IS a plastic step
            h = 1.0e-7
            Dfd = zeros(3, 3)
            for j = 1:3
                δε = j == 3 ? (0.0, 0.0, h/2) : (j == 1 ? (h, 0.0, 0.0) : (0.0, h, 0.0))   # γ = 2εxy
                σp = zeros(4); σm = zeros(4)
                Jexpresso.solid_constitutive_update!(σp, zeros(3, 3), zeros(5), m, (ε0 .+ δε)..., hold, 0.0)
                Jexpresso.solid_constitutive_update!(σm, zeros(3, 3), zeros(5), m, (ε0 .- δε)..., hold, 0.0)
                Dfd[:, j] = ([σp[1], σp[2], σp[4]] .- [σm[1], σm[2], σm[4]])./(2h)
            end
            rel = maximum(abs, D .- Dfd)/maximum(abs, D)
            @printf("   Hiso = %5.1f, Hkin = %5.1f: max |D − D_fd| / max|D| = %.2e\n", Hiso, Hkin, rel)
            @test rel < 1.0e-6
            @test maximum(abs, D .- D') < 1.0e-9*maximum(abs, D)   # associative: symmetric
        end
    end

    #---------------------------------------------------------------------------
    @testset "thick cylinder: first yield, 1D radial reference, limit load" begin
        withenv("JEXPRESSO_PC_OUTFORMAT" => "none", "JEXPRESSO_PC_PLOT" => "false") do
            rm(CSV_PATH; force = true)
            redirect_stdout(devnull) do
                Jexpresso.run_case("Elasticity", "plastic_cylinder2d")
            end
        end
        data, header = readdlm(CSV_PATH, ',', Float64; header = true)
        col(name) = data[:, findfirst(==(name), vec(header))]
        P   = case_fn(:plastic_cylinder_properties)
        fy  = case_fn(:plastic_cylinder_first_yield)
        pL  = case_fn(:plastic_cylinder_limit_pressure)
        t   = col("t")
        δ   = P.δmax .* t
        p   = col("R_inner_x") ./ P.a
        py  = col("R_inner_y") ./ P.a
        epA = col("A_ep_bar")
        its = col("newton_iters")
        ref = case_fn(:plastic_cylinder_radial_reference, δ)

        # (a) elastic range
        ela  = findall(x -> 0 < x < 0.999*fy.δ_e, δ)
        pLam = fy.p_e .* δ[ela] ./ fy.δ_e
        e_el = maximum(abs, p[ela] .- pLam)/fy.p_e
        @printf("   elastic increments: %d, max |p − p_Lamé|/p_e = %.2e, max ēᵖ(A) = %.1e\n",
                length(ela), e_el, maximum(epA[ela]))
        @test length(ela) >= 3
        @test e_el < 1.0e-6
        @test maximum(epA[ela]) == 0.0

        # (b) the whole curve against the radial reference; the bore pressure is
        #     the same seen from either direction (axisymmetry)
        e_ref = maximum(abs, p .- ref.p)/pL
        e_sym = maximum(abs, p .- py)/pL
        @printf("   max |p − p_ref|/p_L = %.2e over %d increments; |p_x − p_y|/p_L = %.1e\n",
                e_ref, length(p) - 1, e_sym)
        @test e_ref < 2.0e-3
        @test e_sym < 1.0e-6

        # (c) the limit load
        @printf("   p(δ = %.3g) = %.5f,  reference %.5f,  exact limit p_L = %.5f  (%.3f%%)\n",
                δ[end], p[end], ref.p[end], pL, 100*(p[end]/pL - 1))
        @test abs(p[end]/pL - 1) < 5.0e-3
        @test abs(ref.p[end]/pL - 1) < 5.0e-3

        # (d) Newton
        @printf("   Newton iterations per increment: mean %.2f, max %d\n",
                sum(its[2:end])/(length(its) - 1), Int(maximum(its)))
        @test maximum(its) <= 8
    end
end

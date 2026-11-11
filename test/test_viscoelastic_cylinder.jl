#---------------------------------------------------------------------------------
# test/test_viscoelastic_cylinder.jl — the quasi-static viscoelastic solver
# (src/kernel/solvers/quasistatic_solid.jl) on the thick-walled cylinder of
# Artioli, Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60 (2017), Sec. 4.2
# (problems/Elasticity/viscoelastic_cylinder2d).
#
#   julia --project=. test/test_viscoelastic_cylinder.jl
#
# Needs the package (`using Jexpresso`): it runs the real case 24 times, with
# the JEXPRESSO_VE_* overrides of its deck, and reads back solid_history.csv.
# About six minutes on four cores; the first run also pays the compilation.
#
# Two references, which is what lets the errors be TOLD APART:
#
#   closed form     u_r(r,t) by the correspondence principle, Heaviside load
#                   (user_analytic.jl: cylinder_ur_exact). Exact in space AND
#                   time.
#   time-discrete   the same exponential recursion the solver runs, on the
#                   Lamé family u_r = A r + B/r, which this geometry never
#                   leaves (user_analytic.jl: cylinder_time_discrete_reference).
#                   Exact in space, discrete in time with the SAME Δt.
#
#   SEM − time-discrete   = spatial error alone
#   SEM − closed form     = spatial + time error; at nop = 4 the time error
#                           dominates by orders of magnitude, so its rate is
#                           the integrator's order.
#
# What is pinned down:
#   1. t = 0⁺ is the elastic Lamé solution.
#   2. The spatial error alone, every step of every run, and its p-convergence.
#   3. The time error converges at second order for the "exponential" (default)
#      and "midpoint" recursions and at first order for backward Euler, both
#      Prony sets — the integrator and its order, as a Δt sweep.
#   4. The paper's Fig. 5.2 (FEAP, Δt = 1, the load entering in the first
#      increment): reproduced to the precision the figure carries.
#---------------------------------------------------------------------------------

using Test
using Printf
using DelimitedFiles

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)                      # the deck's mesh and output paths are relative to the repo root

using Jexpresso

const CSV_PATH = joinpath(ROOT, "output", "Elasticity", "viscoelastic_cylinder2d", "output", "solid_history.csv")

#
# Artioli et al. (2017) Fig. 5.2, the "Ref" curves (an overkill mixed u-p-εᵛ FEM
# run, Δt = 1), at integration steps k = 0..20. Read off the VECTOR paths of the
# arXiv PDF (arXiv:1701.06676v2, p. 22) — marker centres mapped through the
# axis ticks, calibration residual ~1e-5 — not eyeballed from a raster. The u_B
# markers are triangles whose bounding-box centre sits a constant 0.00124
# (panel a) / 0.00021 (panel b) below the data point; that offset, read from
# step 0 where u = 0, has been removed. Ref, Q9 and VEM k = 2 agree to the
# fourth digit in the figure.
#
const FIG52 = Dict(
    "ve1" => (uA = [0.00000, 0.05828, 0.09239, 0.12600, 0.15927, 0.19222, 0.22484, 0.25712, 0.28908, 0.32070, 0.35199,
                    0.38295, 0.41358, 0.44404, 0.47418, 0.50398, 0.53345, 0.56259, 0.59156, 0.62021, 0.64852],
              uB = [0.00000, 0.03460, 0.05199, 0.06887, 0.08560, 0.10215, 0.11838, 0.13460, 0.15050, 0.16639, 0.18196,
                    0.19752, 0.21292, 0.22798, 0.24305, 0.25795, 0.27269, 0.28725, 0.30182, 0.31606, 0.33030]),
    "ve2" => (uA = [0.00000, 0.05030, 0.06801, 0.08115, 0.09094, 0.09821, 0.10363, 0.10763, 0.11064, 0.11285, 0.11453,
                    0.11575, 0.11669, 0.11737, 0.11788, 0.11828, 0.11853, 0.11876, 0.11890, 0.11904, 0.11913],
              uB = [0.00000, 0.03054, 0.03960, 0.04627, 0.05118, 0.05484, 0.05753, 0.05955, 0.06105, 0.06216, 0.06301,
                    0.06364, 0.06409, 0.06443, 0.06469, 0.06489, 0.06503, 0.06511, 0.06520, 0.06525, 0.06531]))

#
# The case's own functions (cylinder_properties, the two references) are defined
# by run_case's `include`, i.e. in a newer world than the one this file was
# compiled in: call them through invokelatest.
#
case_fn(f::Symbol, args...; kw...) = Base.invokelatest(getfield(Jexpresso, f), args...; kw...)

"""
Run the case once and compare its probe history with both references.
"""
function run_ve(; set = "ve1", dt = 1.0, load = "step", integrator = "exponential", nop = 4)
    withenv("JEXPRESSO_VE_SET"        => set,
            "JEXPRESSO_VE_DT"         => string(dt),
            "JEXPRESSO_VE_LOAD"       => load,
            "JEXPRESSO_VE_INTEGRATOR" => integrator,
            "JEXPRESSO_VE_NOP"        => string(nop),
            "JEXPRESSO_VE_OUTFORMAT"  => "none") do

        rm(CSV_PATH; force = true)
        redirect_stdout(devnull) do
            Jexpresso.run_case("Elasticity", "viscoelastic_cylinder2d")
        end
        data, header = readdlm(CSV_PATH, ',', Float64; header = true)
        col(name) = data[:, findfirst(==(name), vec(header))]
        t, uA, uB = col("t"), col("A_ux"), col("B_ux")

        P   = case_fn(:cylinder_properties)
        ref = case_fn(:cylinder_time_discrete_reference, t; integrator = Symbol(integrator),
                      initial_solve = (load == "step"))
        refA = ref.A .* P.a .+ ref.B ./ P.a
        refB = ref.A .* P.b .+ ref.B ./ P.b
        exA  = [case_fn(:cylinder_ur_exact, P.a, ti) for ti in t]
        exB  = [case_fn(:cylinder_ur_exact, P.b, ti) for ti in t]

        # The solver writes error columns only when it compares (a step load):
        # with the load in the first increment the closed form is another problem.
        isig  = findfirst(==("err_sigma_L2_rel"), vec(header))
        scale = maximum(abs, exA)
        return (; t, uA, uB,
                  e_space = max(maximum(abs, uA .- refA), maximum(abs, uB .- refB))/scale,
                  e_total = max(maximum(abs, uA .- exA),  maximum(abs, uB .- exB))/scale,
                  e_sigma0 = isig === nothing ? NaN : data[1, isig])
    end
end

rate(e1, e2, r = 2.0) = log(e1/e2)/log(r)

@testset verbose = true "viscoelastic thick cylinder (Artioli et al. 2017, Sec. 4.2)" begin

    #---------------------------------------------------------------------------
    @testset "t = 0⁺ is the Lamé solution; the spatial error alone is tiny" begin
        r = run_ve()                                     # the deck as shipped
        P = case_fn(:cylinder_properties)
        lame(x) = P.p*P.a^2/(P.b^2 - P.a^2)*(1 + P.ν)/P.E*((1 - 2P.ν)*x + P.b^2/x)
        @printf("   t=0⁺   u_A = %.10f (Lamé %.10f)   u_B = %.10f (Lamé %.10f)\n",
                r.uA[1], lame(P.a), r.uB[1], lame(P.b))
        @printf("   nop=4, Δt=1: spatial error %.2e, total error %.2e (relative to max u_A)\n",
                r.e_space, r.e_total)
        @test abs(r.uA[1] - lame(P.a)) < 1.0e-6*lame(P.a)
        @test abs(r.uB[1] - lame(P.b)) < 1.0e-6*lame(P.b)
        @test r.e_space < 1.0e-6
        @test r.e_sigma0 < 5.0e-4                        # stresses: one order below u
    end

    #---------------------------------------------------------------------------
    @testset "spatial error: p-convergence (Δt = 1, against the time-discrete reference)" begin
        errs = Float64[]
        for nop in 2:6
            push!(errs, run_ve(nop = nop).e_space)
            @printf("   nop = %d   spatial error %.3e\n", nop, errs[end])
        end
        # Spectral, with an odd/even pattern (measured 2.3e-3, 2.5e-6, 6.0e-7,
        # 4.1e-9, 1.0e-9): monotone, over five decades from nop = 2 to 6, and
        # below 1e-6 at the deck's nop = 4.
        @test all(diff(errs) .< 0)
        @test errs[1]/errs[end] > 1.0e5
        @test errs[3] < 1.0e-6
    end

    #---------------------------------------------------------------------------
    @testset "time error: Δt sweep against the closed form" begin
        dts = [1.0, 0.5, 0.25, 0.125]
        for (set, integrator, order) in (("ve1", "exponential",    2),
                                         ("ve2", "exponential",    2),
                                         ("ve1", "midpoint",       2),
                                         ("ve2", "backward_euler", 1))
            runs = [run_ve(set = set, dt = dt, integrator = integrator) for dt in dts]
            e = [r.e_total for r in runs]
            rates = [rate(e[i], e[i+1]) for i = 1:length(e)-1]
            @printf("   %-4s %-15s  e(Δt) = %s   rates %s\n", set, integrator,
                    join((@sprintf("%.3e", x) for x in e), "  "),
                    join((@sprintf("%.3f", x) for x in rates), "  "))
            @test all(r -> r.e_space < 1.0e-6, runs)    # space never pollutes the rate
            @test all(x -> abs(x - order) < 0.1, rates)
        end
    end

    #---------------------------------------------------------------------------
    @testset "the paper's Fig. 5.2 (FEAP convention: load in the first increment, Δt = 1)" begin
        for set in ("ve1", "ve2")
            r   = run_ve(set = set, load = "first_increment")
            fig = FIG52[set]
            dA  = maximum(abs, r.uA .- fig.uA)
            dB  = maximum(abs, r.uB .- fig.uB)
            @printf("   %s: max |SEM − Fig. 5.2| = %.2e at A, %.2e at B   (u_A(20) = %.5f vs %.5f)\n",
                    set, dA, dB, r.uA[end], fig.uA[end])
            # The figure's own resolution: 1/150 of the axis per point, markers
            # placed to ~0.02 pt — i.e. ~1e-4 on panel (a), ~2e-5 on panel (b).
            tol = set == "ve1" ? 2.0e-4 : 5.0e-5
            @test dA < tol
            @test dB < tol
            @test r.e_space < 1.0e-6
        end
    end
end

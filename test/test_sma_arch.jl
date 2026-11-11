#---------------------------------------------------------------------------------
# test/test_sma_arch.jl — the shape-memory-alloy material of the quasi-static
# solid path (src/kernel/solvers/quasistatic_solid_sma.jl) and the SMA arch of
# Artioli, Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60 (2017), Sec. 4.4
# (problems/Elasticity/sma_arch2d).
#
#   julia --project=. test/test_sma_arch.jl
#
# Needs the package (`using Jexpresso`). About five minutes.
#
#   1. THE MATERIAL POINT (no mesh)
#      a. The step is the minimiser of the incremental energy
#         Φ(x) = G‖e − x‖² + ½h‖x‖² + τ‖x‖ + R‖x − xₙ‖ over ‖x‖ ≤ εL: no
#         feasible perturbation, small or large, lowers it — 4000 random states
#         covering the four cases (elastic, flow, saturated, full reverse).
#      b. Uniaxial stress (σyy = 0 by Newton, σzz = 0 by the plane-stress
#         wrapper) against the exact one-dimensional problem — the energy
#         ½E(ε − √(2/3)q)² + ½hq² + τ|q| + R|q − qₙ| minimised directly: the
#         shape-memory loop at T = Mf and the superelastic loop at Mf + 80 K.
#      c. Heating at zero stress after a load: ‖eᵗʳ‖ = min(εL, max(0, (R − τ)/h)),
#         zero (full recovery) once τ ≥ R.
#      d. The plane-stress tangent against central differences of the update;
#         the elastic one is E/(1−ν²)[1 ν 0; ν 1 0; 0 0 (1−ν)/2].
#   2. THE ELASTIC ARCH against its exact solution (Timoshenko & Goodier's
#      curved bar, user_analytic.jl; JEXPRESSO_SMA_VERIFY=elastic):
#      a. the exact displacement satisfies the plane-stress strain–stress
#         relations (finite differences of u against Hooke's law);
#      b. the SEM error converges spectrally, nop = 3 … 6, plane stress, and
#         is ~1e-9 at nop = 4 in plane strain too.
#   3. THE DEVICE, the load–temperature history of Table 3 (plane stress, nop 4):
#      a. no transformation at q = 1.5 N/mm, and there the compliance of the
#         exact elastic arch;
#      b. the clamp reaction balances the load at every increment;
#      c. shape memory: a residual displacement after each unloading, the
#         reversed cycle the mirror image of the first;
#      d. shape recovery: the transformation strain is gone once
#         T ≥ Mf + R/β = 242.4 K, the arch back in its initial configuration.
#
# The paper's Fig. 5.6 is NOT a pass/fail reference: with the stated data the
# computed arch is softer than the paper's (see the README of the case).
#---------------------------------------------------------------------------------

using Test
using Printf
using Random
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using Jexpresso
const J = Jexpresso

include(joinpath(ROOT, "tools", "vtu_reader.jl"))
using .VTUReader

const OUT_DIR  = joinpath(ROOT, "output", "Elasticity", "sma_arch2d", "output")
const CSV_PATH = joinpath(OUT_DIR, "solid_history.csv")

# Functions the case defines are `include`d by run_case into a newer world.
case_fn(f::Symbol, args...; kw...) = Base.invokelatest(getfield(Jexpresso, f), args...; kw...)

function read_csv(path)
    lines = readlines(path)
    hdr   = String.(split(lines[1], ","))
    rows  = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    return Dict(h => [r[i] for r in rows] for (i, h) in enumerate(hdr))
end

function run_arch(; verify = false, nop = 4, hypothesis = "plane_stress", vtk = false)
    redirect_stdout(devnull) do
        withenv("JEXPRESSO_SMA_VERIFY" => (verify ? "elastic" : nothing), "JEXPRESSO_SMA_NOP" => string(nop),
                "JEXPRESSO_SMA_HYPOTHESIS" => hypothesis, "JEXPRESSO_SMA_PLOT" => "false",
                "JEXPRESSO_SMA_OUTFORMAT" => (vtk ? "vtk" : "none"),
                "JEXPRESSO_SMA_STEPS" => "40", "JEXPRESSO_SMA_HEAT_STEPS" => "10",
                "JEXPRESSO_SMA_E" => nothing, "JEXPRESSO_SMA_THICKNESS" => nothing) do
            Jexpresso.run_case("Elasticity", "sma_arch2d")
        end
    end
    return read_csv(CSV_PATH)
end

# the material of the paper
const E, ν = 53000.0, 0.36
material(ps = true) = J.SolidSMA(E, ν, E/(3(1 - 2ν)), E/(2(1 + ν)), 1000.0, 2.1, 223.0, 0.04, 50.0,
                                 ps, t -> 223.0, Ref(223.0))
const M = material()
const R = sqrt(2/3)*M.σy0

rdev(s) = (x = s*randn(); y = s*randn(); (x, y, -x - y, s*randn()))
Φ(x, a, xn, τ) = 0.5*(2M.G + M.h)*J._dd(x, x) - J._dd(a, x) + τ*J._nn(x) + R*J._nn(J._lin(1.0, x, -1.0, xn))

# uniaxial stress, strain-controlled in x: εyy by Newton (σyy = 0), σzz = 0 in the wrapper
function uniaxial(m, path, τ)
    xn = J._SMA0; εzz = 0.0; εyy = 0.0
    σ = zeros(4); D = zeros(3, 3); out = Float64[]
    for εxx in path
        for it = 1:100
            J.sma_point!(σ, D, m, εxx, εyy, 0.0, xn, εzz, τ)
            abs(σ[2]) < 1e-10 && break
            εyy -= σ[2]/D[2,2]
            it == 100 && error("uniaxial: σyy did not converge")
        end
        xn, εzz, _ = J.sma_point!(σ, D, m, εxx, εyy, 0.0, xn, εzz, τ)
        push!(out, σ[1])
    end
    return out
end

# the exact uniaxial response: minimise the 1D energy over q (golden section; convex)
function uniaxial_exact(m, path, τ)
    q = 0.0; out = Float64[]
    for ε in path
        f(qq) = 0.5*m.E*(ε - sqrt(2/3)*qq)^2 + 0.5*m.h*qq^2 + τ*abs(qq) + R*abs(qq - q)
        lo, hi = -m.εL, m.εL
        for _ = 1:300
            a1 = lo + 0.381966011250105*(hi - lo); a2 = lo + 0.618033988749895*(hi - lo)
            f(a1) < f(a2) ? (hi = a2) : (lo = a1)
        end
        q = 0.5*(lo + hi)
        push!(out, m.E*(ε - sqrt(2/3)*q))
    end
    return out
end

# stress-controlled point (σxx = target, σyy = 0): Newton with backtracking
function stress_point!(m, xn, εzz, τ, target, εx, εy)
    σ = zeros(4); D = zeros(3, 3)
    res(ex, ey) = (J.sma_point!(σ, D, m, ex, ey, 0.0, xn, εzz, τ); [σ[1] - target, σ[2]])
    r = res(εx, εy)
    for it = 1:200
        norm(r) < 1e-10 && break
        d = D[1:2, 1:2] \ r
        λ = 1.0
        while true
            rt = res(εx - λ*d[1], εy - λ*d[2])
            if norm(rt) < (1 - 1e-4*λ)*norm(r) || λ < 1e-8
                εx -= λ*d[1]; εy -= λ*d[2]; r = rt; break
            end
            λ /= 2
        end
        it == 200 && error("stress_point!: no convergence")
    end
    x, ez, _ = J.sma_point!(σ, D, m, εx, εy, 0.0, xn, εzz, τ)
    return x, ez, εx, εy
end

@testset verbose = true "SMA: material point and the arch of Artioli et al. (2017), Sec. 4.4" begin

    #---------------------------------------------------------------------------
    @testset "material point: the step minimises the incremental energy" begin
        Random.seed!(3)
        worst = 0.0; counts = Dict{Symbol,Int}(); feasible = true
        for trial = 1:4000
            τ  = rand() < 0.5 ? 0.0 : 200*rand()
            a  = rdev(rand()*400)
            xn = rdev(0.03*rand())
            rn = J._nn(xn); rn > M.εL && (xn = J._sc(M.εL/rn*(rand() < 0.3 ? 1.0 : rand()), xn))
            rand() < 0.1 && (xn = J._SMA0)
            x, case = J.sma_return_map(M, a, xn, τ)
            counts[case] = get(counts, case, 0) + 1
            feasible &= J._nn(x) <= M.εL*(1 + 1e-12) && abs(x[1] + x[2] + x[3]) < 1e-15
            f0 = Φ(x, a, xn, τ)
            for k = 1:100
                y = J._lin(1.0, x, 1.0, rdev(10.0^(-rand()*6)*M.εL))
                J._nn(y) > M.εL && (y = J._sc(M.εL/J._nn(y), y))
                worst = min(worst, (Φ(y, a, xn, τ) - f0)/max(abs(f0), 1.0))
            end
        end
        @printf("    cases %s; min (Φ(y) − Φ(x))/|Φ(x)| = %.2e\n", string(counts), worst)
        @test feasible
        @test all(c -> get(counts, c, 0) > 10, (:elastic, :flow, :saturated, :reverse))
        @test worst >= -1e-13
    end

    @testset "material point: uniaxial loops against the exact 1D problem" begin
        # shape memory at T = Mf: load to 5 %, reverse to −5 %
        p1 = vcat(range(0, 0.05, length = 201), range(0.05, -0.05, length = 401)[2:end])
        e1 = maximum(abs.(uniaxial(M, p1, 0.0) .- uniaxial_exact(M, p1, 0.0)))
        # superelastic at Mf + 80 K: load to 5 % and back
        τh = M.β*80
        p2 = vcat(range(0, 0.05, length = 201), range(0.05, 0.0, length = 201)[2:end])
        s2 = uniaxial(M, p2, τh)
        e2 = maximum(abs.(s2 .- uniaxial_exact(M, p2, τh)))
        @printf("    max |σ − σ_1D|: shape-memory loop %.1e MPa, superelastic loop %.1e MPa\n", e1, e2)
        @test e1 < 1e-4
        @test e2 < 1e-4
        # superelastic: the forward transformation starts at √(3/2)(R + τ) — the
        # last elastic point lies below it, the next above — and unloading to
        # ε = 0 leaves no stress (the transformation fully reversed)
        σAM = sqrt(1.5)*(R + τh)
        k   = findfirst(>(σAM), s2)
        @test s2[k-1] ≈ E*p2[k-1] rtol = 1e-12
        @test s2[k] < E*p2[k]
        @test abs(s2[end]) < 1e-8
    end

    @testset "material point: shape recovery on heating at zero stress" begin
        xn = J._SMA0; εzz = 0.0; εx = 0.0; εy = 0.0
        for σt in vcat(range(0, 120, length = 61), range(120, 0, length = 61)[2:end])
            xn, εzz, εx, εy = stress_point!(M, xn, εzz, 0.0, σt, εx, εy)
        end
        @test J._nn(xn) ≈ M.εL rtol = 1e-12                 # saturated, kept at σ = 0
        @test εx ≈ sqrt(2/3)*M.εL rtol = 1e-9               # the residual strain
        worst = 0.0
        for T in range(223, 303, length = 41)
            τ = J.sma_tau(M, T)
            xn, εzz, εx, εy = stress_point!(M, xn, εzz, τ, 0.0, εx, εy)
            worst = max(worst, abs(J._nn(xn) - min(M.εL, max(0.0, (R - τ)/M.h))))
        end
        @printf("    max |‖eᵗʳ‖ − min(εL, max(0, (R − τ)/h))| while heating: %.1e; final εxx = %.1e\n", worst, εx)
        @test worst < 1e-12
        @test J._nn(xn) == 0.0
        @test abs(εx) < 1e-14
    end

    @testset "material point: tangents" begin
        Random.seed!(7)
        worst = 0.0; n = 0
        for trial = 1:300
            τ  = rand() < 0.5 ? 0.0 : 150*rand()
            xn = rdev(0.02*rand()); J._nn(xn) > M.εL && (xn = J._sc(M.εL/J._nn(xn), xn))
            ε  = (0.02*randn(), 0.02*randn(), 0.01*randn())
            σ0 = zeros(4); D0 = zeros(3, 3)
            _, _, case = J.sma_point!(σ0, D0, M, ε..., xn, 0.0, τ)
            case === :elastic && continue
            n += 1
            Dfd = zeros(3, 3)
            for j = 1:3
                δ = 1e-7; ep = collect(ε); em = collect(ε)
                ep[j] += (j == 3 ? δ/2 : δ); em[j] -= (j == 3 ? δ/2 : δ)
                sp = zeros(4); sm = zeros(4)
                J.sma_point!(sp, zeros(3, 3), M, ep..., xn, 0.0, τ); J.sma_point!(sm, zeros(3, 3), M, em..., xn, 0.0, τ)
                Dfd[:, j] = ([sp[1], sp[2], sp[4]] .- [sm[1], sm[2], sm[4]]) ./ (2δ)
            end
            worst = max(worst, maximum(abs.(D0 .- Dfd))/maximum(abs.(Dfd)))
        end
        Del = zeros(3, 3); J.sma_point!(zeros(4), Del, M, 1e-5, 2e-5, 0.0, J._SMA0, 0.0, 0.0)
        eel = maximum(abs.(Del .- E/(1 - ν^2) .* [1 ν 0; ν 1 0; 0 0 (1 - ν)/2]))/E
        @printf("    plane-stress tangent vs central differences (%d inelastic states): %.1e; elastic: %.1e\n", n, worst, eel)
        @test worst < 1e-5
        @test eel < 1e-14
    end

    #---------------------------------------------------------------------------
    @testset "elastic arch against the exact curved-bar solution" begin
        # load the case once (defines curved_bar_exact)
        run_arch(verify = true, nop = 4)
        # the exact displacement against Hooke's law: strains by central differences
        worst = 0.0
        for (r, θ) in ((3.6, 0.3), (4.0, 1.2), (4.4, 2.0), (3.8, 2.9))
            x, y = r*cos(θ), r*sin(θ); h = 1e-5
            ex = case_fn(:curved_bar_exact, x, y)
            dux = (case_fn(:curved_bar_exact, x + h, y).ux - case_fn(:curved_bar_exact, x - h, y).ux)/(2h)
            duy = (case_fn(:curved_bar_exact, x, y + h).uy - case_fn(:curved_bar_exact, x, y - h).uy)/(2h)
            gxy = (case_fn(:curved_bar_exact, x, y + h).ux - case_fn(:curved_bar_exact, x, y - h).ux)/(2h) +
                  (case_fn(:curved_bar_exact, x + h, y).uy - case_fn(:curved_bar_exact, x - h, y).uy)/(2h)
            hk = ((ex.sxx - ν*ex.syy)/E, (ex.syy - ν*ex.sxx)/E, 2(1 + ν)*ex.sxy/E)
            sc = maximum(abs.(hk))
            worst = max(worst, maximum(abs.((dux, duy, gxy) .- hk))/sc)
        end
        @printf("    exact displacement vs plane-stress Hooke's law: %.1e\n", worst)
        @test worst < 1e-7
        errs = Float64[]
        for nop = 3:6
            H = run_arch(verify = true, nop = nop)
            push!(errs, H["err_u_max"][end])
        end
        Hs = run_arch(verify = true, nop = 4, hypothesis = "plane_strain")
        @printf("    |e_u|max, nop = 3..6 (plane stress): %s;  plane strain, nop = 4: %.1e\n",
                join((@sprintf("%.1e", e) for e in errs), ", "), Hs["err_u_max"][end])
        @test all(errs[k+1] < errs[k]/20 for k = 1:length(errs)-1)
        @test errs[2] < 2e-9                                   # nop = 4
        @test Hs["err_u_max"][end] < 2e-9
    end

    #---------------------------------------------------------------------------
    @testset "the device: Table 3, plane stress, nop = 4" begin
        H  = run_arch(vtk = true)
        t  = H["t"]; uA = H["A_ux"]
        at(tb) = findfirst(x -> abs(x - tb) < 1e-9, t)
        P  = case_fn(:sma_arch_properties)
        q  = [case_fn(:sma_arch_load, s) for s in t]

        # a. elastic start: no transformation at q = 1.5 N/mm, and the exact compliance
        files = [m.captures[1] for m in eachmatch(r"file=\"([^\"]+)\"", read(joinpath(OUT_DIR, "simulation.pvd"), String))]
        @test length(files) == length(t)
        g1 = read_vtk_grid(joinpath(OUT_DIR, files[at(0.025)]))
        @test maximum(g1.data["etr_sat"]) == 0
        comp = case_fn(:curved_bar_exact, 4.5, 0.0; P = 1.0, hypothesis = :plane_stress).ux
        @printf("    q = 1.5: u_A/q = %.6f mm per N/mm, the exact elastic arch %.6f (point clamp)\n", uA[at(0.025)]/1.5, comp)
        @test uA[at(0.025)]/1.5 ≈ comp rtol = 1e-3

        # b. equilibrium: the clamp reaction balances q over the 1 mm end face
        eq = maximum(abs.(H["R_clamp_x"] .+ q .* (P.b - P.a)))/P.qmax
        @printf("    max |R_clamp_x + q|/qmax = %.1e,  max |R_clamp_y|/qmax = %.1e\n", eq, maximum(abs.(H["R_clamp_y"]))/P.qmax)
        @test eq < 1e-6
        @test maximum(abs.(H["R_clamp_y"]))/P.qmax < 1e-6

        # c. shape memory and the mirrored cycle
        @printf("    u_A at t = 1..5: %s mm\n", join((@sprintf("%.4f", uA[at(tb)]) for tb in 1:5), ", "))
        @test uA[at(1)] > uA[at(2)] > 1.0
        @test uA[at(3)] < uA[at(4)] < -1.0
        @test abs(uA[at(3)] + uA[at(1)]) < 1e-2*abs(uA[at(1)])
        @test abs(uA[at(4)] + uA[at(2)]) < 1e-2*abs(uA[at(2)])

        # d. recovery: gone once T ≥ Mf + R/β (t ≥ 4.24), not before
        Trec = P.Mf + R/P.β
        g42 = read_vtk_grid(joinpath(OUT_DIR, files[at(4.2)]))
        g43 = read_vtk_grid(joinpath(OUT_DIR, files[at(4.3)]))
        g5  = read_vtk_grid(joinpath(OUT_DIR, files[end]))
        @printf("    max ‖eᵗʳ‖/εL: t = 4.2 (T = %.0f K) %.3f, t = 4.3 (T = %.0f K) %.1e; recovery from %.1f K; final max |u| = %.1e mm\n",
                case_fn(:sma_arch_temperature, 4.2), maximum(g42.data["etr_sat"]),
                case_fn(:sma_arch_temperature, 4.3), maximum(g43.data["etr_sat"]), Trec,
                maximum(hypot.(g5.data["ux"], g5.data["uy"])))
        @test maximum(g42.data["etr_sat"]) > 0.01
        @test maximum(g43.data["etr_sat"]) == 0
        @test maximum(hypot.(g5.data["ux"], g5.data["uy"])) < 1e-8
        @test abs(uA[end]) < 1e-8
        its = H["newton_iters"][2:end]
        @printf("    Newton: %.2f iterations per increment on average, %d increment(s) cut\n",
                sum(its)/length(its), count(>(0), H["newton_cuts"][2:end]))
    end
end

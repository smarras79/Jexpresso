#---------------------------------------------------------------------------------
# test/test_perforated_plate.jl — the perforated plastic plate of Artioli,
# Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60 (2017), Sec. 4.3
# (problems/Elasticity/perforated_plate2d): J2 perfect plasticity, plane
# strain, the top edge pulled up by δ = 2 mm in 400 increments, nop = 4.
#
#   julia --project=. test/test_perforated_plate.jl
#
# Needs the package (`using Jexpresso`). About five minutes, the first run's
# compilation included. There is no closed form; the references are the
# paper's numbers:
#
#   1. Table 1: the horizontal displacement of A = (R, 0) and the vertical one
#      of B = (0, R) at the end of the loading, against the paper's reference
#      (an overkill mixed u-p-εᵛ FEM solution): within 0.5 %.
#   2. Fig. 5.4, the reaction sum on the top edge (tools/artioli2017_data,
#      read from the PDF's vector graphics):
#      a. elastic steps: the paper's sums are exactly 10x Jexpresso's force per
#         unit thickness — the paper's force is that of a 10 mm thick plate —
#         to 0.2 % (its five methods agree among themselves to 0.1 %). From
#         step 4: the digitized curves carry a ~0.2 kg offset (-0.2 at step
#         0), 0.1 % of the reaction at step 1;
#      b. the plateau (limit load), x 10: within 1.5 % of the paper's quadratic
#         VEM/FEM values at every load step past 150;
#      c. perfect plasticity: from step 200 to 400 the reaction rises < 0.2 %.
#   3. Global equilibrium of the converged states: no body force and a
#      traction-free hole, so the top and bottom reactions cancel, and so do
#      the lateral ones, to the Newton tolerance.
#   4. Newton: no increment needs cutting; on average at most 8 iterations per
#      increment (the paper's Table 2: 6.1–6.7 for its quadratic elements).
#   5. The nodal reaction field ry the VTK files carry, read back with
#      tools/vtu_reader.jl and summed over the top edge, is R_top_y.
#---------------------------------------------------------------------------------

using Test
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using Jexpresso

include(joinpath(ROOT, "tools", "vtu_reader.jl"))
using .VTUReader

const OUT_DIR  = joinpath(ROOT, "output", "Elasticity", "perforated_plate2d", "output")
const CSV_PATH = joinpath(OUT_DIR, "solid_history.csv")

function read_csv(path)
    lines = readlines(path)
    hdr   = String.(split(lines[1], ","))
    rows  = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    return Dict(h => [r[i] for r in rows] for (i, h) in enumerate(hdr))
end

function read_fig5_4()
    lines = filter(l -> !startswith(l, "#"), readlines(joinpath(ROOT, "tools", "artioli2017_data", "fig5_4.csv")))
    hdr   = String.(split(lines[1], ","))
    rows  = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    return Dict(h => [r[i] for r in rows] for (i, h) in enumerate(hdr))
end

const METHODS = ["VEM_Quad", "VEM_Tri", "VEM_Voronoi", "FEM_Q9", "FEM_T6"]

println("Running the perforated plate (400 increments, nop = 4) ...")
withenv("JEXPRESSO_PP_STEPS" => "400", "JEXPRESSO_PP_NOP" => "4",
        "JEXPRESSO_PP_OUTFORMAT" => "vtk", "JEXPRESSO_PP_PLOT" => "false") do
    Jexpresso.run_case("Elasticity", "perforated_plate2d")
end

const H   = read_csv(CSV_PATH)
const REF = read_fig5_4()
const STEP = round.(Int, H["t"] .* 400)

@testset verbose = true "perforated plastic plate (Artioli et al. 2017, Sec. 4.3)" begin

    @testset "Table 1: u_A and v_B at the end of the loading" begin
        @test H["t"][end] ≈ 1.0
        uA = H["A_ux"][end]; vB = H["B_uy"][end]
        @printf("    u_A = %.4f mm (Ref 2.741, %+.2f%%), v_B = %.4f mm (Ref 1.859, %+.2f%%)\n",
                uA, 100*(uA/2.741 - 1), vB, 100*(vB/1.859 - 1))
        @test abs(uA/2.741 - 1) < 0.005
        @test abs(vB/1.859 - 1) < 0.005
    end

    @testset "Fig. 5.4: reaction sum on the top edge" begin
        R = H["R_top_y"]
        # a. elastic steps
        ratios = Float64[]
        for s = 4:40, m in METHODS
            k = findfirst(==(s), STEP); j = findfirst(==(s), round.(Int, REF["step"]))
            push!(ratios, REF[m][j]/R[k])
        end
        @printf("    elastic steps 4-40: paper / Jexpresso = %.4f (%.4f - %.4f)\n",
                sum(ratios)/length(ratios), minimum(ratios), maximum(ratios))
        @test all(r -> abs(r/10 - 1) < 0.002, ratios)
        # b. plateau, on the paper's scale
        worst = 0.0
        for s = 150:400
            k = findfirst(==(s), STEP); j = findfirst(==(s), round.(Int, REF["step"]))
            avg = sum(REF[m][j] for m in METHODS)/length(METHODS)
            worst = max(worst, abs(10*R[k]/avg - 1))
        end
        @printf("    steps 150-400: |10 R / paper mean - 1| <= %.2f%%; end: 10 R = %.0f kg, paper %.0f - %.0f kg\n",
                100*worst, 10*R[end], minimum(REF[m][end] for m in METHODS), maximum(REF[m][end] for m in METHODS))
        @test worst < 0.015
        # c. perfect plasticity: the limit load
        k200 = findfirst(==(200), STEP)
        @test 0 <= R[end]/R[k200] - 1 < 0.002
    end

    @testset "global equilibrium of every converged increment" begin
        ty = maximum(abs.(H["R_top_y"] .+ H["R_bottom_y"]) ./ max.(abs.(H["R_top_y"]), 1.0))
        tx = maximum(abs.(H["R_left_x"] .+ H["R_right_x"]) ./ max.(abs.(H["R_right_x"]), 1.0))
        @printf("    max |R_top + R_bottom| / |R_top| = %.2e,  max |R_left + R_right| / |R_right| = %.2e\n", ty, tx)
        @test ty < 1e-7
        @test tx < 1e-7
    end

    @testset "Newton" begin
        its = H["newton_iters"][2:end]; cuts = H["newton_cuts"][2:end]
        @printf("    %.2f iterations per increment on average, at most %d; %d increment(s) cut\n",
                sum(its)/length(its), maximum(its), count(>(0), cuts))
        @test all(==(0), cuts)
        @test sum(its)/length(its) <= 8
    end

    @testset "the nodal reaction field ry in the VTK output" begin
        pvd   = read(joinpath(OUT_DIR, "simulation.pvd"), String)
        files = [m.captures[1] for m in eachmatch(r"file=\"([^\"]+)\"", pvd)]
        @test length(files) == 101                 # t = 0, 0.01, ..., 1
        g   = read_vtk_grid(joinpath(OUT_DIR, files[end]))
        Hy  = maximum(g.y)
        top = findall(y -> abs(y - Hy) <= 1e-8*Hy, g.y)
        @test g.time ≈ 1.0
        @test sum(g.data["ry"][top]) ≈ H["R_top_y"][end] rtol = 1e-10
        # u_y is prescribed on the top and the bottom edge only: ry vanishes in between
        @test all(iszero, g.data["ry"][(g.y .> 1e-8*Hy) .& (g.y .< Hy*(1 - 1e-8))])
        @test !isempty(g.cells) && all(c -> length(c) == 4, g.cells)
    end
end

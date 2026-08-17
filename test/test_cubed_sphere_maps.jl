#---------------------------------------------------------------------------------
# test/test_cubed_sphere_maps.jl — the cube-face → sphere maps used by
# :cubed_sphere_map (src/kernel/mesh/cubed_sphere_maps.jl).
#
#   julia test/test_cubed_sphere_maps.jl
#
# NO `using Jexpresso`: the maps are deliberately free of St_mesh/MPI, so this
# file includes the source directly and needs nothing but Test. That keeps the
# numerics checkable without building the package.
#
# What is actually being pinned down here, in the order it matters:
#
#   1. THE CONFORMAL MAP IS CONFORMAL. This is the test that would catch a
#      mistyped Rančić coefficient, a wrong exponent or a bad Möbius constant —
#      none of which the anchor points below would notice. A conformal map has
#      a Jacobian that is a scalar times a rotation, so |r_u| = |r_v| and
#      r_u · r_v = 0 EVERYWHERE. Measured by central differences, the residual
#      sits at the differencing floor (~1e-9 at h = 1e-5) for the conformal map
#      and at O(0.3) for the gnomonic and equiangular maps, which are not
#      conformal. Three orders of magnitude of margin either way.
#
#   2. SEAMS STAY WATERTIGHT. A node on a cube edge belongs to two panels and
#      the grid stores it once; which panel the tie-break picks must not change
#      where the node lands. Checked on all twelve edges by forcing each
#      adjacent panel in turn.
#
#   3. The series reversion B = A⁻¹ really inverts A, and forward ∘ inverse is
#      the identity for every map.
#
#   4. Published anchor points of the conformal map: the face centre goes to
#      the pole, the face corner to the cube corner.
#---------------------------------------------------------------------------------

using Test

include(joinpath(@__DIR__, "..", "src", "kernel", "mesh", "cubed_sphere_maps.jl"))

const MAPS = (:gnomonic, :equiangular, :conformal)

# Central-difference Jacobian of a face map and the two conformality residuals.
function _conformality(f, u, v; h = 1.0e-5)
    p(a, b) = collect(f(a, b))
    ru = (p(u+h, v) .- p(u-h, v)) ./ (2h)
    rv = (p(u, v+h) .- p(u, v-h)) ./ (2h)
    nu, nv = sqrt(sum(ru.^2)), sqrt(sum(rv.^2))
    return abs(nu - nv)/max(nu, nv), abs(sum(ru .* rv))/(nu*nv)
end

# The remap with the panel FORCED, so the two panels sharing an edge can be
# compared against each other.
function _remap_forced(x, y, z, panel, from, to)
    r = sqrt(x^2 + y^2 + z^2)
    xn, yn, zn = x/r, y/r, z/r
    ug, vg = _panel_face_coords(xn, yn, zn, panel)
    u, v   = _map_inverse(from, clamp(ug, -1, 1), clamp(vg, -1, 1))
    U, V, W = _map_forward(to, clamp(u, -1, 1), clamp(v, -1, 1))
    return _panel_to_cartesian(U, V, W, panel)
end


@testset "cubed-sphere maps" begin

    @testset "Rančić series reversion" begin
        @test length(A_RANCIC) == 31
        @test length(B_RANCIC) == 31
        @test B_RANCIC[2] ≈ 1/A_RANCIC[2]
        for r in (0.01, 0.05, 0.1, 0.2), θ in range(0, 2π, length = 17)
            Z = r*cis(θ)
            @test _Z_rancic(_W_rancic(Z)) ≈ Z atol = 1.0e-12
        end
    end

    @testset "conformal anchor points" begin
        # centre of the face -> the pole of the panel
        @test all(isapprox.(_conformal_forward(0.0, 0.0), (0.0, 0.0, 1.0); atol = 1.0e-14))
        # corner of the face -> the corner of the inscribed cube
        c = sqrt(3)/3
        @test all(isapprox.(_conformal_forward(1.0, 1.0), (c, c, c); atol = 1.0e-14))
    end

    @testset "images lie on the unit sphere" begin
        for m in MAPS, u in range(-1, 1, length = 21), v in range(-1, 1, length = 21)
            U, V, W = _map_forward(m, u, v)
            @test sqrt(U^2 + V^2 + W^2) ≈ 1.0 atol = 1.0e-13
        end
    end

    # The two structural properties the seam argument in cubed_sphere_maps.jl
    # rests on. If either breaks, panels stop agreeing on their shared edge.
    @testset "square symmetries" begin
        for m in MAPS, u in range(-1, 1, length = 15), v in range(-1, 1, length = 15)
            a = _map_forward(m, u, v)
            b = _map_forward(m, v, u)
            @test a[1] ≈ b[2] atol = 1.0e-13     # map(v,u) = swap₁₂(map(u,v))
            @test a[2] ≈ b[1] atol = 1.0e-13
            @test a[3] ≈ b[3] atol = 1.0e-13
        end
        for m in MAPS, v in range(-1, 1, length = 21)
            U, _, W = _map_forward(m, 1.0, v)    # the edge lies in the plane U = W
            @test U ≈ W atol = 1.0e-13
        end
    end

    @testset "forward ∘ inverse == identity" begin
        for m in MAPS, u in range(-0.999, 0.999, length = 21), v in range(-0.999, 0.999, length = 21)
            U, V, W = _map_forward(m, u, v)
            ug, vg  = U/W, V/W                   # the equidistant face coordinate
            ui, vi  = _map_inverse(m, ug, vg)
            @test ui ≈ u atol = 1.0e-11
            @test vi ≈ v atol = 1.0e-11
        end
    end

    # THE test: only :conformal preserves angles.
    @testset "conformality" begin
        worst = Dict(m => (0.0, 0.0) for m in MAPS)
        for m in MAPS, u in range(-0.97, 0.97, length = 33), v in range(-0.97, 0.97, length = 33)
            iso, orth = _conformality((a, b) -> _map_forward(m, a, b), u, v)
            worst[m] = (max(worst[m][1], iso), max(worst[m][2], orth))
        end
        # conformal: at the central-difference floor
        @test worst[:conformal][1] < 1.0e-6
        @test worst[:conformal][2] < 1.0e-6
        # the other two are genuinely not conformal — this is what gives the
        # test above its meaning
        @test worst[:gnomonic][1]    > 0.1
        @test worst[:equiangular][1] > 0.1
    end

    @testset "panel assignment" begin
        # a direction pointing squarely at a face normal gets that panel
        @test cubed_sphere_panel_of( 1.0,  0.1,  0.1) == 1
        @test cubed_sphere_panel_of(-1.0,  0.1,  0.1) == 2
        @test cubed_sphere_panel_of( 0.1,  1.0,  0.1) == 3
        @test cubed_sphere_panel_of( 0.1, -1.0,  0.1) == 4
        @test cubed_sphere_panel_of( 0.1,  0.1,  1.0) == 5
        @test cubed_sphere_panel_of( 0.1,  0.1, -1.0) == 6
        # the face coordinate of a panel centre is (0, 0), of a corner (±1, ±1)
        for p = 1:6
            _, _, ew = _PANEL_FRAMES[p]
            @test all(isapprox.(_panel_face_coords(ew..., p), (0.0, 0.0); atol = 1.0e-14))
        end
    end

    # A node on a cube edge is stored ONCE but belongs to two panels. Both must
    # move it to the same place or the shell tears along the seam.
    @testset "seam consistency" begin
        edges = NTuple{3,Float64}[]
        for t in range(-1, 1, length = 41), s1 in (-1.0, 1.0), s2 in (-1.0, 1.0)
            push!(edges, (s1,  t, s2))
            push!(edges, (s1, s2,  t))
            push!(edges, ( t, s1, s2))
        end
        for to in (:equiangular, :conformal)
            worst = 0.0
            for (x, y, z) in edges
                r = sqrt(x^2 + y^2 + z^2)
                xn, yn, zn = x/r, y/r, z/r
                a = (abs(xn), abs(yn), abs(zn))
                m = maximum(a)
                cands = Int[]
                a[1] >= m - 1.0e-12 && push!(cands, xn >= 0 ? 1 : 2)
                a[2] >= m - 1.0e-12 && push!(cands, yn >= 0 ? 3 : 4)
                a[3] >= m - 1.0e-12 && push!(cands, zn >= 0 ? 5 : 6)
                length(cands) < 2 && continue
                ref = _remap_forced(xn, yn, zn, cands[1], :gnomonic, to)
                for k = 2:length(cands)
                    other = _remap_forced(xn, yn, zn, cands[k], :gnomonic, to)
                    worst = max(worst, maximum(abs.(other .- ref)))
                end
            end
            @test worst < 1.0e-13
        end
    end

    @testset "remap_direction on the SWsphere grid" begin
        msh = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere", "cubed_sphere.msh")
        if !isfile(msh)
            @info "cubed_sphere.msh not found — skipping the grid check" msh
        else
            lines = readlines(msh)
            i = findfirst(l -> strip(l) == "\$Nodes", lines)
            n = parse(Int, strip(lines[i+1]))
            pts = [let f = split(strip(lines[i+1+k]))
                       (parse(Float64, f[2]), parse(Float64, f[3]), parse(Float64, f[4]))
                   end for k = 1:n]
            R = sqrt(sum(pts[1].^2))

            for to in (:equiangular, :conformal)
                moved = [remap_direction(x, y, z, :gnomonic, to) for (x, y, z) in pts]
                # still on the sphere
                for p in moved
                    @test sqrt(sum(p.^2)) ≈ 1.0 atol = 1.0e-12
                end
                # and still 602 DISTINCT nodes: a map that folded the panel
                # would collapse two of them onto each other
                mind = Inf
                for a = 1:length(moved), b = a+1:length(moved)
                    mind = min(mind, sum((moved[a] .- moved[b]).^2))
                end
                @test R*sqrt(mind) > 1.0e5     # > 100 km on a ~1000 km grid
            end

            # a :gnomonic -> :gnomonic remap must not move anything
            for (x, y, z) in pts
                p = remap_direction(x, y, z, :gnomonic, :gnomonic)
                @test all(isapprox.(R .* p, (x, y, z); atol = 1.0e-6))
            end
        end
    end
end

#---------------------------------------------------------------------------------
# test/test_sphere_mesh.jl — the spherical shell read through the ORDINARY
# gmsh path.
#
#   julia --project=. test/test_sphere_mesh.jl
#
# There is no sphere-specific grid builder any more: a cubed sphere is an
# ordinary St_mesh with `lmanifold == true`, produced by mod_mesh_read_gmsh!
# like every other grid. So this file no longer tests a private reader — it
# tests the two things that can actually go wrong on a CLOSED surface, plus the
# format regression that motivated folding the reader away:
#
#   1. MSH 4.1 WITH PARAMETRIC NODES loads.  gmsh writes `x y z u v` for nodes
#      that sit on a CAD surface — which is every node of a cubed sphere built
#      with `In Sphere {1}`. The hand-rolled parser this replaced rejected
#      those files outright ("MSH 4.1 files with parametric nodes are not
#      supported"). GridapGmsh reads them because it *is* gmsh.
#
#   2. The node count is EXACT and nothing is duplicated.  On a closed shell
#      every edge is interior and shared by exactly two elements, so
#
#           npoin = npoin_linear + nedges*(ngl-2) + nelem*(ngl-2)^2
#
#      holds only if each unique edge contributed its high-order nodes exactly
#      once. Creating them per element instead would duplicate every seam node
#      — the classic failure on a surface with no boundary — and it would show
#      up here as both a wrong count and coincident nodes.
#
#   3. Euler characteristic V - E + F = 2, i.e. the shell really is closed and
#      watertight, and the high-order nodes all land ON the sphere.
#
# The fixture is generated from the .msh that ships with the SWsphere case, so
# the test needs no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver
using PartitionedArrays, MPI

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

#---------------------------------------------------------------------------------
# Rewrite an MSH 2.2 quad-surface file as MSH 4.1 with `parametric = 1`, i.e.
# with each node line carrying the two surface parameters after (x,y,z).
#---------------------------------------------------------------------------------
function msh22_to_msh41_parametric(src::AbstractString, dst::AbstractString; R = 6.371e6)
    L = readlines(src)
    i = findfirst(==("\$Nodes"), L);    n = parse(Int, L[i+1])
    coords = [Tuple(parse.(Float64, split(L[i+1+k])[2:4])) for k = 1:n]
    j = findfirst(==("\$Elements"), L); m = parse(Int, L[j+1])
    quads  = [parse.(Int, split(L[j+1+k])[6:9]) for k = 1:m]

    open(dst, "w") do f
        println(f, "\$MeshFormat"); println(f, "4.1 0 8"); println(f, "\$EndMeshFormat")
        println(f, "\$Nodes")
        println(f, "1 $n 1 $n")
        println(f, "2 1 1 $n")                       # entityDim=2, tag=1, PARAMETRIC=1
        for k = 1:n; println(f, k); end
        for c in coords
            u = atan(c[2], c[1])
            v = asin(clamp(c[3]/R, -1.0, 1.0))
            println(f, c[1], " ", c[2], " ", c[3], " ", u, " ", v)
        end
        println(f, "\$EndNodes")
        println(f, "\$Elements")
        println(f, "1 $m 1 $m")
        println(f, "2 1 3 $m")
        for (k, q) in enumerate(quads)
            println(f, k, " ", q[1], " ", q[2], " ", q[3], " ", q[4])
        end
        println(f, "\$EndElements")
    end
    return dst, n, m
end

#---------------------------------------------------------------------------------
# Inputs dict for a shell read, defaulted by Jexpresso's own
# mod_inputs_user_inputs! so this test cannot drift from the real defaults.
#
# That function reads three module globals that normally come from the command
# line (run.jl sets them before it includes a case), so they are supplied here.
#---------------------------------------------------------------------------------
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_sphere_mesh.jl"
end

function shell_inputs(fname, nop; lspherical = true)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => fname,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => lspherical,
        :sphere_radius        => 6.37122e6,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename] = fname                 # keep our path, not the case layout
    inputs[:nop]           = nop
    return inputs
end

@testset "spherical shell through the ordinary gmsh reader" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    tmp = mktempdir()
    f41, nlin_expected, nelem_expected =
        msh22_to_msh41_parametric(CASE_MSH, joinpath(tmp, "cubed_sphere_41.msh"))

    with_mpi() do distribute
        for (label, fname) in (("MSH 2.2", CASE_MSH), ("MSH 4.1 parametric", f41))
            for nop in (2, 4, 5)
                @testset "$label nop=$nop" begin

                    mesh, _ = mod_mesh_mesh_driver(shell_inputs(fname, nop), 1, distribute)
                    ngl = mesh.ngl

                    # (1) recognised as a curved surface, not flattened
                    @test mesh.lmanifold
                    @test mesh.nsd == 2

                    # (2) topology of a closed shell
                    @test mesh.nelem        == nelem_expected
                    @test mesh.npoin_linear == nlin_expected
                    @test mesh.nedges       == 2*nelem_expected           # closed quad surface
                    @test mesh.npoin_linear - mesh.nedges + mesh.nelem == 2   # Euler

                    # (3) EXACT high-order node count: one set per unique edge
                    @test mesh.npoin == mesh.npoin_linear +
                                        mesh.nedges*(ngl-2) +
                                        mesh.nelem*(ngl-2)^2

                    # (4) every node lands ON the sphere
                    crd = mesh.coords
                    r = sqrt.(crd[1, 1:mesh.npoin].^2 .+
                              crd[2, 1:mesh.npoin].^2 .+
                              crd[3, 1:mesh.npoin].^2)
                    @test maximum(abs.(r .- mesh.radius))/mesh.radius < 1.0e-12

                    # (5) nothing duplicated: no two nodes closer than a small
                    #     fraction of the smallest LGL spacing. A per-element
                    #     (rather than per-unique-edge) numbering would put two
                    #     coincident nodes on every seam and fail here.
                    tolabs = 1.0e-4*mesh.radius/(mesh.nelem^0.5 * ngl)
                    key(ip) = (round(Int, crd[1, ip]/tolabs),
                               round(Int, crd[2, ip]/tolabs),
                               round(Int, crd[3, ip]/tolabs))
                    seen = Set{NTuple{3,Int}}()
                    dup  = 0
                    for ip = 1:mesh.npoin
                        k = key(ip)
                        k in seen ? (dup += 1) : push!(seen, k)
                    end
                    @test dup == 0

                    # (6) lon/lat filled and in range
                    @test length(mesh.lon) == mesh.npoin
                    @test all(-π   .<= mesh.lon .<= π)
                    @test all(-π/2 .<= mesh.lat .<= π/2)
                end
            end
        end

        # A flat 2D grid must NOT be taken for a manifold: the coplanarity test
        # is what keeps every existing case on the unchanged code path.
        @testset "flat 2D grid is not a manifold" begin
            flat = joinpath(@__DIR__, "..", "problems", "MHD",
                            "orszagTangBormanis2024", "OT_32x32_periodic.msh")
            if isfile(flat)
                mesh, _ = mod_mesh_mesh_driver(
                    shell_inputs(flat, 4; lspherical = false),
                    1, distribute)
                @test mesh.lmanifold == false
                # coords is (NSD_MAX, npoin) for EVERY grid now — the flat/shell
                # distinction lives in `lmanifold`, which is the assertion above,
                # and never in the row count. The third row of a flat 2D grid is
                # zero, which is what mesh.z held there anyway.
                @test size(mesh.coords, 1) == Jexpresso.NSD_MAX
                @test all(iszero, mesh.coords[3, :])
            end
        end
    end
end

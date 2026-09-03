#---------------------------------------------------------------------------------
# test/test_sphere_orientation.jl — a cubed-sphere .msh with a panel wound
# INWARD must be read and re-oriented, not refused.
#
#   julia --project=. test/test_sphere_orientation.jl
#
# WHAT IT GUARDS. On the shell the surface Jacobian is the signed triple
# product x̂·(a_ξ × a_η), so an element wound clockwise as seen from outside has
# a negative Jacobian and build_sphere_metrics stops:
#
#   ERROR sphere_metrics.jl: non-positive surface Jacobian in element 1 at node (1,1).
#
# gmsh orients each panel of a cubed sphere on its own, and a 32×32 grid written
# by gmsh 4.1 from the shipped cubed_sphere.geo arrived with its whole -z panel
# wound inward (the other five outward), which is exactly that error on element
# 1. orient_shell_elements_outward! in mesh.jl now transposes the (ξ,η)
# numbering of every inward element at read time.
#
# The fixture is made here, from the shipped 10×10 MSH 2.2 grid: every quad
# whose centroid lies on the -z panel is written with its vertices in reverse
# order (100 of the 600), into a temporary file. That file is read through the
# ordinary path, and then:
#
#   1. the reader reports 100 flipped elements and every element of the mesh
#      that comes out has a_ξ × a_η pointing outward;
#   2. build_sphere_metrics does not throw, every surface Jacobian is positive,
#      the mass matrix sums to 4πR², and check_sphere_metrics passes all its
#      tests — i.e. the re-oriented grid is as good as the original;
#   3. the original grid is untouched by the pass (0 flips).
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, check_sphere_metrics
using PartitionedArrays, MPI
using Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_sphere_orientation.jl"
end

function shell_inputs(msh, nop)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => msh,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 6.37122e6,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)
    inputs[:gmsh_filename] = msh
    inputs[:nop]           = nop
    return inputs
end

#
# Write a copy of an MSH 2.2 file with the quads of one cube panel reversed.
# `panel` is (axis, sign), e.g. (3, -1) for the -z panel. Returns the number of
# quads reversed.
#
function write_inverted_panel(src, dst; panel = (3, -1))
    lines = readlines(src)
    inod  = findfirst(==("\$Nodes"), lines)
    nnod  = parse(Int, lines[inod+1])
    xyz   = Dict{Int, NTuple{3,Float64}}()
    for k = 1:nnod
        t = split(lines[inod+1+k])
        xyz[parse(Int, t[1])] = (parse(Float64, t[2]), parse(Float64, t[3]), parse(Float64, t[4]))
    end
    iel  = findfirst(==("\$Elements"), lines)
    nel  = parse(Int, lines[iel+1])
    nrev = 0
    for k = 1:nel
        t = split(lines[iel+1+k])
        parse(Int, t[2]) == 3 || continue              # a 4-node quad
        ids = parse.(Int, t[end-3:end])
        c   = sum(xyz[i][panel[1]] for i in ids)/4
        cabs = ntuple(d -> abs(sum(xyz[i][d] for i in ids)/4), 3)
        # on the requested panel: that axis dominates the centroid, with that sign
        (argmax(cabs) == panel[1] && sign(c) == panel[2]) || continue
        lines[iel+1+k] = join(vcat(t[1:end-4], string.(reverse(ids))), " ")
        nrev += 1
    end
    open(dst, "w") do io
        for l in lines
            println(io, l)
        end
    end
    return nrev
end

# a_ξ × a_η · p₁ on the linear corners of every element, as the reader tests it
function count_inward(mesh)
    ngl = Int(mesh.ngl); n = 0
    for iel = 1:Int(mesh.nelem)
        i1 = mesh.connijk[iel,1,1]; i2 = mesh.connijk[iel,ngl,1]; i4 = mesh.connijk[iel,1,ngl]
        a = (mesh.x[i2]-mesh.x[i1], mesh.y[i2]-mesh.y[i1], mesh.z[i2]-mesh.z[i1])
        b = (mesh.x[i4]-mesh.x[i1], mesh.y[i4]-mesh.y[i1], mesh.z[i4]-mesh.z[i1])
        s = (a[2]*b[3]-a[3]*b[2])*mesh.x[i1] + (a[3]*b[1]-a[1]*b[3])*mesh.y[i1] + (a[1]*b[2]-a[2]*b[1])*mesh.z[i1]
        s < 0 && (n += 1)
    end
    return n
end


@testset "spherical shell: inward-wound elements are re-oriented at read time" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    tmpdir = mktempdir()
    bad    = joinpath(tmpdir, "cubed_sphere_inverted.msh")
    nrev   = write_inverted_panel(CASE_MSH, bad; panel = (3, -1))
    @info "  reversed $nrev quads (the -z panel) into $bad"
    @test nrev == 100

    with_mpi() do distribute

        nop = 4
        R   = 6.37122e6

        #--- the untouched grid: nothing to flip
        mesh0, _ = mod_mesh_mesh_driver(shell_inputs(CASE_MSH, nop), 1, distribute)
        @test count_inward(mesh0) == 0
        @test Jexpresso.orient_shell_elements_outward!(mesh0) == 0

        #--- the inverted grid: read, and every element comes out outward
        mesh, _ = mod_mesh_mesh_driver(shell_inputs(bad, nop), 1, distribute)
        @test Int(mesh.nelem) == 600
        @test Int(mesh.npoin) == Int(mesh0.npoin)
        @test count_inward(mesh) == 0
        # a second pass finds nothing left to do
        @test Jexpresso.orient_shell_elements_outward!(mesh) == 0

        #--- the metrics build and are as good as on the original grid
        metrics = build_sphere_metrics(mesh, shell_inputs(bad, nop); verbose = false)
        @test minimum(metrics.Je) > 0
        @test abs(sum(metrics.M) - 4π*R^2) < 1.0e-10*4π*R^2
        @test check_sphere_metrics(mesh, metrics; verbose = false)

        metrics0 = build_sphere_metrics(mesh0, shell_inputs(CASE_MSH, nop); verbose = false)
        @test abs(sum(metrics.M) - sum(metrics0.M)) < 1.0e-10*sum(metrics0.M)
        # the transposition permutes the nodes of an element; the element areas
        # are the same set of numbers
        @test isapprox(sort(vec(sum(metrics.Je, dims = (2, 3)))),
                       sort(vec(sum(metrics0.Je, dims = (2, 3)))); rtol = 1.0e-12)
    end
end

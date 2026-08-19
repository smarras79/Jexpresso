#---------------------------------------------------------------------------------
# test/test_mesh_coords_sync.jl — mesh.coords and mesh.x/y/z must agree.
#
#   julia --project=. test/test_mesh_coords_sync.jl
#
# WHY THIS EXISTS
#
# mesh.coords (NSD_MAX, npoin) is the canonical coordinate storage. mesh.x, mesh.y
# and mesh.z are DUPLICATES of its rows, kept in step by hand: eight explicit
# `mesh.coords[r, :] = mesh.<r>[:]` assignments scattered through mesh.jl, against
# forty-one places that assign mesh.x/y/z. Nothing enforces the correspondence —
# a code path that updates one and forgets the other simply produces a mesh whose
# two halves disagree, and which half a given routine reads is a matter of which
# was written first.
#
# That is not a hypothetical either. mesh.jl has at least one site that replaces
# mesh.x/mesh.y wholesale as npoin grows (the Laguerre semi-infinite extension)
# with no coords assignment anywhere after it, and another where the sync line is
# commented out.
#
# This file pins the invariant so it cannot rot further, and so the staged removal
# of mesh.x/y/z has something to hold onto: every step of that migration has to
# leave this test green.
#
# NOTE ON WHAT IT STILL PROVES, now that stage 2 has landed. mesh.x/y/z are views
# of coords' rows, so the equality assertions below are TAUTOLOGICAL — there is no
# longer any storage that could disagree, which was the point. What they still
# catch is a regression that reintroduces separate storage, and what the rest of
# the file checks is not tautological at all: that coords is (NSD_MAX, npoin) on
# every grid shape, that the row count does not depend on nsd or on lmanifold, and
# that row 3 of a flat grid is zero. Those are real properties of the mesh builder
# and they are what the 1D/Laguerre case actually exercises.
#
# AND THE VALUES. Shape is not enough: a 1D build that reallocated coords to zeros
# produced an array of exactly the right shape holding no grid at all, and every 1D
# case plotted as one vertical line at x = 0. The two 1D testsets therefore assert
# what the coordinates ARE -- element vertices on the linear grid, LGL nodes
# increasing across each element, the span equal to [xmin, xmax].
#
# It is deliberately about the INVARIANT, not about any one case. Add a case to
# CASES and it is covered.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, NSD_MAX
using PartitionedArrays, MPI
using Printf

@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SoliWaveIsland"
    user_input_file            = "test/test_mesh_coords_sync.jl"
end

const ROOT = normpath(joinpath(@__DIR__, ".."))

#
# (label, case directory, mesh file relative to it, nsd)
#
# Kept to grids that ship inside the repository, so the test needs no gmsh binary
# and no JexpressoMeshes checkout.
#
const CASES = [
    ("2D flat  (SoliWaveIsland)", "problems/ShallowWater/SoliWaveIsland", "SoliWaveIsland.msh", 2),
    ("2D flat  (CompEuler/theta)", "problems/CompEuler/theta",            "hexa_TFI_10x10.msh", 2),
    ("2D shell (SWsphere)",        "problems/ShallowWater/SWsphere",      "cubed_sphere.msh",   2),
]

#
# A natively-built 1D grid WITH Laguerre semi-infinite elements. This is the path
# that motivated the test: mod_mesh_build_mesh! replaces mesh.x/mesh.y wholesale
# as npoin grows to take in the semi-infinite layers (mesh.jl, the
# `mesh.x = x_new` block), and there is no coords assignment after it.
#
const LAGUERRE_1D = Dict{Symbol,Any}(
    :lread_gmsh         => false,
    :nsd                => 1,
    :nop                => 6,
    :nop_laguerre       => 24,
    :llaguerre_1d_right => true,
    :llaguerre_1d_left  => true,
    :laguerre_beta      => 1.0,
    :xmin               => -2.5,
    :xmax               =>  2.5,
    :nelx               =>  50,
)

#
# A natively-built 1D grid with NO Laguerre layers -- the plain sod1d/wav1d shape.
# Ten elements on [0, 1] at nop = 4, so the element vertices are exactly the
# multiples of 0.1 and every interior LGL node sits strictly between two of them.
#
const PLAIN_1D = Dict{Symbol,Any}(
    :lread_gmsh => false,
    :nsd        => 1,
    :nop        => 4,
    :xmin       =>  0.0,
    :xmax       =>  1.0,
    :nelx       =>  10,
)

function build_mesh(casedir, msh; nop = 4, extra = Dict{Symbol,Any}())
    inputs = Dict{Symbol,Any}(
        :lread_gmsh          => msh !== nothing,
        :nop                 => nop,
        :interpolation_nodes => "lgl",
        :backend             => Jexpresso.CPU(),
        :lgrid_only          => true,
    )
    msh === nothing || (inputs[:gmsh_filename] = joinpath(ROOT, casedir, msh))
    for (k, v) in extra; inputs[k] = v; end
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)
    msh === nothing || (inputs[:gmsh_filename] = joinpath(ROOT, casedir, msh))
    inputs[:nop]           = nop
    #
    # mod_mesh_build_mesh! (the NATIVE 1D builder) reads the module-level global
    # `inputs` rather than its argument — run.jl sets it before it ever gets
    # called, so nothing in a normal run notices. A test that drives the mesh
    # driver directly has to set it too, or the 1D path dies with
    # `UndefVarError: inputs not defined in Jexpresso`.
    #
    @eval Jexpresso global inputs = $inputs
    for (k, v) in extra; inputs[k] = v; end
    mesh, _ = with_mpi() do distribute
        mod_mesh_mesh_driver(inputs, 1, distribute)
    end
    return mesh
end

@testset verbose = true "mesh.coords rows vs the x/y/z views" begin

    for (label, casedir, msh, nsd) in CASES
        path = joinpath(ROOT, casedir, msh)
        if !isfile(path)
            @warn "skipping (mesh not in the repository)" label path
            continue
        end
        @testset "$label" begin
            mesh = build_mesh(casedir, msh)
            npoin = Int(mesh.npoin)

            # Stage 0 of the migration: three rows, always, whatever nsd is.
            @test size(mesh.coords, 1) == NSD_MAX
            @test size(mesh.coords, 2) == npoin

            # the duplicates must be the same length as the thing they duplicate
            @test length(view(mesh.coords, 1, :)) == npoin
            @test length(view(mesh.coords, 2, :)) == npoin
            @test length(view(mesh.coords, 3, :)) == npoin

            # ...and hold the same numbers, exactly. Not approximately: these are
            # copies of each other, so any difference at all is a missed sync.
            for (r, v, nm) in ((1, view(mesh.coords, 1, :), "x"), (2, view(mesh.coords, 2, :), "y"), (3, view(mesh.coords, 3, :), "z"))
                d = maximum(abs, mesh.coords[r, :] .- v)
                d == 0 || @info "coords row $r disagrees with mesh.$nm" maxdiff=d
                @test mesh.coords[r, :] == v
            end
        end
    end

    @testset "1D + Laguerre (natively built, no mesh file)" begin
        mesh = build_mesh("", nothing; nop = 6, extra = LAGUERRE_1D)
        npoin = Int(mesh.npoin)
        @test size(mesh.coords, 1) == NSD_MAX
        @test size(mesh.coords, 2) == npoin
        @test length(view(mesh.coords, 1, :)) == npoin
        for (r, v, nm) in ((1, view(mesh.coords, 1, :), "x"), (2, view(mesh.coords, 2, :), "y"), (3, view(mesh.coords, 3, :), "z"))
            length(v) == npoin || @info "mesh.$nm has the wrong length" length(v) npoin
            @test length(v) == npoin
            d = maximum(abs, mesh.coords[r, 1:min(end,length(v))] .- v[1:min(end,npoin)])
            d == 0 || @info "coords row $r disagrees with mesh.$nm" maxdiff=d
            @test mesh.coords[r, :] == v
        end

        #
        # ...and the numbers themselves. See the note on the plain-1D testset
        # below: shape alone does not tell you the grid survived the build.
        #
        x  = mesh.coords[1, 1:npoin]
        xi = mesh.coords[1, 1:Int(mesh.npoin_original)]   # the LGL part
        @test allunique(x)                                # no two nodes on top of each other
        @test minimum(xi) ≈ LAGUERRE_1D[:xmin]
        @test maximum(xi) ≈ LAGUERRE_1D[:xmax]
        # the semi-infinite layers reach out past both ends
        @test minimum(x) < LAGUERRE_1D[:xmin]
        @test maximum(x) > LAGUERRE_1D[:xmax]
    end

    #
    # THE VALUES, not just the shape.
    #
    # Everything above this point is about lengths and row counts, and a mesh can
    # pass all of it while holding no grid at all: add_high_order_nodes_1D_native_mesh!
    # used to REALLOCATE coords to zeros (harmless when mesh.x was separate storage
    # and coords got refilled from it afterwards, fatal once coords became the
    # storage), and the result was an npoin-column array of the right shape with
    # every 1D node sitting on x = 0. Shape tests saw nothing; every 1D case
    # plotted as a single vertical line. So: check the coordinates are the grid.
    #
    @testset "1D native grid: coords are the grid, not zeros" begin
        mesh  = build_mesh("", nothing; nop = 4, extra = PLAIN_1D)
        npoin = Int(mesh.npoin)
        nelem = Int(mesh.nelem)
        ngl   = Int(mesh.ngl)
        xmin  = PLAIN_1D[:xmin]
        xmax  = PLAIN_1D[:xmax]
        Δx    = (xmax - xmin)/nelem
        x     = mesh.coords[1, 1:npoin]

        # exactly one node is at the origin: the left endpoint. This single
        # assertion is the one that fails when the grid collapses.
        @test count(iszero, x) == 1
        @test allunique(x)
        @test minimum(x) ≈ xmin
        @test maximum(x) ≈ xmax

        # the element vertices are still the linear grid, in order
        for i = 1:Int(mesh.npoin_linear)
            @test x[i] ≈ xmin + (i - 1)*Δx
        end

        # and each element's LGL nodes increase across it, endpoint to endpoint
        for iel = 1:nelem
            xe = [mesh.coords[1, mesh.connijk[iel, i, 1, 1]] for i = 1:ngl]
            @test issorted(xe)
            @test allunique(xe)
            @test xe[1]   ≈ xmin + (iel - 1)*Δx
            @test xe[end] ≈ xmin + iel*Δx
        end

        # rows 2 and 3 stay zero on a 1D grid
        @test all(iszero, mesh.coords[2, 1:npoin])
        @test all(iszero, mesh.coords[3, 1:npoin])
    end

    #
    # A flat grid's third row is zero — the same thing mesh.z held. This is the
    # property that lets row 3 exist unconditionally without changing any answer.
    #
    @testset "flat grid: row 3 is zero" begin
        path = joinpath(ROOT, "problems/ShallowWater/SoliWaveIsland/SoliWaveIsland.msh")
        if isfile(path)
            mesh = build_mesh("problems/ShallowWater/SoliWaveIsland", "SoliWaveIsland.msh")
            @test all(iszero, mesh.coords[3, :])
            @test mesh.lmanifold == false
        end
    end
end

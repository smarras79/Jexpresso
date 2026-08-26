#=============================================================================
 test_stretching.jl -- the vertical stretch must reach BOTH copies of the
 vertical coordinate, and must actually stretch.

 WHY THIS EXISTS. The mesh stores node positions twice: as mesh.x/y/z and as
 the mesh.coords matrix. stretch_mesh_3D! wrote only mesh.coords[3,:], so a
 deck with :lstretch => true produced

   * stretched  coords[3,:]  -- what build_column_topology reads (columns.jl)
                               and what one of the two 3D VTK writers emits
   * uniform    mesh.z       -- what the 3D METRIC TERMS are built from
                               (metric_terms.jl:361/485/545), i.e. what the
                               discretisation actually used

with nothing anywhere reporting a problem. The run completed, on a grid nobody
asked for. Asserting "the stretch happened" alone would NOT have caught it --
it did happen, in one of the two arrays. So this checks agreement as well.

     julia --project=<env> test/mesh/test_stretching.jl
=============================================================================#

using Test, MPI, Printf
MPI.Initialized() || MPI.Init()

# stretching.jl calls get_mpi_comm() for rank-0 printing only. Supply it rather
# than loading the whole package: this test is about the arithmetic and the two
# arrays, and pulling in Jexpresso would make it a two-minute test of nothing
# extra.
get_mpi_comm() = MPI.COMM_WORLD

include(joinpath(@__DIR__, "..", "..", "src", "kernel", "mesh", "stretching.jl"))

# The smallest thing the stretcher will accept: it reads zmax, coords, z, and
# for the two-block types connijk/ngl/nelem.
mutable struct MiniMesh
    npoin::Int; nelem::Int; ngl::Int
    zmax::Float64; ymax::Float64
    coords::Matrix{Float64}
    x::Vector{Float64}; y::Vector{Float64}; z::Vector{Float64}
    connijk::Array{Int,4}
end

# One column of NZ elements, with the LGL node distribution INSIDE each element.
# Uniform levels would misstate the one quantity the stretch formula keys on:
# sigma_1, the smallest positive z, which for nop=4 is 0.3453 of the half
# element and not h/4.
const LGL5 = [-1.0, -sqrt(3/7), 0.0, sqrt(3/7), 1.0]     # nop = 4
function build(nz_el, ngl, ztop)
    @assert ngl == 5 "LGL5 is the nop=4 node set"
    h    = ztop / nz_el
    zs   = Float64[]
    for e = 1:nz_el, k = 1:ngl
        e > 1 && k == 1 && continue                      # shared face node
        push!(zs, (e - 1) * h + (LGL5[k] + 1) * h / 2)
    end
    nlev  = length(zs)
    npoin = nlev
    coords = zeros(3, npoin); coords[3, :] .= zs
    conn = zeros(Int, nz_el, ngl, ngl, ngl)
    for e = 1:nz_el, k = 1:ngl
        conn[e, 1, 1, k] = (e - 1) * (ngl - 1) + k
    end
    MiniMesh(npoin, nz_el, ngl, ztop, ztop, coords,
             zeros(npoin), zeros(npoin), copy(zs), conn)
end

# The deck's own settings, verbatim.
inputs = Dict{Symbol,Any}(
    :stretch_factor       => 1.15,
    :stretch_type         => "fixed_first_twoblocks_strong",
    :first_zelement_size  => 10.0,
    :zlevel_transition    => 2000.0)

const NZ, NGL, ZTOP = 60, 5, 5000.0
m = build(NZ, NGL, ZTOP)
z_before = copy(m.coords[3, :])

stretch_mesh_3D!(m, inputs, m.npoin)

@testset "vertical stretching" begin

    @testset "it stretches, monotonically" begin
        d = diff(m.coords[3, :])
        @printf("  spacing: first %.2f m, last %.2f m  (uniform grid was %.2f m)\n",
                d[1], d[end], minimum(diff(z_before)))
        @test m.coords[3, :] != z_before          # something happened
        @test all(>(0), d)                        # monotone: no folded cells
        # NOT asserted: the DIRECTION. With these deck values this map coarsens
        # towards the ground rather than refining, and whether that is the
        # formula's fault or the deck's is the user's call, not this test's.
        # Printed above so a change in it is visible.
    end

    @testset "BOTH copies of the vertical coordinate agree" begin
        # The bug: coords[3,:] stretched, mesh.z left uniform. The metrics are
        # built from mesh.z, so this is the assertion that matters.
        worst = maximum(abs, m.z .- m.coords[3, :])
        @printf("  max |mesh.z - coords[3,:]| = %.3e\n", worst)
        @test worst == 0.0
        @test m.z != z_before                     # z really moved, not just coords
    end

    @testset "the domain is preserved" begin
        @test m.coords[3, 1]   ≈ 0.0   atol = 1e-9
        @test m.coords[3, end] ≈ ZTOP  atol = 1e-6
    end
end

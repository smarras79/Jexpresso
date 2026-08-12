#---------------------------------------------------------------------------------
# test/test_sphere_visc.jl — the ARTIFICIAL DIFFUSION δν∇²q on the spherical
# shell (src/kernel/operators/sphere_rhs.jl).
#
#   julia --project=. test/test_sphere_visc.jl
#
# The term is the SURFACE (Laplace-Beltrami) Laplacian, assembled in weak form
# from the manifold metrics. The sharpest possible test of it is that it has the
# right eigenvalues: every homogeneous harmonic polynomial of degree l,
# restricted to a sphere of radius R, is a spherical harmonic, and
#
#   Δₛ Y_l = -l(l+1)/R² Y_l    exactly.
#
# So the assembled operator is applied to z, xy, xyz and Re((x+iy)⁴) — degrees
# 1, 2, 3 and 4 — and each one has to come back scaled by its own -l(l+1)/R².
# Getting the metric terms, the index convention of dψ, or the ξ/η pairing wrong
# breaks this immediately; a term that merely "looks like diffusion" and damps
# things would still pass a blow-up test but not this one.
#
# Then the three structural properties that make it a stabiliser rather than
# just a term: it annihilates constants, it is symmetric, and it is negative
# semi-definite (qᵀMΔq = -∫|∇ₛq|² ≤ 0, so the term can only remove energy).
#
# The fixture is the .msh that ships with the SWsphere case, so the test needs
# no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_sphere_viscosity
using PartitionedArrays, MPI
using Random, Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_sphere_visc.jl"
end

function shell_inputs(nop; kwargs...)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 6.37122e6,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename] = CASE_MSH
    inputs[:nop]           = nop
    for (k, v) in kwargs
        inputs[k] = v
    end
    return inputs
end

#
# Δₛq, as the solver sees it: the viscous element kernel accumulated over the
# grid and divided by the (diagonal) mass matrix, which is exactly what
# _sphere_rhs_kernel! does with the term.
#
function surface_laplacian(q, mesh, metrics, sp; ieq = 1)

    RHS = zeros(Float64, Int(mesh.npoin), Int(sp.neqs))
    ngl = Int(mesh.ngl)

    for iel = 1:Int(mesh.nelem)
        Jexpresso._sphere_visc_el!(RHS, q, mesh.connijk, iel, ngl,
                                   metrics.Je,
                                   metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                   metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                   metrics.dψ, metrics.ω,
                                   sp.μ, sp.gξ, sp.gη, Int(sp.neqs))
    end

    # the un-inverted form too: it is the symmetric one, and the quadratic form
    # qᵀ(∫∇ψ·∇q) is what the definiteness test needs
    return RHS[:, ieq] .* metrics.Minv, RHS[:, ieq]
end

l2(v, M) = sqrt(sum(M[ip]*v[ip]^2 for ip = 1:length(v)))


@testset "spherical shell: surface Laplacian of the artificial diffusion" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    with_mpi() do distribute

        nop    = 5
        inputs = shell_inputs(nop; :lvisc => true, :μ => 1.0, :ivisc_equations => [1])

        mesh, _  = mod_mesh_mesh_driver(inputs, 1, distribute)
        metrics  = build_sphere_metrics(mesh, inputs; verbose = false)
        sp       = build_sphere_params(mesh, metrics, inputs; neqs = 1)

        @test sp.lvisc
        @test sp.μ[1] == 1.0

        npoin = Int(mesh.npoin)
        R     = Float64(mesh.radius)
        x = mesh.coords[1, 1:npoin] ./ R
        y = mesh.coords[2, 1:npoin] ./ R
        z = mesh.coords[3, 1:npoin] ./ R
        M = metrics.M             # the L² norm below is mass-weighted

        #-------------------------------------------------------------------------
        # (1) the eigenvalues. Δₛ Y_l = -l(l+1)/R² Y_l for a harmonic polynomial
        #     of degree l restricted to the sphere.
        #-------------------------------------------------------------------------
        harmonics = (
            (1, z),                                  #  z
            (2, x .* y),                             #  xy
            (3, x .* y .* z),                        #  xyz
            (4, x.^4 .- 6 .* x.^2 .* y.^2 .+ y.^4),  #  Re((x+iy)⁴)
        )

        for (l, Y) in harmonics
            q = zeros(Float64, npoin, 1)
            q[:, 1] .= Y

            Δq, _ = surface_laplacian(q, mesh, metrics, sp)

            exact = @. -l*(l+1)/R^2 * Y
            err   = l2(Δq .- exact, M) / l2(exact, M)

            @info @sprintf("  l = %d : relative L² error of Δₛ = %.3e", l, err)
            @test err < 1.0e-4
        end

        #-------------------------------------------------------------------------
        # (2) constants are annihilated exactly (Σ_k dψ[k,i] = 0 per element,
        #     and the shell has no boundary for a term to leak out of)
        #-------------------------------------------------------------------------
        qc = ones(Float64, npoin, 1)
        Δc, _ = surface_laplacian(qc, mesh, metrics, sp)
        @test maximum(abs, Δc)*R^2 < 1.0e-10

        #-------------------------------------------------------------------------
        # (3) symmetry and negative semi-definiteness. RHS = -∫∇ₛψ·(ν∇ₛq), so
        #     q₁ᵀRHS(q₂) = -∫∇q₁·∇q₂ is symmetric, and q ᵀRHS(q) = -∫|∇q|² ≤ 0.
        #     This is why the term can only ever REMOVE energy.
        #-------------------------------------------------------------------------
        rng = MersenneTwister(20250812)
        q1  = reshape(randn(rng, npoin), npoin, 1)
        q2  = reshape(randn(rng, npoin), npoin, 1)

        _, L1 = surface_laplacian(q1, mesh, metrics, sp)
        _, L2 = surface_laplacian(q2, mesh, metrics, sp)

        a12 = sum(q1[:,1] .* L2)
        a21 = sum(q2[:,1] .* L1)
        @test abs(a12 - a21) <= 1.0e-10*max(abs(a12), abs(a21))

        @test sum(q1[:,1] .* L1) < 0.0
        @test sum(q2[:,1] .* L2) < 0.0

        #-------------------------------------------------------------------------
        # (4) the deck switches. :lvisc => false must zero ν whatever :μ says,
        #     and :ivisc_equations must select the equations it names.
        #-------------------------------------------------------------------------
        @test all(iszero, build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => false, :μ => 1.0e5), 4))

        μ = build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :μ => 1.0e5,
                             :ivisc_equations => [2,3,4]), 4)
        @test μ == [0.0, 1.0e5, 1.0e5, 1.0e5]

        # the documented default is the momentum equations, not all of them
        μd = build_sphere_viscosity(Dict{Symbol,Any}(:lvisc => true, :μ => 7.0), 4)
        @test μd == [0.0, 7.0, 7.0, 7.0]

        # a per-equation vector is taken per equation
        μv = build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :μ => [0.0, 1.0, 2.0, 3.0],
                             :ivisc_equations => 2:4), 4)
        @test μv == [0.0, 1.0, 2.0, 3.0]

        @test_throws ErrorException build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :μ => 1.0,
                             :ivisc_equations => [5]), 4)
        @test_throws ErrorException build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :μ => -1.0), 4)
    end
end

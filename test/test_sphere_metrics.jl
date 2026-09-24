#---------------------------------------------------------------------------------
# test/test_sphere_metrics.jl — the metric terms of the 2-D manifold
# (src/kernel/mesh/sphere_metrics.jl), in both of the forms it can build them.
#
#   julia --project=. test/test_sphere_metrics.jl
#
# Kopriva's curl-invariant form (J. Sci. Comput. 26(3), 301-327, 2006, Eq. (15);
# Sec. 3.2.3 / Eq. (23) of Kelly, Alves, Eckermann et al., JCP 552, 2026,
# 114683) DEGENERATES on a 2-D manifold: it fixes the two in-surface metric
# terms only up to the direction v in which the surface is extended off itself,
#
#   J a¹ = a_η × v ,   J a² = v × a_ξ ,   J = v·(a_ξ × a_η) ,
#
# and v = n̂ recovers the textbook cross-product formulas. So BOTH forms this
# file exercises are curl-invariant, and the tests below are written to hold
# that distinction straight. In particular:
#
#   * the curved-surface metric identity
#         ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂) = 0
#     is at ROUND-OFF FOR BOTH, because J n̂ is DEFINED as the curl. It is a
#     check on sphere_metrics.jl, not a ranking — asserting otherwise (by
#     comparing it against the same identity computed with J n̂ = a_ξ × a_η)
#     measures the definition of the third term and nothing else;
#
#   * what DOES separate them is the strong-form annihilation of a rigid
#     rotation, ∇ₛ·(Ω × x) = 0, which is exact only for v = n̂ — the manifold
#     analogue of free-stream preservation for the chain-rule RHS.
#
# Everything else — aⁱ·a_j = δⁱⱼ, tangency, orientation, the area — has to keep
# holding under both, and is checked for both.
#
# The fixture is the .msh that ships with the SWsphere case, so the test needs
# no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, check_sphere_metrics,
                 sphere_metrics_form
using PartitionedArrays, MPI
using Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_sphere_metrics.jl"
end

function shell_inputs(nop, form)
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
        :sphere_metrics       => form,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename]  = CASE_MSH
    inputs[:nop]            = nop
    inputs[:sphere_metrics] = form
    return inputs
end

#
# max over the grid of the relative residual of
#
#   ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂) = 0 ,
#
# normalised by the size of the curvature term itself, (2/R)J.
#
function curvature_identity_residual(mesh, metrics)

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    dψ    = metrics.dψ
    H     = 2.0/mesh.radius
    worst = 0.0

    for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        for comp = 1:3
            d = 0.0
            for k = 1:ngl
                c1 = comp == 1 ? metrics.dξdx[iel,k,j] : comp == 2 ? metrics.dξdy[iel,k,j] : metrics.dξdz[iel,k,j]
                c2 = comp == 1 ? metrics.dηdx[iel,i,k] : comp == 2 ? metrics.dηdy[iel,i,k] : metrics.dηdz[iel,i,k]
                d += dψ[k,i]*metrics.Je[iel,k,j]*c1 + dψ[k,j]*metrics.Je[iel,i,k]*c2
            end
            Jn = comp == 1 ? metrics.Jnx[iel,i,j] : comp == 2 ? metrics.Jny[iel,i,j] : metrics.Jnz[iel,i,j]
            worst = max(worst, abs(d + H*Jn)/(H*metrics.Je[iel,i,j]))
        end
    end

    return worst
end

#
# aⁱ·a_j - δⁱⱼ, and aⁱ·x̂ (tangency to the TRUE sphere, not to the discrete one),
# recomputed from the coordinates.
#
function basis_residuals(mesh, metrics)

    crd   = mesh.coords
    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    dψ    = metrics.dψ

    worst_id  = 0.0
    worst_rad = 0.0

    for iel = 1:nelem, j = 1:ngl, i = 1:ngl

        a1x = a1y = a1z = 0.0
        a2x = a2y = a2z = 0.0
        for k = 1:ngl
            ipk = mesh.connijk[iel, k, j]
            a1x += dψ[k,i]*crd[1,ipk]; a1y += dψ[k,i]*crd[2,ipk]; a1z += dψ[k,i]*crd[3,ipk]
            ipl = mesh.connijk[iel, i, k]
            a2x += dψ[k,j]*crd[1,ipl]; a2y += dψ[k,j]*crd[2,ipl]; a2z += dψ[k,j]*crd[3,ipl]
        end

        c1x, c1y, c1z = metrics.dξdx[iel,i,j], metrics.dξdy[iel,i,j], metrics.dξdz[iel,i,j]
        c2x, c2y, c2z = metrics.dηdx[iel,i,j], metrics.dηdy[iel,i,j], metrics.dηdz[iel,i,j]

        worst_id = max(worst_id,
                       abs(c1x*a1x + c1y*a1y + c1z*a1z - 1.0),
                       abs(c2x*a2x + c2y*a2y + c2z*a2z - 1.0),
                       abs(c1x*a2x + c1y*a2y + c1z*a2z),
                       abs(c2x*a1x + c2y*a1y + c2z*a1z))

        ip = mesh.connijk[iel,i,j]
        r  = sqrt(crd[1,ip]^2 + crd[2,ip]^2 + crd[3,ip]^2)
        # scaled by R so the number is a sine of an angle rather than a 1/length
        worst_rad = max(worst_rad,
                        abs(c1x*crd[1,ip] + c1y*crd[2,ip] + c1z*crd[3,ip])*mesh.radius/r,
                        abs(c2x*crd[1,ip] + c2y*crd[2,ip] + c2z*crd[3,ip])*mesh.radius/r)
    end

    return worst_id, worst_rad
end


#
# strong-form ∇ₛ·(Ω × x) for a 40 m/s solid-body rotation — exactly zero in the
# continuum on any sphere, and D_ξ(Ω × x) = Ω × a_ξ holds exactly for the nodal
# interpolant, so whatever comes out is pure metric error.
#
function rigid_rotation_divergence(mesh, metrics)

    ngl = Int(mesh.ngl); crd = mesh.coords; dψ = metrics.dψ
    Ω = (0.3, -0.7, 0.5) ./ sqrt(0.83) .* (40.0/mesh.radius)
    worst = 0.0

    for iel = 1:Int(mesh.nelem), j = 1:ngl, i = 1:ngl
        d = 0.0
        for c = 1:3
            dξ = dη = 0.0
            for k = 1:ngl
                ipk = mesh.connijk[iel,k,j]; ipl = mesh.connijk[iel,i,k]
                dξ += dψ[k,i]*(Ω[mod1(c+1,3)]*crd[mod1(c+2,3),ipk] - Ω[mod1(c+2,3)]*crd[mod1(c+1,3),ipk])
                dη += dψ[k,j]*(Ω[mod1(c+1,3)]*crd[mod1(c+2,3),ipl] - Ω[mod1(c+2,3)]*crd[mod1(c+1,3),ipl])
            end
            c1 = c == 1 ? metrics.dξdx[iel,i,j] : c == 2 ? metrics.dξdy[iel,i,j] : metrics.dξdz[iel,i,j]
            c2 = c == 1 ? metrics.dηdx[iel,i,j] : c == 2 ? metrics.dηdy[iel,i,j] : metrics.dηdz[iel,i,j]
            d += dξ*c1 + dη*c2
        end
        worst = max(worst, abs(d))
    end
    return worst
end


@testset "spherical shell: the two extension directions of the metric terms" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    with_mpi() do distribute

        @test sphere_metrics_form(Dict{Symbol,Any}())                      == :cross_product
        @test sphere_metrics_form(Dict(:sphere_metrics => "cp"))           == :cross_product
        @test sphere_metrics_form(Dict(:sphere_metrics => :normal))        == :cross_product
        @test sphere_metrics_form(Dict(:sphere_metrics => :radial))        == :radial
        @test sphere_metrics_form(Dict(:sphere_metrics => "CI"))           == :radial   # alias
        @test_throws ErrorException sphere_metrics_form(Dict(:sphere_metrics => :nonsense))

        for nop in (3, 5, 7)

            @testset "nop = $nop" begin

                res = Dict{Symbol,Float64}()
                rot = Dict{Symbol,Float64}()

                for form in (:cross_product, :radial)

                    inputs  = shell_inputs(nop, form)
                    mesh, _ = mod_mesh_mesh_driver(inputs, 1, distribute)
                    metrics = build_sphere_metrics(mesh, inputs; verbose = false)

                    @test metrics.form == form

                    #-------------------------------------------------------------
                    # (1) everything check_sphere_metrics asserts, at its own
                    #     form-dependent tolerances
                    #-------------------------------------------------------------
                    @test check_sphere_metrics(mesh, metrics; verbose = false)

                    #-------------------------------------------------------------
                    # (2) aⁱ·a_j = δⁱⱼ — exact under BOTH forms. This is the
                    #     definition of the contravariant basis and neither form
                    #     is allowed to trade it away.
                    #-------------------------------------------------------------
                    worst_id, worst_rad = basis_residuals(mesh, metrics)
                    @test worst_id < 1.0e-11

                    #-------------------------------------------------------------
                    # (3) tangency to the TRUE sphere, aⁱ·x̂ = 0. Exact for the
                    #     curl-invariant form, whose J aⁱ are cross products with
                    #     x̂ itself; only O(hᴺ) for the cross-product form, which
                    #     is tangent to the polynomial interpolant instead.
                    #-------------------------------------------------------------
                    if form === :radial
                        @test worst_rad < 1.0e-8
                    end

                    #-------------------------------------------------------------
                    # (4) the area, Σ M[ip] = 4πR²: the mass matrix must not have
                    #     moved when the metrics changed form.
                    #-------------------------------------------------------------
                    aerr = abs(sum(metrics.M) - 4π*mesh.radius^2)/(4π*mesh.radius^2)
                    @test aerr < 1.0e-6

                    #-------------------------------------------------------------
                    # (5) the curved-surface metric identity. Round-off for BOTH
                    #     forms: J n̂ is DEFINED as the curl, so this closes by
                    #     algebra whatever v is. It ranks nothing — it checks
                    #     that (†) and (‡) were coded consistently.
                    #-------------------------------------------------------------
                    res[form] = curvature_identity_residual(mesh, metrics)
                    @test res[form] < 1.0e-11

                    #-------------------------------------------------------------
                    # (6) the rigid rotation — the property that DOES separate
                    #     the two, and the strong-form analogue of free-stream
                    #     preservation on a manifold.
                    #-------------------------------------------------------------
                    rot[form] = rigid_rotation_divergence(mesh, metrics)

                    @info @sprintf("  nop=%d %-14s : identity %.3e, rigid-rot %.3e s⁻¹, aⁱ·a_j-δ %.3e, aⁱ·x̂ %.3e, area %.3e",
                                   nop, string(form), res[form], rot[form], worst_id, worst_rad, aerr)
                end

                # v = n̂ annihilates a rigid rotation exactly; v = x̂ does not,
                # but only at the level of the truncation error of the sphere.
                @test rot[:cross_product] < 1.0e-16
                @test rot[:radial] > rot[:cross_product]
                @test rot[:radial] < 1.0e-8
            end
        end
    end
end


#---------------------------------------------------------------------------------
# THE 120° CORNER, from the metrics' side.
#
# This is the regression that a corner-stretched :conformal map slipped past: the
# map's own unit tests (orthogonality, watertight seams, forward ∘ inverse) all
# passed, and the metrics it produced were still wrong by O(1) — M6 = 0.41,
# unmoved by nop from 3 to 7 or by n from 5 to 40, because forcing a 90° corner
# with a non-zero Jacobian makes the corner a cone point no polynomial can fit.
#
# There is no map that gives a 90° corner AND a non-singular Jacobian (see
# cubed_sphere_maps.jl), so :conformal must be REFUSED on a grid with a node on a
# cube corner rather than regularised. Assert the refusal, and assert that asking
# for the map the grid already carries is still a no-op that builds clean metrics.
#---------------------------------------------------------------------------------
@testset "spherical shell: :conformal is refused on a grid with corner nodes" begin

    with_mpi() do distribute

        inputs = shell_inputs(4, :radial)
        inputs[:cubed_sphere_map] = :conformal
        @test_throws ErrorException mod_mesh_mesh_driver(inputs, 1, distribute)

        # the deprecated alias must be refused identically, not slip through
        inputs = shell_inputs(4, :radial)
        inputs[:cubed_sphere_map] = :conformal_exact
        @test_throws ErrorException mod_mesh_mesh_driver(inputs, 1, distribute)

        # and the map the grid DOES carry stays a no-op with clean metrics
        inputs = shell_inputs(4, :radial)
        inputs[:cubed_sphere_map] = :equiangular
        mesh, _ = mod_mesh_mesh_driver(inputs, 1, distribute)
        metrics = build_sphere_metrics(mesh, inputs; verbose = false)
        @test check_sphere_metrics(mesh, metrics; verbose = false)
    end
end

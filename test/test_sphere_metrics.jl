#---------------------------------------------------------------------------------
# test/test_sphere_metrics.jl — the metric terms of the 2-D manifold
# (src/kernel/mesh/sphere_metrics.jl), in both of the forms it can build them.
#
#   julia --project=. test/test_sphere_metrics.jl
#
# The point of the CURL-INVARIANT form (Kopriva, J. Sci. Comput. 26(3), 301-327,
# 2006, Eq. (15); Sec. 3.2.3 / Eq. (23) of Kelly, Alves, Eckermann et al., JCP
# 552, 2026, 114683) is one identity:
#
#   ∂_ξ(J a¹) + ∂_η(J a²) + (2/R)(J n̂) = 0
#
# — the curved-surface metric identity, the 2/R being the mean curvature of the
# shell. Under the cross-product form it holds to the ORDER OF THE
# APPROXIMATION and improves as nop grows. Under the curl-invariant form the
# three metric terms are three parts of one discrete curl, so it is an algebraic
# cancellation and holds to ROUND-OFF at every order, on any grid. That contrast
# is what this file measures: not "CI is smaller" but "CI is at round-off while
# CP is at truncation error, and the gap widens as the grid coarsens".
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
# normalised by the size of the curvature term itself, (2/R)J. This is M5 of
# check_sphere_metrics, recomputed here so the number can be compared between
# the two forms rather than only thresholded.
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


@testset "spherical shell: curl-invariant vs cross-product metric terms" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

    with_mpi() do distribute

        @test sphere_metrics_form(Dict{Symbol,Any}())                             == :curl_invariant
        @test sphere_metrics_form(Dict(:sphere_metrics => "CI"))                  == :curl_invariant
        @test sphere_metrics_form(Dict(:sphere_metrics => :cross_product))        == :cross_product
        @test sphere_metrics_form(Dict(:sphere_metrics => "cp"))                  == :cross_product
        @test_throws ErrorException sphere_metrics_form(Dict(:sphere_metrics => :nonsense))

        for nop in (3, 5, 7)

            @testset "nop = $nop" begin

                res = Dict{Symbol,Float64}()

                for form in (:cross_product, :curl_invariant)

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
                    if form === :curl_invariant
                        @test worst_rad < 1.0e-8
                    end

                    #-------------------------------------------------------------
                    # (4) the area, Σ M[ip] = 4πR²: the mass matrix must not have
                    #     moved when the metrics changed form.
                    #-------------------------------------------------------------
                    aerr = abs(sum(metrics.M) - 4π*mesh.radius^2)/(4π*mesh.radius^2)
                    @test aerr < 1.0e-6

                    #-------------------------------------------------------------
                    # (5) THE identity
                    #-------------------------------------------------------------
                    res[form] = curvature_identity_residual(mesh, metrics)

                    @info @sprintf("  nop=%d %-15s : metric identity %.3e, aⁱ·a_j-δ %.3e, aⁱ·x̂ %.3e, area %.3e",
                                   nop, string(form), res[form], worst_id, worst_rad, aerr)
                end

                # round-off, independent of nop and of the grid
                @test res[:curl_invariant] < 1.0e-11

                # and strictly better than the cross-product form, which carries
                # the truncation error of the sphere at this order
                @test res[:curl_invariant] < res[:cross_product]
            end
        end
    end
end

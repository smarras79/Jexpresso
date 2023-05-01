function build_rhs(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis1, basis2, ω1, ω2, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, Δt, T)

    F      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    G      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    rhs_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem)

    #B.C.
    Fuser, Guser = user_flux(T, SD, qp, mesh)

    for iel=1:mesh.nelem_semi_inf
        for i=1:mesh.ngl
            for j=1:mesh.ngr
                ip = mesh.connijk[i,j,iel]

                F[i,j,iel] = Fuser[ip]
                G[i,j,iel] = Guser[ip]
            end
        end
    end

   # for ieq = 1:neqs
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngr

                    dFdξ = 0.0
                    dFdη = 0.0
                    dGdξ = 0.0
                    dGdη = 0.0
                    for k = 1:mesh.ngl
                        dFdξ = dFdξ + basis1.dψ[k, i]*F[k,j,iel]
                        dFdη = dFdη + basis2.dψ[k, j]*F[i,k,iel]

                        dGdξ = dGdξ + basis1.dψ[k, i]*G[k,j,iel]
                        dGdη = dGdη + basis2.dψ[k, j]*G[i,k,iel]
                    end
                    dFdx = dFdξ*metrics.dξdx[i,j,iel] + dFdη*metrics.dηdx[i,j,iel]
                    dGdy = dGdξ*metrics.dξdy[i,j,iel] + dGdη*metrics.dηdy[i,j,iel]
                    rhs_el[i, j, iel] -= ω1[i]*ω2[j]*metrics.Je[i,j,iel]*(dFdx + dGdy)
                end
            end
        end
    #end
    #show(stdout, "text/plain", el_matrices.D)

    #Build rhs_el(diffusion)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqs, basis, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)


    #DSS(rhs_el)
    RHS = DSS_rhs_laguerre(SD, rhs_el + rhs_diff_el, mesh.connijk, neqs, T)
    divive_by_mass_matrix!(RHS, M, QT,neqs)

    
    return RHS

end

function build_rhs_diff(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, basis1, basis2, ω1, ω2, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem)

    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)

    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem_semi_inf

        for j=1:mesh.ngr, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            qnel[i,j,iel,1] = qp[m,1]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngr
            ωJkl = ω[k]1*ω2[l]*metrics.Je[k, l, iel]

            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + basis1.dψ[i,k]*qnel[i,l,iel]
                dqdη = dqdη + basis2.dψ[i,l]*qnel[k,i,iel]
            end
            dqdx = dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel]
            dqdy = dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel]

            ∇ξ∇q_kl = metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy
            ∇η∇q_kl = metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy

            for i = 1:mesh.ngl
                hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                dhdξ_ik, dhdη_il = basis1.dψ[i,k], basis2.dψ[i,l]

                rhsdiffξ_el[i,l,iel] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                rhsdiffη_el[k,i,iel] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
            end
        end
    end

    return (rhsdiffξ_el*νx + rhsdiffη_el*νy)

end

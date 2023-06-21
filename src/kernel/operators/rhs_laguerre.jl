function build_rhs(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, neqs, basis1, basis2, ω1, ω2, mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, M, De, Le, time, inputs, Δt, deps, T; qnm1=zeros(1,1), qnm2=zeros(1,1))

    F      = zeros(mesh.ngl, mesh.ngr, mesh.nelem)
    G      = zeros(mesh.ngl, mesh.ngr, mesh.nelem)
    rhs_el = zeros(mesh.ngl, mesh.ngr, mesh.nelem)

    #B.C.

    for iel=1:mesh.nelem_semi_inf
        for i=1:mesh.ngl
            for j=1:mesh.ngr
                ip = mesh.connijk_lag[i,j,iel]

                F[i,j,iel,1:neqs], G[i,j,iel,1:neqs] = user_flux(T, SD, qp[ip,1:neqs], mesh; neqs=neqs,ip)
            end
        end
    end

   # for ieq = 1:neqs
        for iel=1:mesh.nelem_semi_inf
            for i=1:mesh.ngl
                for j=1:mesh.ngr

                    dFdξ = 0.0
                    dFdη = 0.0
                    dGdξ = 0.0
                    dGdη = 0.0
                    for k = 1:mesh.ngl
                        dFdξ = dFdξ + basis1.dψ[k, i]*F[k,j,iel]

                        dGdξ = dGdξ + basis1.dψ[k, i]*G[k,j,iel]
                    end
                    for k = 1:mesh.ngr
                        dFdη = dFdη + basis2.dψ[k, j]*F[i,k,iel]

                        dGdη = dGdη + basis2.dψ[k, j]*G[i,k,iel]
                    end
                    dFdx = dFdξ*metrics2.dξdx[i,j,iel] + dFdη*metrics2.dηdx[i,j,iel]
                    dGdy = dGdξ*metrics2.dξdy[i,j,iel] + dGdη*metrics2.dηdy[i,j,iel]
                    rhs_el[i, j, iel] -= ω1[i]*ω2[j]*metrics2.Je[i,j,iel]*(dFdx + dGdy)
                end
            end
        end
    #end
    #show(stdout, "text/plain", el_matrices.D)

    #Build rhs_el(diffusion)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqs, basis1, basis2, ω1, ω2, inputs[:νx], inputs[:νy], mesh, metrics1, metrics2, T)


    #DSS(rhs_el)
    RHS = DSS_rhs_laguerre(SD, rhs_el + rhs_diff_el, mesh, neqs, T)
    divive_by_mass_matrix!(RHS, M, QT,neqs)

    
    return RHS

end

function build_rhs_diff(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, basis1, basis2, ω1, ω2, νx, νy, mesh::St_mesh, metrics1::St_metrics, metrics2::St_metrics, T)

    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngr,mesh.nelem)

    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngr,mesh.nelem)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngr,mesh.nelem)

    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem_semi_inf

        for j=1:mesh.ngr, i=1:mesh.ngl
            m = mesh.connijk_lag[i,j,iel]
            qnel[i,j,iel,1] = qp[m,1]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngr
            ωJkl = ω1[k]*ω2[l]*metrics2.Je[k, l, iel]

            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + basis1.dψ[i,k]*qnel[i,l,iel]
            end
            for i = 1:mesh.ngr
                dqdη = dqdη + basis2.dψ[i,l]*qnel[k,i,iel]
            end
            dqdx = dqdξ*metrics2.dξdx[k,l,iel] + dqdη*metrics2.dηdx[k,l,iel]
            dqdy = dqdξ*metrics2.dξdy[k,l,iel] + dqdη*metrics2.dηdy[k,l,iel]
            dqdx *= νx
            dqdy *= νy
            ∇ξ∇q_kl = metrics2.dξdx[k,l,iel]*dqdx + metrics2.dξdy[k,l,iel]*dqdy
            ∇η∇q_kl = metrics2.dηdx[k,l,iel]*dqdx + metrics2.dηdy[k,l,iel]*dqdy

            for i = 1:mesh.ngl
                hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                dhdξ_ik = basis1.dψ[i,k]

                rhsdiffξ_el[i,l,iel] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
            end
            for i = 1:mesh.ngr
                hll,     hkk     =  basis2.ψ[l,l],  basis1.ψ[k,k]
                dhdη_il = basis2.dψ[i,l]

                rhsdiffη_el[k,i,iel] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
            end
        end
    end

    return (rhsdiffξ_el*νx + rhsdiffη_el*νy)

end

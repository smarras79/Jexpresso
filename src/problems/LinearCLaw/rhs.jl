include("../AbstractProblems.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/mesh/metric_terms.jl")

function build_rhs(SD::NSD_2D, QT, AP::LinearCLaw, neqs, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, T)
    
    F    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    G    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    
    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqs)
    
    dFdx = dGdy = zeros(neqs)
    dFdξ = dFdη = zeros(neqs)
    dGdξ = dGdη = zeros(neqs)
    
    c = 1.0
    
    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                
                p = qp.qn[ip,1]
                u = qp.qn[ip,2]
                v = qp.qn[ip,3]
                
                F[i,j,iel,1] = c^2*u
                F[i,j,iel,2] = p
                F[i,j,iel,3] = 0

                G[i,j,iel,1] = c^2*v
                G[i,j,iel,2] = 0
                G[i,j,iel,3] = p
                
            end
        end
        
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                
                dFdξ = dFdξ = zeros(T, neqs)
                dGdξ = dGdη = zeros(T, neqs)
                for k = 1:mesh.ngl
                    dFdξ = dFdξ[1:neqs] .+ dψ[k,i]*F[k,j,iel,1:neqs]
                    dFdη = dFdη[1:neqs] .+ dψ[k,j]*F[i,k,iel,1:neqs]
                    
                    dGdξ = dGdξ[1:neqs] .+ dψ[k,i]*G[k,j,iel,1:neqs]
                    dGdη = dGdη[1:neqs] .+ dψ[k,j]*G[i,k,iel,1:neqs]
                end
                dFdx = dFdξ[1:neqs]*metrics.dξdx[i,j,iel] .+dFdη[1:neqs]*metrics.dηdx[i,j,iel]
                dGdy = dGdξ[1:neqs]*metrics.dξdy[i,j,iel] .+ dGdη[1:neqs]*metrics.dηdy[i,j,iel]
                
                rhs_el[i,j,iel,1] = -ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[1] + dGdy[1])
                rhs_el[i,j,iel,2] = -ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[2] + dGdy[2])
                rhs_el[i,j,iel,3] = -ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[3] + dGdy[3])
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)

    return rhs_el
end

function build_rhs_diff(SD::NSD_2D, QT, AP::LinearCLaw, neqs, qp, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

    N = mesh.ngl - 1
    
    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem,3)
    
    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            
            qnel[i,j,iel,1] = qp.qn[m,1]
            qnel[i,j,iel,2] = qp.qn[m,2]
            qnel[i,j,iel,3] = qp.qn[m,3]
        end
        
        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
            
            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + dψ[i,k]*qnel[i,l,iel,1]
                dqdη = dqdη + dψ[i,l]*qnel[k,i,iel,1]
            end
            dqdx = dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel]
            dqdy = dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel]
            
            ∇ξ∇q_kl = metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy
            ∇η∇q_kl = metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy
            
            for i = 1:mesh.ngl
                
                hll,     hkk     =  ψ[l,l],  ψ[k,k]
                dhdξ_ik, dhdη_il = dψ[i,k], dψ[i,l]
                
                rhsdiffξ_el[i,l,iel] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                rhsdiffη_el[k,i,iel] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
            end
        end
    end

    return (rhsdiffξ_el*νx + rhsdiffη_el*νy)
    
end


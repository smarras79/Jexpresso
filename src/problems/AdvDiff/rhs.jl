include("../AbstractProblems.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/mesh/metric_terms.jl")
include("./user_flux.jl")
include("./user_source.jl")

function build_rhs(SD::NSD_1D, QT, AP::AdvDiff, neqns, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, T)
    
    Fuser  = user_flux(SD, mesh.npoin, qp, T)
    F      = zeros(mesh.ngl, mesh.nelem)
    rhs_el = zeros(mesh.ngl, mesh.nelem)
    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.connijk[i,iel]            
            F[i, iel] = Fuser[ip]
        end
    end

    for iel=1:mesh.nelem
        dξdx = 2.0/mesh.Δx[iel]
        for i=1:mesh.ngl
            
            dFdξ = 0.0
            for k = 1:mesh.ngl
                dFdξ = dFdξ + dψ[k, i]*F[k,iel]
            end
            dFdx = dFdξ*dξdx
            
            rhs_el[i, iel] = -ω[i]*dFdx
        end
    end
    
    return rhs_el
end


function build_rhs(SD::NSD_2D, QT, AP::AdvDiff, neqns, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, T)

    Fuser, Guser = user_flux(SD, qp, mesh.npoin, T)
    F      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    G      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    rhs_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                
                F[i,j,iel] = Fuser[ip] #qp[ip,1]*qp[ip,2]
                G[i,j,iel] = Guser[ip] #qp[ip,1]*qp[ip,3]
                
            end
        end
    end
    
   # for ieq = 1:neqns
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    
                    dFdξ = dFdη = 0.0
                    dGdξ = dGdη = 0.0
                    for k = 1:mesh.ngl
                        dFdξ = dFdξ + dψ[k, i]*F[k,j,iel]
                        dFdη = dFdη + dψ[k, j]*F[i,k,iel]

                        dGdξ = dGdξ + dψ[k, i]*G[k,j,iel]
                        dGdη = dGdη + dψ[k, j]*G[i,k,iel]
                    end
                    dFdx = dFdξ*metrics.dξdx[i,j,iel] + dFdη*metrics.dηdx[i,j,iel]
                    dGdy = dGdξ*metrics.dξdy[i,j,iel] + dGdη*metrics.dηdy[i,j,iel]
                    
                    rhs_el[i, j, iel] = -ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx + dGdy)
                end
            end
        end
    #end
    #show(stdout, "text/plain", el_matrices.D)

    return rhs_el
end

function build_rhs_diff(SD::NSD_1D, QT, AP::AdvDiff, nvars, qp, ψ, dψ, ω, νx, _, mesh::St_mesh, metrics::St_metrics, T)

    N           = mesh.ngl - 1
    qnel        = zeros(mesh.ngl, mesh.nelem)
    rhsdiffξ_el = zeros(mesh.ngl, mesh.nelem)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2.0
        dξdx = 2.0/mesh.Δx[iel]
        
        for i=1:mesh.ngl
            qnel[i,iel,1] = qp[mesh.connijk[i,iel], 1]
        end
        
        for k = 1:mesh.ngl
            ωJk = ω[k]*Jac
            
            dqdξ = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + dψ[i,k]*qnel[i,iel]
            end
            dqdx = dqdξ*dξdx            
            ∇ξ∇q = dξdx*dqdx
            
            for i = 1:mesh.ngl
                hll     =  ψ[k,k]
                dhdξ_ik = dψ[i,k]
                
                rhsdiffξ_el[i, iel] -= ωJk * dψ[i,k] * ψ[k,k]*∇ξ∇q
            end
        end
    end
    #rhsdiffξ_el = zeros(mesh.ngl, mesh.nelem)

    return rhsdiffξ_el
end

function build_rhs_diff(SD::NSD_2D, QT, AP::AdvDiff, nvars, qp, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

    N = mesh.ngl - 1
    
    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    
    rhsdiffξ_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    rhsdiffη_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    
    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]            
            qnel[i,j,iel,1] = qp[m,1]
        end
        
        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
            
            dqdξ = 0.0
            dqdη = 0.0
            for i = 1:mesh.ngl
                dqdξ = dqdξ + dψ[i,k]*qnel[i,l,iel]
                dqdη = dqdη + dψ[i,l]*qnel[k,i,iel]
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

function build_rhs_diff_MxV(SD::NSD_2D, QT, AP::AdvDiff, nvars, qn, L, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

    nothing
    
end


function build_rhs(SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, D, f)

    #
    # Linear RHS in flux form: f = u*q
    #  
    rhs = zeros(mesh.npoin)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            fe[i] = f[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                rhs[I] = rhs[I] - D[i,j,iel]*fe[j]
            end
        end
    end
    # M⁻¹*rhs where M is diagonal
    rhs .= rhs./M
    
    return rhs
end

function build_rhs(SD::NSD_1D, QT::Exact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, el_mat, f)

    #
    # Linear RHS in flux form: f = u*q
    #
    
    rhs = zeros(mesh.npoin)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            fe[i] = f[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                rhs[I] = rhs[I] - el_mat.D[i,j,iel]*fe[j]
            end
        end
    end
    
    # M⁻¹*rhs where M is not full
    rhs = M\rhs
    
    return rhs
end

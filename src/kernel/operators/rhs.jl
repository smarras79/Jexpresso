include("../abstractTypes.jl")
include("../mesh/mesh.jl")
include("../mesh/metric_terms.jl")
include("../../problems/AbstractProblems.jl")
#include("../../io/print_matrix.jl")

#---------------------------------------------------------------------------
# Fetch problem name to access the user_rhs functions
#---------------------------------------------------------------------------
problem_name = ARGS[1]
user_flux_dir = string("../../problems/", problem_name, "/user_flux.jl")
user_source_dir = string("../../problems/", problem_name, "/user_source.jl")
include(user_flux_dir)
include(user_source_dir)
#---------------------------------------------------------------------------
function rhs!(du, u, params, time)

    #SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, De, u)
    T       = params.T
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    De      = params.De
    Le      = params.Le

    RHS = build_rhs(SD, QT, PT, u, neqns, basis.ψ, basis.dψ, ω, mesh, metrics, M, De, Le, time, inputs, T)    
    du = copy(RHS)
    
    return du #This is already DSSed
end



function build_rhs(SD::NSD_1D, QT::Inexact, PT::AdvDiff, qp::Array, neqns, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, T)

    Fuser = user_flux(T, SD, qp, mesh)
    
    #
    # Linear RHS in flux form: f = u*u
    #  
    RHS = zeros(mesh.npoin)
    qe  = zeros(mesh.ngl)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            qe[i] = qp[I,1]
            fe[i] = Fuser[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                RHS[I] = RHS[I] - De[i,j,iel]*fe[j] + inputs[:νx]*Le[i,j,iel]*qe[j]
            end
        end
    end
    
    # M⁻¹*rhs where M is diagonal
    RHS .= RHS./M

    #apply_periodicity!(SD, ~, qp, mesh, inputs,  QT, ~, ~, ~, ω, time, ~, ~)
    
    return RHS
    
end

function build_rhs(SD::NSD_1D, QT::Exact, PT::AdvDiff, qp::Array, neqns, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, T) nothing end


function build_rhs(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, neqns, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, T)
    
    F      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    G      = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    rhs_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem)
    
    #B.C.
    apply_boundary_conditions!(SD, rhs_el, qp, mesh, inputs, QT, metrics, ψ, dψ, ω, time, neqns)   
    Fuser, Guser = user_flux(T, SD, qp, mesh)
    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                
                F[i,j,iel] = 0.8*qp[ip]#Fuser[ip]
                G[i,j,iel] = 0.8*qp[ip]#Guser[ip]
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

    #Build rhs_el(diffusion)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqns, ψ, dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    
    
    #DSS(rhs_el)
    RHS = DSS_rhs(SD, rhs_el + rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin,neqns, mesh.nop, T)
    divive_by_mass_matrix!(RHS, M, QT)

    
    return RHS
    
end

function build_rhs(SD::NSD_2D, QT::Exact, PT::AdvDiff, qp::Array, neqns, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, T) nothing end


function build_rhs(SD::NSD_2D, QT::Inexact, PT::LinearCLaw, qp::Array, neqns, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics, M, De, Le, time, inputs, T)    

    F    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqns)
    G    = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqns)

    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqns)
    #B.C.
    apply_boundary_conditions!(SD, rhs_el, qp, mesh, inputs, QT, metrics, ψ, dψ, ω, time, neqns)
    Fuser, Guser = user_flux(T, SD, qp, mesh)
    dFdx = dGdy = zeros(neqns)
    dFdξ = dFdη = zeros(neqns)
    dGdξ = dGdη = zeros(neqns)

    c = 1.0
    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]

                p = qp[ip,1]
                u = qp[ip,2]
                v = qp[ip,3]

                F[i,j,iel,1] = Fuser[ip,1]
                F[i,j,iel,2] = Fuser[ip,2]
                F[i,j,iel,3] = Fuser[ip,3]

                G[i,j,iel,1] = Guser[ip,1]
                G[i,j,iel,2] = Guser[ip,2]
                G[i,j,iel,3] = Guser[ip,3]
            end
        end

        for i=1:mesh.ngl
            for j=1:mesh.ngl

                dFdξ = dFdξ = zeros(T, neqns)
                dGdξ = dGdη = zeros(T, neqns)
                for k = 1:mesh.ngl
                    dFdξ[1:neqns] = dFdξ[1:neqns] .+ dψ[k,i]*F[k,j,iel,1:neqns]
                    dFdη[1:neqns] = dFdη[1:neqns] .+ dψ[k,j]*F[i,k,iel,1:neqns]

                    dGdξ[1:neqns] = dGdξ[1:neqns] .+ dψ[k,i]*G[k,j,iel,1:neqns]
                    dGdη[1:neqns] = dGdη[1:neqns] .+ dψ[k,j]*G[i,k,iel,1:neqns]
                end
                dFdx = dFdξ[1:neqns]*metrics.dξdx[i,j,iel] .+ dFdη[1:neqns]*metrics.dηdx[i,j,iel]
                 dGdy = dGdξ[1:neqns]*metrics.dξdy[i,j,iel] .+ dGdη[1:neqns]*metrics.dηdy[i,j,iel]

                rhs_el[i,j,iel,1:neqns] = -ω[i]*ω[j]*metrics.Je[i,j,iel]*(dFdx[1:neqns] + dGdy[1:neqns])
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, qp,  neqns, ψ, dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)

    #DSS(rhs_el)
    RHS = DSS_rhs(SD, rhs_el + rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, neqns, mesh.nop, T)
    divive_by_mass_matrix!(RHS, M, QT)

    return rhs_el
end

function build_rhs_diff(SD::NSD_1D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

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
      
    return rhsdiffξ_el*νx
end

function build_rhs_diff(SD::NSD_2D, QT::Inexact, PT::AdvDiff, qp::Array, nvars, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

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

function build_rhs_diff(SD::NSD_2D, QT, PT::LinearCLaw, qp, neqns, ψ, dψ, ω, νx, νy, mesh::St_mesh, metrics::St_metrics, T)

    N = mesh.ngl - 1

    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem, neqns)

    rhsdiffξ_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqns)
    rhsdiffη_el = zeros(mesh.ngl, mesh.ngl, mesh.nelem, neqns)

    #
    # Add diffusion ν∫∇ψ⋅∇q (ν = const for now)
    #
    for iel=1:mesh.nelem

        for j=1:mesh.ngl, i=1:mesh.ngl
            m = mesh.connijk[i,j,iel]
            qnel[i,j,iel,1:neqns] = qp[m,1:neqns]
        end

        for k = 1:mesh.ngl, l = 1:mesh.ngl
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]

            for ieq = 1:neqns
                dqdξ = 0.0
                dqdη = 0.0
                for i = 1:mesh.ngl
                    dqdξ = dqdξ + dψ[i,k]*qnel[i,l,iel,ieq]
                    dqdη = dqdη + dψ[i,l]*qnel[k,i,iel,ieq]
                end
                dqdx = dqdξ*metrics.dξdx[k,l,iel] + dqdη*metrics.dηdx[k,l,iel]
                dqdy = dqdξ*metrics.dξdy[k,l,iel] + dqdη*metrics.dηdy[k,l,iel]

                ∇ξ∇q_kl = metrics.dξdx[k,l,iel]*dqdx + metrics.dξdy[k,l,iel]*dqdy
                ∇η∇q_kl = metrics.dηdx[k,l,iel]*dqdx + metrics.dηdy[k,l,iel]*dqdy

                for i = 1:mesh.ngl

                    hll,     hkk     =  ψ[l,l],  ψ[k,k]
                    dhdξ_ik, dhdη_il = dψ[i,k], dψ[i,l]

                    rhsdiffξ_el[i,l,iel, ieq] -= ωJkl*dhdξ_ik*hll*∇ξ∇q_kl
                    rhsdiffη_el[k,i,iel, ieq] -= ωJkl*hkk*dhdη_il*∇η∇q_kl
                end
            end
        end
     end

    return (rhsdiffξ_el*νx + rhsdiffη_el*νy)

end

function build_rhs_source(SD::NSD_2D,
                          QT::Inexact,
                          q::Array,
                          mesh::St_mesh,
                          M::AbstractArray, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M.*S    
end

function build_rhs_source(SD::NSD_2D,
                          QT::Exact,
                          q::Array,
                          mesh::St_mesh,
                          M::Matrix, #M is sparse for exact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M*S   
end


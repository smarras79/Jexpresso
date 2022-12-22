using Test

include("../mesh/mesh.jl")
include("../mesh/metric_terms.jl")
include("../basis/basis_structs.jl")
include("../Problems/AbstractProblems.jl")

abstract type AbstractIntegrationType end
struct Inexact <: AbstractIntegrationType end

abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

function build_rhs(SD::NSD_2D, QT::Inexact, AP::Adv2D, q, ψ, dψdξ, ω, mesh::St_mesh, metrics::St_metrics, M, f)

    N   = mesh.nop
    ngl = N+1
    
    rhs_el = zeros(ngl*ngl,mesh.nelem)
    for iel=1:mesh.nelem
        for i=1:N+1
            for j=1:N+1
                ip_el = i + 1 + j*(N + 1)
                dqdξ = 0
                dqdη = 0
                for k = 1:N+1
                    dqdξ = dqdξ + dψdξ[k,i]*q[k,j,iel]
                    dqdη = dqdη + dψdη[k,j]*q[i,k,iel]
                end
                dqdx = dqdξ*dξdx[i,j] + dqdη*dηdx[i,j]
                dqdy = dqdξ*dξdy[i,j] + dqdη*dηdy[i,j]

                rhe_el[ip_el] = ω[i]*ω[j]*
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
    
end



function build_rhs(SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, el_mat, f)

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
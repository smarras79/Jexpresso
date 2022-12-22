using Test

include("../mesh/mesh.jl")
include("../mesh/metric_terms.jl")
include("../basis/basis_structs.jl")

abstract type AbstractIntegrationType end
struct Inexact <: AbstractIntegrationType end

abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

abstract type AbstractProblem end
struct Wave1D <: AbstractProblem end
struct AD1D <: AbstractProblem end
struct NS1D <: AbstractProblem end
struct Burgers1D <: AbstractProblem end
struct Adv2D <: AbstractProblem end
struct Heat2D <: AbstractProblem end

function rhs(SD::NSD_2D, AP::Adv2D, QT::Inexact, q, ψ, dψdξ, ω, mesh, N, Q, T)

    rhs = zeros(mesh.npoin);    
    for iel=1:mesh.nelem
        for i=1:N+1
            for j=1:N+1
                ip = i + 1 + j*(N + 1)
                dqdξ = 0
                dqdη = 0
                for k = 1:N+1
                    dqdξ = dqdξ + dψdξ[k,i]*q[k,j,iel]
                    dqdη = dqdη + dψdη[k,j]*q[i,k,iel]
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
    
end

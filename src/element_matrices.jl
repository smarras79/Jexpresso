using Test
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using SparseArrays

include("./mesh/mod_mesh.jl")
include("./basis/basis_structs.jl")

abstract type AbstractIntegrationType end
struct Exact <: AbstractIntegrationType end
struct Inexact <: AbstractIntegrationType end


abstract type AbstractMassType end
mutable struct St_ElMat_Exact{TFloat} <: AbstractMassType
    M::Array{TFloat, 3} #Full mass matrix with exact integration
    D::Array{TFloat, 3} #Full differentiation matrix with exact integration
end

mutable struct St_ElMat_Inexact{TFloat} <: AbstractMassType
    M::Array{TFloat, 2} #Diagonal mass matrix for inexact integration
    D::Array{TFloat, 3} #Sparse differentiation matrix also for inexact int.
end
    
#Exact mass matrix

function build_element_matrices!(TP::Exact, ψ, dψdξ, ω, nelem, N, Q, T)

    el_matrices = St_ElMat_Exact{T}(zeros(N+1, N+1, nelem),
                                    zeros(N+1, N+1, nelem))
    
    for iel=1:nelem
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    el_matrices.M[i,j,iel] = el_matrices.M[i,j,iel] + (1.0/8.0)*ω[k]*ψ[i,k]*ψ[j,k]
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] + ω[k]*ψ[i,k]*dψdξ[j,k]
                end
            end
        end
    end
    #@show el_matrices.D
    
    return el_matrices
end

function build_element_matrices!(TP::Inexact, ψ, dψdξ, ω, nelem, N, Q, T)
    
    el_matrices = St_ElMat_Inexact{T}(zeros(N+1, nelem),
                                      zeros(N+1, N+1, nelem))

    for iel=1:nelem
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] + ω[k]*ψ[i,k]*dψdξ[j,k] #Sparse
                    if (i == j)
                        el_matrices.M[i,iel]   = el_matrices.M[i,iel] + 0.5*ω[k]*ψ[i,k]*ψ[j,k] #Store only the diagonal elements
                    end
                end
            end
        end
    end
    #@show el_matrices.D
        
    return el_matrices
    
end


function DSSmatrix!(M::Matrix, Me::Matrix, conn, nelem, npoin, N, T)

    M = zeros(npoin+1, npoin+1)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            for j=1:N+1
                J = conn[j,iel]                
                M[I,J] = M[I,J] + Me[i,j,iel]                
            end
        end
    end
    
end


function DSSarray!(A::Array, Ae::Array, conn, nelem, npoin, N, T)

    A = zeros(npoin+1)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            A[I] = A[I] + Ae[i,iel]
        end
    end

end


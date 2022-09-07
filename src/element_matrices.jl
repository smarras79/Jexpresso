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
    
function build_element_matrices!(QT::Exact, ψ, dψdξ, ω, mesh, N, Q, T)

    el_matrices = St_ElMat_Exact{T}(zeros(N+1, N+1, mesh.nelem),
                                    zeros(N+1, N+1, mesh.nelem))
    
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    el_matrices.M[i,j,iel] = el_matrices.M[i,j,iel] + Jac*ω[k]*ψ[i,k]*ψ[j,k]
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +      ω[k]*ψ[i,k]*dψdξ[j,k]
                end
            end
        end
    end
    #@show el_matrices.D
    
    return el_matrices
end

function build_element_matrices!(QT::Inexact, ψ, dψdξ, ω, mesh, N, Q, T)
    
    el_matrices = St_ElMat_Inexact{T}(zeros(N+1,      mesh.nelem),
                                      zeros(N+1, N+1, mesh.nelem))

    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1                   
                    if (i == j)
                        el_matrices.M[i,iel] = el_matrices.M[i,iel] + Jac*ω[k]*ψ[i,k]*ψ[j,k] #Store only the diagonal elements
                    end
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +      ω[k]*ψ[i,k]*dψdξ[j,k] #Sparse
                end
            end
        end
    end
    #@show el_matrices.D
        
    return el_matrices
    
end


function DSS(QT::Exact, Me::AbstractArray, conn, nelem, npoin, N, T)

    M    = zeros(npoin, npoin)
    Minv = zeros(npoin, npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            for j=1:N+1
                J = conn[j,iel]                
                M[I,J] = M[I,J] + Me[i,j,iel]                
            end
        end
    end
    Minv = inv(M)
    
    return M , Minv
end


function DSS(QT::Inexact, Ae::AbstractArray, conn, nelem, npoin, N, T)

    A = zeros(npoin)
    Ainv = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            A[I] = A[I] + Ae[i,iel]
        end
    end
    Ainv = 1.0./A
    
    return A, Ainv
end


function DSSarray(Ae::AbstractArray, conn, nelem, npoin, N, T)

    A = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            A[I] = A[I] + Ae[i,iel]
        end
    end

    return A
end


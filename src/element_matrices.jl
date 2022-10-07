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

abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

abstract type AbstractMassType end
mutable struct St_ElMat{TFloat} <: AbstractMassType
    M::Array{TFloat}
    D::Array{TFloat}
end


function build_element_matrices!(SD::NSD_1D, QT::Exact, ψ, dψdξ, ω, mesh, N, Q, T)

    el_matrices = St_ElMat{T}(zeros(N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem))
    
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for iq=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    el_matrices.M[i,j,iel] = el_matrices.M[i,j,iel] + Jac*ω[iq]*ψ[i,iq]*ψ[j,iq]
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +      ω[iq]*ψ[i,iq]*dψdξ[j,iq]
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
end

function build_element_matrices!(SD::NSD_1D, QT::Inexact, ψ, dψdξ, ω, mesh, N, Q, T)
    
    el_matrices = St_ElMat{T}(zeros(N+1,      mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem))

    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for iq=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    if (i == j)
                        el_matrices.M[i,iel] = el_matrices.M[i,iel] + Jac*ω[iq]*ψ[i,iq]*ψ[j,iq] #Store only the diagonal elements
                    end
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +      ω[iq]*ψ[i,iq]*dψdξ[j,iq] #Sparse
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
    
end


function build_element_matrices!(SD::NSD_2D, QT::Inexact, ψ, dψdξ, ω, mesh, N, Q, T)
    
    el_matrices = St_ElMat{T}(zeros(N+1, N+1, N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, N+1, N+1, mesh.nelem))

    
    for iel=1:mesh.nelem
        
        for k = 1:Q+1
            for l = 1:Q+1
                
                Jkl = 1
                for i = 1:N+1
                    for j = 1:N+1
                        ψJK = ψ[i,k]*ψ[j,l]

                        for m = 1:N+1
                            for n = 1:N+1
                                ψIK = ψ[m,k]*ψ[n,l]                                
                                el_matrices.M[i,j,m,n,iel] = el_matrices.M[i,j,m,n,iel] +  ω[k,l]*Jkl*ψIK*ψJK #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
    
end


function DSS(SD::NSD_1D, QT::Exact, Me::AbstractArray, conn, nelem, npoin, N, T)

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


function DSS(SD::NSD_1D, QT::Inexact, Ae::AbstractArray, conn, nelem, npoin, N, T)

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


function DSSarray(SD::NSD_1D, Ae::AbstractArray, conn, nelem, npoin, N, T)

    A = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            A[I] = A[I] + Ae[i,iel]
        end
    end

    return A
end


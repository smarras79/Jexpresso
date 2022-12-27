using Test
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using SparseArrays

include("../AbstractTypes.jl")
include("../mesh/mesh.jl")
include("../mesh/metric_terms.jl")
include("../basis/basis_structs.jl")


abstract type AbstractMassType end
mutable struct St_ElMat{TFloat} <: AbstractMassType
    M::Array{TFloat} #Mass
    D::Array{TFloat} #Differentiation
    L::Array{TFloat} #Laplacian
end

#
# Element matrices:
#
function build_element_matrices!(SD::NSD_1D, QT::Exact, ψ, dψdξ, ω, mesh, N, Q, T)

    el_matrices = St_ElMat{T}(zeros(N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem),
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
                              zeros(N+1, N+1, mesh.nelem),
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


function build_element_matrices!(SD::NSD_2D, QT::Inexact, MT::Monolithic, ψ, dψdξ, ω, mesh, metrics, N, Q, T)
    
    el_matrices = St_ElMat{T}(zeros(N+1, N+1, N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, N+1, N+1, mesh.nelem))

    
    for iel=1:mesh.nelem
        
        for k = 1:Q+1
            for l = 1:Q+1

                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[k, l, iel]
                
                for i = 1:N+1
                    for j = 1:N+1
                        ψJK = ψ[i,k]*ψ[j,l]

                        for m = 1:N+1
                            for n = 1:N+1
                                ψIK = ψ[m,k]*ψ[n,l]                                
                                el_matrices.M[i,j,m,n,iel] = el_matrices.M[i,j,m,n,iel] + ωkl*Jkle*ψIK*ψJK #Sparse
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

# Mass
function build_mass_matrix!(SD::NSD_2D, QT::Inexact, MT::TensorProduct, ψ, ω, mesh, metrics, N, Q, T)
    
    M = zeros(N+1, N+1, N+1, N+1, mesh.nelem)
    
    for iel=1:mesh.nelem
        
        for k = 1:Q+1
            for l = 1:Q+1

                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[k, l, iel]
                
                for i = 1:N+1
                    for j = 1:N+1
                        ψJK = ψ[i,k]*ψ[j,l]

                        for m = 1:N+1
                            for n = 1:N+1
                                ψIK = ψ[m,k]*ψ[n,l]                                
                                M[i,j,m,n,iel] = M[i,j,m,n,iel] + ωkl*Jkle*ψIK*ψJK #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return M
end

function build_mass_matrix!(SD::NSD_2D, QT::Inexact, MT::Monolithic, ψ, ω, mesh, metrics, N, Q, T)

    M = zeros((N+1)^2, mesh.nelem)
    @info size(M)

    for iel=1:mesh.nelem
        
        for i = 1:N+1
            for j = 1:N+1

                ωij  = ω[i]*ω[j]
                Jije = metrics.Je[i, j, iel]

                m = i + (j - 1)*(N + 1)
                n = m
                
                M[m,iel] = M[m,iel] + ωij*Jije #Sparse
            end
        end
    end
    #show(stdout, "text/plain", M)
    
    return M
end

# Laplace
function build_laplace_matrix!(SD::NSD_2D, QT::Inexact, dψdξ, ω, mesh, metrics, N, Q, T)
    
    L = zeros(N+1, N+1, N+1, N+1, mesh.nelem)
    
    for iel=1:mesh.nelem
        
        for k = 1:Q+1
            for l = 1:Q+1
                
                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[k, l, iel]
                
                for i = 1:N+1
                    for j = 1:N+1
                        
                        dψJKdx = dψdξ[i,k]*ψ[j,l]*metrics.dξdx[k,l,iel] + ψ[i,k]*dψdη[j,l]*metrics.dηdx[k,l,iel]
                        dψJKdy = dψdξ[i,k]*ψ[j,l]*metrics.dξdy[k,l,iel] + ψ[i,k]*dψdη[j,l]*metrics.dηdy[k,l,iel]
                        
                        for m = 1:N+1
                            for n = 1:N+1
                                
                                dψIKdx = dψdξ[m,k]*ψ[n,l]*metrics.dξdx[k,l,iel] + ψ[m,k]*dψdη[n,l]*metrics.dηdx[k,l,iel]
                                dψIKdy = dψdξ[m,k]*ψ[n,l]*metrics.dξdy[k,l,iel] + ψ[m,k]*dψdη[n,l]*metrics.dηdy[k,l,iel]
                                
                                L[i,j,m,n,iel] = L[i,j,m,n,iel] + ωkl*Jkle*(dψIKdx*dψJKdx + dψIKdy*dψJKdy)
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.L)
    
    return L
end

#
# DSS
#
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


function DSS(SD::NSD_2D, QT::Inexact, Ae::AbstractArray, conn, nelem, npoin, N, T)

    A  = zeros(npoin)
    MN = (N+1)^2
    
    for iel=1:nelem
        for i=1:MN
            I = conn[i,iel]
            A[I] = A[I] + Ae[i,iel]
        end
    end
        
    return A
end

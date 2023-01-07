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
function build_mass_matrix(SD::NSD_2D, QT::Inexact, MT::TensorProduct, ψ, ω, mesh, metrics, N, Q, T)

    MN = N + 1
    QN = Q + 1
    
    M = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    
    for iel=1:mesh.nelem
        
        for l = 1:QN
            for k = 1:QN
                
                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[k, l, iel]
                
                for j = 1:MN
                    for i = 1:MN
                        I = i + (j - 1)*(N + 1)
                        ψJK = ψ[i,k]*ψ[j,l]
                        for n = 1:MN
                            for m = 1:MN
                                J = m + (n - 1)*(N + 1)
                                ψIK = ψ[m,k]*ψ[n,l]
                                M[I,J,iel] = M[I,J,iel] + ωkl*Jkle*ψIK*ψJK #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", M)
    
    return M
end

function build_mass_matrix(SD::NSD_2D, QT::Inexact, MT::Monolithic, ψ, ω, mesh, metrics, N, Q, T)

    MN = (N+1)^2
    QN = MN
    M = zeros((N+1)^2, mesh.nelem)
    
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
function build_laplace_matrix(SD::NSD_2D, QT::Inexact, MT::TensorProduct, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    L = zeros((N+1)^2, (N+1)^2, mesh.nelem)

    show(stdout, "text/plain", ψ)
    println("\n")
    show(stdout, "text/plain", metrics.dξdx[:,:,1])
    
    println("\n")
    show(stdout, "text/plain", metrics.dξdy[:,:,1])
    
    println("\n")
    show(stdout, "text/plain", metrics.dηdx[:,:,1])
    
    println("\n")
    show(stdout, "text/plain", metrics.dηdy[:,:,1])

    #L = zeros((N+1), (N+1), N+1, N+1, mesh.nelem)
    for iel=1:1#mesh.nelem
        for l = 1:QN, k = 1:QN          
            ωkl  = ω[k]*ω[l]
            Jkle = metrics.Je[k, l, iel]
            for j = 1:MN, i = 1:MN     
                J = i + (j - 1)*(N + 1)
                
                hjl = ψ[j,l]
                hik = ψ[i,k]

                dhik_dξ = dψ[i,k]
                dhjl_dη = dψ[j,l]

                dψJK_dx = dhik_dξ*hjl*metrics.dξdx[k,l,iel] + hik*dhjl_dη*metrics.dηdx[k,i,iel]
                dψJK_dy = dhik_dξ*hjl*metrics.dξdy[k,j,iel] + hik*dhjl_dη*metrics.dηdy[k,l,iel]
                
                for n = 1:N+1, m = 1:N+1
                    I = m + (n - 1)*(N + 1)
                   
                    hnl, hmk        =  ψ[n,l],  ψ[m,k]
                    dhmk_dξ,dhnl_dη = dψ[m,k], dψ[n,l]
                    
                   
                    dψIK_dx = dhmk_dξ*hnl*metrics.dξdx[j,k,iel] + hmk*dhnl_dη*metrics.dηdx[j,k,iel]
                    dψIK_dy = dhmk_dξ*hnl*metrics.dξdy[k,j,iel] + hmk*dhnl_dη*metrics.dηdy[k,i,iel]

                    
                    #L[m,n,i,j,iel] += (dψIK_dx*dψJK_dx + dψIK_dy*dψJK_dy) #ωkl*Jkle
                    #L[I,J,iel] += (dψIK_dx*dψJK_dx) # + dψIK_dy*dψJK_dy)
                    L[I,J,iel] +=  dψIK_dy*dψJK_dy
                    if(iel < 2)
                        @info  L[I,J,iel]
                    end
                    #if (I == 1 && J == 1)
                    #    @info m, n, i, j, l, k
                    #    #@info  hnl, hmk
                    #    #@info dhmk_dξ,dhnl_dη
                    #    #@info dhik_dξ*hjl*metrics.dξdy[j,k,iel]
                    #    #@info hik*dhjl_dη*metrics.dηdy[j,k,iel]
                    #    @info dψIK_dy
                    #    @info dψJK_dy
                    #    @info  L[I,J,iel]
                    #    println("\n")
                    #end
                    
                end
            end
        end
    end
    
    #@info size(L)
    #show(stdout, "text/plain", L[:,:,1])
    
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


function DSS(SD::NSD_2D, QT::Inexact, Ae::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    A  = zeros(npoin)
    #=
    MN = (N+1)^2
    for iel=1:nelem
        for i=1:MN
            @info I = conn[i,iel]
    A[I] = A[I] + Ae[i,iel]
    end
    end
    =#
    for iel=1:nelem
        for i=1:N+1
            for j=1:N+1
                m = i + (j - 1)*(N + 1)
                
                I = conn[i,j,iel]
                #I = conn[m,iel] #this doesn't work if ngl>4! it returns a ZERO. CHECK and debug!
                if (I == 0)
                    error( "ELEMENT_MATRICES.jl ZEROOOOOO")
                end
                A[I] = A[I] + Ae[m,iel]
            end
        end
    end
    #show(stdout, "text/plain", M)
    return A
end

function DSSijk_mass(SD::NSD_2D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    M  = zeros(npoin)
    for iel=1:nelem
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[i,j,iel]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[m,n,iel]
                        
                        #M[IP,JP] = M[IP,JP] + Mel[I,J,iel] #if exact
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", M)
    return M
end

function DSSijk_laplace(SD::NSD_2D, QT::Inexact, Lel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    L  = zeros(npoin, npoin)
    
    for iel=1:nelem
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[i,j,iel]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[m,n,iel]
                        
                        L[IP,JP] = L[IP,JP] + Lel[I,J,iel] #if exact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", L)
    return L
end



function DSSijk_rhs(SD::NSD_2D, QT::Inexact, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)   
    
    V  = zeros(npoin)
    for iel = 1:nelem
        for j = 1:N+1
            for i = 1:N+1
                I = conn[i,j,iel]
                
                if (I == 0)
                    error( "ELEMENT_MATRICES.jl ZEROOOOOO")
                end
                V[I] = V[I] + Vel[i,j,iel]
            end
        end
    end
    #show(stdout, "text/plain", V)
    return V
end

function DSS_rhs(SD::NSD_2D, QT::Inexact, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)   
    
    V  = zeros(npoin)
    for iel = 1:nelem
        for i = 1:N+1
            for j = 1:N+1
                
                I   = conn[i,j,iel]
                Iel = i + (j - 1)*(N + 1)
                
                V[I] = V[I] + Vel[Iel,iel]
            end
        end
    end
    #show(stdout, "text/plain", V)
    return V
end

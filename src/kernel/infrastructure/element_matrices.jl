using Test
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using SparseArrays

include("../abstractTypes.jl")
include("../mesh/mesh.jl")
include("../mesh/metric_terms.jl")
include("../bases/basis_structs.jl")

abstract type AbstractMassType end
mutable struct St_ElMat{TFloat} <: AbstractMassType
    M::Array{TFloat} #Mass
    D::Array{TFloat} #Differentiation
    L::Array{TFloat} #Laplacian
end

#
# Element matrices:
#
function build_element_matrices(SD::NSD_1D, QT::Exact, ψ, dψdξ, ω, mesh, N, Q, T)

    el_matrices = St_ElMat{T}(zeros(N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem))
    
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for iq=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    el_matrices.M[i,j,iel] = el_matrices.M[i,j,iel] + Jac*ω[iq]*ψ[i,iq]*ψ[j,iq]
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +     ω[iq]*ψ[i,iq]*dψdξ[j,iq]
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
end

function build_element_matrices(SD::NSD_1D, QT::Inexact, ψ, dψdξ, ω, mesh, N, Q, T)
    
    el_matrices = St_ElMat{T}(zeros(N+1,      mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem),
                              zeros(N+1, N+1, mesh.nelem))

    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for i=1:N+1
            el_matrices.M[i,iel] += Jac*ω[i] #Store only the diagonal elements
            for iq=1:Q+1, j=1:N+1
                el_matrices.D[i,j,iel] += ω[iq]*ψ[i,iq]*dψdξ[j,iq] #Sparse
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return el_matrices
    
end

function build_element_matrices(SD::NSD_2D, QT, ψ, dψdξ, ω, mesh, N, Q, T)
    nothing
end

function build_element_matrices(SD::NSD_3D, QT, ψ, dψdξ, ω, mesh, N, Q, T)
    nothing
end


function build_differentiation_matrix(SD::NSD_1D, ψ, dψdξ, ω, mesh, N, Q, T)
    
    Del = zeros(N+1, N+1, mesh.nelem)

    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for i=1:N+1
            for iq=1:Q+1, j=1:N+1
                Del[i,j,iel] += ω[iq]*ψ[i,iq]*dψdξ[j,iq] #Sparse
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    
    return Del
    
end

function build_differentiation_matrix(SD::NSD_2D, ψ, dψdξ, ω, mesh, N, Q, T)
    nothing
end

function build_differentiation_matrix(SD::NSD_3D, ψ, dψdξ, ω, mesh, N, Q, T)
    nothing
end


#
# Element mass matrix
#
function build_mass_matrix(SD::NSD_1D, QT::Inexact, ψ, ω, mesh, metrics, N, Q, T)
        
    Me = zeros(T, N+1, mesh.nelem)
    
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for i=1:N+1
            Me[i,iel] += Jac*ω[i]
        end
    end
    #show(stdout, "text/plain", Me)
    
    return Me
    
end

function build_mass_matrix(SD::NSD_1D, QT::Exact, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    M = zeros((N+1), (N+1), mesh.nelem)
    
    for iel=1:mesh.nelem
        for k = 1:QN
            ωk = ω[k]*mesh.Δx[iel]
            for i = 1:MN
                ψik = ψ[i,k]
                for j = 1:MN
                    ψjk = ψ[j,k]
                    M[i,j,iel] += ωk*ψik*ψjk
                end
            end
        end
    end
    #show(stdout, "text/plain", M)
    
    return M
end

function build_mass_matrix(SD::NSD_2D, QT, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    Me = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    
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
                                Me[I,J,iel] += ωkl*Jkle*ψIK*ψJK #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", Me)
    
    return Me
end

#
# Element Laplace matrix
#
function build_laplace_matrix(SD::NSD_1D, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    L = zeros(MN, MN, mesh.nelem)
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2.0
        dξdx = 2.0/mesh.Δx[iel]
        
        for k = 1:QN
            ωJk = ω[k]*Jac
            
            for i = 1:MN, j = 1:MN

                dψik_dx = dψ[i,k]
                dψjk_dx = dψ[j,k]
                
                L[i,j,iel] -= ωJk*(dψik_dx*dψjk_dx)*dξdx
            end
        end
    end
    #@info size(L)
    #show(stdout, "text/plain", L)
    
    return L
end
#
function build_laplace_matrix(SD::NSD_2D, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    L = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    for iel=1:mesh.nelem
        for l = 1:QN, k = 1:QN
            if (metrics.Je[k, l, iel] < 0)
                @info " NEGATIVE JACOBIAN:" metrics.Je[k, l, iel] iel                
            end
            
            ωJkl = ω[k]*ω[l]
            #ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
            for j = 1:MN, i = 1:MN     
                J = i + (j - 1)*(N + 1)
                #J = mesh.connijk[i,j,iel]
                
                hjl = ψ[j,l]
                hik = ψ[i,k]

                dhik_dξ = dψ[i,k]
                dhjl_dη = dψ[j,l]
                
                dψJK_dx = dhik_dξ*hjl*metrics.dξdx[k,l,iel] + hik*dhjl_dη*metrics.dηdx[k,l,iel]
                dψJK_dy = dhik_dξ*hjl*metrics.dξdy[k,l,iel] + hik*dhjl_dη*metrics.dηdy[k,l,iel]
                
                for n = 1:N+1, m = 1:N+1
                    I = m + (n - 1)*(N + 1)
                    #I = mesh.connijk[m,n,iel]
                   
                    hnl, hmk        =  ψ[n,l],  ψ[m,k]
                    dhmk_dξ,dhnl_dη = dψ[m,k], dψ[n,l]
                    
                    dψIK_dx = dhmk_dξ*hnl*metrics.dξdx[k,l,iel] + hmk*dhnl_dη*metrics.dηdx[k,l,iel]
                    dψIK_dy = dhmk_dξ*hnl*metrics.dξdy[k,l,iel] + hmk*dhnl_dη*metrics.dηdy[k,l,iel]
                    
                    L[I,J, iel] -= ωJkl*(dψIK_dx*dψJK_dx + dψIK_dy*dψJK_dy)
                end
            end
        end
    end
    #@info size(L)
    #show(stdout, "text/plain", L)
    
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


function DSS(SD::NSD_2D, QT::Inexact, Ae::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    A  = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            for j=1:N+1
                m = i + (j - 1)*(N + 1)
                
                I = conn[i,j,iel]
                                
                A[I] = A[I] + Ae[m,iel]
            end
        end
    end
    #show(stdout, "text/plain", M)
    return A
end


function DSS_mass(SD::NSD_2D, QT::Exact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    M  = zeros(npoin, npoin)
    for iel=1:nelem
        
        #show(stdout, "text/plain", 36.0*Mel[:,:,iel])
        
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[i,j,iel]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[m,n,iel]
                        
                        M[IP,JP] = M[IP,JP] + Mel[I,J,iel] #if exact
                    end
                end
            end
        end
        #println("\n")
        #show(stdout, "text/plain", M[:,:, iel])
    end
    return M
end

function DSS_mass(SD::NSD_2D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    M  = zeros(npoin)
    for iel=1:nelem

        #show(stdout, "text/plain", 36.0*Mel[:,:,iel])

        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[i,j,iel]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[m,n,iel]
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end
        #println("\n")
        #show(stdout, "text/plain", M[:,:, iel])
    end
    return M
end

function DSS_mass(SD::NSD_1D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    
    M = zeros(npoin)
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            M[I] = M[I] + Mel[i,iel]
        end
    end
    #show(stdout, "text/plain", M)
    return M
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
                        
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", M)
    return M
end

function DSS_laplace(SD::NSD_1D, Lel::AbstractArray, mesh::St_mesh, T)

    L = zeros(mesh.npoin, mesh.npoin)    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.connijk[i,iel]
            for j=1:mesh.ngl
                J = mesh.connijl[j,iel]
                L[I,J] = L[I,J] + Le[i,j,iel]                
            end
        end
    end
        
    return L
end


function DSS_laplace(SD::NSD_2D, Lel::AbstractArray, mesh::St_mesh, T)
    
    L  = zeros(mesh.npoin, mesh.npoin)
    for iel=1:mesh.nelem
        for j = 1:mesh.ngl
            for i = 1:mesh.ngl
                J = i + (j - 1)*mesh.ngl
                JP = mesh.connijk[i,j,iel]
                for n = 1:mesh.ngl
                    for m = 1:mesh.ngl
                        I = m + (n - 1)*mesh.ngl
                        IP = mesh.connijk[m,n,iel]
                        
                        L[IP,JP] = L[IP,JP] + Lel[I,J,iel] #if exact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", L)
    return L
end

function DSS(SD::NSD_1D, QT::Inexact, Ae::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

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


function DSS_rhs(SD::NSD_1D, Ve::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    V = zeros(npoin)
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            V[I] = V[I] + Ve[i,iel]
        end
    end
    Ainv = 1.0./A
    
    return A, Ainv
end

function DSS_rhs(SD::NSD_2D, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)   
    
    V  = zeros(T, npoin)
    for iel = 1:nelem
        for j = 1:N+1
            for i = 1:N+1
                I = conn[i,j,iel]
                
                V[I] = V[I] + Vel[i,j,iel]
            end
        end
    end
    #show(stdout, "text/plain", V)
    return V
end



function divive_by_mass_matrix!(RHS::AbstractArray, M::AbstractArray, QT::Exact)
    RHS = M\RHS #M is not iagonal
end

function divive_by_mass_matrix!(RHS::AbstractArray, M::AbstractArray, QT::Inexact) 
    RHS .= RHS./M #M is diagonal (stored as a vector)
end


using Test
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using SparseArrays

include("../abstractTypes.jl")
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
                    el_matrices.D[i,j,iel] = el_matrices.D[i,j,iel] +     ω[iq]*ψ[i,iq]*dψdξ[j,iq]
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

#
# Element mass matrix
# 
function build_mass_matrix!(SD::NSD_2D, MT::TensorProduct, ψ, ω, mesh, metrics, N, Q, T)
    
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

function build_mass_matrix_absorbing!(SD::NSD_2D, MT::TensorProduct, ψ, ω, ψGR, ωGR, mesh, metrics, N, Q, NGR, QGR, T;dir = "x", side = "min")

    MN = N + 1
    QN = Q + 1
    MNGR = NGR + 1
    QNGR = QGR + 1
    
    if (dir == "x")
        MN1 = MNGR
        MN2 = MN
        QN1 = QNGR
        QN2 = QN
        ψ1 = ψGR
        ψ2 = ψ
        ω1 = ωGR
        ω2 = ω
        if (side == "min")
            nelem = size(mesh.xmin_faces,2)
        else
            nelem = size(mesh.xmax_faces,2)
        end 
    else
        MN2 = MNGR
        MN1 = MN
        QN2 = QNGR
        QN1 = QN
        ψ2 = ψGR
        ψ1 = ψ
        ω2 = ωGR
        ω1 = ω
        if (side == "min")
            nelem = size(mesh.ymin_faces,2)
        else
            nelem = size(mesh.ymax_faces,2)
        end
    end

    M = zeros((N+1)^2, (N+1)^2, nelem)

    for iel=1:nelem

        for l = 1:QN2
            for k = 1:QN1

                ωkl  = ω1[k]*ω2[l]
                Jkle = metrics.Je[k, l, iel]

                for j = 1:MN2
                    for i = 1:MN1
                        I = i + (j - 1)*(N + 1)
                        ψJK = ψ1[i,k]*ψ2[j,l]
                        for n = 1:MN2
                            for m = 1:MN1
                                J = m + (n - 1)*(N + 1)
                                ψIK = ψ1[m,k]*ψ2[n,l]
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

function expand_mass_matrix_laguerre!(SD::NSD_2D, QT::Inexact, mesh, M1,M2,M3,M4,M5,inputs)
    if(M2 == 0.0)
        npoin1 = 0
    else
        npoin1 = size(M2,1)
    end
    if(M3 == 0.0)
        npoin2 = 0
    else
        npoin2 = size(M3,1)
    end
    if(M4 == 0.0)
        npoin3 = 0
    else
        npoin3 = size(M4,1)
    end
    if(M5 == 0.0)
        npoin4 = 0
    else
        npoin4 = size(M5,1)
    end
    npoint = mesh.npoin+npoin1+npoin2+npoin3+npoin4
    M = zeros(npoint)
    M[1:mesh.npoin] .= M1[:]
    if (npoin1 > 0.1)
       #add contribution to in-domain (boundary) nodes
       for e =1:size(mesh.xmin_faces,2)
           i=1
           for j=1:mesh.ngl
               J = i + (j - 1)*(mesh.ngl)
               for m=1:mesh.ngl
                   I = i + (m - 1)*(mesh.ngl)
                   ip = mesh.xmin_faces[m,e]
                   M[ip] += M2[I,J,e]
               end
           end
       end
       #append new nodes
       for e =1:size(mesh.xmin_faces,2)
           for i=2:mesh.ngr
               for j=1:mesh.ngl
                   J = i + (j - 1)*(mesh.ngl)
                   for n=2:mesh.ngr
                       for j=1:mesh.ngl
                           I = n + (m - 1)*(mesh.ngl)
                           ip = mesh.npoin + mesh.conn_xminlag[i,j,e]
                           M[ip] += M2[I,J,e]
                       end
                   end
               end
           end
       end
   end
   
   if (npoin2 > 0.1)
       #add contribution to in-domain (boundary) nodes
       for e =1:size(mesh.xmax_faces,2)
           i=1
           for j=1:mesh.ngl
               J = i + (j - 1)*(mesh.ngl)
               for m=1:mesh.ngl
                   I = i + (m - 1)*(mesh.ngl)
                   ip = mesh.xmax_faces[m,e]
                   M[ip] += M3[I,J,e]
               end
           end
       end
       #append new nodes
       for e =1:size(mesh.xmax_faces,2)
           for i=2:mesh.ngr
               for j=1:mesh.ngl
                   J = i + (j - 1)*(mesh.ngl)
                   for n=2:mesh.ngr
                       for j=1:mesh.ngl
                           I = n + (m - 1)*(mesh.ngl)
                           ip = mesh.npoin + npoin1 + mesh.conn_xmaxlag[i,j,e]
                           M[ip] += M3[I,J,e]
                       end
                   end
               end
           end
       end
   end
   if (npoin3 > 0.1)
       #add contribution to in-domain (boundary) nodes
       for e =1:size(mesh.ymin_faces,2)
           j=1
           for i=1:mesh.ngl
               J = i + (j - 1)*(mesh.ngl)
               for m=1:mesh.ngl
                   I = m + (i - 1)*(mesh.ngl)
                   ip = mesh.ymin_faces[m,e]
                   M[ip] += M4[I,J,e]
               end
           end
       end
       #append new nodes
       for e =1:size(mesh.ymin_faces,2)
           for i=1:mesh.ngl
               for j=2:mesh.ngr
                   J = i + (j - 1)*(mesh.ngl)
                   for n=1:mesh.ngl
                       for j=2:mesh.ngr
                           I = n + (m - 1)*(mesh.ngl)
                           ip = mesh.npoin + npoin1 + npoin2 + mesh.conn_yminlag[i,j,e]
                           M[ip] += M4[I,J,e]
                       end
                   end
               end
           end
       end
   end
   
   if (npoin4 > 0.1)
       #add contribution to in-domain (boundary) nodes
       for e =1:size(mesh.ymax_faces,2)
           j=1
           for i=1:mesh.ngl
               J = i + (j - 1)*(mesh.ngl)
               for m=1:mesh.ngl
                   I = m + (i - 1)*(mesh.ngl)
                   ip = mesh.ymax_faces[m,e]
                   M[ip] += M5[I,J,e]
               end
           end
       end
       #append new nodes
       for e =1:size(mesh.ymax_faces,2)
           for i=1:mesh.ngl
               for j=2:mesh.ngr
                   J = i + (j - 1)*(mesh.ngl)
                   for n=1:mesh.ngl
                       for j=2:mesh.ngr
                           I = n + (m - 1)*(mesh.ngl)
                           ip = mesh.npoin + npoin1 + npoin2 + npoin3 + mesh.conn_ymaxlag[i,j,e]
                           M[ip] += M5[I,J,e]
                       end
                   end
               end
           end
       end
   end 
               
   return M
end   
            

function build_mass_matrix!(SD::NSD_2D, QT::Inexact, MT::Monolithic, ψ, ω, mesh, metrics, N, Q, T)

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



#
# Element Laplace matrix
#
function build_laplace_matrix(SD::NSD_2D, MT::TensorProduct, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    L = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    for iel=1:mesh.nelem
        for l = 1:QN, k = 1:QN          
            ωJkl = ω[k]*ω[l]*metrics.Je[k, l, iel]
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
                    
                    #L[m,n,i,j,iel] += (dψIK_dx*dψJK_dx + dψIK_dy*dψJK_dy) 
                    L[I,J, iel] += ωJkl*(dψIK_dx*dψJK_dx + dψIK_dy*dψJK_dy)
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


function DSSijk_mass(SD::NSD_2D, QT::Exact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
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

function DSSijk_laplace(SD::NSD_2D, Lel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
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

function DSSijk_rhs(SD::NSD_2D, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)   
    
    V  = zeros(npoin)
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
    RHS = M\RHS #M may is sparse but not necessarily dsiagonal
end

function divive_by_mass_matrix!(RHS::AbstractArray, M::AbstractArray, QT::Inexact) 
    RHS .= RHS./M #M is diagonal (stored as a vector)
end

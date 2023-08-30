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

function build_mass_matrix!(Me, SD::NSD_2D, QT, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    for iel=1:mesh.nelem
        
        for l = 1:Q+1
            for k = 1:Q+1
                
                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[iel, k, l]
                
                for j = 1:N+1
                    for i = 1:N+1
                        I = i + (j - 1)*(N + 1)
                        ψJK = ψ[i,k]*ψ[j,l]
                        for n = 1:N+1
                            for m = 1:N+1
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
  
end
function build_mass_matrix(SD::NSD_2D, QT, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    
    Me = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    
    for iel=1:mesh.nelem
        
        for l = 1:Q+1
            for k = 1:Q+1
                
                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[iel, k, l]
                
                for j = 1:N+1
                    for i = 1:N+1
                        I = i + (j - 1)*(N + 1)
                        ψJK = ψ[i,k]*ψ[j,l]
                        for n = 1:N+1
                            for m = 1:N+1
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
    
    Le = zeros((N+1)^2, (N+1)^2, mesh.nelem)
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                
                for j = 1:N+1
                    for i = 1:N+1
                        J = i + (j - 1)*(N + 1)
                        
                        dψJK_dx = dψ[i,k]*ψ[j,l]*metrics.dξdx[iel,k,l] + ψ[i,k]*dψ[j,l]*metrics.dηdx[iel,k,l]
                        dψJK_dy = dψ[i,k]*ψ[j,l]*metrics.dξdy[iel,k,l] + ψ[i,k]*dψ[j,l]*metrics.dηdy[iel,k,l]
                        
                        for n = 1:N+1
                            for m = 1:N+1
                                I = m + (n - 1)*(N + 1)
                                
                                dψIK_dx = dψ[m,k]*ψ[n,l]*metrics.dξdx[iel,k,l] + ψ[m,k]*dψ[n,l]*metrics.dηdx[iel,k,l]
                                dψIK_dy = dψ[m,k]*ψ[n,l]*metrics.dξdy[iel,k,l] + ψ[m,k]*dψ[n,l]*metrics.dηdy[iel,k,l]
                                
                                Le[I,J, iel] += ω[k]*ω[l]*(dψIK_dx*dψJK_dx + dψIK_dy*dψJK_dy)
                            end
                        end
                    end
                end
            end
        end
    end
    
    #@info size(L)
    #show(stdout, "text/plain", L)
    
    return -Le
end

#
# DSS
#
function DSS(SD::NSD_1D, QT::Exact, Me::AbstractArray, conn, nelem, npoin, N, T)

    M    = zeros(npoin, npoin)
    Minv = zeros(npoin, npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[iel,i]
            for j=1:N+1
                J = conn[iel,j]
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
            I = conn[iel,i]
            A[I] = A[I] + Ae[iel,i]
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
                
                I = conn[iel,i,j]
                                
                A[I] = A[I] + Ae[iel,m]
            end
        end
    end
    #show(stdout, "text/plain", M)
    return A
end


function DSS_mass(SD::NSD_2D, QT::Exact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    M  = zeros(T, npoin, npoin)
    for iel=1:nelem
        
        #show(stdout, "text/plain", 36.0*Mel[:,:,iel])
        
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[iel,i,j]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[iel,m,n]
                        
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

function DSS_mass!(M, SD::NSD_2D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    for iel=1:nelem
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[iel,i,j]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[iel,m,n]
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end    
    end
end


function DSS_mass(SD::NSD_2D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)
    
    M  = zeros(T, npoin)
    for iel=1:nelem
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[iel,i,j]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[iel,m,n]
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
    M = zeros(T, npoin)
    for iel=1:nelem
        for i=1:N+1
            I = conn[iel,i]
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
                JP = conn[iel,i,j]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[iel,m,n]
                        
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", M)
    return M
end

function DSS_generic_matrix(SD::NSD_1D, Lel::AbstractArray, mesh::St_mesh, T)

    L = zeros(mesh.npoin, mesh.npoin)    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.connijk[iel,i]
            for j=1:mesh.ngl
                J = mesh.connijl[iel,j]
                L[I,J] = L[I,J] + Le[i,j,iel]                
            end
        end
    end
        
    return L
end


function DSS_generic_matrix(SD::NSD_2D, Lel::AbstractArray, mesh::St_mesh, T)
    
    L  = zeros(mesh.npoin, mesh.npoin)
    for iel=1:mesh.nelem
        for j = 1:mesh.ngl
            for i = 1:mesh.ngl
                J = i + (j - 1)*mesh.ngl
                JP = mesh.connijk[iel,i,j]
                for n = 1:mesh.ngl
                    for m = 1:mesh.ngl
                        I = m + (n - 1)*mesh.ngl
                        IP = mesh.connijk[iel,m,n]
                        
                        L[IP,JP] = L[IP,JP] + Lel[I,J,iel] #if exact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", L)
    return L
end


function DSS_laplace!(L, Lel::AbstractArray, mesh::St_mesh, T, SD::NSD_2D)
    
    for iel=1:mesh.nelem
        for j = 1:mesh.ngl
            for i = 1:mesh.ngl
                J = i + (j - 1)*mesh.ngl
                JP = mesh.connijk[iel,i,j]
                for n = 1:mesh.ngl
                    for m = 1:mesh.ngl
                        I = m + (n - 1)*mesh.ngl
                        IP = mesh.connijk[iel,m,n]
                        
                        L[IP,JP] = L[IP,JP] + Lel[I,J,iel] #if exact
                    end
                end
            end
        end
    end    
    #show(stdout, "text/plain", L)
end


function DSS(SD::NSD_1D, QT::Inexact, Ae::AbstractArray, conn::AbstractArray, nelem, npoin, N, T)

    A = zeros(npoin)
    Ainv = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[iel,i]
            A[I] = A[I] + Ae[i,iel]
        end
    end
    Ainv = 1.0./A
    
    return A, Ainv
end


function DSS_rhs(SD::NSD_1D, Ve::AbstractArray, conn::AbstractArray, nelem, npoin, neqs, N, T)

    V = zeros(npoin,neqs)
    for iel=1:nelem
        for i=1:N+1
            I = conn[iel,i]
            V[I,:] = V[I,:] + Ve[iel,i,:]
        end
    end
    
    return V
end


function DSS_rhs!(SD::NSD_1D, V::SubArray{Float64}, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, neqs, N, T)   
    
    for iel = 1:nelem
        for i = 1:N+1
            I = conn[iel,i]
            
            V[I,:] += Vel[iel,i,:]
        end
    end
    #show(stdout, "text/plain", V)
end



function newDSS_rhs!(SD::NSD_2D, du::AbstractArray, Vel::AbstractArray, conn::AbstractArray, nelem, npoin, neqs, N, T)   

    for ieq = 1:neqs
        for iel = 1:nelem
            for j = 1:N+1
                for i = 1:N+1
                    I1d = (ieq - 1)*npoin + conn[iel,i,j]
                    
                    du[I1d] += Vel[iel,i,j,ieq]
                end
            end
        end
    end
    #show(stdout, "text/plain", V)
end

function DSS_rhs!(SD::NSD_2D, V, Vel, mesh, nelem, ngl, neqs)

    for ieq = 1:neqs
        for iel = 1:nelem
            for j = 1:ngl
                for i = 1:ngl
                    I = mesh.connijk[iel,i,j]
                    V[I,ieq] += Vel[iel,i,j,ieq]
                end
            end
        end
    end
    #show(stdout, "text/plain", V)
end



function divive_by_mass_matrix!(RHS::AbstractArray, M::AbstractArray, QT::Exact)
    RHS = M\RHS #M is not iagonal
end

function divive_by_mass_matrix!(RHS::AbstractArray, M::AbstractArray, QT::Inexact, neqs)

    for i = 1:neqs
        for j = 1:length(M)
            RHS[j, i] /= M[j]
        end
    end
    
end

function matrix_wrapper(SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
                        ldss_laplace=false, ldss_differentiation=false)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false
    
    Me = zeros(TFloat, (N+1)^2, (N+1)^2, mesh.nelem)
    build_mass_matrix!(Me, SD, QT, basis.ψ, ω, mesh, metrics, N, Q, TFloat)
    M  = zeros(TFloat, mesh.npoin)
    DSS_mass!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat)
    
    Le = zeros(TFloat, 1, 1)
    L  = zeros(TFloat, 1,1)
    if lbuild_laplace_matrix
        Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics, N, Q, TFloat)
        if ldss_laplace
            L  = DSS_generic_matrix(SD, Le, mesh, TFloat)
        end
    end
     
    De = zeros(TFloat, 1, 1)
    D  = zeros(TFloat, 1,1)
    if lbuild_differentiation_matrix
        De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh,  N, Q, TFloat)
        if ldss_differentiation
            D  = DSS_generic_matrix(SD, De, mesh, TFloat)
        end
    end
    
    return (; Me, De, Le, M, D, L)
end


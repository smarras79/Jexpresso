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

function build_differentiation_matrix_Laguerre!(De, SD::NSD_2D, QT, ψ, ψ1, dψ, dψ1, ω, ω1, mesh, metrics, N, Q, T)


    for iel=1:mesh.nelem_semi_inf

        for l = 1:mesh.ngr
            for k = 1:mesh.ngl

                ωkl  = ω[k]*ω1[l]
                Jkle = metrics.Je[iel, k, l]

                for j = 1:mesh.ngr
                    for i = 1:mesh.ngl
                        J = i + (j - 1)*(mesh.ngl)
                        ψJK = ψ[i,k]*ψ1[j,l]
                        for n = 1:mesh.ngr
                            for m = 1:mesh.ngl
                                I = m + (n - 1)*(mesh.ngl)
                                dψIK_dx = dψ[i,k]*ψ1[j,l]*metrics.dξdx[iel,k,l] + ψ[i,k]*dψ1[j,l]*metrics.dηdx[iel,k,l]
                                dψIK_dy = dψ[i,k]*ψ1[j,l]*metrics.dξdy[iel,k,l] + ψ[i,k]*dψ1[j,l]*metrics.dηdy[iel,k,l]
				De[I,J,iel] += ωkl*Jkle*ψJK*(dψIK_dx+dψIK_dy) #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", Me)

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
function build_mass_matrix!(Me, SD::NSD_1D, QT::Inexact, ψ, ω, mesh, metrics, N, Q, T)
    
    for iel=1:mesh.nelem
        Jac = mesh.Δx[iel]/2
        
        for i=1:N+1
            Me[i,iel] += Jac*ω[i]
        end
    end
end

function build_mass_matrix!(Me, SD::NSD_2D, QT::Inexact, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    for iel=1:mesh.nelem
        
        for l = 1:Q+1
            for k = 1:Q+1
                
                ωkl  = ω[k]*ω[l]
                Jkle = metrics.Je[iel, k, l]
                ωJ   = ωkl*Jkle
                
                for j = 1:N+1
                    for i = 1:N+1
                        I = i + (j - 1)*(N + 1)
                        ψikjl = ψ[i,k]*ψ[j,l]
                        
                        for n = 1:N+1
                            for m = 1:N+1
                                J = m + (n - 1)*(N + 1)
                                ψmknl = ψ[m,k]*ψ[n,l]
                                
                                Me[I,J,iel] += ωJ * ψikjl * ψmknl #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    
end

function build_mass_matrix!(Me, SD::NSD_3D, QT::Inexact, ψ, ω, mesh, metrics, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    for iel=1:mesh.nelem
        
        for o = 1:Q+1
            for n = 1:Q+1
                for m = 1:Q+1
                    
                    ωmno  = ω[m]*ω[n]*ω[o]
                    Jmnoe = metrics.Je[iel, m, n, o]
                    ωJ    = ωmno*Jmnoe

                    for k = 1:N+1
                        for j = 1:N+1
                            for i = 1:N+1
                                
                                I = i + (j - 1)*(N + 1) + (k - 1)*(N + 1)*(N + 1)
                                
                                ψijk = ψ[i,m]*ψ[j,n]*ψ[k,o]

                                for r = 1:N+1
                                    for q = 1:N+1
                                        for p = 1:N+1

                                            J = p + (q - 1)*(N + 1) + (r - 1)*(N + 1)*(N + 1)

                                            ψpqr = ψ[p,m]*ψ[q,n]*ψ[r,o]
                                            Me[I, J, iel] += ωJ * ψijk * ψpqr #Sparse
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
end

@kernel function build_mass_matrix_1d_gpu!(Me, ψ, ω, Je, Q)
    ie = @index(Group, Linear)
    i = @index(Local, Linear)
    
    Me[i,ie] += Je[ie, i, 1] * ω[i]
end

@kernel function build_mass_matrix_2d_gpu!(Me, ψ, ω, Je, N, Q)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    I = i_x + (i_y - 1)*(N+1)
    for l=1:Q+1
        for k=1:Q+1
            ωkl = ω[k]*ω[l]
            Jkle = Je[ie,k,l]
            ψJK = ψ[i_x,k]*ψ[i_y,l]
            for n=1:N+1
                for m=1:N+1
                    J = m + (n-1)*(N+1)
                    ψIK = ψ[m,k]*ψ[n,l]
                    Me[I,J,ie] += ωkl*Jkle*ψIK*ψJK
                end
            end
        end
    end
end

@kernel function build_mass_matrix_3d_gpu!(Me, ψ, ω, Je, N, Q)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    i_z = il[3]
    I = i_x + (i_y - 1)*(N+1) + (i_z - 1)*(N+1)*(N+1)
    for o=1:Q+1
        for n=1:Q+1
            for m=1:Q+1 
                ωmno = ω[m]*ω[n]*ω[o]
                Jmnoe = Je[ie,m,n,o]
                ψJK = ψ[i_x,m]*ψ[i_y,n]*ψ[i_z,o]
                for r=1:N+1
                    for q=1:N+1
                        for p=1:N+1
                            J = p + (q-1)*(N+1) + (r-1)*(N+1)*(N+1)
                            ψIK = ψ[p,m]*ψ[q,n]*ψ[r,o]
                            Me[I,J,ie] += ωmno*Jmnoe*ψIK*ψJK
                        end
                    end
                end
            end
        end
    end
end


function build_mass_matrix_Laguerre!(Me, SD::NSD_1D, QT, ψ, ω, mesh, metrics, N, Q, T)

    for iel=1:mesh.nelem_semi_inf
        for i=1:mesh.ngr
            Jac = metrics.Je[iel,i,1]
            Me[i,iel] += Jac*ω[i]
        end
    end

end

function build_mass_matrix_Laguerre!(Me, SD::NSD_2D, QT, ψ, ψ1, ω, ω1, mesh, metrics, N, Q, T)


    for iel=1:mesh.nelem_semi_inf

        for l = 1:mesh.ngr
            for k = 1:mesh.ngl

                ωkl  = ω[k]*ω1[l]
                Jkle = metrics.Je[iel, k, l]

                for j = 1:mesh.ngr
                    for i = 1:mesh.ngl
                        J = i + (j - 1)*(mesh.ngl)
                        ψJK = ψ[i,k]*ψ1[j,l]
                        for n = 1:mesh.ngr
                            for m = 1:mesh.ngl
                                I = m + (n - 1)*(mesh.ngl)
                                ψIK = ψ[m,k]*ψ1[n,l]
                                Me[I,J,iel] += ωkl*Jkle*ψIK*ψJK #Sparse
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", Me)

end

@kernel function build_mass_matrix_Laguerre_2d_gpu!(Me, ψ, ψ1, ω, ω1, Je, ngl, ngr)


    iel = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    #ωkl  = ω[i_x]*ω1[i_y]
    #Jkle = Je[iel, k, l]
    I = i_x + (i_y - 1)*ngl
    for l = 1:ngr
        for k = 1:ngl
            ωkl = ω[k]*ω1[l]
            Jkle = Je[iel,k,l]
            ψJK = ψ[i_x,k]*ψ1[i_y,l]
            for n = 1:ngr
                for m = 1:ngl
                    J = m + (n - 1)*(ngl)
                    ψIK = ψ[m,k]*ψ1[n,l]
                    Me[I,J,iel] += ωkl*Jkle*ψIK*ψJK #Sparse
                end
            end
        end
    end
    #show(stdout, "text/plain", Me)

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


function build_laplace_matrix(SD::NSD_2D, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    Le = zeros((N+1),(N+1))
    
    for i=1:N+1
        for j=1:N+1
            for k=1:Q+1
                sum = ω[k]*dψ[i,k]*dψ[j,k]
                Le[i,j] = Le[i,j] + sum
            end
        end
    end 
    
    #@info size(L)
    #show(stdout, "text/plain", L)
    
    return Le
end


@kernel function build_laplace_matrix_gpu!(Le, dψ, ω, Q)
    idx = @index(Global, NTuple)
    i = idx[1]
    j = idx[2]

    for k=1:Q+1
        sum = ω[k]*dψ[i,k]*dψ[j,k]
        Le[i,j] = Le[i,j] + sum
    end

end
#
# DSS
#
function DSS_mass!(M, SD::NSD_2D, QT::Exact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T; llump=false)

    if llump == true
        for iel=1:nelem
            for j = 1:N+1
                for i = 1:N+1
                    J = i + (j - 1)*(N + 1)
                    JP = conn[iel,i,j]
                    for n = 1:N+1
                        for m = 1:N+1
                            I = m + (n - 1)*(N + 1)
                            IP = conn[iel,m,n]
                            M[IP] += Mel[I,J,iel] #if inexact
                        end
                    end
                end
            end    
        end
        
    else
        
        for iel=1:nelem    
            for j = 1:N+1
                for i = 1:N+1
                    J = i + (j - 1)*(N + 1)
                    JP = conn[iel,i,j]
                    for n = 1:N+1
                        for m = 1:N+1
                            I = m + (n - 1)*(N + 1)
                            IP = conn[iel,m,n]
                            M[IP,JP] += Mel[I,J,iel] #if exact
                        end
                    end
                end
            end
        end
    end
end

function DSS_mass_Laguerre!(M, SD::NSD_1D, Mel::AbstractArray, Mel_lag::AbstractArray, mesh, N, T; llump=false)

    for iel=1:mesh.nelem
        for i = 1:mesh.ngl
            IP = mesh.connijk[iel,i,1]
            M[IP] = M[IP] + Mel[i,iel] #if inexact
        end
    end
    #@info M[mesh.npoin_linear] 
    for iel=1:mesh.nelem_semi_inf
        for i = 1:mesh.ngr
            IP = mesh.connijk_lag[iel,i,1]
            M[IP] = M[IP] + Mel_lag[i,iel] #if inexact
        end
    end
    #@info M[mesh.npoin_linear]
end

function DSS_mass_Laguerre!(M, SD::NSD_2D, Mel::AbstractArray, Mel_lag::AbstractArray, mesh, N, T; llump=false)

    for iel=1:mesh.nelem

        #show(stdout, "text/plain", 36.0*Mel[:,:,iel])

        for j = 1:mesh.ngl
            for i = 1:mesh.ngl
                J = i + (j - 1)*(mesh.ngl)
                JP = mesh.connijk[iel,i,j]
                for n = 1:mesh.ngl
                    for m = 1:mesh.ngl
                        I = m + (n - 1)*(mesh.ngl)
                        IP = mesh.connijk[iel,m,n]
                        #@info I,J,iel,Mel[I,J,iel]
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                        #@info IP, M[IP]
                    end
                end
            end
        end
        #println("\n")
        #show(stdout, "text/plain", M[:,:, iel])
    end

    for iel=1:mesh.nelem_semi_inf

        #show(stdout, "text/plain", 36.0*Mel[:,:,iel])

        for j = 1:mesh.ngr
            for i = 1:mesh.ngl
                J = i + (j - 1)*(mesh.ngl)
                JP = mesh.connijk_lag[iel,i,j]
                for n = 1:mesh.ngr
                    for m = 1:mesh.ngl
                        I = m + (n - 1)*(mesh.ngl)
                        IP = mesh.connijk_lag[iel,m,n]
                        M[IP] = M[IP] + Mel_lag[I,J,iel] #if inexact
                    end
                end
            end
        end
        
        #println("\n")
        #show(stdout, "text/plain", M[:,:, iel])
    end
end

@kernel function DSS_mass_Laguerre_gpu_2D!(M, Mel_lag, connijk_lag, ngl, ngr)

    iel = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    J = i_x + (i_y - 1)*(ngl)
    for n = 1:ngr
        for m = 1:ngl
            I = m + (n - 1)*(ngl)
            IP = connijk_lag[iel,m,n]
            KernelAbstractions.@atomic M[IP] += Mel_lag[I,J,iel] #if inexact
        end
    end
end



function DSS_mass!(M, SD::NSD_1D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T; llump=false)
    
    for iel=1:nelem
        for i = 1:N+1
            IP = conn[iel,i,1]
            M[IP] = M[IP] + Mel[i,iel] #if inexact
        end
    end
    
end


function DSS_mass!(M, SD::NSD_2D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T; llump=false)
    
    for iel=1:nelem
        for j = 1:N+1
            for i = 1:N+1
                J = i + (j - 1)*(N + 1)
                JP = conn[iel,i,j]
                for n = 1:N+1
                    for m = 1:N+1
                        I = m + (n - 1)*(N + 1)
                        IP = conn[iel,m,n]
                        
                        M[IP] += Mel[I,J,iel] #if inexact
                    end
                end
            end
        end    
    end
end

function DSS_mass!(M, SD::NSD_3D, QT::Inexact, Mel::AbstractArray, conn::AbstractArray, nelem, npoin, N, T; llump=false)
    
    for iel=1:nelem
        
        for k = 1:N+1
            for j = 1:N+1
                for i = 1:N+1
                    J = i + (j - 1)*(N + 1) + (k - 1)*(N + 1)*(N + 1)
                    JP = conn[iel,i,j,k]

                    for o = 1:N+1
                        for n = 1:N+1
                            for m = 1:N+1
                                I = m + (n - 1)*(N + 1) + (o - 1)*(N + 1)*(N + 1)
                                IP = conn[iel,m,n,o]
                                
                                M[IP] += Mel[I,J,iel] #if inexact
                            end
                        end
                    end
                end    
            end
        end     
    end
    
end

@kernel function DSS_Mass_gpu_1D!(M, Mel, conn)
    ie = @index(Group, Linear)
    i = @index(Local, Linear)

    IP = conn[ie, i, 1]
    M[IP] = M[IP] + Mel[i,ie]
end

@kernel function DSS_Mass_gpu_2D!(M, Mel, conn, nelem, npoin, N)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    J = i_x + (i_y - 1)*(N+1)
    JP = conn[ie,i_x,i_y]

    for n = 1:N+1
        for m = 1:N+1
            I = m + (n-1)*(N+1)
            IP = conn[ie,m,n]
            KernelAbstractions.@atomic M[IP] += Mel[I,J,ie]
        end
    end
end

@kernel function DSS_Mass_gpu_3D!(M, Mel, conn, nelem, npoin, N)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    i_z = il[3]
    J = i_x + (i_y - 1)*(N+1) + (i_z - 1)*(N+1)*(N+1)
    JP = conn[ie,i_x,i_y,i_z]

    for n = 1:N+1
        for m = 1:N+1
            for k = 1:N+1
                I = k + (m-1)*(N+1) + (n-1)*(N+1)*(N+1)
                IP = conn[ie,k,m,n]
                KernelAbstractions.@atomic M[IP] += Mel[I,J,ie]
            end
        end
    end
end

function DSS_laplace!(L, Lel::AbstractArray, mesh::St_mesh, T, ::NSD_2D)
    
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

@kernel function DSS_laplace_gpu!(L, Lel, connijk, ωx, ωy, nx, ny, dξdx, dydη, dηdy, dxdξ)
    ie = @index(Group, Linear)
    idx = @index(Local, NTuple)
    i = idx[1]
    j = idx[2]
    
    ip = connijk[ie, i, j]
    for k=1:nx
        jp = connijk[ie, k, j]
        KernelAbstractions.@atomic L[ip, jp] += dξdx[ie, i, k] * Lel[i, k] * dydη[ie, i, k] * ωx[j]
    end 
    
    for l=1:ny
        jp = connijk[ie,i,l]
        KernelAbstractions.@atomic L[ip, jp] += dηdy[ie, i, l] * Lel[j, l]*dxdξ[ie, i, l] * ωy[i]
    end 
end

@kernel function DSS_laplace_gpu_lag!(L, Lel, connijk, ωx, ωy, nx, ny, dηdx_lag, dydη, dηdy, dxdη_lag)
    ie = @index(Group, Linear)
    idx = @index(Local, NTuple)
    i = idx[1]
    j = idx[2]

    ip = connijk[ie, j, i]
    for k=1:ny
        jp = connijk[ie, j, k]
        KernelAbstractions.@atomic L[ip, jp] += dηdx_lag[ie, j, k] * Lel[i, k] * dydη[ie, j, j] * ωx[j]
    end

    for l=1:nx
        jp = connijk[ie,l,i]
        KernelAbstractions.@atomic L[ip, jp] += dηdy[ie, j, l] * Lel[j, l]*dxdη_lag[ie, l, i] * ωy[i]
    end
end

function DSS_laplace!(L, SD::NSD_2D, Lel::AbstractArray, ω, mesh, metrics, N, T; llump=false)

    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[iel,i,j]
                for k =1:mesh.ngl
                    jp = mesh.connijk[iel,k,j]
                    L[ip,jp] += metrics.dξdx[iel,i,k]*Lel[i,k]*ω[j]*metrics.dydη[iel,i,k]
                end

                for l = 1:mesh.ngl
                    jp = mesh.connijk[iel,i,l]
                    L[ip,jp] += metrics.dηdy[iel,i,l]*Lel[j,l]*ω[i]*metrics.dxdξ[iel,i,l]
                end
            end
        end
    end
end

function DSS_laplace_Laguerre!(L, SD::NSD_2D, Lel::AbstractArray, Lel_lag::AbstractArray, ω, ω_lag, mesh, metrics, metrics_lag, N, T; llump=false)

    for iel=1:mesh.nelem

        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[iel,i,j]
                for k =1:mesh.ngl
                    jp = mesh.connijk[iel,k,j]
                    L[ip,jp] += metrics.dξdx[iel,i,k]*Lel[i,k]*ω[j]*metrics.dydη[iel,i,k]
                end
                
                for l = 1:mesh.ngl
                    jp = mesh.connijk[iel,i,l]
                    L[ip,jp] += metrics.dηdy[iel,i,l]*Lel[j,l]*ω[i]*metrics.dxdξ[iel,i,l]
                end
            end
        end 
    end

    for iel=1:mesh.nelem_semi_inf

        for i=1:mesh.ngr
            for j=1:mesh.ngl
                ip = mesh.connijk_lag[iel,j,i]
                for k =1:mesh.ngr
                    jp = mesh.connijk_lag[iel,j,k]
                    L[ip,jp] += metrics_lag.dηdx[iel,j,k]*Lel_lag[i,k]*ω[j]*metrics.dydη[iel,j,j]
                end
                ### this only works for a standard grid fix this in the future 
                for l = 1:mesh.ngl
                    jp = mesh.connijk_lag[iel,l,i]
                    L[ip,jp] += metrics.dηdy[iel,j,l]*Lel[j,l]*ω_lag[i]*metrics_lag.dxdη[iel,l,i]
                end
            end
        end
    end

end



function DSS_rhs!(RHS, rhs_el, mesh, nelem, ngl, neqs, ::NSD_1D, ::FD)
    nothing
end
function DSS_rhs!(RHS, rhs_el, mesh, nelem, ngl, neqs, ::NSD_2D, ::FD)
    nothing
end

function DSS_rhs!(RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_1D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem
            for i = 1:ngl
                I = connijk[iel,i,1]
                RHS[I,ieq] += rhs_el[iel,i,ieq]
            end
        end
    end
    
end

function DSS_rhs!(RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_2D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem
            for j = 1:ngl
                for i = 1:ngl
                    #I = Ref{Int64}(mesh.connijk[iel,i,j])
                    #RHS[I[],ieq] += Ref{Float64}(rhs_el[iel,i,j,ieq])[]
                    I = connijk[iel,i,j]
                    RHS[I,ieq] += rhs_el[iel,i,j,ieq]
                end
            end
        end
    end
    #show(stdout, "text/plain", V)
end


function DSS_rhs!(RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_3D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem
            for k = 1:ngl
                for j = 1:ngl
                    for i = 1:ngl
                        I = connijk[iel,i,j,k]
                        RHS[I,ieq] += rhs_el[iel,i,j,k,ieq]
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", V)
end

function DSS_rhs_laguerre!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, ::NSD_1D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem_semi_inf
            for i = 1:ngr
                I = connijk_lag[iel,i,1]

                RHS[I,ieq] += rhs_el[iel,i,ieq]
            end
        end
    end
end

function DSS_rhs_laguerre!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, ::NSD_2D, ::ContGal)

    for ieq = 1:neqs
        for iel = 1:nelem_semi_inf
            for j = 1:ngr
                for i = 1:ngl
                    I = connijk_lag[iel,i,j]
                    
                    RHS[I,ieq] += rhs_el[iel,i,j,ieq]
                end
            end
        end
    end
end


function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin, ::FD)
    nothing
end



function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, ::FD)
    nothing
end

function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin, ::ContGal)
    
    RHSaux .= RHS
    for ip=1:npoin
        a = zero(eltype(RHS))
        for jp = 1:npoin
            a += Minv[ip,jp]*RHSaux[jp]
        end
        RHS[ip] = a
    end
    
end

function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, ::ContGal)
    
    for ip=1:npoin
        RHS[ip] = Minv[ip]*RHS[ip]
    end
    
end

function matrix_wrapper(::FD, SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
                        ldss_laplace=false, ldss_differentiation=false)

    if typeof(SD) == NSD_1D
        Me = zeros(TFloat, 1, 1)
    elseif typeof(SD) == NSD_2D
        Me = zeros(TFloat, 1, 1, 1)
    elseif typeof(SD) == NSD_3D
        Me = zeros(TFloat, 1, 1, 1, 1)
    end
    
    if (QT == Exact() && inputs[:llump] == false)
        M    = zeros(TFloat, 1, 1)
        Minv = zeros(TFloat, 1, 1)
    else
        M    = zeros(TFloat, 1)
        Minv = zeros(TFloat, 1)
    end
    
    Le = zeros(TFloat, 1, 1)
    L  = zeros(TFloat, 1,1)
    
    De = zeros(TFloat, 1, 1)
    D  = zeros(TFloat, 1,1)
    
    return (; Me, De, Le, M, Minv, D, L)
    
end

function DSS_global_RHS!(RHS, ip2gip, gip2owner, parts, npoin, gnpoin, neqs)
    for i = 1:neqs
        DSS_global_mass!(@view(RHS[:,i]), ip2gip, gip2owner, parts, npoin, gnpoin)
    end
end

function DSS_global_mass!(M, ip2gip, gip2owner, parts, npoin, gnpoin)
    # @info ip2gip
    row_partition = map(parts) do part
        row_partition = LocalIndices(gnpoin,part,ip2gip,gip2owner)
        # gM = M
        row_partition
    end
    pM = pvector(values->@view(M[:]), row_partition)
    # map(parts,local_values(pM)) do part,values
    #     # if part == 1
    #         @info values
    # #     end
    # end

    # map(partition(pM),row_partition) do values, indices
    #     local_index_to_owner = local_to_owner(indices)
    #     @info local_index_to_owner
    #     # for lid in 1:length(local_index_to_owner)
    #     #     owner = local_index_to_owner[lid]
    #     #     @test values[lid] == 10*owner
    #     # end
    # end


    assemble!(pM) |> wait
    consistent!(pM) |> wait
    M = map(local_values(pM)) do values
        M = values
        M
    end
end

function matrix_wrapper(::ContGal, SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
        ldss_laplace=false, ldss_differentiation=false, backend = CPU(), interp)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false

    if typeof(SD) == NSD_1D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, Int64(mesh.nelem))
    elseif typeof(SD) == NSD_2D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, (N+1)^2, Int64(mesh.nelem))
    elseif typeof(SD) == NSD_3D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^3, (N+1)^3, Int64(mesh.nelem))
    end
    if (backend == CPU())
        build_mass_matrix!(Me, SD, QT, basis.ψ, ω, mesh, metrics, N, Q, TFloat)
    elseif (SD == NSD_1D())
        k = build_mass_matrix_1d_gpu!(backend, (N+1))
        k(Me, basis.ψ, ω, metrics.Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
    elseif (SD == NSD_2D())
        k= build_mass_matrix_2d_gpu!(backend,(N+1,N+1))
        k(Me, basis.ψ, ω, metrics.Je, N, Q;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
    elseif (SD == NSD_3D())
        k= build_mass_matrix_3d_gpu!(backend,(N+1,N+1,N+1))
        k(Me, basis.ψ, ω, metrics.Je, N, Q;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl,mesh.ngl))
    end
    if (QT == Exact() && inputs[:llump] == false)
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
    else
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    end
    if (backend == CPU() || SD == NSD_1D())
        @info "start DSS_mass"
        @time DSS_nc_gather_mass!(M, mesh, SD, QT, Me, mesh.connijk, mesh.poin_in_edge,
                                mesh.non_conforming_facets, mesh.non_conforming_facets_parents_ghost,
                                mesh.ip2gip, mesh.gip2ip, mesh.pgip_ghost, mesh.pgip_owner, N, interp)
        M1 = copy(M)
        DSS_mass!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
        M2 = copy(M)
        @info "end DSS_mass"
        # @info [ [mesh.x[i], mesh.y[i], v] for (i, v) in enumerate(M) ]
    elseif (SD == NSD_2D())
        k = DSS_Mass_gpu_2D!(backend,(N+1,N+1))
        connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1, N+1)
        KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
        k(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
    elseif (SD == NSD_3D())
        k = DSS_Mass_gpu_3D!(backend,(N+1,N+1,N+1))
        connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1, N+1,N+1)
        KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
        k(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl,mesh.ngl))
    end
    DSS_global_mass!(M, mesh.ip2gip, mesh.gip2owner, mesh.parts, mesh.npoin, mesh.gnpoin)
    M3 = copy(M)
    @time DSS_nc_scatter_mass!(M, SD, QT, Me, mesh.connijk, mesh.poin_in_edge, mesh.non_conforming_facets,
    mesh.non_conforming_facets_children_ghost, mesh.ip2gip, mesh.gip2ip, mesh.cgip_ghost, mesh.cgip_owner, N, interp)
    M4 = copy(M)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    # @info rank, [(x, y, m_i, m1_i, m2_i, m3_i, m4_i-m3_i) for (x, y, m_i, m1_i, m2_i, m3_i, m4_i) in zip(mesh.x, mesh.y, M, M1, M2, M3, M4)]
    # @info rank, [(mesh.x[mesh.connijk[iel,i,j]], mesh.y[mesh.connijk[iel,i,j]], M[mesh.connijk[iel,i,j]], M1[mesh.connijk[iel,i,j]]) for j in 1: mesh.ngl, i in 1: mesh.ngl, iel in 1: mesh.nelem]
    # @mystop("my stop at DSS_global_mass!")

    mass_inverse!(Minv, M, QT)
    Le = KernelAbstractions.zeros(backend,TFloat, 1, 1)
    L  = KernelAbstractions.zeros(backend, TFloat, 1,1)
    if lbuild_laplace_matrix
        if (backend == CPU())
            Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics, N, Q, TFloat)
            L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))            
            if ldss_laplace
                DSS_laplace_Laguerre!(L, SD, Le, ω, mesh, metrics, N, TFloat; llump=inputs[:llump])
            end
        else
            Le = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.ngl), Int64(mesh.ngl))

            k = build_laplace_matrix_gpu!(backend)
            k(Le, basis.dψ, ω, TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
            KernelAbstractions.synchronize(backend)
            L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
            KernelAbstractions.synchronize(backend)
            k = DSS_laplace_gpu!(backend)
            k(L, Le, connijk, ω, ω, mesh.ngl, mesh.ngl, metrics.dξdx, metrics.dydη, metrics.dηdy, metrics.dxdξ;
              ndrange = (mesh.nelem*mesh.ngl, mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
            KernelAbstractions.synchronize(backend)
        end
    end
    
    De = KernelAbstractions.zeros(backend, TFloat, 1, 1)
    D  = KernelAbstractions.zeros(backend, TFloat, 1,1)
    if lbuild_differentiation_matrix
        De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh,  N, Q, TFloat)
        if ldss_differentiation
            D  = DSS_generic_matrix(SD, De, mesh, TFloat)
        end
    end
    
    return (; Me, De, Le, M, Minv, D, L)
end


function mass_inverse!(Minv, M::AbstractMatrix, QT, ::FD)
    nothing
end

function mass_inverse!(Minv, M::AbstractMatrix, QT, ::ContGal)
    Minv .= inv(M)
end

function matrix_wrapper_laguerre(::FD, SD, QT, basis, ω, mesh, metrics, N, Q, TFloat; ldss_laplace=false, ldss_differentiation=false)
    
    return 0
end

function matrix_wrapper_laguerre(::ContGal, SD, QT, basis, ω, mesh, metrics, N, Q, TFloat; ldss_laplace=false, ldss_differentiation=false, backend = CPU())

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false
    
    if (ldss_differentiation) lbuild_differentiation_matrix = true end
    if (ldss_laplace) lbuild_laplace_matrix = true end

    if typeof(SD) == NSD_1D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, Int64(mesh.nelem))
    elseif typeof(SD) == NSD_2D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, (N+1)^2, Int64(mesh.nelem))
    end
    if (backend == CPU())
        build_mass_matrix!(Me, SD, QT, basis[1].ψ, ω[1], mesh, metrics[1], N, Q, TFloat)
    elseif (SD == NSD_1D())
        k = build_mass_matrix_1d_gpu!(backend, (N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
    elseif (SD == NSD_2D())
        k = build_mass_matrix_2d_gpu!(backend, (N+1, N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, N, Q;ndrange = (mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
    else
        nothing
    end


    if typeof(SD) == NSD_1D
        M_lag = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.ngr*mesh.ngr), Int64(mesh.nelem_semi_inf))
    elseif typeof(SD) == NSD_2D
        M_lag = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.ngl*mesh.ngr), Int64(mesh.ngl*mesh.ngr), Int64(mesh.nelem_semi_inf))
    end
    if (QT == Exact() && inputs[:llump] == false)
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
    else
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    end
    if (backend == CPU())
        if typeof(SD) == NSD_1D
            build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[2].ψ, ω[2], mesh, metrics[2], N, Q, TFloat)
        elseif typeof(SD) == NSD_2D
            build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[1].ψ, basis[2].ψ, ω[1], ω[2], mesh, metrics[2], N, Q, TFloat)
        end
        
        @time DSS_mass_Laguerre!(M, SD, Me, M_lag, mesh, N, TFloat; llump=inputs[:llump])
    else
        if (typeof(SD) == NSD_1D)
            k = build_mass_matrix_1d_gpu!(backend)
            k(M_lag, basis[2].ψ, ω[2], metrics[2].Je, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngr), workgroupsize = (mesh.ngr))
            KernelAbstractions.synchronize(backend)

            connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1)
            KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
            k1 = DSS_Mass_gpu_1D!(backend,(N+1))
            k1(M, Me, connijk;ndrange =(mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))

            KernelAbstractions.synchronize(backend)

            connijk_lag = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem_semi_inf), Int64(mesh.ngr))
            KernelAbstractions.copyto!(backend, connijk_lag, mesh.connijk_lag)
            k2 = DSS_Mass_gpu_1D!(backend)
            k2(M, M_lag, connijk_lag;ndrange = (mesh.nelem_semi_inf*mesh.ngr), workgroupsize = (mesh.ngr))

            KernelAbstractions.synchronize(backend)
        elseif (typeof(SD) == NSD_2D)
            k = build_mass_matrix_Laguerre_2d_gpu!(backend)
            k(M_lag, basis[1].ψ, basis[2].ψ, ω[1], ω[2], metrics[2].Je, mesh.ngl, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngl,mesh.ngr), workgroupsize = (mesh.ngl,mesh.ngr))
            KernelAbstractions.synchronize(backend)
            
            connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1, N+1)
            KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
            k1 = DSS_Mass_gpu_2D!(backend,(N+1,N+1))
            k1(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
            
            KernelAbstractions.synchronize(backend)

            connijk_lag = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem_semi_inf), Int64(mesh.ngl), Int64(mesh.ngr))
            KernelAbstractions.copyto!(backend, connijk_lag, mesh.connijk_lag)
            k2 = DSS_mass_Laguerre_gpu_2D!(backend)
            k2(M, M_lag, connijk_lag, mesh.ngl, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngl,mesh.ngr), workgroupsize = (mesh.ngl,mesh.ngr))
        
            KernelAbstractions.synchronize(backend)
        end
    end

    mass_inverse!(Minv, M, QT)

    Le = KernelAbstractions.zeros(backend, TFloat, 1, 1)
    L  = KernelAbstractions.zeros(backend, TFloat, 1,1)
    Le_Lag = KernelAbstractions.zeros(backend, TFloat, 1,1)
    if lbuild_laplace_matrix
        if (backend == CPU())
            L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
            Le = build_laplace_matrix(SD, basis[1].ψ, basis[1].dψ, ω[1], mesh, metrics[1], N, Q, TFloat)
            Le_lag = build_laplace_matrix(SD, basis[2].ψ, basis[2].dψ, ω[2], mesh, metrics[2], mesh.ngr-1, mesh.ngr-1, TFloat)
            L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin),Int64(mesh.npoin))
        
            if ldss_laplace
                DSS_laplace_Laguerre!(L, SD, Le, Le_lag, ω[1], ω[2], mesh, metrics[1], metrics[2], N, TFloat; llump=inputs[:llump])
            end
        else
            Le = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
            Le_lag = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.ngr), Int64(mesh.ngr))
            
            k = build_laplace_matrix_gpu!(backend)
            k(Le, basis[1].dψ, ω[1], TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
            KernelAbstractions.synchronize(backend)
            k(Le_lag, basis[2].dψ, ω[2], TInt(mesh.ngr-1); ndrange = (mesh.ngr, mesh.ngr))
            L = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
            KernelAbstractions.synchronize(backend)
            k = DSS_laplace_gpu!(backend)
            k(L, Le, connijk, ω[1], ω[1], mesh.ngl, mesh.ngl, metrics[1].dξdx, metrics[1].dydη, metrics[1].dηdy, metrics[1].dxdξ; 
              ndrange = (mesh.nelem*mesh.ngl, mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
            KernelAbstractions.synchronize(backend)
            k = DSS_laplace_gpu_lag!(backend)
            k(L, Le_lag, connijk_lag, ω[1], ω[2], mesh.ngl, mesh.ngr, metrics[2].dηdx, metrics[1].dydη, metrics[1].dηdy, metrics[2].dxdη; 
              ndrange = (mesh.nelem_semi_inf*mesh.ngr, mesh.ngl), workgroupsize = (mesh.ngr, mesh.ngl))
            KernelAbstractions.synchronize(backend)
        end
    end

    De = KernelAbstractions.zeros(backend, TFloat, 1, 1)
    D  = KernelAbstractions.zeros(backend, TFloat, 1, 1)
    if lbuild_differentiation_matrix
        De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh,  N, Q, TFloat)
        if ldss_differentiation
            D  = DSS_generic_matrix(SD, De, mesh, TFloat)
        end
    end

    return (; Me, De, Le, M, Minv, D, L)
end

function mass_inverse!(Minv, M::AbstractVector, QT)
    Minv .= TFloat(1.0)./M
end

@kernel function diagm_gpu!(Minv_d, Minv)
    ip = @index(Global, Linear)

    Minv_d[ip,ip] = Minv[ip]
end

@kernel function add_to_diag!(M, val)
    ip = @index(Global, Linear)

    M[ip,ip] += val
end

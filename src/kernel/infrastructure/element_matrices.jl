using SparseArrays
using Base.Threads

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
function build_mass_matrix!(Me, SD::NSD_1D, QT::Inexact, ψ, ω, nelem, Je, Δx, N, Q, T)
    
    for iel=1:nelem
        #Jac = Δx[iel]/2
        
        for i=1:N+1
            Me[i,iel] += Je[iel,i]*ω[i]
        end
    end
end

function build_mass_matrix!(Me, SD::NSD_2D, QT::Inexact, ψ, ω, nelem, Je, Δx, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    @inbounds for iel=1:nelem
        
        for l = 1:Q+1
            for k = 1:Q+1
                
                ωkl  = ω[k]*ω[l]
                Jkle = Je[iel, k, l]
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

function build_mass_matrix!(Me, SD::NSD_3D, QT::Inexact, ψ, ω, nelem, Je, Δx, N, Q, T)
    
    MN = N + 1
    QN = Q + 1
    @inbounds for iel=1:nelem
        
        for o = 1:Q+1
            for n = 1:Q+1
                for m = 1:Q+1
                    
                    ωmno  = ω[m]*ω[n]*ω[o]
                    Jmnoe = Je[iel, m, n, o]
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


function build_mass_matrix_Laguerre!(Me, SD::NSD_1D, QT, ψ, ω, mesh, metrics, Δx, N, Q, T)

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

function build_laplace_matrix(SD::NSD_2D, ψ, dψ, ω, nelem, mesh, metrics, N, Q, T)
    #=
    Inexact integration (Q = N): nodal basis ψ[i,k] = δ_{ik}.
    The quadrature sums over (k,l) collapse analytically into four terms
    (A, B, C, D), reducing complexity from O(N^6) to O(N^4) per element.
    ψ is no longer needed.
    =#
    Np1 = N + 1
    Le = zeros(T, nelem, Np1, Np1, Np1, Np1)

    # Preallocate temporaries outside element loop → zero heap allocs inside
    G11 = zeros(T, Np1, Np1)   # ωk ωl J (ξx² + ξy²)
    G12 = zeros(T, Np1, Np1)   # ωk ωl J (ξx ηx + ξy ηy)
    G22 = zeros(T, Np1, Np1)   # ωk ωl J (ηx² + ηy²)
    K11 = zeros(T, Np1, Np1, Np1)  # K11[l,m,i] = Σ_k G11[k,l] dψ[m,k] dψ[i,k]
    K22 = zeros(T, Np1, Np1, Np1)  # K22[m,n,j] = Σ_l G22[m,l] dψ[n,l] dψ[j,l]

    for iel = 1:nelem

        # ── Step 1: Weighted contravariant metric at each quadrature point  O(N²) ──
        @inbounds for l = 1:Np1, k = 1:Np1
            wJ  = ω[k] * ω[l] * metrics.Je[iel,k,l]
            ξx  = metrics.dξdx[iel,k,l];  ξy = metrics.dξdy[iel,k,l]
            ηx  = metrics.dηdx[iel,k,l];  ηy = metrics.dηdy[iel,k,l]
            G11[k,l] = wJ * (ξx*ξx + ξy*ξy)
            G12[k,l] = wJ * (ξx*ηx + ξy*ηy)
            G22[k,l] = wJ * (ηx*ηx + ηy*ηy)
        end

        # ── Step 2a: K11[l,m,i] = Σ_k G11[k,l] dψ[m,k] dψ[i,k]  O(N⁴) ──
        # Feeds Term A:  Le[m, n, i, n] += K11[n, m, i]   (δ_{n,j} collapses j→n)
        fill!(K11, zero(T))
        @inbounds for l = 1:Np1, m = 1:Np1, i = 1:Np1
            s = zero(T)
            for k = 1:Np1
                s += G11[k,l] * dψ[m,k] * dψ[i,k]
            end
            K11[l,m,i] = s
        end

        # ── Step 2b: K22[m,n,j] = Σ_l G22[m,l] dψ[n,l] dψ[j,l]  O(N⁴) ──
        # Feeds Term D:  Le[m, n, m, j] += K22[m, n, j]   (δ_{m,i} collapses i→m)
        fill!(K22, zero(T))
        @inbounds for m = 1:Np1, n = 1:Np1, j = 1:Np1
            s = zero(T)
            for l = 1:Np1
                s += G22[m,l] * dψ[n,l] * dψ[j,l]
            end
            K22[m,n,j] = s
        end

        # ── Step 3: Assemble Le ──

        # Term A  (δ_{n,j} already applied → only j=n entries, O(N³))
        @inbounds for m = 1:Np1, n = 1:Np1, i = 1:Np1
            Le[iel,m,n,i,n] += K11[n,m,i]
        end

        # Term D  (δ_{m,i} already applied → only i=m entries, O(N³))
        @inbounds for m = 1:Np1, n = 1:Np1, j = 1:Np1
            Le[iel,m,n,m,j] += K22[m,n,j]
        end

        # Terms B + C  (cross-metric coupling, no Kronecker constraint, O(N⁴))
        # B: G12[i,n] dψ[m,i] dψ[j,n]   (k→i, l→n via δ_{ik}δ_{nl})
        # C: G12[m,j] dψ[i,m] dψ[n,j]   (k→m, l→j via δ_{mk}δ_{jl})
        @inbounds for n = 1:Np1, m = 1:Np1
            for j = 1:Np1, i = 1:Np1
                Le[iel,m,n,i,j] += G12[i,n]*dψ[m,i]*dψ[j,n] +
                                    G12[m,j]*dψ[i,m]*dψ[n,j]
            end
        end

    end
    return Le
end

function build_laplace_matrix_lag(SD::NSD_2D, ψ, dψ, ω, mesh, metrics, N, Q, T)
    
    Le = zeros((N+1),(N+1))
    
    for i=1:N+1
        for j=1:N+1
            for k=1:Q+1
                sum = ω[k]*dψ[i,k]*dψ[j,k]
                Le[i,j] = Le[i,j] + sum
            end
        end
    end 
    
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

# =============================================================================
#  Fused Laplace assembly + DSS — 2D SEM, inexact GLL integration
#
#  Original pipeline:
#    build_laplace_matrix  →  Le[nelem, Np1, Np1, Np1, Np1]   (45.8 MiB)
#    DSS_laplace!          →  L[ndof, ndof]  via scatter-add
#
#  This file replaces both with a single fused kernel.
#  Le is never allocated.  Each of the four inexact-integration terms is
#  scattered directly into L[IP,JP] via mesh.connijk as it is computed.
#
#  Memory per element (N=6, Float64):
#    G11,G12,G22   3 × 7²  × 8 =  1.2 KB
#    K11,K22       2 × 7³  × 8 =  5.4 KB
#    Le_loc                    =  0 KB  ← eliminated
#    ─────────────────────────────────
#    Total working set         ≈  6.6 KB  (fits in L1)
#
#  Complexity: O(N⁴) per element  (vs O(N⁶) in the unfused original)
# =============================================================================


# -----------------------------------------------------------------------------
"""
    build_sparsity_pattern(mesh, N, ndof, T)

One pass over mesh.connijk to collect all (IP, JP) pairs for every element.
Returns a structurally correct SparseMatrixCSC with nzval = 0.
The COO vectors are freed by GC immediately after sparse() returns.
"""
function build_sparsity_pattern(mesh, N, ::Type{T}) where {T}
    Np1    = N + 1
    nelem  = mesh.nelem
    ndof   = mesh.npoin                   # total unique DOFs
    nnz_ub = nelem * Np1^4               # upper bound; duplicates merged below

    Iv = Vector{Int32}(undef, nnz_ub)
    Jv = Vector{Int32}(undef, nnz_ub)

    ptr = 0
    @inbounds for iel = 1:nelem
        for j = 1:Np1, i = 1:Np1
            JP = mesh.connijk[iel, i, j]
            for n = 1:Np1, m = 1:Np1
                ptr += 1
                Iv[ptr] = mesh.connijk[iel, m, n]
                Jv[ptr] = JP
            end
        end
    end

    # sparse() merges duplicate (IP,JP) entries, setting nzval = 0 everywhere.
    # rowval / colptr are now fixed for the lifetime of the solve.
    return sparse(Iv, Jv, zeros(T, nnz_ub), ndof, ndof)
end


# -----------------------------------------------------------------------------
"""
    assemble_and_DSS_laplace!(L, dψ, ω, mesh, metrics, N)

Fused element assembly + DSS scatter.

Inexact GLL integration identity  ψᵢ(ξₖ) = δᵢₖ  reduces the O(N⁶) double
quadrature to four closed-form terms.  Each term is scattered immediately
into L[IP,JP] via mesh.connijk — Le (full or local) is never allocated.

Four terms after index collapse (row DOF = (m,n), col DOF = (i,j)):

  A  j = n  (δ_{nj})   Σ_k G11[k,n]  dψ[m,k] dψ[i,k]    ← K11[n,m,i]
  B                     G12[i,n]      dψ[m,i] dψ[j,n]
  C                     G12[m,j]      dψ[i,m] dψ[n,j]
  D  i = m  (δ_{mi})   Σ_l G22[m,l]  dψ[n,l] dψ[j,l]    ← K22[m,n,j]

L may be any AbstractMatrix (SparseMatrixCSC recommended, or dense for testing).
"""
function assemble_and_DSS_laplace!(L::AbstractMatrix{T},
                                    dψ, ω,
                                    mesh, metrics, N) where {T}
    Np1  = N + 1

    # ── Per-element temporaries — allocated once, reused every iteration ──────
    G11  = zeros(T, Np1, Np1)          # ωₖωₗ J (ξₓ² + ξᵧ²)
    G12  = zeros(T, Np1, Np1)          # ωₖωₗ J (ξₓηₓ + ξᵧηᵧ)
    G22  = zeros(T, Np1, Np1)          # ωₖωₗ J (ηₓ² + ηᵧ²)
    K11  = zeros(T, Np1, Np1, Np1)     # K11[l,m,i]  = Σ_k G11[k,l] dψ[m,k]dψ[i,k]
    K22  = zeros(T, Np1, Np1, Np1)     # K22[m,n,j]  = Σ_l G22[m,l] dψ[n,l]dψ[j,l]
    # Total: 6.6 KB for N=6 → fits in L1 cache

    for iel = 1:mesh.nelem

        # ── Step 1: Weighted contravariant metric tensors  O(N²) ──────────────
        @inbounds for l = 1:Np1, k = 1:Np1
            wJ        = ω[k] * ω[l] * metrics.Je[iel, k, l]
            ξx, ξy    = metrics.dξdx[iel,k,l], metrics.dξdy[iel,k,l]
            ηx, ηy    = metrics.dηdx[iel,k,l], metrics.dηdy[iel,k,l]
            G11[k,l]  = wJ * (ξx*ξx + ξy*ξy)
            G12[k,l]  = wJ * (ξx*ηx + ξy*ηy)
            G22[k,l]  = wJ * (ηx*ηx + ηy*ηy)
        end

        # ── Step 2a: K11[l,m,i] = Σ_k G11[k,l] dψ[m,k] dψ[i,k]  O(N⁴) ──────
        fill!(K11, zero(T))
        @inbounds for l = 1:Np1, m = 1:Np1, i = 1:Np1
            s = zero(T)
            for k = 1:Np1
                s += G11[k,l] * dψ[m,k] * dψ[i,k]
            end
            K11[l,m,i] = s
        end

        # ── Step 2b: K22[m,n,j] = Σ_l G22[m,l] dψ[n,l] dψ[j,l]  O(N⁴) ──────
        fill!(K22, zero(T))
        @inbounds for m = 1:Np1, n = 1:Np1, j = 1:Np1
            s = zero(T)
            for l = 1:Np1
                s += G22[m,l] * dψ[n,l] * dψ[j,l]
            end
            K22[m,n,j] = s
        end

        # ── Step 3: Scatter all four terms into L[IP,JP] ─────────────────────
        #
        # This IS the DSS scatter — mesh.connijk maps (iel,m,n) → global IP.
        # No Le_loc needed: each term goes straight to the global matrix.

        # ── Term A:  row=(m,n), col=(i,n)  [j locked to n by δ_{nj}]  O(N³) ──
        @inbounds for n = 1:Np1
            for m = 1:Np1
                IP = mesh.connijk[iel, m, n]
                for i = 1:Np1
                    JP = mesh.connijk[iel, i, n]    # j = n
                    L[IP, JP] += K11[n, m, i]
                end
            end
        end

        # ── Term D:  row=(m,n), col=(m,j)  [i locked to m by δ_{mi}]  O(N³) ──
        @inbounds for m = 1:Np1
            for n = 1:Np1
                IP = mesh.connijk[iel, m, n]
                for j = 1:Np1
                    JP = mesh.connijk[iel, m, j]    # i = m
                    L[IP, JP] += K22[m, n, j]
                end
            end
        end

        # ── Terms B+C:  no Kronecker constraint, full O(N⁴) ──────────────────
        @inbounds for n = 1:Np1, m = 1:Np1
            IP = mesh.connijk[iel, m, n]
            for j = 1:Np1, i = 1:Np1
                JP = mesh.connijk[iel, i, j]
                L[IP, JP] += G12[i,n]*dψ[m,i]*dψ[j,n] +
                             G12[m,j]*dψ[i,m]*dψ[n,j]
            end
        end

    end  # iel

    return L
end


# -----------------------------------------------------------------------------
"""
    build_laplace_matrix(SD, dψ, ω, mesh, metrics, N, T) → SparseMatrixCSC

Public entry point.  Two phases:
  1. build_sparsity_pattern  — allocates L with correct CSC structure, nzval=0
  2. assemble_and_DSS_laplace! — fills L.nzval by streaming element contributions

Replaces the original build_laplace_matrix + DSS_laplace! pair.
Le is never allocated.
"""
function build_laplace_matrix(SD::NSD_2D, dψ, ω,
                               mesh, metrics, N,
                               ::Type{T} = Float64) where {T}
    L = build_sparsity_pattern(mesh, N, T)
    assemble_and_DSS_laplace!(L, dψ, ω, mesh, metrics, N)
    return L
end


# -----------------------------------------------------------------------------
# COMPATIBILITY SHIM (optional)
# If other code still calls DSS_laplace! with a pre-built Le, this wrapper
# preserves that interface.  Remove once all callers are migrated.
#
function DSS_laplace!(L, SD::NSD_2D, Lel::AbstractArray, mesh, N; llump=false)
    Np1 = N + 1
    @inbounds for iel = 1:mesh.nelem
        for j = 1:Np1, i = 1:Np1
            JP = mesh.connijk[iel, i, j]
            for n = 1:Np1, m = 1:Np1
                IP = mesh.connijk[iel, m, n]
                L[IP, JP] += Lel[iel, m, n, i, j]
            end
        end
    end
end


# Alternative version with pre-computed element matrices for better performance
function DSS_laplace_sparse(mesh, Lel)

    # Pre-allocate arrays for triplet format
    # Estimate size: ngl^2 entries per element * number of elements
    max_entries = mesh.ngl^2 * mesh.ngl^2 * mesh.nelem
    
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    
    # Reserve space to avoid frequent reallocations
    sizehint!(I_vec, max_entries)
    sizehint!(J_vec, max_entries)
    sizehint!(V_vec, max_entries)
    
    # Assembly loop
    for iel = 1:mesh.nelem
        for j = 1:mesh.ngl, i = 1:mesh.ngl
            JP = mesh.connijk[iel, i, j]
            
            for n = 1:mesh.ngl, m = 1:mesh.ngl
                IP = mesh.connijk[iel, m, n]
                
                val = Lel[iel, m, n, i, j]
                if abs(val) > eps(Float64)  # Skip near-zero entries
                    push!(I_vec, IP)
                    push!(J_vec, JP)
                    push!(V_vec, val)
                end
            end
        end
    end
    
    # Create sparse matrix and sum duplicate entries automatically
    return sparse(I_vec, J_vec, V_vec)
end


function DSS_laplace_sparse_threaded(mesh, Lel)
    #
    # CSC aasembly
    #
    
    # Thread-local storage for triplets
    thread_triplets = [Tuple{Int, Int, Float64}[] for _ in 1:nthreads()]
    
    # Thread-parallel assembly
    @threads for iel = 1:mesh.nelem
        tid = threadid()
        local_triplets = thread_triplets[tid]
        
        for j = 1:mesh.ngl, i = 1:mesh.ngl
            JP = mesh.connijk[iel, i, j]
            
            for n = 1:mesh.ngl, m = 1:mesh.ngl
                IP = mesh.connijk[iel, m, n]
                
                val = Lel[iel, m, n, i, j]
                if abs(val) > eps(Float64)
                    push!(local_triplets, (IP, JP, val))
                end
            end
        end
    end
    
    # Combine thread-local results
    all_triplets = vcat(thread_triplets...)
    
    # Extract vectors and create sparse matrix
    I_vec = [t[1] for t in all_triplets]
    J_vec = [t[2] for t in all_triplets]
    V_vec = [t[3] for t in all_triplets]
    
    return sparse(I_vec, J_vec, V_vec)
end


# Utility function to convert to different sparse formats if needed
function convert_sparse_format(A::SparseMatrixCSC; format=:CSR)
    if format == :CSR
        # Julia's SparseMatrixCSC is essentially CSC format
        # For true CSR, you'd need to transpose and use rowvals/nzval
        return A'  # This gives CSR-like access pattern
    elseif format == :COO
        I, J, V = findnz(A)
        return (I, J, V)
    else
        return A  # Default CSC format
    end
end

# Example usage:
# L = assemble_diffusion_matrix_sparse(mesh, metrics, Lel, ω)
# 
# # For iterative solvers, you might want:
# using LinearAlgebra
# F = lu(L)  # Direct solver factorization
# 
# # Or for iterative methods:
# using IterativeSolvers
# x = cg(L, b)  # Conjugate gradient
####

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

function DSS_global_normals!(nx, ny, nz, mesh, SD::NSD_1D) nothing end

function DSS_global_normals!(nx, ny, nz, mesh, SD::NSD_2D)
   normals = zeros(Float64, mesh.npoin, 2)

    @inbounds for iedge = 1:mesh.nedges_bdy

        poin_edge = @view mesh.poin_in_bdy_edge[iedge, :]
        for i = 1:mesh.ngl
            ip = poin_edge[i]
            normals[ip, 1] += nx[iedge, i]
            normals[ip, 2] += ny[iedge, i]
            #end
        end
    end

    pM = setup_assembler(mesh.SD, normals, mesh.ip2gip, mesh.gip2owner)
    if pM != nothing
        assemble_mpi!(@view(normals[:,:]),pM)
    end
    @inbounds for iedge = 1:mesh.nedges_bdy

        poin_edge = @view mesh.poin_in_bdy_edge[iedge, :]
        for i = 1:mesh.ngl
            ip = poin_edge[i]
            mag = sqrt(normals[ip, 1]^2+ normals[ip, 2]^2)
            normx=0
            normy=0
            if (mag > 0)
                normx = normals[ip, 1]/mag
                normy = normals[ip, 2]/mag
            end
            if mesh.bdy_edge_type[iedge] != "periodicx" && (abs(nx[iedge, i] - normx) < 0.25)
                nx[iedge, i] = normx
            end
            if mesh.bdy_edge_type[iedge] != "periodicy" && (abs(ny[iedge, i] - normy) < 0.25)
                ny[iedge, i] = normy
            end

            #makesure vectors are unit vectors 
            mag = sqrt(nx[iedge, i]^2 + ny[iedge, i]^2)
            if (mag > 0)
                nx[iedge, i] = nx[iedge, i]/mag
                ny[iedge, i] = ny[iedge, i]/mag
            end
        end
    end

end
 
function DSS_global_normals!(nx, ny, nz, mesh, SD::NSD_3D)

    normals = zeros(Float64, mesh.npoin, 3)

    @inbounds for iface = 1:mesh.nfaces_bdy

        poin_face = @view mesh.poin_in_bdy_face[iface, :, :]
        for j = 1:mesh.ngl, i = 1:mesh.ngl
            ip              = poin_face[i, j]
            normals[ip, 1] += nx[iface, i, j]
            normals[ip, 2] += ny[iface, i, j]
            normals[ip, 3] += nz[iface, i, j]
        end
    end

    pM = setup_assembler(mesh.SD, normals, mesh.ip2gip, mesh.gip2owner)
    if pM != nothing
        assemble_mpi!(@view(normals[:,:]),pM)
    end
    @inbounds for iface = 1:mesh.nfaces_bdy

        poin_face = @view mesh.poin_in_bdy_face[iface, :, :]
        for j = 1:mesh.ngl, i = 1:mesh.ngl
            ip = poin_face[i, j]
            
            mag = sqrt(normals[ip, 1]^2+ normals[ip, 2]^2 + normals[ip, 3]^2)
            normx=0
            normy=0
            normz=0
            if (mag > 0)
            
                normx = normals[ip, 1]/mag
                normy = normals[ip, 2]/mag
                normz = normals[ip, 3]/mag
                
            end
            if mesh.bdy_face_type[iface] != "periodicx" && (abs(nx[iface, i, j] - normx) < 0.25)
                nx[iface, i, j] = normx
            end
            if mesh.bdy_face_type[iface] != "periodicy" && (abs(ny[iface, i, j] - normy) < 0.25)
                ny[iface, i, j] = normy
            end
            if mesh.bdy_face_type[iface] != "periodicz" && (abs(nz[iface, i, j] - normz) < 0.25)
                nz[iface, i, j] = normz
            end

            #makesure vectors are unit vectors 
            mag = sqrt(nx[iface, i, j]^2 + ny[iface, i, j]^2 + nz[iface, i, j]^2)
            if (mag > 0) 
                nx[iface, i, j] = nx[iface, i, j]/mag
                ny[iface, i, j] = ny[iface, i, j]/mag
                nz[iface, i, j] = nz[iface, i, j]/mag
            end
        end
    end
    
end



function DSS_global_RHS!(RHS, g_dss_cache, neqs)

    if g_dss_cache === nothing return end
    
    assemble_mpi!(@view(RHS[:,:]),g_dss_cache)
    
end

function DSS_global_RHS_pvector!(RHS, g_dss_cache, neqs)
    for i = 1:neqs
       DSS_global_RHS_v0!(@view(RHS[:,i]), g_dss_cache)
    end
end

function DSS_global_RHS_v0!(M, pM)
    # # @info ip2gip

    # pM = pvector(values->@view(M[:]), row_partition)
    sizeM = length(M)
    # pM = map(parts, local_values(pM)) do part, localpM
    #     @info part, length(localpM), sizeM
    #     localpM = copy(M)
    # end

    map( partition(pM)) do values
        for i = 1:sizeM
            values[i] = M[i]
        end
    end


    assemble!(pM) |> wait
    consistent!(pM) |> wait
    map(local_values(pM)) do values
        for i = 1:sizeM
            M[i] = values[i]
        end
    end
end


function DSS_global_mass!(SD, M, ip2gip, gip2owner, parts, npoin, gnpoin)

    if SD == NSD_1D()
        return nothing
    end
    
    #check_memory(" in sem_setup before setup_assembler.")
    g_dss_cache = setup_assembler(SD, M, ip2gip, gip2owner)
    #check_memory(" in sem_setup after setup_assembler.")
    
    assemble_mpi!(M,g_dss_cache)

    return g_dss_cache
    
end

function DSS_global_mass_pvector!(SD, M, ip2gip, gip2owner, parts, npoin, gnpoin)
    # @info ip2gip
    row_partition = map(parts) do part
        row_partition = LocalIndices(gnpoin,part,ip2gip,gip2owner)
        # gM = M
        row_partition
    end
    g_dss_cache = pvector(values->@view(M[:]), row_partition)

    assemble!(g_dss_cache) |> wait
    consistent!(g_dss_cache) |> wait
    M = map(local_values(g_dss_cache)) do values
        M = values
        M
    end
    return g_dss_cache
end



function matrix_wrapper(::ContGal, SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
                        ldss_laplace=false, ldss_differentiation=false, backend = CPU(), interp)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false
    if (ldss_differentiation) lbuild_differentiation_matrix = true end
    if (ldss_laplace) lbuild_laplace_matrix = true end

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if typeof(SD) == NSD_1D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, Int64(mesh.nelem))
    elseif typeof(SD) == NSD_2D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^2, (N+1)^2, Int64(mesh.nelem))
    elseif typeof(SD) == NSD_3D
        Me = KernelAbstractions.zeros(backend, TFloat, (N+1)^3, (N+1)^3, Int64(mesh.nelem))
    end
    
    if (backend == CPU())
        build_mass_matrix!(Me, SD, QT, basis.ψ, ω, mesh.nelem, metrics.Je, mesh.Δx, N, Q, TFloat)
    else
        if (SD == NSD_1D())
            k = build_mass_matrix_1d_gpu!(backend, (N+1))
            k(Me, basis.ψ, ω, metrics.Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
        elseif (SD == NSD_2D())
            k= build_mass_matrix_2d_gpu!(backend,(N+1,N+1))
            k(Me, basis.ψ, ω, metrics.Je, N, Q;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
        elseif (SD == NSD_3D())
            k= build_mass_matrix_3d_gpu!(backend,(N+1,N+1,N+1))
            k(Me, basis.ψ, ω, metrics.Je, N, Q;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl,mesh.ngl))
        end
    end
    
    if (QT == Exact() && inputs[:llump] == false)
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
    else
        M    = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
        Minv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    end
    
    if backend == CPU()
        #-------------------------------
        # backend -> CPU
        #-------------------------------
        if (inputs[:ladapt] == true)
            DSS_nc_gather_mass!(M, mesh, SD, QT, Me, mesh.connijk, mesh.poin_in_edge,
                                    mesh.non_conforming_facets, mesh.non_conforming_facets_parents_ghost,
                                    mesh.ip2gip, mesh.gip2ip, mesh.pgip_ghost, mesh.pgip_owner, N, interp)

        end

        if (rank == 0) println(" # DSS_mass .........................") end
        DSS_mass!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
        if (rank == 0) println(" # DSS_mass ......................... DONE") end
        
    else
        #-------------------------------
        # backend -> GPU
        #-------------------------------
        if SD == NSD_1D()
            DSS_mass!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
        elseif SD == NSD_2D()
            k = DSS_Mass_gpu_2D!(backend,(N+1,N+1))
            connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1, N+1)
            KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
            k(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
        elseif SD == NSD_3D()
            k = DSS_Mass_gpu_3D!(backend,(N+1,N+1,N+1))
            connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem), N+1, N+1,N+1)
            KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
            k(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl,mesh.ngl))
        end
    end
    
    g_dss_cache = DSS_global_mass!(SD, M, mesh.ip2gip, mesh.gip2owner, mesh.parts, mesh.npoin, mesh.gnpoin)

    if (rank == 0) println(" DSS_global_normals ......") end
    DSS_global_normals!(metrics.nx, metrics.ny, metrics.nz, mesh, SD)
    if (rank == 0) println(" DSS_global_normals ...... DONE") end

    if (inputs[:ladapt] == true)
        DSS_nc_scatter_mass!(M, SD, QT, Me, mesh.connijk, mesh.poin_in_edge, mesh.non_conforming_facets,
                                   mesh.non_conforming_facets_children_ghost, mesh.ip2gip, mesh.gip2ip, mesh.cgip_ghost, mesh.cgip_owner, N, interp)
    end

    if (inputs[:bdy_fluxes])
        if SD == NSD_3D()
            M_surf = build_surface_mass_matrix(mesh.nfaces_bdy, mesh.npoin, ω, basis.ψ, mesh.ngl, metrics.Jef, mesh.poin_in_bdy_face, TFloat, mesh.Δx, inputs)
            assemble_mpi!(M_surf,g_dss_cache)
            M_surf_inv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
            mass_inverse!(M_surf_inv, M_surf, QT)
            M_edge_inv = KernelAbstractions.zeros(backend, TFloat, 1)
            
        else
            M_surf_inv = KernelAbstractions.zeros(backend, TFloat, 1)
            M_edge = build_segment_mass_matrix(mesh.nedges_bdy, mesh.npoin, ω, basis.ψ, mesh.ngl, metrics.Jef, mesh.poin_in_bdy_edge, TFloat, mesh.Δx, inputs)
            assemble_mpi!(M_edge,g_dss_cache)
            M_edge_inv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
            mass_inverse!(M_edge_inv, M_edge, QT)
        end
    else
        M_surf_inv = KernelAbstractions.zeros(backend, TFloat, 1)
        M_edge_inv = KernelAbstractions.zeros(backend, TFloat, 1)
    end

    if (rank == 0) println(" mass_inverse ......") end
    mass_inverse!(Minv, M, QT)
    if (rank == 0) println(" mass_inverse ...... DONE") end
    
    Le = KernelAbstractions.zeros(backend,TFloat, 1, 1)
    L  = KernelAbstractions.zeros(backend, TFloat, 1,1)
    
    if lbuild_laplace_matrix
        if (backend == CPU())
            #
            # CPU
            #
            #
            if (rank == 0) println(" build_laplace_matrix (fused DSS) .......... ") end
            L = build_laplace_matrix(SD,
                                     basis.dψ,
                                     ω,
                                     mesh,
                                     metrics,
                                     N,
                                     TFloat)
            
            if (rank == 0) println(" build_laplace_matrix (fused DSS) ......... DONE") end
            
            #=
            if (rank == 0) println(" build_laplace_matrix ...................... ") end
            Le = @time build_laplace_matrix(SD,
            basis.ψ, basis.dψ,
            ω, mesh.nelem,
            mesh,
            metrics,
            N, Q, TFloat)
            if (rank == 0) println(" build_laplace_matrix ..................... DONE") end
            
            if (inputs[:lsparse])
            if (rank == 0) println(" DSS sparse format") end
            L = DSS_laplace_sparse(mesh, Le)
            if (rank == 0) println(" DSS sparse format .................... DONE") end
            else
                if (rank == 0) println(" DSS non-sparse format") end
                L = KernelAbstractions.zeros(backend,
                                             TFloat,
                                             Int64(mesh.npoin),
                                             Int64(mesh.npoin))

                DSS_laplace!(L, SD,
                             Le, ω,
                             mesh, metrics,
                             N, TFloat;
                             llump=inputs[:llump])
                if (rank == 0) println(" DSS non-sparse format .................... DONE") end
            end
            =#
        else
            #
            # GPU
            #
            Le = KernelAbstractions.zeros(backend,
                                          TFloat,
                                          Int64(mesh.ngl),
                                          Int64(mesh.ngl))

            k = build_laplace_matrix_gpu!(backend)

            k(Le, basis.dψ, ω, TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))

            KernelAbstractions.synchronize(backend)
            
            L = KernelAbstractions.zeros(backend,
                                         TFloat,
                                         Int64(mesh.npoin),
                                         Int64(mesh.npoin))
            
            KernelAbstractions.synchronize(backend)
            k = DSS_laplace_gpu!(backend)
            
            k(L, Le,
              connijk,
              ω, ω,
              mesh.ngl,
              mesh.ngl,
              metrics.dξdx,
              metrics.dydη,
              metrics.dηdy,
              metrics.dxdξ;
              ndrange = (mesh.nelem*mesh.ngl,
                         mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
            
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
    
    return (; Me, De, M, Minv, g_dss_cache, D, L, M_surf_inv, M_edge_inv)
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

function matrix_wrapper_laguerre(::ContGal, SD, QT, basis, ω, mesh, metrics, N, Q, TFloat;
                                 ldss_laplace=false, ldss_differentiation=false, backend = CPU(), interp)

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
        build_mass_matrix!(Me, SD, QT, basis[1].ψ, ω[1], mesh.nelem, metrics[1].Je, mesh.Δx, N, Q, TFloat)
    elseif (SD == NSD_1D())
        k = build_mass_matrix_1d_gpu!(backend, (N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
    elseif (SD == NSD_2D())
        k = build_mass_matrix_2d_gpu!(backend, (N+1, N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, N, Q;ndrange = (mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
    else
        nothing
    end
    
    if (inputs[:bdy_fluxes])
        if SD == NSD_3D()
            M_surf = build_surface_mass_matrix(mesh.nfaces_bdy, mesh.npoin, ω, basis.ψ, mesh.ngl, metrics.Jef, mesh.poin_in_bdy_face, TFloat, mesh.Δx, inputs)
            assemble_mpi!(M_surf,g_dss_cache)
            M_surf_inv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
            mass_inverse!(M_surf_inv, M_surf, QT)
            M_edge_inv = KernelAbstractions.zeros(backend, TFloat, 1)
            
        else
            M_surf_inv = KernelAbstractions.zeros(backend, TFloat, 1)
            M_edge = build_segment_mass_matrix(mesh.nedges_bdy, mesh.npoin, ω, basis.ψ, mesh.ngl, metrics.Jef, mesh.poin_in_bdy_edge, TFloat, mesh.Δx, inputs)
            assemble_mpi!(M_edge,g_dss_cache)
            M_edge_inv = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
            mass_inverse!(M_edge_inv, M_edge, QT)
        end
    else
        M_surf_inv = KernelAbstractions.zeros(backend, TFloat, 1)
        M_edge_inv = KernelAbstractions.zeros(backend, TFloat, 1)
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
            build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[2].ψ, ω[2], mesh, metrics[2], mesh.Δx, N, Q, TFloat)
        elseif typeof(SD) == NSD_2D
            build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[1].ψ, basis[2].ψ, ω[1], ω[2], mesh, metrics[2], N, Q, TFloat)
        end
        
        DSS_mass_Laguerre!(M, SD, Me, M_lag, mesh, N, TFloat; llump=inputs[:llump])
        g_dss_cache = DSS_global_mass!(SD, M, mesh.ip2gip, mesh.gip2owner, mesh.parts, mesh.npoin, mesh.gnpoin)
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
            Le = build_laplace_matrix_lag(SD,
                                          basis[1].ψ,
                                          basis[1].dψ,
                                          ω[1],
                                          mesh,
                                          metrics[1],
                                          N, Q, TFloat)
            
            Le_lag = build_laplace_matrix_lag(SD,
                                              basis[2].ψ,
                                              basis[2].dψ,
                                              ω[2],
                                              mesh,
                                              metrics[2],
                                              mesh.ngr-1,
                                              mesh.ngr-1,
                                              TFloat)
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

    return (; Me, De, Le, M, Minv, M_edge_inv, M_surf_inv, g_dss_cache, D, L)
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

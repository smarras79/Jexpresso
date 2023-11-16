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
                        M[IP] = M[IP] + Mel[I,J,iel] #if inexact
                    end
                end
            end
        end    
    end
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


function DSS_rhs!(RHS, rhs_el, mesh, nelem, ngl, neqs, SD::NSD_1D)

    for ieq = 1:neqs
        for iel = 1:nelem
            for i = 1:ngl
                I = mesh.connijk[iel,i,1]
                RHS[I,ieq] += rhs_el[iel,i,1,ieq]
            end
        end
    end
    
end

function DSS_rhs!(RHS, rhs_el, mesh, nelem, ngl, neqs, SD::NSD_2D)

    for ieq = 1:neqs
        for iel = 1:nelem
            for j = 1:ngl
                for i = 1:ngl
                    I = mesh.connijk[iel,i,j]
                    RHS[I,ieq] += rhs_el[iel,i,j,ieq]
                end
            end
        end
    end
    #show(stdout, "text/plain", V)
end

function DSS_rhs_laguerre!(RHS, rhs_el, mesh, nelem, ngl, neqs, SD::NSD_1D)

  for ieq = 1:neqs
    for iel = 1:mesh.nelem_semi_inf
        for i = 1:mesh.ngr
            I = mesh.connijk_lag[iel,i,1]

            RHS[I,ieq] += rhs_el[iel,i,1,ieq]
        end
    end
  end
end

function DSS_rhs_laguerre!(RHS, rhs_el, mesh, nelem, ngl, neqs, SD::NSD_2D)

  for ieq = 1:neqs
    for iel = 1:mesh.nelem_semi_inf
        for j = 1:mesh.ngr
            for i = 1:mesh.ngl
                I = mesh.connijk_lag[iel,i,j]
                
                 RHS[I,ieq] += rhs_el[iel,i,j,ieq]
            end
        end
    end
  end
end

#function divide_by_mass_matrix!(RHS::AbstractArray, RHSaux, Minv, neqs, npoin, ::Exact)
function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin)
    
    RHSaux .= RHS
    for ip=1:npoin
        a = zero(eltype(RHS))
        for jp = 1:npoin
            a += Minv[ip,jp]*RHSaux[jp]
        end
        RHS[ip] = a
    end
    
end

function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin)

    for ip=1:npoin
        RHS[ip] = Minv[ip]*RHS[ip]
    end
end

function matrix_wrapper(SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
                        ldss_laplace=false, ldss_differentiation=false)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false

    if typeof(SD) == NSD_1D
        Me = zeros(TFloat, (N+1)^2, mesh.nelem)
    elseif typeof(SD) == NSD_2D
        Me = zeros(TFloat, (N+1)^2, (N+1)^2, mesh.nelem)
    elseif typeof(SD) == NSD_3D
        Me = zeros(TFloat, (N+1)^2, (N+1)^2, (N+1)^2, mesh.nelem)
    end
    build_mass_matrix!(Me, SD, QT, basis.ψ, ω, mesh, metrics, N, Q, TFloat)
    
    if (QT == Exact() && inputs[:llump] == false)
        M    = zeros(TFloat, mesh.npoin, mesh.npoin)
        Minv = zeros(TFloat, mesh.npoin, mesh.npoin)
    else
        M    = zeros(TFloat, mesh.npoin)
        Minv = zeros(TFloat, mesh.npoin)
    end
    DSS_mass!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
    mass_inverse!(Minv, M, QT)
    
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
    
    return (; Me, De, Le, M, Minv, D, L)
end

function mass_inverse!(Minv, M::AbstractMatrix, QT)
    Minv .= inv(M)
end

function matrix_wrapper_laguerre(SD, QT, basis, ω, mesh, metrics, N, Q, TFloat; ldss_laplace=false, ldss_differentiation=false)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false

    if typeof(SD) == NSD_1D
        Me = zeros(TFloat, (N+1)^2, mesh.nelem)
    elseif typeof(SD) == NSD_2D
        Me = zeros(TFloat, (N+1)^2, (N+1)^2, mesh.nelem)
    end
    
    build_mass_matrix!(Me, SD, QT, basis[1].ψ, ω[1], mesh, metrics[1], N, Q, TFloat)
    if typeof(SD) == NSD_1D
      M_lag = zeros(TFloat, mesh.ngr*mesh.ngr, mesh.nelem_semi_inf)
    elseif typeof(SD) == NSD_2D
      M_lag = zeros(TFloat, mesh.ngl*mesh.ngr, mesh.ngl*mesh.ngr, mesh.nelem_semi_inf)
    end
    if (QT == Exact() && inputs[:llump] == false)
        M    = zeros(TFloat, mesh.npoin, mesh.npoin)
        Minv = zeros(TFloat, mesh.npoin, mesh.npoin)
    else
        M    = zeros(TFloat, mesh.npoin)
        Minv = zeros(TFloat, mesh.npoin)
    end
    if typeof(SD) == NSD_1D
        build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[2].ψ, ω[2], mesh, metrics[2], N, Q, TFloat)
    elseif typeof(SD) == NSD_2D
        build_mass_matrix_Laguerre!(M_lag, SD, QT, basis[1].ψ, basis[2].ψ, ω[1], ω[2], mesh, metrics[2], N, Q, TFloat)
    end
    #@info maximum(Me), minimum(Me), maximum(M_lag),minimum(M_lag) 
    DSS_mass_Laguerre!(M, SD, Me, M_lag, mesh, N, TFloat; llump=inputs[:llump])
    mass_inverse!(Minv, M, QT)
   # @info maximum(M),minimum(M),maximum(Minv),minimum(Minv)
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

    return (; Me, De, Le, M, Minv, D, L)
end

function mass_inverse!(Minv, M::AbstractVector, QT)
    Minv .= 1.0./M
end

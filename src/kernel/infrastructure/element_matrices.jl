using SparseArrays
using JACC
using Base.Threads
using Atomix


abstract type AbstractMassType end
mutable struct St_ElMat{TFloat} <: AbstractMassType
    M::Array{TFloat} #Mass
    D::Array{TFloat} #Differentiation
    L::Array{TFloat} #Laplacian
end

function DSS_rhs!(RHS, rhs_el, connijk, nelem, ngl, neqs, SD, method)
    total_iterations = neqs * nelem * ngl
    if(SD == NSD_1D())
        @info "Were here"
        JACC.parallel_for(total_iterations, DSS_rhs_jacc!, RHS, rhs_el, connijk, nelem, ngl, neqs, SD, method)
        @info "And here"
    end
end

function DSS_rhs_laguerre!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, SD, method)
    DSS_rhs_laguerre_jacc!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, SD, method)
end

function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin, method)
    divide_by_mass_matrix_jacc_all_parallel!(RHS, RHSaux, Minv, neqs, npoin, method)
end

function divide_by_mass_matrix!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, method)
    divide_by_mass_matrix_jacc!(RHS, RHSaux, Minv, neqs, npoin, method)
end

function mass_inverse!(Minv, M::AbstractVector, QT)
    Minv .= TFloat(1.0)./M
end

function build_mass_matrix_1d_jacc!(N, Me, ψ, ω, Je, Q)
    ie = (N - 1) ÷ size(Je, 2) + 1
    i = (N - 1) % size(Je, 2) + 1
    
    Me[i, ie] += Je[ie, i, 1] * ω[i]

    #@info M[i, ie], Je[ie,i, 1], ω[i] # Debugging information
end


#JACCIFY
function build_mass_matrix_2d_jacc!(idx_linear, Me, ψ, ω, Je, N, Q)
    # For JACC.parallel_for, idx_linear is the single linear index passed to each thread
    # We need to decompose this into our 5D indices manually
    
    # Total dimensions for decomposition
    total_I = (N + 1)^2
    total_J = (N + 1)^2
    total_ie = size(Me, 3)
    total_k = Q + 1 
    total_l = Q + 1
    
    # Check if we're within the valid range first
    total_iterations = total_I * total_J * total_ie * total_k * total_l
    if idx_linear > total_iterations
        return
    end
    
    # Decompose the linear index into 5D indices
    idx = idx_linear - 1  # Convert to 0-based for easier math
    
    l = mod(idx, total_l) + 1
    idx = div(idx, total_l)
    
    k = mod(idx, total_k) + 1
    idx = div(idx, total_k)
    
    ie = mod(idx, total_ie) + 1
    idx = div(idx, total_ie)
    
    J = mod(idx, total_J) + 1
    idx = div(idx, total_J)
    
    I = idx + 1

    # Decompose the linear indices 'I' and 'J' into 2D basis function indices.
    # Note: Using 0-based indexing for the decomposition, then adding 1 for Julia indexing
    i_y = div(I - 1, N + 1) + 1
    i_x = mod(I - 1, N + 1) + 1
    n = div(J - 1, N + 1) + 1
    m = mod(J - 1, N + 1) + 1
    
    #@info i_x i_y m n k l
    # Comprehensive bounds checking
    #=if (i_x < 1 || i_x > N+1 || 
        i_y < 1 || i_y > N+1 || 
        m < 1 || m > N+1 || 
        n < 1 || n > N+1 ||
        k < 1 || k > Q || 
        l < 1 || l > Q ||
        ie < 1 || ie > size(Me, 3) ||
        I < 1 || I > total_I ||
        J < 1 || J > total_J)
        return
    end=#

    # Additional check: ensure ψ array has correct dimensions
    if (size(ψ, 1) < max(i_x, m) || size(ψ, 2) < max(k, l))
        return
    end

    # Calculate the single contribution to the integral from quadrature point (k,l).
    # This computes: ω_k * ω_l * J_e(k,l) * ψ_{i_x}(k) * ψ_{i_y}(l) * ψ_m(k) * ψ_n(l)
    term = (ω[k] * ω[l]) * Je[ie, k, l] * (ψ[i_x, k] * ψ[i_y, l]) * (ψ[m, k] * ψ[n, l])

    # Atomically add this thread's contribution to the correct matrix entry.
    
    if(J == 5 && I == 5 && ie == 1)
        @info term i_x i_y m n k l
    end
    
    # This operation is essential to safely update the same memory location from multiple threads.
    Atomix.@atomic Me[I, J, ie] += term

end




#=    ie = @index(Group, Linear)
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
end=#

#JACCIFY
function build_mass_matrix_Laguerre_2d_gpu!(Me, ψ, ψ1, ω, ω1, Je, ngl, ngr)


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


function build_laplace_matrix_jacc_parallel_all!(idx_linear, Le, dψ, ω, Q)
    #this calculation is done twice, maybe pass it in?
    n_nodes = size(dψ, 1)

   # JACC.parallel_for(1:(n_nodes * n_nodes * (Q + 1))) do idx_linear
        i = (idx_linear - 1) % n_nodes + 1
        temp_idx = (idx_linear - 1) ÷ n_nodes
        j = temp_idx % n_nodes + 1
        k = temp_idx ÷ n_nodes + 1

        term = ω[k] * dψ[i, k] * dψ[j, k]

        Le[i,j] = Le[i,j] + term
        @info "Le[i,j]", i, j, term, Le[i,j]
    #end
end

#JACCIFY
function DSS_mass_Laguerre_gpu_2D!(M, Mel_lag, connijk_lag, ngl, ngr)

    iel = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    J = i_x + (i_y - 1)*(ngl)
    for n = 1:ngr
        for m = 1:ngl
            I = m + (n - 1)*(ngl)
            IP = connijk_lag[iel,m,n]
           # KernelAbstractions.@atomic M[IP] += Mel_lag[I,J,iel] #if inexact
            M[IP] += Mel_lag[I,J,iel] #if exact
        end
    end
end

#NEED TO STILL MAKE THIS A JACC function
function DSS_Mass_jacc_1D!(M, SD, QT, Me, conn, nelem, npoin, N, TFloat; llump=false)    
    n_elements = size(conn, 1)
    n_local_nodes = size(conn, 2)
    
    for idx_linear in 1:(n_elements * n_local_nodes)
        ie = (idx_linear - 1) ÷ n_local_nodes + 1
        i = (idx_linear - 1) % n_local_nodes + 1
        IP = conn[ie, i, 1]
        M[IP] += Me[i, ie]
    end
end

#JACCIFY
#=function DSS_Mass_gpu_2D_kernel!(ie, M, Mel, conn, N)
    for i_y = 1:N+1
        for i_x = 1:N+1
            J = i_x + (i_y - 1)*(N+1)
            JP = conn[ie, i_x, i_y]
            
            for n = 1:N+1
                for m = 1:N+1
                    I = m + (n-1)*(N+1)
                    IP = conn[ie, m, n]
                    Atomix.@atomic M[IP] += Mel[I, J, ie]
                end
            end
        end
    end
end=#


function DSS_Mass_jacc_2D!(M, Mel, conn, nelem, npoin, N)
    n_elements = size(conn, 1)
    n_local_nodes_per_dim = N + 1
    total_local_nodes = n_local_nodes_per_dim * n_local_nodes_per_dim
    
    @info Mel[:,:,1]

    for idx_linear in 1:(n_elements * total_local_nodes)
        ie = (idx_linear - 1) ÷ total_local_nodes + 1
        local_idx = (idx_linear - 1) % total_local_nodes + 1
        
        i_x = (local_idx - 1) % n_local_nodes_per_dim + 1
        i_y = (local_idx - 1) ÷ n_local_nodes_per_dim + 1
        
        J = i_x + (i_y - 1) * (N + 1)
        JP = conn[ie, i_x, i_y]
        
        for n = 1:N+1
            for m = 1:N+1
                I = m + (n - 1) * (N + 1)
                IP = conn[ie, m, n]
                M[IP] += Mel[I, J, ie]
            end
        end
    end
end






function DSS_laplace_jacc!(idx_linear, L, Lel, connijk, ωx, ωy, nx, ny, dξdx, dydη, dηdy, dxdξ)
    # If something goes wrong, check this file and try to see if
    # the parallel fors in this method need to be seperated
   
    n_elements = size(connijk, 1)
    n_local_i = size(connijk, 2)
    n_local_j = size(connijk, 3)

    #ACC.parallel_for(1:(n_elements * n_local_i * n_local_j * nx)) do idx_linear
        ie = (idx_linear - 1) ÷ (n_local_i * n_local_j * nx) + 1
        rem1 = (idx_linear - 1) % (n_local_i * n_local_j * nx)
        i = (rem1 ÷ (n_local_j * nx)) + 1
        rem2 = rem1 % (n_local_j * nx)
        j = (rem2 ÷ nx) + 1
        k = (rem2 % nx) + 1

        ip = connijk[ie, i, j]
        jp = connijk[ie, k, j]
        term = dξdx[ie, i, k] * Lel[i, k] * dydη[ie, i, k] * ωx[j]

        L[ip, jp] += term
      # JACC.atomic_add!(L, ip, jp, term) MADE THIS NOT ATOMIC
    #end
    
    #JACC.parallel_for(1:(n_elements * n_local_i * n_local_j * ny)) do idx_linear
        ie = (idx_linear - 1) ÷ (n_local_i * n_local_j * ny) + 1
        rem1 = (idx_linear - 1) % (n_local_i * n_local_j * ny)
        i = (rem1 ÷ (n_local_j * ny)) + 1
        rem2 = rem1 % (n_local_j * ny)
        j = (rem2 ÷ ny) + 1
        l = (rem2 % ny) + 1

        ip = connijk[ie,i,j]
        jp = connijk[ie,i,l]
        term = dηdy[ie, i, l] * Lel[j, l] * dxdξ[ie, i, l] * ωy[i]

        L[ip, jp] += term
      # JACC.atomic_add!(L, ip, jp, term) MADE THIS NOT ATOMIC
   # end
end

function DSS_laplace_jacc_lag!(L, Lel, connijk, ωx, ωy, nx, ny, dηdx_lag, dydη, dηdy, dxdη_lag)
    n_elements = size(connijk, 1)
    n_local_i = size(connijk, 2)
    n_local_j = size(connijk, 3)

    JACC.parallel_for(1:(n_elements * n_local_i * n_local_j * ny)) do idx_linear
        ie = (idx_linear - 1) ÷ (n_local_i * n_local_j * ny) + 1
        rem1 = (idx_linear - 1) % (n_local_i * n_local_j * ny)
        j = (rem1 ÷ (n_local_i * ny)) + 1
        rem2 = rem1 % (n_local_i * ny)
        i = (rem2 ÷ ny) + 1
        k = (rem2 % ny) + 1

        ip = connijk[ie, j, i]
        jp = connijk[ie, j, k]
        term = dηdx_lag[ie, j, k] * Lel[i, k] * dydη[ie, j, j] * ωx[j]

        JACC.atomic_add!(L, ip, jp, term)
    end
    
    JACC.parallel_for(1:(n_elements * n_local_i * n_local_j * nx)) do idx_linear
        ie = (idx_linear - 1) ÷ (n_local_i * n_local_j * nx) + 1
        rem1 = (idx_linear - 1) % (n_local_i * n_local_j * nx)
        j = (rem1 ÷ (n_local_i * nx)) + 1
        rem2 = rem1 % (n_local_i * nx)
        i = (rem2 ÷ nx) + 1
        l = (rem2 % nx) + 1

        ip = connijk[ie,j,i]
        jp = connijk[ie,l,i]
        term = dηdy[ie, j, l] * Lel[j, l] * dxdη_lag[ie, l, i] * ωy[i]
        
        JACC.atomic_add!(L, ip, jp, term)
    end
end

function DSS_global_RHS!(RHS, pM, neqs)

    if pM == nothing return end
    
    assemble_mpi!(@view(RHS[:,:]),pM)
    
end


function DSS_rhs_jacc!(RHS, rhs_el, mesh, nelem, ngl, neqs, ::NSD_1D, ::FD)
    nothing
end
function DSS_rhs_jacc!(RHS, rhs_el, mesh, nelem, ngl, neqs, ::NSD_2D, ::FD)
    nothing
end

function DSS_rhs_jacc!(N, RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_1D, ::ContGal)

    ieq = (N - 1) ÷ (nelem * ngl) + 1
    rem1 = (N - 1) % (nelem * ngl)
    iel = rem1 ÷ ngl + 1
    i = rem1 % ngl + 1

    I = connijk[iel, i, 1]

    RHS[I, ieq] += rhs_el[iel, i, ieq]
end



function DSS_rhs_laguerre_jacc!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, ::NSD_1D, ::ContGal)
    total_tasks = neqs * nelem_semi_inf * ngr

    JACC.parallel_for(1:total_tasks) do idx_linear
        ieq = (idx_linear - 1) ÷ (nelem_semi_inf * ngr) + 1
        rem = (idx_linear - 1) % (nelem_semi_inf * ngr)
        iel = rem ÷ ngr + 1
        i = rem % ngr + 1

        I = connijk_lag[iel, i, 1]

        RHS[I,ieq] += rhs_el[iel,i,ieq]
    end
end


function divide_by_mass_matrix_jacc_all_parallel!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin, ::ContGal)

    total_elements = npoin * npoin
    JACC.parallel_for(1:total_elements) do idx_linear
        ip = (idx_linear - 1) % npoin + 1
        jp = (idx_linear - 1) ÷ npoin + 1
        
         partial_product = Minv[ip, jp] * RHSaux[jp]
        
        JACC.atomic_add!(RHS, ip, partial_product)
    end
end

#=
function divide_by_mass_matrix_jacc!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, ::ContGal)
    # The loop iterates over each point, 'ip'
    JACC.parallel_for(1:npoin) do ip
        # Each thread handles a unique 'ip' and performs the element-wise multiplication
        RHS[ip] = Minv[ip] * RHS[ip]
    end
end
=#

function divide_by_mass_matrix_jacc!(RHS, RHSaux, Minv::AbstractVector, neqs, npoin, method::ContGal)
    for ip = 1:npoin
        RHS[ip] = Minv[ip] * RHS[ip]
    end
end

# I believe this is not used anywhere

function matrix_wrapper(::FD, SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
                        ldss_laplace=false, ldss_differentiation=false)

    if typeof(SD) == NSD_1D
        Me_base = zeros(TFloat, 1, 1)
        Me = JACC.array(Me_base)
    elseif typeof(SD) == NSD_2D
        Me_base = zeros(TFloat, 1, 1)
        Me = JACC.array(Me_base)   
     elseif typeof(SD) == NSD_3D
        Me_base = zeros(TFloat, 1, 1)
        Me = JACC.array(Me_base)
    end
    
    if (QT == Exact() && inputs[:llump] == false)
        M_base = zeros(TFloat, 1, 1)
        M = JACC.array(M_base)      
        Minv_base = zeros(TFloat, 1, 1)
        Minv = JACC.array(Minv_base)
    else
        M_base = zeros(TFloat, 1, 1)
        M = JACC.array(M_base)      
        Minv_base = zeros(TFloat, 1, 1)
        Minv = JACC.array(Minv_base)
    end
    
    Le_base = zeros(TFloat, 1, 1)
    Le = JACC.array(Le_base)
    L_base  = zeros(TFloat, 1,1)
    L  = JACC.array(L_base)
    
    De_base = zeros(TFloat, 1, 1)
    De = JACC.array(De_base)
    D_base  = zeros(TFloat, 1,1)
    D  = JACC.array(D_base)
   
    @info "We are here", maximum(M), maximum(Minv), maximum(Me), maximum(De), maximum(D), maximum(Le), maximum(L)
    return (; Me, De, Le, M, Minv, D, L, M_surf_inv=Minv, M_edge_inv=Minv)
    
end


function matrix_wrapper(::ContGal, SD, QT, basis::St_Lagrange, ω, mesh, metrics, N, Q, TFloat;
    ldss_laplace=false, ldss_differentiation=false, backend = CPU(), interp=nothing)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false
    if (ldss_differentiation) lbuild_differentiation_matrix = true end
    if (ldss_laplace) lbuild_laplace_matrix = true end

    if typeof(SD) == NSD_1D
        Me_base = zeros(TFloat, (N+1)^2, Int64(mesh.nelem))
        Me = JACC.array(Me_base)
        #build_mass_matrix_1d_jacc!(Me, basis.ψ, ω, metrics.Je, Q)
        total_iterations = size(metrics.Je, 1)* size(metrics.Je, 2) 

        JACC.parallel_for(total_iterations, build_mass_matrix_1d_jacc!, Me, basis.ψ, ω, metrics.Je, Q)
    
    elseif typeof(SD) == NSD_2D
        Me_base = zeros(TFloat, (N+1)^2, (N+1)^2, Int64(mesh.nelem))
        Me = JACC.array(Me_base)

        #Make this a parallel for call
        num_elements = size(Me, 3)
        total_iterations = num_elements * (N + 1)^4 * (Q + 1)^2

        #JACC.parallel_for(N, build_mass_matrix_2d_jacc,Me, basis.ψ, ω, metrics.Je, N, Q)
        JACC.parallel_for(total_iterations, build_mass_matrix_2d_jacc!, Me, basis.ψ, ω, metrics.Je, N, Q)

        #=else
        if (SD == NSD_1D())
            k = build_mass_matrix_1d_gpu!((N+1))
            k(Me, basis.ψ, ω, metrics.Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
        end=#
    end
    if (QT == Exact() && inputs[:llump] == false)
        M_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        M = JACC.array(M_base)
        Minv_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        Minv = JACC.array(Minv_base)
    else
        M_base = zeros(TFloat, Int64(mesh.npoin))
        M = JACC.array(M_base)
        Minv_base = zeros(TFloat, Int64(mesh.npoin))
        Minv = JACC.array(Minv_base)
    end
   
    if SD == NSD_1D()
        @info "Building 1D Mass Matrix!!!!!!!!!!"
        DSS_Mass_jacc_1D!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
    elseif SD == NSD_2D()
        @info "Building 2D Mass Matrix!!!!!!!!!!"
        #total_iterations = mesh.nelem * (N + 1)^4
        # N is the polynomial order, nelem is the number of elements
        # Use the actual connectivity array dimensions
        @info N
     #   JACC.parallel_for(mesh.nelem, DSS_Mass_gpu_2D_kernel!, M, Me, mesh.connijk, N)

        # Launch the JACC kernel - you can still pass N but the function will use actual dimensions
    #    JACC.parallel_for(total_iterations, DSS_Mass_jacc_2D!, M, Me, mesh.connijk, mesh.nelem, N)
      #  JACC.parallel_for(total_iterations, DSS_Mass_jacc_2D!, M, Me, mesh.connijk, mesh.nelem, mesh.npoin, N)
    
        DSS_Mass_jacc_2D!(M, Me, mesh.connijk, mesh.nelem, mesh.npoin, N)



        #DSS_Mass_gpu_2D!(M, Me, mesh.connijk, mesh.nelem, mesh.npoin, N)
        #synchronize()
    end
    @info "Done Building 2D Mass Matrix!!!!!!!!!!"

    mass_inverse!(Minv, M, QT)

    Le_base = zeros(TFloat, 1, 1)
    Le = JACC.array(Le_base)
    L_base  = zeros(TFloat, 1,1)
    L = JACC.array(L_base)
    @info "BEFORE"
    if true # this true needs to be changed to lbuild_laplace_matrix
        @info "AFTER"
        Le_base = zeros(TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
        Le = JACC.array(Le_base)
        # I think the issue is here


        n_nodes = size(basis.dψ, 1)
        total_iterations= (n_nodes * n_nodes * (Q + 1))

        @info N
        JACC.parallel_for(total_iterations, build_laplace_matrix_jacc_parallel_all!, Le, basis.dψ, ω, Q)
        #k = build_laplace_matrix_jacc_parallel_all!(backend)
        #k(Le, basis.dψ, ω, TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
        #synchronize()
        L_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        L = JACC.array(L_base)
        #synchronize()
        
        n_elements = size(mesh.connijk, 1)
        n_local_i = size(mesh.connijk, 2)
        n_local_j = size(mesh.connijk, 3)

        total_iterations = (n_elements * n_local_i * n_local_j * mesh.ngl)
        JACC.parallel_for(total_iterations, DSS_laplace_jacc!, L, Le, mesh.connijk, ω, ω, mesh.ngl, mesh.ngl, metrics.dξdx, metrics.dydη, metrics.dηdy, metrics.dxdξ;
        #=ndrange = (mesh.nelem*mesh.ngl,
        mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl)=#)
       
        
        #k = DSS_laplace_gpu!(backend)
        #k(L, Le, connijk, ω, ω, mesh.ngl, mesh.ngl, metrics.dξdx, metrics.dydη, metrics.dηdy, metrics.dxdξ;
        #ndrange = (mesh.nelem*mesh.ngl,
        #mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
        
       # synchronize()
    end
    
    De_base = zeros(TFloat, 1, 1)
    De = JACC.array(De_base)
    D_base  = zeros(TFloat, 1,1)
    D = JACC.array(D_base)
    if lbuild_differentiation_matrix
        De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh,  N, Q, TFloat)
        if ldss_differentiation
            D  = DSS_generic_matrix(SD, De, mesh, TFloat)
        end
    end
    # Should be this
    # ("We are here 2", 126419.75308725929, 0.00040000000000316936, 126419.75308725929, 0.0, 0.0, 0.0, 0.0)

    @info "We are here 2", maximum(M), maximum(Minv), maximum(Me), maximum(De), maximum(D), maximum(Le), maximum(L)
    return (; Me, De, Le, M, Minv, D, L, M_surf_inv=Minv, M_edge_inv=Minv)
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

function matrix_wrapper_laguerre(::ContGal, SD, QT, basis, ω, mesh, metrics, N, Q, TFloat; ldss_laplace=false, ldss_differentiation=false, backend = CPU(), interp=nothing)

    lbuild_differentiation_matrix = false
    lbuild_laplace_matrix = false    
    if (ldss_differentiation) lbuild_differentiation_matrix = true end
    if (ldss_laplace) lbuild_laplace_matrix = true end

    if typeof(SD) == NSD_1D
        Me_base = zeros(TFloat, (N+1)^2, Int64(mesh.nelem))
        Me = JACC.array(Me_base)
    elseif typeof(SD) == NSD_2D
        Me_base = zeros(TFloat, (N+1)^2, (N+1)^2, Int64(mesh.nelem))
        Me = JACC.array(Me_base)
    end
    #=if (backend == CPU())
        build_mass_matrix!(Me, SD, QT, basis[1].ψ, ω[1], mesh.nelem, metrics[1].Je, mesh.Δx, N, Q, TFloat)
    =#
    if (SD == NSD_1D())
        k = build_mass_matrix_1d_gpu!(backend, (N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, Q; ndrange = (mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))
    elseif (SD == NSD_2D())
        k = build_mass_matrix_2d_gpu!(backend, (N+1, N+1))
        k(Me, basis[1].ψ, ω[1], metrics[1].Je, N, Q;ndrange = (mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
    else
        nothing
    end


    if typeof(SD) == NSD_1D
        M_lag_base = zeros(TFloat, Int64(mesh.ngr*mesh.ngr), Int64(mesh.nelem_semi_inf))
        M_lag = JACC.array(M_lag_base)
    elseif typeof(SD) == NSD_2D
        M_lag_base = zeros(TFloat, Int64(mesh.ngl*mesh.ngr), Int64(mesh.ngl*mesh.ngr), Int64(mesh.nelem_semi_inf))
        M_lag = JACC.array(M_lag_base)
    end
    if (QT == Exact() && inputs[:llump] == false)
        M_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        M = JACC.array(M_base)
        Minv_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        Minv = JACC.array(Minv_base)
    else
        M_base = zeros(TFloat, Int64(mesh.npoin))
        M = JACC.array(M_base)
        Minv_base = zeros(TFloat, Int64(mesh.npoin))
        Minv = JACC.array(Minv_base)
    end
        if (typeof(SD) == NSD_1D)
            k = build_mass_matrix_1d_gpu!(backend)
            k(M_lag, basis[2].ψ, ω[2], metrics[2].Je, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngr), workgroupsize = (mesh.ngr))
            
            connijk_base = zeros(TInt, Int64(mesh.nelem), N+1)
            connijk = JACC.array(connijk_base)
            k1 = DSS_Mass_gpu_1D!(backend,(N+1))
            k1(M, Me, connijk;ndrange =(mesh.nelem*mesh.ngl), workgroupsize = (mesh.ngl))

            JACC.synchronize(backend)

            connijk_lag_base = zeros(TInt, Int64(mesh.nelem_semi_inf), Int64(mesh.ngr))
            connijk_lag = JACC.array(connijk_lag_base)
            k2 = DSS_Mass_gpu_1D!(backend)
            k2(M, M_lag, connijk_lag;ndrange = (mesh.nelem_semi_inf*mesh.ngr), workgroupsize = (mesh.ngr))

            JACC.synchronize(backend)

        elseif (typeof(SD) == NSD_2D)
            k = build_mass_matrix_Laguerre_2d_gpu!(backend)
            k(M_lag, basis[1].ψ, basis[2].ψ, ω[1], ω[2], metrics[2].Je, mesh.ngl, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngl,mesh.ngr), workgroupsize = (mesh.ngl,mesh.ngr))
            
            connijk_base = zeros(TInt, Int64(mesh.nelem), N+1, N+1)
            connijk= JACC.array(connijk_base)
            k1 = DSS_Mass_gpu_2D!(backend,(N+1,N+1))
            k1(M,Me,connijk,mesh.nelem, mesh.npoin, N;ndrange =(mesh.nelem*mesh.ngl,mesh.ngl), workgroupsize = (mesh.ngl,mesh.ngl))
            
            JACC.synchronize(backend)

            connijk_lag_base = zeros(TInt, Int64(mesh.nelem_semi_inf), Int64(mesh.ngl), Int64(mesh.ngr))
            connijk_lag = JACC.array(connijk_lag_base)
            k2 = DSS_mass_Laguerre_gpu_2D!(backend)
            k2(M, M_lag, connijk_lag, mesh.ngl, mesh.ngr;ndrange = (mesh.nelem_semi_inf*mesh.ngl,mesh.ngr), workgroupsize = (mesh.ngl,mesh.ngr))
            
            JACC.synchronize(backend)
        end
    

    mass_inverse!(Minv, M, QT)

    Le_base = zeros(TFloat, 1, 1)
    Le = JACC.array(Le_base)
    L_base  = zeros(TFloat, 1, 1)
    L = JACC.array(L_base)
    Le_lag_base = zeros(TFloat, 1, 1)
    Le_lag = JACC.array(Le_lag_base)
    
    if lbuild_laplace_matrix
            Le_base = zeros(TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
            Le = JACC.array(Le_base)
            Le_lag_base = zeros(TFloat, Int64(mesh.ngr), Int64(mesh.ngr))
            Le_lag = JACC.array(Le_lag_base)

            n_nodes = size(dψ, 1)
            total_iterations = (n_nodes * n_nodes * (Q + 1))

            parallel_for(total_iterations, build_laplace_matrix_jacc_parallel_all!, Le, basis[1].dψ, ω[1], TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
            #k = build_laplace_matrix_jacc_parallel_all!(backend)
            #k(Le, basis[1].dψ, ω[1], TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
            JACC.synchronize(backend)
            
            parallel_for(N, build_laplace_matrix_jacc_parallel_all!, Le_lag, basis[2].dψ, ω[2], TInt(mesh.ngr-1); ndrange = (mesh.ngr, mesh.ngr))
            
            #k(Le_lag, basis[2].dψ, ω[2], TInt(mesh.ngr-1); ndrange = (mesh.ngr, mesh.ngr))
            L_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
            L = JACC.array(L_base)
            JACC.synchronize(backend)
            
            
            
            
            k = DSS_laplace_gpu!(backend)
            k(L, Le, connijk, ω[1], ω[1], mesh.ngl, mesh.ngl, metrics[1].dξdx, metrics[1].dydη, metrics[1].dηdy, metrics[1].dxdξ; 
            ndrange = (mesh.nelem*mesh.ngl, mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
            KernelAbstractions.synchronize(backend)
                      if lbuild_laplace_matrix
                    
                        Le_base = zeros(TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                        Le = JACC.array(Le_base)
                        Le_lag_base = zeros(TFloat, Int64(mesh.ngr), Int64(mesh.ngr))
                        Le_lag = JACC.array(Le_lag_base)
            
                        k = build_laplace_matrix_gpu!(backend)
                        k(Le, basis[1].dψ, ω[1], TInt(mesh.ngl-1); ndrange=(mesh.ngl,mesh.ngl))
                        JACC.synchronize(backend)
                        event = k(Le_lag, basis[2].dψ, ω[2], TInt(mesh.ngr-1); ndrange = (mesh.ngr, mesh.ngr))
                        L_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
                        L = JACC.array(L_base)
                        wait(event)
                        JACC.synchronize(backend)
                        k = DSS_laplace_gpu!(backend)
                        k(L, Le, connijk, ω[1], ω[1], mesh.ngl, mesh.ngl, metrics[1].dξdx, metrics[1].dydη, metrics[1].dηdy, metrics[1].dxdξ; 
                          ndrange = (mesh.nelem*mesh.ngl, mesh.ngl), workgroupsize = (mesh.ngl, mesh.ngl))
                          JACC.synchronize(backend)
                        k = DSS_laplace_gpu_lag!(backend)
                        k(L, Le_lag, connijk_lag, ω[1], ω[2], mesh.ngl, mesh.ngr, metrics[2].dηdx, metrics[1].dydη, metrics[1].dηdy, metrics[2].dxdη; 
                          ndrange = (mesh.nelem_semi_inf*mesh.ngr, mesh.ngl), workgroupsize = (mesh.ngr, mesh.ngl))
                          JACC.synchronize(backend)
                    end
                k = DSS_laplace_gpu_lag!(backend)
            k(L, Le_lag, connijk_lag, ω[1], ω[2], mesh.ngl, mesh.ngr, metrics[2].dηdx, metrics[1].dydη, metrics[1].dηdy, metrics[2].dxdη; 
              ndrange = (mesh.nelem_semi_inf*mesh.ngr, mesh.ngl), workgroupsize = (mesh.ngr, mesh.ngl))
              JACC.synchronize(backend)
        end
    

    De_base = zeros(TFloat, 1, 1)
    De = JACC.array(De_base)
    D_base = zeros(TFloat, 1, 1)
    D = JACC.array(D_base)
    if lbuild_differentition_matrix
        De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh,  N, Q, TFloat)
        if ldss_differentiation
            D  = DSS_generic_matrix(SD, De, mesh, TFloat)
        end
    end

    @info "We are here 3", maximum(M), maximum(Minv), maximum(Me), maximum(De), maximum(D), maximum(Le), maximum(L)

    return (; Me, De, Le, M, Minv, D, L, M_surf_inv=Minv, M_edge_inv=Minv)
end


@kernel function diagm_gpu!(Minv_d, Minv)
    ip = @index(Global, Linear)

    Minv_d[ip,ip] = Minv[ip]
end

@kernel function add_to_diag!(M, val)
    ip = @index(Global, Linear)

    M[ip,ip] += val
end

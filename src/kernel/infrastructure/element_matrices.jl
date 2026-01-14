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
    if(SD == NSD_1D())
        total_iterations = neqs * nelem * ngl
       # @info "Were here"
        JACC.parallel_for(total_iterations, DSS_rhs_jacc!, RHS, rhs_el, connijk, nelem, ngl, neqs, SD, method)
       # @info "And here
    elseif(SD == NSD_2D())
        total_iterations = neqs * nelem * ngl * ngl
        #@info "Were here 2"
        JACC.parallel_for(total_iterations, DSS_rhs_jacc!, RHS, rhs_el, connijk, nelem, ngl, neqs, SD, method)
       # @info "And here 2"
    else
        total_iterations = ngl * ngl * ngl * nelem * neqs
        JACC.parallel_for(total_iterations, DSS_rhs_jacc!, RHS, rhs_el, connijk, nelem, ngl, neqs, SD, method)

        #error("DSS_rhs! not implemented for 3D")
    end
end

function DSS_global_mass!(SD, M, ip2gip, gip2owner, parts, npoin, gnpoin)

    if SD == NSD_1D()
        return nothing
    end
    
    pM = setup_assembler(SD, M, ip2gip, gip2owner)
    
    @time assemble_mpi!(M,pM)

    return pM
    
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

function build_mass_matrix_3d_jacc!(idx, Me, ψ, ω, Je, N, Q, nelem)
    # --- 1. Define Dimensionality Constants ---
    n_basis_1d = N + 1
    n_quad_1d = Q + 1
    n_basis_3d = n_basis_1d^3
    n_quad_3d = n_quad_1d^3

    # --- 2. Decompose the Global Linear Index 'idx' ---
    idx0 = idx - 1

    # Extract the flat indices for the quadrature point, column (J), row (I), and element (ie)
    quad_idx_flat = idx0 % n_quad_3d
    remainder = idx0 ÷ n_quad_3d

    J_flat = remainder % n_basis_3d
    remainder = remainder ÷ n_basis_3d

    I_flat = remainder % n_basis_3d
    ie = (remainder ÷ n_basis_3d) + 1

    # --- 3. Decompose Flat Indices into 3D Coordinates ---
    # Decompose quadrature index into (m, n, o)
    o_idx = quad_idx_flat % n_quad_1d; temp_quad = quad_idx_flat ÷ n_quad_1d
    n_idx = temp_quad % n_quad_1d; m_idx = temp_quad ÷ n_quad_1d
    m = m_idx + 1; n = n_idx + 1; o = o_idx + 1

    # Decompose column index J into (p, q, r)
    p_idx = J_flat % n_basis_1d; temp_J = J_flat ÷ n_basis_1d
    q_idx = temp_J % n_basis_1d; r_idx = temp_J ÷ n_basis_1d
    p = p_idx + 1; q = q_idx + 1; r = r_idx + 1

    # Decompose row index I into (i_x, i_y, i_z)
    i_x_idx = I_flat % n_basis_1d; temp_I = I_flat ÷ n_basis_1d
    i_y_idx = temp_I % n_basis_1d; i_z_idx = temp_I ÷ n_basis_1d
    i_x = i_x_idx + 1; i_y = i_y_idx + 1; i_z = i_z_idx + 1

    # --- 4. Calculate Contribution and Atomically Add ---
    ωmno = ω[m] * ω[n] * ω[o]
    Jmnoe = Je[ie, m, n, o]
    
    # Basis function for column J at point (m,n,o)
    ψIK = ψ[p, m] * ψ[q, n] * ψ[r, o]
    # Basis function for row I at point (m,n,o)
    ψJK = ψ[i_x, m] * ψ[i_y, n] * ψ[i_z, o]
    
    term = ωmno * Jmnoe * ψIK * ψJK

    # Get the final 1-based flat indices for writing to Me
    I = I_flat + 1
    J = J_flat + 1
    
    # Atomically add the result. This is crucial because many threads
    # (for different quadrature points) will be writing to the same Me[I, J, ie].
    @Atomix.atomic Me[I, J, ie] += term
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

    @info "Inside Laguerre Mass Matrix Build"
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

        Atomix.@atomic Le[i,j] = Le[i,j] + term
        @info "Le[i,j]", i, j, term, Le[i,j]
    #end
end

#JACCIFY
function DSS_mass_Laguerre_gpu_2D!(M, Mel_lag, connijk_lag, ngl, ngr)
    @info "Inside Laguerre Mass Matrix Assembly"
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

#=
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
=#

function DSS_Mass_jacc_2D!(idx, M, Mel, conn, nelem, npoin, N)
    # --- 1. Calculate Constants ---
    n_local_nodes_per_dim = N + 1
    total_local_nodes = n_local_nodes_per_dim^2

    # --- 2. Decompose Global Linear Index ---
    # The single 'idx' now represents a unique (element, local_row, local_col) tuple.
    
    # Get the element index 'ie'.
    ie = (idx - 1) ÷ (total_local_nodes^2) + 1

    # Get the flattened local row 'I' and column 'J' indices.
    local_flat_index = (idx - 1) % (total_local_nodes^2)
    J = local_flat_index ÷ total_local_nodes + 1
    I = local_flat_index % total_local_nodes + 1

    # --- 3. Get Global Index and Perform Update ---
    
    # We need the 2D local coordinates (m, n) of the 'row' index 'I'
    # to find its corresponding global index in the connectivity array.
    m = (I - 1) % n_local_nodes_per_dim + 1
    n = (I - 1) ÷ n_local_nodes_per_dim + 1

    # Get the global index 'IP' for the destination vector M.
    IP = conn[ie, m, n]

    # Atomically add the single contribution Mel[I, J, ie].
    @Atomix.atomic M[IP] += Mel[I, J, ie]
end


function DSS_Mass_jacc_3D!(idx, M, Mel, conn, nelem, N)
    # --- 1. Define Dimensionality Constants ---
    n_basis_1d = N + 1
    n_basis_3d = n_basis_1d^3

    # --- 2. Decompose the Global Linear Index 'idx' ---
    # The single index now represents a unique (element, row, column) tuple.
    idx0 = idx - 1

    # Extract the flat column index (J), row index (I), and element index (ie)
    J_flat = idx0 % n_basis_3d
    remainder = idx0 ÷ n_basis_3d

    I_flat = remainder % n_basis_3d
    ie = (remainder ÷ n_basis_3d) + 1

    # --- 3. Decompose Flat Row Index to Get Global Index ---
    # Get the 3D local coordinates (k, m, n) of the 'row' index 'I'
    # to find its corresponding global index in the connectivity array.
    k_idx = I_flat % n_basis_1d
    temp = I_flat ÷ n_basis_1d
    m_idx = temp % n_basis_1d
    n_idx = temp ÷ n_basis_1d
    k = k_idx + 1; m = m_idx + 1; n = n_idx + 1

    # Get the global 'row' index 'IP' from the connectivity map.
    IP = conn[ie, k, m, n]

    # --- 4. Perform the Atomic Update ---
    # Get the final 1-based flat indices for accessing the element matrix.
    I = I_flat + 1
    J = J_flat + 1
    
    # Atomically add the contribution.
    @Atomix.atomic M[IP] += Mel[I, J, ie]
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

        Atomix.@atomic L[ip, jp] += term
      # JACC.atomic_add!(L, ip, jp, term) MADE THIS NOT ATOMIC
   # end
end

function DSS_laplace_jacc_lag!(L, Lel, connijk, ωx, ωy, nx, ny, dηdx_lag, dydη, dηdy, dxdη_lag)
    @info "Inside Laguerre Laplace Assembly"
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
    
    pM = setup_assembler(SD, M, ip2gip, gip2owner)
    
    @time assemble_mpi!(M,pM)

    return pM
    
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

    Atomix.@atomic RHS[I, ieq] += rhs_el[iel, i, ieq]
end

#N = total_iterations = neqs * nelem * ngl * ngl
 
function DSS_rhs_jacc!(N, RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_2D, ::ContGal)
    #@info "Inside 2D RHS JACC"
    ieq = (N - 1) ÷ (nelem * ngl*ngl) + 1
    rem1 = (N - 1) % (nelem * ngl*ngl)
    iel = rem1 ÷ (ngl*ngl) + 1
    rem2 = rem1 % (ngl*ngl)
    i = rem2 ÷ ngl + 1
    j = rem2 % ngl + 1

    I = connijk[iel, i, j]

    Atomix.@atomic RHS[I, ieq] += rhs_el[iel, i, j, ieq]
end


function DSS_rhs_jacc!(N, RHS, rhs_el, connijk, nelem, ngl, neqs, ::NSD_3D, ::ContGal)
    idx0 = N - 1

    # Total number of nodes per element
    ngl3 = ngl * ngl * ngl
    
    # Get equation index
    ieq = (idx0 ÷ (nelem * ngl3)) + 1
    rem1 = idx0 % (nelem * ngl3)
    
    # Get element index
    iel = (rem1 ÷ ngl3) + 1
    rem2 = rem1 % ngl3
    
    # Decompose the local 3D node index into i, j, k
    i_idx = rem2 % ngl
    rem3 = rem2 ÷ ngl
    j_idx = rem3 % ngl
    k_idx = rem3 ÷ ngl
    
    # Convert from 0-based to 1-based indices
    i = i_idx + 1
    j = j_idx + 1
    k = k_idx + 1

    # --- Perform the DSS Operation ---
    # 1. Find the global node index
    I = connijk[iel, i, j, k]

    # 2. Atomically add the local contribution to the global RHS
    @Atomix.atomic RHS[I, ieq] += rhs_el[iel, i, j, k, ieq]
end


function DSS_rhs_laguerre_jacc!(RHS, rhs_el, connijk_lag, nelem_semi_inf, ngl, ngr, neqs, ::NSD_1D, ::ContGal)
    total_tasks = neqs * nelem_semi_inf * ngr

    JACC.parallel_for(1:total_tasks) do idx_linear
        ieq = (idx_linear - 1) ÷ (nelem_semi_inf * ngr) + 1
        rem = (idx_linear - 1) % (nelem_semi_inf * ngr)
        iel = rem ÷ ngr + 1
        i = rem % ngr + 1

        I = connijk_lag[iel, i, 1]

        Atomix.@atomic RHS[I,ieq] += rhs_el[iel,i,ieq]
    end
end


function divide_by_mass_matrix_jacc_all_parallel!(RHS, RHSaux, Minv::AbstractMatrix, neqs, npoin, ::ContGal)
    @info "Inside divide by mass matrix all parallel"
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

        num_elements = size(Me, 3)
        total_iterations = num_elements * (N + 1)^4 * (Q + 1)^2

        JACC.parallel_for(total_iterations, build_mass_matrix_2d_jacc!, Me, basis.ψ, ω, metrics.Je, N, Q)
    elseif typeof(SD) == NSD_3D
        #Me = JACC.zeros(Float64, n_basis_3d, n_basis_3d, nelem)
        Me_base = zeros(Float64, (N + 1)^3, (N + 1)^3, mesh.nelem)
        Me = JACC.array(Me_base)
        
        # 1. Calculate the total number of iterations
        n_basis_1d = N + 1
        n_quad_1d = Q + 1

        n_basis_3d = n_basis_1d^3
        n_quad_3d = n_quad_1d^3

        # The total number of threads is the product of all loop dimensions
        total_iterations = mesh.nelem * n_basis_3d * n_basis_3d * n_quad_3d

        # 2. Launch the single, fully flattened parallel kernel
        JACC.parallel_for(
            total_iterations,
            build_mass_matrix_3d_jacc!,
            Me, basis.ψ, ω, metrics.Je, N, Q, mesh.nelem
        )
        #JACC parallel for for build mass matrix 3D
    
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
        DSS_Mass_jacc_1D!(M, SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, N, TFloat; llump=inputs[:llump])
    elseif SD == NSD_2D()
        n_local_nodes_per_dim = N + 1
        total_local_nodes = n_local_nodes_per_dim^2
        total_iterations = mesh.nelem * total_local_nodes * total_local_nodes
        JACC.parallel_for(total_iterations, DSS_Mass_jacc_2D!, M, Me, mesh.connijk, mesh.nelem, mesh.npoin, N)
        

    elseif SD == NSD_3D()
        n_basis_1d = N + 1
        n_basis_3d = n_basis_1d^3

        # Total iterations = (num_elements) * (nodes_per_element) * (nodes_per_element)
        total_iterations = mesh.nelem * n_basis_3d * n_basis_3d

        JACC.parallel_for(
            total_iterations,
            DSS_Mass_jacc_3D!,
            M, Me, mesh.connijk, mesh.nelem, N
        )

    end

    pM = DSS_global_mass!(SD, M, mesh.ip2gip, mesh.gip2owner, mesh.parts, mesh.npoin, mesh.gnpoin)

    mass_inverse!(Minv, M, QT)

    Le_base = zeros(TFloat, 1, 1)
    Le = JACC.array(Le_base)
    L_base  = zeros(TFloat, 1,1)
    L = JACC.array(L_base)
    if lbuild_differentiation_matrix # this true needs to be changed to lbuild_laplace_matrix
        Le_base = zeros(TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
        Le = JACC.array(Le_base)

        n_nodes = size(basis.dψ, 1)
        total_iterations= (n_nodes * n_nodes * (Q + 1))

        @info N
        JACC.parallel_for(total_iterations, build_laplace_matrix_jacc_parallel_all!, Le, basis.dψ, ω, Q)
        
        L_base = zeros(TFloat, Int64(mesh.npoin), Int64(mesh.npoin))
        L = JACC.array(L_base)
        
        n_elements = size(mesh.connijk, 1)
        n_local_i = size(mesh.connijk, 2)
        n_local_j = size(mesh.connijk, 3)

        total_iterations = (n_elements * n_local_i * n_local_j * mesh.ngl)
        JACC.parallel_for(total_iterations, DSS_laplace_jacc!, L, Le, mesh.connijk, ω, ω, mesh.ngl, mesh.ngl, metrics.dξdx, metrics.dydη, metrics.dηdy, metrics.dxdξ;)
       
        
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
    return (; Me, De, Le, M, Minv,pM, D, L, M_surf_inv=Minv, M_edge_inv=Minv)
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

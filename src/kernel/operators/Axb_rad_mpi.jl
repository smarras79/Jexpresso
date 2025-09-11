using Krylov

mutable struct MatvecCache
    # number of points on each process
    proc_npoin :: Vector{Int}
    # Local 2 global indices for each process
    proc_ip2gip :: Vector{Vector{Int}}
    # Matvec contributions in local indexing
    proc_vec    :: Vector{Vector{Float64}}
end

mutable struct MatvecCacheSparse
    # number of points on each process
    proc_npoin  :: Vector{Int}
    # Local 2 global indices for each process
    proc_ip2gip :: Vector{Vector{Int}}
    # Matvec sparse contributions
    proc_vec    :: Vector{Vector{Float64}}
    # Indices of non-zero contributions by process
    proc_sp_ip  :: Vector{Vector{Int}}
    # Indices of non-zero contribution of local process
    sp_ip       :: Vector{Int}
    # Values of non-zero contributions of local process
    sp_val      :: Vector{Float64}
    # Number of non-zero contributions from each process
    nnzeros     :: Vector{Int}
end

function setup_Matvec_assembler(index_a,gip2owner)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = MPI.Allreduce(maximum(index_a), MPI.MAX, comm)

    npoin = size(index_a,1)
    proc_npoin = MPI.Allgather(npoin,comm)

    proc_ip2gip = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]
    
    @info "sizes for proc_ip2gip", size(proc_ip2gip[1]), size(proc_ip2gip[2]), rank
    
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if i != rank
            push!(requests, MPI.Isend(index_a, i, 1, comm))
        end
        if i != rank
            push!(requests, MPI.Irecv!(proc_ip2gip[i+1], i, 1, comm))
        end
    end

    MPI.Waitall!(requests)

    proc_ip2gip[rank+1] .= index_a

    proc_vec = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    cache = MatvecCache(proc_npoin, proc_ip2gip, proc_vec)

    return cache
end

function setup_MatvecSparse_assembler(index_a,gip2owner)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = MPI.Allreduce(maximum(index_a), MPI.MAX, comm)

    npoin = size(index_a,1)
    proc_npoin = MPI.Allgather(npoin,comm)

    proc_ip2gip = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    @info "sizes for proc_ip2gip", size(proc_ip2gip[1]), size(proc_ip2gip[2]), rank

    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if i != rank
            push!(requests, MPI.Isend(index_a, i, 1, comm))
        end
        if i != rank
            push!(requests, MPI.Irecv!(proc_ip2gip[i+1], i, 1, comm))
        end
    end

    MPI.Waitall!(requests)

    proc_ip2gip[rank+1] .= index_a

    proc_vec = [Vector{Float64}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    proc_sp_ip = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    sp_ip  = zeros(Int, npoin)

    sp_val = zeros(Float64, npoin)

    nnzeros = zeros(Int, rank_sz)

    cache = MatvecCacheSparse(proc_npoin, proc_ip2gip, proc_vec, proc_sp_ip, sp_ip, sp_val, nnzeros)

    return cache
end

function assemble_mpi_matvec!(y_local, y, cache::MatvecCache)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    #=iter = 1
    for ip=1:cache.proc_npoin[rank+1]
        val = y_local[ip]

        if (abs(val) > 1e-14)
            iter += 1
        end
    end
    
    @info rank, iter-1
    =#
    requests = MPI.Request[]
    @inbounds for i in 0:rank_sz-1
        fill!(cache.proc_vec[i+1],zero(Float64))
        if i != rank
            push!(requests, MPI.Isend(y_local, i, 1, comm))
        end
        if i != rank
            push!(requests, MPI.Irecv!(cache.proc_vec[i+1], i, 1, comm))
        end
    end

    MPI.Waitall!(requests)
    
    @inbounds cache.proc_vec[rank+1] .= y_local

    @inbounds for i in 0:rank_sz-1
        npoin = cache.proc_npoin[i+1]
        for j=1:npoin
            gip = cache.proc_ip2gip[i+1][j]
            y[gip] += cache.proc_vec[i+1][j]
        end
    end

end

function assemble_mpi_matvec_sparse!(y_local, y, cache::MatvecCacheSparse)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    npoin = size(y_local,1)
    ## Find sparse indices and corresponding values
    iter = 1
    @inbounds for ip=1:npoin
        val = y_local[ip]
        
        if (abs(val) > 1e-7)
            cache.sp_ip[iter] = ip
            cache.sp_val[iter] = val
            iter += 1
        end
    end

    nnzero = iter - 1
    #=requests_nzeros = MPI.Request[]
    
    @inbounds for i in 0:rank_sz-1
        if i != rank
            push!(requests_nzeros, MPI.Isend(nnzero, i, 1, comm))
        end
        if i != rank
            push!(requests_nzeros, MPI.Irecv!(@view(cache.nnzeros[i+1]), i, 1, comm))
        end
    end
    
    MPI.Waitall!(requests_nzeros)

    cache.nnzeros[rank+1] = nnzero
    =#
    @time cache.nnzeros = MPI.Allgather(nnzero,comm)

    requests_indices = MPI.Request[]


    @inbounds for i in 0:rank_sz-1
        fill!(cache.proc_sp_ip[i+1],zero(Int))
        idx_end = cache.nnzeros[i+1]
        #send_buffer = @view(cache.sp_ip[1:nnzero])
        #recv_buffer = @view(cache.proc_sp_ip[i+1][1:idx_end])
        resize!(cache.proc_sp_ip[i+1],idx_end)
        if i != rank
            push!(requests_indices, MPI.Isend(cache.sp_ip[1:nnzero], i, 1, comm))
        end
        if i != rank
            push!(requests_indices, MPI.Irecv!(cache.proc_sp_ip[i+1], i, 1, comm))
        end
    end

    #if (rank == 0)
    #    @info cache.proc_sp_ip[2]
    #end
    

    @inbounds cache.proc_sp_ip[rank+1][1:nnzero] .= cache.sp_ip[1:nnzero]

    requests = MPI.Request[]
    @inbounds for i in 0:rank_sz-1
        fill!(cache.proc_vec[i+1],zero(Float64))
        idx_end = cache.nnzeros[i+1]
        #send_buffer = @view(cache.sp_val[1:nnzero])
        #recv_buffer = @view(cache.proc_vec[i+1][1:idx_end])
        resize!(cache.proc_vec[i+1],idx_end)
        if i != rank
            push!(requests, MPI.Isend(cache.sp_val[1:nnzero], i, 1, comm))
        end
        if i != rank
            push!(requests, MPI.Irecv!(cache.proc_vec[i+1], i, 1, comm))
        end
    end

    MPI.Waitall!(requests)
    MPI.Waitall!(requests_indices)
    #if (rank == 0)
    #    @info cache.proc_proc_vec[2]
    #end
    @inbounds cache.proc_vec[rank+1][1:nnzero] .= cache.sp_val[1:nnzero]
    

    @inbounds for i in 0:rank_sz-1
        idx_end = cache.nnzeros[i+1]
        for j=1:idx_end
            ip = cache.proc_sp_ip[i+1][j]
            gip = cache.proc_ip2gip[i+1][ip]
            y[gip] += cache.proc_vec[i+1][j]
        end
    end

end

function create_parallel_linear_operator(A_local::AbstractMatrix{T},
                                        ip2gip::AbstractVector{Int},
                                        gip2owner::AbstractVector{Int},
                                        npoin::Int,
                                        gnpoin::Int,
                                        comm::MPI.Comm = MPI.COMM_WORLD) where T

    # Forward operator: y = A*x
    nprocs = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    tolerance = 1e-14
    y_local = zeros(Float64,gnpoin)
    y_temp = zeros(Float64,npoin)
    A_rows = A_local.rowval
    A_ranges = zeros(Int64, npoin,2)
    x_local = zeros(Float64,npoin)
    local_gips = zeros(Int64,npoin)
    local_values = zeros(Float64,npoin)
    for ip=1:npoin
        A_ranges[ip,1], A_ranges[ip,2] = first(nzrange(A_local, ip)), last(nzrange(A_local, ip))
    end
    
    pM = setup_MatvecSparse_assembler(ip2gip,gip2owner)
    
    function matvec!(y, x)
        fill!(y_local,zero(Float64))
        #Perform sparse A*x and store partial result
        #y_local = A_local*x
        for ip=1:npoin
            gip = ip2gip[ip]
            x_local[ip] = x[gip]
        end
        y_temp .= A_local * x_local
        if (nprocs == 0)
            for ip=1:npoin
                gip = ip2gip[ip]
                y_local[gip] = y_temp[ip]
            end
        else
            assemble_mpi_matvec_sparse!(y_temp, y_local, pM)
        end
        copyto!(y, y_local)  # y_local now contains the reduced result
        return y

        #=    #Naive version build global size resulting vector
        fill!(y_local,zero(Float64))
        #Perform sparse A*x and store partial result
        #y_local = A_local*x
        for ip=1:npoin
            gip = ip2gip[ip]
            x_local[ip] = x[gip]
        end
        y_temp .= A_local * x_local
        for ip=1:npoin
            gip = ip2gip[ip]
            y_local[gip] = y_temp[ip]
        end
        # Reduce-sum across all processes (each process contributes to global columns)
        MPI.Allreduce!(y_local, +, comm)
        # Copy to output vector
        copyto!(y, y_local)  # y_local now contains the reduced result
        return y=#
    end

    # Transpose operator: y = A'*x
    function rmatvec!(y, x)
        fill!(y_local,zero(Float64))
        for ip=1:npoin
            gip = ip2gip[ip]
            x_local[ip] = x[gip]
        end
        y_temp .= A_local' * x_local
        if (nprocs == 0)
            for ip=1:npoin
                gip = ip2gip[ip]
                y_local[gip] = y_temp[ip]
            end
        else
            assemble_mpi_matvec_sparse!(y_temp, y_local, pM)
        end
        copyto!(y, y_local)  # y_local now contains the reduced result
        return y

        #=fill!(y_local,zero(Float64))
        for ip=1:npoin
            gip = ip2gip[ip]
            x_local[ip] = x[gip]
        end
        y_temp .= A_local' * x_local

        for ip=1:npoin
            gip = ip2gip[ip]
            y_local[gip] = y_temp[ip]
        end
        # Reduce-sum across all processes (each process contributes to global columns)
        MPI.Allreduce!(y_local, +, comm)
        # Copy to output vector
        copyto!(y, y_local)  # y_local now contains the reduced result
        return y=#

    end

    # Create LinearOperator
    return LinearOperator{T}(gnpoin, gnpoin, false, false, matvec!, rmatvec!, rmatvec!)
end

function solve_parallel_lsqr(ip2gip, gip2owner, A_local, b, gnpoin, npoin; tol::Float64 = 1e-7)
   
     # Create parallel linear operator (better than AbstractMatrix)
    A_parallel = create_parallel_linear_operator(A_local, ip2gip, gip2owner, npoin, gnpoin)

    # Gather the full RHS vector with proper Allgatherv
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    nprocs = MPI.Comm_size(MPI.COMM_WORLD)

    #build contributions to global RHS
    b_global = zeros(Float64,gnpoin)
    for ip=1:npoin
        gip = ip2gip[ip]
        b_global[gip] = b[ip]
    end
    #Use allreduce to build global rhs vector
    MPI.Allreduce!(b_global, +, MPI.COMM_WORLD)

    # Solve using LSQR
    if rank == 0
        println("Starting parallel LSQR solve...")
        println("Global problem size: $(gnpoin) × $(gnpoin)")
        println("Number of processes: $(MPI.Comm_size(MPI.COMM_WORLD))")
    end
    maxiter = gnpoin
    @info tol, maxiter
    x, stats = Krylov.lsqr(A_parallel, b_global;
                   atol = tol,
                   rtol = tol,
                   btol = tol,
                   etol = tol,
                   axtol = tol,
                   itmax = maxiter,
                   verbose = (rank == 0) ? 1 : 0)  # Only print on rank 0

    #return to local indexing
    x_local = zeros(Float64,npoin)
    for ip=1:npoin
        gip = ip2gip[ip]
        x_local[ip] = x[gip]
    end

    @info stats
    return x_local
end


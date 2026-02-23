using Krylov

mutable struct MatvecCacheSparse
    # number of points on each process
    proc_npoin      :: Vector{Int}
    # Local 2 global indices for each process
    proc_ip2gip     :: Vector{Vector{Int}}
    # Matvec sparse contributions
    proc_vec        :: Vector{Vector{Float64}}
    # Indices of non-zero contributions by process
    proc_sp_ip      :: Vector{Vector{Int}}
    # Indices of non-zero contribution of local process
    sp_ip           :: Vector{Int}
    # Values of non-zero contributions of local process
    sp_val          :: Vector{Float64}
    # Number of non-zero contributions from each process
    nnzeros         :: Vector{Int}
    #Preallocated request arrays
    request         :: MPI.MultiRequest
    request_indices :: MPI.MultiRequest
end

function setup_MatvecSparse_assembler(index_a,gip2owner)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = MPI.Allreduce(maximum(index_a), MPI.MAX, comm)

    npoin = size(index_a,1)
    proc_npoin = MPI.Allgather(npoin,comm)

    proc_ip2gip = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    nreq = 0
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if i != rank
            push!(requests, MPI.Isend(index_a, i, 1, comm))
            nreq += 1
        end
        if i != rank
            push!(requests, MPI.Irecv!(proc_ip2gip[i+1], i, 1, comm))
            nreq += 1
        end
    end

    MPI.Waitall!(requests)

    proc_ip2gip[rank+1] .= index_a

    proc_vec = [Vector{Float64}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    proc_sp_ip = [Vector{Int}(undef, proc_npoin[i+1]) for i in 0:rank_sz-1]

    sp_ip  = zeros(Int, npoin)

    sp_val = zeros(Float64, npoin)

    nnzeros = zeros(Int, rank_sz)

    cache = MatvecCacheSparse(proc_npoin, proc_ip2gip, proc_vec, proc_sp_ip, sp_ip, sp_val, nnzeros, MPI.MultiRequest(nreq), MPI.MultiRequest(nreq))

    return cache
end

function assemble_mpi_matvec_sparse!(y_local, y, cache::MatvecCacheSparse)
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    npoin  = size(y_local, 1)

    # Find nonzero indices and values — no threshold, use exact zero only
    iter = 1
    @inbounds for ip = 1:npoin
        val = y_local[ip]
        if val != zero(Float64)
            cache.sp_ip[iter]  = ip
            cache.sp_val[iter] = val
            iter += 1
        end
    end
    nnzero = iter - 1

    # ── Round 1: exchange nnzero counts ──────────────────────────────────────
    req_idx = 1
    @inbounds for i in 0:rank_sz-1
        i == rank && continue
        MPI.Isend(nnzero, comm, cache.request[req_idx]; dest=i, tag=0)
        req_idx += 1
        MPI.Irecv!(@view(cache.nnzeros[i+1:i+1]), comm, cache.request[req_idx]; source=i, tag=0)
        req_idx += 1
    end
    cache.nnzeros[rank+1] = nnzero
    MPI.Waitall(cache.request)

    # ── Resize receive buffers BEFORE posting receives ────────────────────────
    @inbounds for i in 0:rank_sz-1
        i == rank && continue
        n = cache.nnzeros[i+1]
        if length(cache.proc_sp_ip[i+1]) != n
            resize!(cache.proc_sp_ip[i+1], n)
        end
        if length(cache.proc_vec[i+1]) != n
            resize!(cache.proc_vec[i+1], n)
        end
    end

    # ── Round 2: exchange indices ─────────────────────────────────────────────
    req_idx = 1
    send_ip_buf = @view(cache.sp_ip[1:nnzero])
    @inbounds for i in 0:rank_sz-1
        i == rank && continue
        if nnzero > 0
            MPI.Isend(send_ip_buf, comm, cache.request_indices[req_idx]; dest=i, tag=1)
        else
            MPI.Isend(Int[], comm, cache.request_indices[req_idx]; dest=i, tag=1)
        end
        req_idx += 1
        n = cache.nnzeros[i+1]
        if n > 0
            MPI.Irecv!(cache.proc_sp_ip[i+1], comm, cache.request_indices[req_idx]; source=i, tag=1)
        else
            MPI.Irecv!(Int[], comm, cache.request_indices[req_idx]; source=i, tag=1)
        end
        req_idx += 1
    end
    cache.proc_sp_ip[rank+1][1:nnzero] .= send_ip_buf

    # ── Round 3: exchange values ──────────────────────────────────────────────
    req_idx = 1
    send_val_buf = @view(cache.sp_val[1:nnzero])
    @inbounds for i in 0:rank_sz-1
        i == rank && continue
        if nnzero > 0
            MPI.Isend(send_val_buf, comm, cache.request[req_idx]; dest=i, tag=2)
        else
            MPI.Isend(Float64[], comm, cache.request[req_idx]; dest=i, tag=2)
        end
        req_idx += 1
        n = cache.nnzeros[i+1]
        if n > 0
            MPI.Irecv!(cache.proc_vec[i+1], comm, cache.request[req_idx]; source=i, tag=2)
        else
            MPI.Irecv!(Float64[], comm, cache.request[req_idx]; source=i, tag=2)
        end
        req_idx += 1
    end
    cache.proc_vec[rank+1][1:nnzero] .= send_val_buf
    MPI.Waitall(cache.request_indices)
    MPI.Waitall(cache.request)

    # ── Accumulate into global vector ─────────────────────────────────────────
    fill!(y, zero(Float64))
    @inbounds for i in 0:rank_sz-1
        n = cache.nnzeros[i+1]
        for j = 1:n
            ip  = cache.proc_sp_ip[i+1][j]
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
                                        npoin_g::Int,
                                        g_ip2gip::AbstractVector{Int},
                                        g_gip2ip,
                                        comm::MPI.Comm = MPI.COMM_WORLD) where T

    # Forward operator: y = A*x
    nprocs = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    tolerance = 1e-14
    y_local = zeros(Float64,gnpoin)
    y_temp = zeros(Float64,npoin_g)
    A_rows = A_local.rowval
    A_ranges = zeros(Int64, npoin_g,2)
    x_local = zeros(Float64,npoin_g)
    local_gips = zeros(Int64,npoin)
    local_values = zeros(Float64,npoin)

    ip2gip_g = zeros(Int,npoin_g)
    ip2gip_g[1:npoin] .= ip2gip[:]
    #Set up ghost_parent ip2gip if necessary
    if (npoin_g > npoin)
        for ip=npoin+1:npoin_g
            i = ip - npoin
            ip2gip_g[ip] = g_ip2gip[i]
        end
    end
    for ip=1:npoin_g
        A_ranges[ip,1], A_ranges[ip,2] = first(nzrange(A_local, ip)), last(nzrange(A_local, ip))
    end
    
    pM = setup_MatvecSparse_assembler(ip2gip_g,gip2owner)
    
    function matvec!(y, x)
        @inbounds for ip = 1:npoin_g
            x_local[ip] = x[ip2gip_g[ip]]
        end
        mul!(y_temp, A_local, x_local)

        if nprocs <= 4
            fill!(y_local, zero(T))
            @inbounds for ip = 1:npoin_g
                y_local[ip2gip_g[ip]] += y_temp[ip]
            end
            MPI.Allreduce!(y_local, +, comm)
            copyto!(y, y_local)
        else
            assemble_mpi_matvec_sparse!(y_temp, y_local, pM)
            copyto!(y, y_local)
        end
        return y
    end

    function rmatvec!(y, x)
        @inbounds for ip = 1:npoin_g
            x_local[ip] = x[ip2gip_g[ip]]
        end
        mul!(y_temp, A_local', x_local)

        if nprocs <= 4
            fill!(y_local, zero(T))
            @inbounds for ip = 1:npoin_g
                y_local[ip2gip_g[ip]] += y_temp[ip]
            end
            MPI.Allreduce!(y_local, +, comm)
            copyto!(y, y_local)
        else
            assemble_mpi_matvec_sparse!(y_temp, y_local, pM)
            copyto!(y, y_local)
        end
        return y
    end

    # Create LinearOperator
    return LinearOperator{T}(gnpoin, gnpoin, false, false, matvec!, rmatvec!, rmatvec!)
end

function solve_parallel_lsqr(ip2gip, gip2owner, A_local, b, gnpoin, npoin, pM; npoin_g=0, g_ip2gip=Int[], g_gip2ip=Int[], tol::Float64 = 1e-10)
   
     # Create parallel linear operator (better than AbstractMatrix)
    A_parallel = create_parallel_linear_operator(A_local, ip2gip, gip2owner, npoin, gnpoin, npoin_g, g_ip2gip, g_gip2ip)

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
                   itmax = maxiter,
                   verbose = (rank == 0) ? 1 : 0)  # Only print on rank 0

    #return to local indexing
    x_local = zeros(Float64,npoin_g)
    for ip =1:npoin
        gip = ip2gip[ip]
        x_local[ip] = x[gip]
    end

    @info stats
    return x_local
end


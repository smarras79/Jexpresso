
mutable struct CyclingReverseDict
    mapping::Dict{Int, Vector{Int}}
    counters::Dict{Int, Int}
    repeated_keys::Vector{Int}
end

function CyclingReverseDict(a::Vector)
    mapping = Dict{Int, Vector{Int}}()
    repeated_keys = Int[]
    for (i, val) in enumerate(a)
        if haskey(mapping, val)
            if length(mapping[val]) == 1
                push!(repeated_keys, val)
            end
            push!(mapping[val], i)
        else
            mapping[val] = [i]
        end
    end
    CyclingReverseDict(mapping, Dict{Int, Int}(), repeated_keys)
end

function get_repeated_keys(crd::CyclingReverseDict)
    return crd.repeated_keys
end

function get_vals(crd::CyclingReverseDict, key::Int; all = false, first = false, order = true)
    indices = crd.mapping[key]
    if all == true
        return indices
    end
    if first == true
        return indices[1]
    end
    if order == true
        # Get current counter for this key (default to 0)
        counter = get(crd.counters, key, 0)
        # Get the index to return
        result = indices[counter % length(indices) + 1]
        # Increment counter
        crd.counters[key] = counter + 1
        return result
    end
end

# for non-periodic only
mutable struct AssemblerCache
    # Index communication buffers
    recv_idx_buffers::Vector{Vector{Int}}
    # combined_recv_idx::Vector{Int}

    # Send-back buffers
    recvback_idx_buffers::Vector{Vector{Int}}

    # auxiliary
    send_i::Vector{Vector{Int}} 
    send_data_buffers::Vector{Vector{Float64}}
    recv_data_buffers::Vector{Vector{Float64}}
    send_data_sizes::Vector{Int}
    recv_data_sizes::Vector{Int}
    # i_local::Dict{Int, Vector{Int}}

    # Preallocated requests
    requests::MPI.MultiRequest
    requests_back::MPI.MultiRequest
end

function setup_assembler(SD, a, index_a, owner_a)

    if SD == NSD_1D() return nothing end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = maximum(index_a)

    m = size(a, 2)


    send_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    send_i = [Int[] for i in 0:rank_sz-1]
    # send list remote
    for (i, idx) in enumerate(index_a)
        owner = owner_a[i]
        if owner != rank
            buf_idx = get!(send_idx, owner, Int[])
            push!(buf_idx, idx)
            push!(send_i[owner+1], i)
        end
    end
    # send list local (for periodic)
    a_g2l_idx      = CyclingReverseDict(index_a)
    a_g2l_repeated = get_repeated_keys(a_g2l_idx)
    for idx in a_g2l_repeated
        local_idx = get_vals(a_g2l_idx, idx; all = true)[2:end]
        for i in local_idx
            owner = owner_a[i]
            if owner == rank
                buf_idx = get!(send_idx, owner, Int[])
                push!(buf_idx, idx)
                push!(send_i[owner+1], i)
            end
        end
    end


    send_idx_sizes = [length(send_idx[i]) for i in 0:rank_sz-1]
    recv_idx_sizes = MPI.Alltoall(MPI.UBuffer(send_idx_sizes, 1), comm)
    MPI.Barrier(comm)

    # Prepare buffers for sending and receiving data
    send_idx_buffers = [send_idx[i] for i in 0:rank_sz-1]
    recv_idx_buffers = [Vector{Int}(undef, recv_idx_sizes[i+1]) for i in 0:rank_sz-1]

    # Communicate data
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if send_idx_sizes[i+1] > 0
            push!(requests, MPI.Isend(send_idx_buffers[i+1], i, 1, comm))
        end
        if recv_idx_sizes[i+1] > 0
            push!(requests, MPI.Irecv!(recv_idx_buffers[i+1], i, 1, comm))
        end
    end

    # Wait for all communication to complete
    MPI.Waitall(requests)


    # send data back to original ranks
    sendback_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    for i in 0:rank_sz-1
        if recv_idx_sizes[i+1] > 0
            buf_idx = get!(sendback_idx, i, Int[])
            for idx in recv_idx_buffers[i+1]
                push!(buf_idx, idx)
            end
        end
    end
    sendback_idx_sizes = recv_idx_sizes
    recvback_idx_sizes = send_idx_sizes
    MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    sendback_idx_buffers = [sendback_idx[i] for i in 0:rank_sz-1]
    recvback_idx_buffers = [Vector{Int}(undef, length(send_idx[i])) for i in 0:rank_sz-1]


    # Communicate back data
    requests_back = MPI.Request[]
    for i in 0:rank_sz-1
        if sendback_idx_sizes[i+1] > 0
            push!(requests_back, MPI.Isend(sendback_idx_buffers[i+1], i, 3, comm))
        end
        if recvback_idx_sizes[i+1] > 0
            push!(requests_back, MPI.Irecv!(recvback_idx_buffers[i+1], i, 3, comm))
        end
    end


    # Wait for all communication to complete
    MPI.Waitall(requests_back)



    # needed_indices = Set{Int}()
    
    # # 2. Add all indices from recv_idx_buffers
    # for i in 0:rank_sz-1
    #     if recv_idx_sizes[i+1] > 0
    #         for idx in recv_idx_buffers[i+1]
    #             push!(needed_indices, idx)
    #         end
    #     end
    # end
    
    # # 3. Add all indices from recvback_idx_buffers
    # for i in 0:rank_sz-1
    #     if recvback_idx_sizes[i+1] > 0
    #         for idx in recvback_idx_buffers[i+1]
    #             push!(needed_indices, idx)
    #         end
    #     end
    # end


    # change recv_idx_buffers and recvback_idx_buffers to local_idx
    for rk in 1:rank_sz
        recv_idx_buffers_rk     = recv_idx_buffers[rk]
        recvback_idx_buffers_rk = recvback_idx_buffers[rk]
        for (i,idx) in enumerate(recv_idx_buffers_rk)
            recv_idx_buffers_rk[i] = get_vals(a_g2l_idx, idx; first = true)
        end
        if rk-1 == rank
            for (i,idx) in enumerate(recvback_idx_buffers_rk)
                recvback_idx_buffers_rk[i] = send_i[rk][i]
            end
        else
            for (i,idx) in enumerate(recvback_idx_buffers_rk)
                recvback_idx_buffers_rk[i] = get_vals(a_g2l_idx, idx)
            end
        end
    end

    send_data_sizes = [send_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    recv_data_sizes = [recv_idx_sizes[i+1] * m for i in 0:rank_sz-1]

    send_data_buffers = [zeros(Float64, send_idx_sizes[i+1] * m) for i in 0:rank_sz-1]
    recv_data_buffers = [zeros(Float64, recv_idx_sizes[i+1] * m) for i in 0:rank_sz-1]

    # Preallocate requests for assemble function
    n_req = sum(send_data_sizes .> 0) + sum(recv_data_sizes .> 0)
    n_req_back = sum(recv_data_sizes .> 0) + sum(send_data_sizes .> 0)

    return AssemblerCache( recv_idx_buffers, recvback_idx_buffers,
            send_i,send_data_buffers,recv_data_buffers, send_data_sizes, recv_data_sizes,
            MPI.MultiRequest(n_req),
            MPI.MultiRequest(n_req_back))
end


function assemble_mpi!(a, cache::AssemblerCache)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    T = eltype(a)

    is1D = ndims(a) == 1
    n = size(a, 1)
    m = is1D ? 1 : size(a, 2)

    for i in 0:rank_sz-1
        fill!(cache.send_data_buffers[i+1], zero(T))
    end
    @inbounds for owner = 0:rank_sz-1
        buf_data = cache.send_data_buffers[owner+1]
        send_i_local = cache.send_i[owner+1]

        for (i,idx) in enumerate(send_i_local)
            for j = 1:m
                buf_data[(i-1)*m + j] = a[idx,j]
            end
        end
    end
    # MPI.Barrier(comm)


    for i in 0:rank_sz-1
        fill!(cache.recv_data_buffers[i+1], zero(T))
    end


    # Communicate data
    req_idx = 1
    @inbounds for i in 0:rank_sz-1
        if cache.send_data_sizes[i+1] > 0
            MPI.Isend(cache.send_data_buffers[i+1], comm, cache.requests[req_idx]; dest=i, tag=0)
            req_idx += 1
        end
        if cache.recv_data_sizes[i+1] > 0
            MPI.Irecv!(cache.recv_data_buffers[i+1], comm, cache.requests[req_idx]; source=i, tag=0)
            req_idx += 1
        end
    end

    # Wait for all communication to complete
    MPI.Waitall(cache.requests)

    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buffer = cache.recv_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recv_idx_buffers[rk+1])
                for j = 1:m
                    a[local_idx, j] += buffer[(i-1)*m+j]
                end
            end
        end
    end

    # send data back to original ranks
    sendback_data_buffers = cache.recv_data_buffers
    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buf_data = sendback_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recv_idx_buffers[rk+1])
                for j = 1:m
                    buf_data[(i-1)*m+j] = a[local_idx, j]
                end
            end
        end
    end
    sendback_data_sizes = cache.recv_data_sizes
    recvback_data_sizes = cache.send_data_sizes
    # MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    recvback_data_buffers = cache.send_data_buffers


    # Communicate back data
    req_idx = 1
    @inbounds for i in 0:rank_sz-1
        if sendback_data_sizes[i+1] > 0
            MPI.Isend(sendback_data_buffers[i+1], comm, cache.requests_back[req_idx]; dest=i, tag=2)
            req_idx += 1
        end
        if recvback_data_sizes[i+1] > 0
            MPI.Irecv!(recvback_data_buffers[i+1], comm, cache.requests_back[req_idx]; source=i, tag=2)
            req_idx += 1
        end
    end


    # Wait for all communication to complete
    MPI.Waitall(cache.requests_back)


    @inbounds for rk = 0:rank_sz-1
        if recvback_data_sizes[rk+1] > 0
            buffer = recvback_data_buffers[rk+1]
            for (i, local_idx) in enumerate(cache.recvback_idx_buffers[rk+1])
                for j = 1:m
                    a[local_idx, j] = buffer[(i-1)*m+j]
                end
            end
        end
    end
end


mutable struct SendReceiveCache{T}
    # MPI info
    comm::MPI.Comm
    size::Int
    rank::Int
    
    # Send-side buffers
    send_counts::Vector{Int}
    send_targets::Vector{Int}
    send_buffers::Vector{Vector{T}}
    send_indices::Vector{Int}
    
    # Receive-side buffers
    recv_sizes::Vector{Int}
    recv_buffers::Vector{Vector{T}}
    
    # Combined output buffers
    combined_recv_data::Vector{T}
    original_senders::Vector{Int}
    
    # MPI request handling (pre-allocated to exact size)
    requests::MPI.MultiRequest
    
    function SendReceiveCache{T}(comm::MPI.Comm, data2send::AbstractArray{T}, send_targets::Vector{Int}) where T
        rank_sz = MPI.Comm_size(comm)
        rank    = MPI.Comm_rank(comm)
    
        # Validate inputs
        # length(data2send) == length(send_targets) || 
        #     error("data2send and send_targets must have the same rank_sz")
        if ndims(data2send) == 1
            m = 1
        else
            m = size(data2send, 2)
        end
    
        # Pre-allocate and organize data more efficiently
        send_counts = zeros(Int, rank_sz)
        # Count items per destination in single pass
        @inbounds for target in send_targets
            (0 ≤ target < rank_sz) || error("send_targets contains invalid rank")
                send_counts[target + 1] += 1
        end
    
        # Pre-allocate send buffers with exact sizes
        send_buffers = [Vector{Float64}(undef, send_counts[i]) for i in 1:rank_sz]
        send_indices = ones(Int, rank_sz)  # Track insertion positions
        # data = ones(T, length(send_targets),m)  # Track insertion positions
        
        # Fill send buffers in single pass
        @inbounds for (i, target) in enumerate(send_targets)
                idx = target + 1
                send_buffers[idx][send_indices[idx]] = data2send[i]
                send_indices[idx] += 1
        end
    
        # Communicate sizes
        recv_sizes = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)
        
        # Pre-allocate receive buffers
        recv_buffers = [Vector{T}(undef, recv_sizes[i]) for i in 1:rank_sz]
        
        # Pre-calculate total receive rank_sz to avoid dynamic resizing
        total_recv_size = sum(recv_sizes)
        combined_recv_data = zeros(T, total_recv_size)
        original_senders = Vector{Int}(undef, total_recv_size)
    
        # Communicate data with pre-allocated request vector
        req_idx = 0
        requests = MPI.Request[]
        @inbounds for i in 0:rank_sz-1
            if send_counts[i+1] > 0
                req_idx += 1
                push!(requests, MPI.Isend(send_buffers[i+1], i, 4, comm))
            end
            if recv_sizes[i+1] > 0
                req_idx += 1
                push!(requests, MPI.Irecv!(recv_buffers[i+1], i, 4, comm))
            end
        end
    
        # Wait for communication (only for actual requests)
        MPI.Waitall(requests)
    
        # Efficiently combine received data without intermediate allocations
        write_idx = 1
        @inbounds for i in 0:rank_sz-1
            recv_count = recv_sizes[i+1]
            if recv_count > 0
                # Unpack the flat buffer (row-major layout) into 2D array
                copyto!(combined_recv_data, write_idx, recv_buffers[i+1], 1, recv_count)
                fill!(view(original_senders, write_idx:write_idx+recv_count-1), i)
                write_idx += recv_count
            end
        end
        # @show rank, send_buffers,combined_recv_data
        # fill!(combined_recv_data,0.0)
        # for j = 1:m
        #     combined_recv_data[:,j] = send_and_receive(@view(data[:,j]), send_targets, comm)[1]
        #     @show rank, j, combined_recv_data[:,j]
        # end
        new{T}(comm, rank_sz, rank, send_counts, send_targets,
               send_buffers, send_indices,
               recv_sizes, recv_buffers, combined_recv_data, original_senders,
               MPI.MultiRequest(req_idx))
    end
end

# Convenience constructor with type inference
function SendReceiveCache(comm::MPI.Comm, data2send::AbstractArray{T}, send_targets::Vector{Int}) where T
    return SendReceiveCache{T}(comm, data2send, send_targets)
end


function send_and_receive!(cache::SendReceiveCache{T}, data2send::AbstractArray) where T
    

    # Reset send_indices for reuse
    fill!(cache.send_indices, 1)
    # Fill send buffers in single pass
    @inbounds for i in 1:length(cache.send_targets)
        target = cache.send_targets[i]
        idx = target + 1
        cache.send_buffers[idx][cache.send_indices[idx]] = data2send[i]
        cache.send_indices[idx] += 1
    end
    # Communicate data
    req_idx = 0
    @inbounds for i in 0:cache.size-1
        if cache.send_counts[i+1] > 0
            req_idx += 1
            MPI.Isend(cache.send_buffers[i+1], cache.comm, cache.requests[req_idx]; dest=i, tag=3)
        end
        if cache.recv_sizes[i+1] > 0
            req_idx += 1
            MPI.Irecv!(cache.recv_buffers[i+1], cache.comm, cache.requests[req_idx]; source=i, tag=3)
        end
    end
    # Wait for all communication
    MPI.Waitall(cache.requests)

    # Combine received data (original_senders already pre-filled!)
    write_idx = 1
    @inbounds for i in 0:cache.size-1
        recv_count = cache.recv_sizes[i+1]
        if recv_count > 0
            copyto!(cache.combined_recv_data, write_idx, cache.recv_buffers[i+1], 1, recv_count)
            write_idx += recv_count
        end
    end
    
    return nothing
end

function send_and_receive!(combined_recv_data, data2send, send_targets, comm)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    T = eltype(data2send)

    # Validate inputs
    length(data2send) == length(send_targets) || 
        error("data2send and send_targets must have the same size")

    # Pre-allocate and organize data more efficiently
    send_counts = zeros(Int, size)
    
    # Count items per destination in single pass
    @inbounds for target in send_targets
        (0 ≤ target < size) || error("send_targets contains invalid rank")
        send_counts[target + 1] += 1
    end

    # Pre-allocate send buffers with exact sizes
    send_buffers = [Vector{T}(undef, send_counts[i]) for i in 1:size]
    send_indices = ones(Int, size)  # Track insertion positions
    
    # Fill send buffers in single pass
    @inbounds for (data, target) in zip(data2send, send_targets)
        idx = target + 1
        send_buffers[idx][send_indices[idx]] = data
        send_indices[idx] += 1
    end

    # Communicate sizes
    recv_sizes = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)
    
    # Pre-allocate receive buffers
    recv_buffers = [Vector{T}(undef, recv_sizes[i]) for i in 1:size]
    
    # Pre-calculate total receive size to avoid dynamic resizing
    total_recv_size = sum(recv_sizes)
    fill!(combined_recv_data, zero(T))
    original_senders = Vector{Int}(undef, total_recv_size)

    # Communicate data with pre-allocated request vector
    max_requests = 2 * count(x -> x > 0, send_counts) + 2 * count(x -> x > 0, recv_sizes)
    requests = Vector{MPI.Request}(undef, max_requests)
    req_idx = 0

    @inbounds for i in 0:size-1
        if send_counts[i+1] > 0
            req_idx += 1
            requests[req_idx] = MPI.Isend(send_buffers[i+1], i, 4, comm)
        end
        if recv_sizes[i+1] > 0
            req_idx += 1
            requests[req_idx] = MPI.Irecv!(recv_buffers[i+1], i, 4, comm)
        end
    end

    # Wait for communication (only for actual requests)
    if req_idx > 0
        resize!(requests, req_idx)
        MPI.Waitall!(requests)
    end

    # Efficiently combine received data without intermediate allocations
    write_idx = 1
    @inbounds for i in 0:size-1
        recv_count = recv_sizes[i+1]
        if recv_count > 0
            # Copy data directly into pre-allocated arrays
            copyto!(combined_recv_data, write_idx, recv_buffers[i+1], 1, recv_count)
            fill!(view(original_senders, write_idx:write_idx+recv_count-1), i)
            write_idx += recv_count
        end
    end

    return combined_recv_data, original_senders
end

function send_and_receive(data2send, send_targets, comm)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    T = eltype(data2send)

    # Validate inputs
    length(data2send) == length(send_targets) || 
        error("data2send and send_targets must have the same size")

    # Pre-allocate and organize data more efficiently
    send_counts = zeros(Int, size)
    
    # Count items per destination in single pass
    @inbounds for target in send_targets
        (0 ≤ target < size) || error("send_targets contains invalid rank")
        send_counts[target + 1] += 1
    end

    # Pre-allocate send buffers with exact sizes
    send_buffers = [Vector{T}(undef, send_counts[i]) for i in 1:size]
    send_indices = ones(Int, size)  # Track insertion positions
    
    # Fill send buffers in single pass
    @inbounds for (data, target) in zip(data2send, send_targets)
        idx = target + 1
        send_buffers[idx][send_indices[idx]] = data
        send_indices[idx] += 1
    end

    # Communicate sizes
    recv_sizes = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)
    
    # Pre-allocate receive buffers
    recv_buffers = [Vector{T}(undef, recv_sizes[i]) for i in 1:size]
    
    # Pre-calculate total receive size to avoid dynamic resizing
    total_recv_size = sum(recv_sizes)
    combined_recv_data = Vector{T}(undef, total_recv_size)
    original_senders = Vector{Int}(undef, total_recv_size)

    # Communicate data with pre-allocated request vector
    max_requests = 2 * count(x -> x > 0, send_counts) + 2 * count(x -> x > 0, recv_sizes)
    requests = Vector{MPI.Request}(undef, max_requests)
    req_idx = 0

    @inbounds for i in 0:size-1
        if send_counts[i+1] > 0
            req_idx += 1
            requests[req_idx] = MPI.Isend(send_buffers[i+1], i, 5, comm)
        end
        if recv_sizes[i+1] > 0
            req_idx += 1
            requests[req_idx] = MPI.Irecv!(recv_buffers[i+1], i, 5, comm)
        end
    end

    # Wait for communication (only for actual requests)
    if req_idx > 0
        resize!(requests, req_idx)
        MPI.Waitall!(requests)
    end

    # Efficiently combine received data without intermediate allocations
    write_idx = 1
    @inbounds for i in 0:size-1
        recv_count = recv_sizes[i+1]
        if recv_count > 0
            # Copy data directly into pre-allocated arrays
            copyto!(combined_recv_data, write_idx, recv_buffers[i+1], 1, recv_count)
            fill!(view(original_senders, write_idx:write_idx+recv_count-1), i)
            write_idx += recv_count
        end
    end

    return combined_recv_data, original_senders
end


function get_ghost_ips(gelm_ghost, gfacets_ghost, gfacets_owner, conn, pelm2elm, ip2gip, ngl, ghost_p_or_c, comm, SD::NSD_2D)
    rank = MPI.Comm_rank(comm)
    gelm_recv, original_senders = send_and_receive(gelm_ghost, gfacets_owner, comm)
    gfacets_recv = send_and_receive(gfacets_ghost, gfacets_owner, comm)[1]
    lcells = [pelm2elm[x] for x in gelm_recv]
    lfacets = gfacets_recv
    nlfacets    = size(lfacets,1)
    ips_send    = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl)
    ips_targets = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl)
    IP          = KernelAbstractions.zeros(CPU(), TInt, ngl)
    cnt = 1
    for (i, (lfacet, lcell)) in enumerate(zip(lfacets, lcells))

        if (lfacet == 1)
            m = ngl
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, 1, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, ngl, n]
            end
        elseif (lfacet == 2)
            m = 1
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, ngl, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, 1, n]
            end
        elseif (lfacet == 3)
            m = 1:ngl
            n = 1
            if ghost_p_or_c == 1
                IP .= conn[lcell, m, ngl]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, m, 1]
            end
        elseif (lfacet == 4)
            m = 1:ngl
            n = ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, m, 1]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, m, ngl]
            end
        end
        for ip in IP
            ips_send[cnt] = ip2gip[ip]
            ips_targets[cnt] = original_senders[i]
            cnt += 1
        end
    end
    ips_recv, ips_owner  = send_and_receive(ips_send, ips_targets, comm)

    return ips_recv, ips_owner
end


function get_ghost_ips(gelm_ghost, gfacets_ghost, gfacets_owner, conn, pelm2elm, ip2gip, ngl, ghost_p_or_c, comm, SD::NSD_3D)
    rank = MPI.Comm_rank(comm)
    gelm_recv, original_senders = send_and_receive(gelm_ghost, gfacets_owner, comm)
    gfacets_recv = send_and_receive(gfacets_ghost, gfacets_owner, comm)[1]
    lcells = [pelm2elm[x] for x in gelm_recv]
    lfacets = gfacets_recv
    nlfacets    = size(lfacets,1)
    ips_send    = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl * ngl)
    ips_targets = KernelAbstractions.zeros(CPU(), TInt, nlfacets * ngl * ngl)
    IP          = KernelAbstractions.zeros(CPU(), TInt, ngl, ngl)
    cnt = 1
    for (i, (lfacet, lcell)) in enumerate(zip(lfacets, lcells))

        if (lfacet == 1) #front
            l = 1:ngl
            m = 1
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, l, 1, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, l, ngl, n]
            end
        elseif (lfacet == 2) #back
            l = 1:ngl
            m = ngl
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, l, ngl, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, l, 1, n]
            end
        elseif (lfacet == 3) #bottom
            l = 1:ngl
            m = 1:ngl
            n = 1
            if ghost_p_or_c == 1
                IP .= conn[lcell, l, m, 1]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, l, m, ngl]
            end
        elseif (lfacet == 4) #top
            l = 1:ngl
            m = 1:ngl
            n = ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, l, m, ngl]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, l, m, 1]
            end
        elseif (lfacet == 5) #right
            l = ngl
            m = 1:ngl
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, ngl, m, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, 1, m, n]
            end
        elseif (lfacet == 6) #left
            l = 1
            m = 1:ngl
            n = 1:ngl
            if ghost_p_or_c == 1
                IP .= conn[lcell, 1, m, n]
            elseif ghost_p_or_c == 2
                IP .= conn[lcell, ngl, m, n]
            end
        end
        for ip in IP
            # @info rank, ip, cnt,i, lfacet, lcell
            ips_send[cnt] = ip2gip[ip]
            ips_targets[cnt] = original_senders[i]
            cnt += 1
        end
    end
    ips_recv, ips_owner  = send_and_receive(ips_send, ips_targets, comm)

    return ips_recv, ips_owner
end
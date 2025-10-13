    
using Random
using MPI
using BenchmarkTools
using SparseArrays
if !MPI.Initialized()
    MPI.Init()
end

include("Geom.jl")

mutable struct AssemblerCache_v2
    global_max_index::Int
    index_a::Vector{Int}
    owner_a::Vector{Int}

    # Index communication buffers
    recv_idx_buffers::Vector{Vector{Int}}
    # combined_recv_idx::Vector{Int}

    # Send-back buffers
    recvback_idx_buffers::Vector{Vector{Int}}
    # combined_recv_back_idx::Vector{Int}

    # Compact global index mapping
    global_to_compact::Dict{Int, Int}     # Maps global index -> compact index
    # local_to_compact::Dict{Int, Int}      # Maps local index -> compact index (precomputed)
    global_is::Vector{Int}
    gcompact_is::Vector{Int}
    local_is::Vector{Int}
    lcompact_is::Vector{Int}

    sum_array::Vector{Float64}

    # form sparse vector for repeated global idx
    # I::Vector{Int}                       # global index that owned by current rank
    # J::Vector{Int}                       # 
    # V::Vector{Float64}
    # n_row::Int
    # n_col::Int

    # auxiliary
    send_i::Vector{Vector{Int}} 
    send_data_buffers::Vector{Vector{Float64}}
    recv_data_buffers::Vector{Vector{Float64}}
    send_data_sizes::Vector{Int}
    recv_data_sizes::Vector{Int}
    # i_local::Dict{Int, Vector{Int}}

end

function setup_assembler_v2(SD, a, index_a, owner_a)

    # if SD == NSD_1D() return nothing end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)

    global_max_index = maximum(index_a)

    m = size(a, 2)

    # i_local = Dict{Int,  Vector{Int}}()
    # for (i, idx) in enumerate(index_a)
    #     owner = owner_a[i]
    #     if owner == rank
    #         idx_i = get!(i_local, idx, Int[])
    #         push!(idx_i, i)
    #     end
    # end
    # filter!(p -> !isempty(p.second), i_local)

    send_idx = Dict(i => Int[] for i in 0:rank_sz-1)
    send_i = [Int[] for i in 0:rank_sz-1]
    for (i, idx) in enumerate(index_a)
        owner = owner_a[i]
        if owner != rank
            buf_idx = get!(send_idx, owner, Int[])
            push!(buf_idx, idx)
            push!(send_i[owner+1], i)
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
    MPI.Waitall!(requests)

    # Combine received data into a single vector
    # combined_recv_idx = Int[]
    # for i in 0:rank_sz-1
    #     if recv_idx_sizes[i+1] > 0
    #         append!(combined_recv_idx, recv_idx_buffers[i+1])
    #     end
    # end

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
    MPI.Waitall!(requests_back)



    # Combine received data into a single vector
    # combined_recv_back_idx = Int[]
    # for i in 0:rank_sz-1
    #     if recvback_idx_sizes[i+1] > 0
    #         append!(combined_recv_back_idx, recvback_idx_buffers[i+1])
    #     end
    # end

    needed_indices = Set{Int}()
    
    # 1. Add repeated indices from local index_a (owned by this rank)
    local_index_counts = Dict{Int, Int}()
    for (i, idx) in enumerate(index_a)
        if owner_a[i] == rank
            local_index_counts[idx] = get(local_index_counts, idx, 0) + 1
        end
    end
    for (idx, count) in local_index_counts
        if count > 1
            push!(needed_indices, idx)
        end
    end
    
    # 2. Add all indices from recv_idx_buffers
    for i in 0:rank_sz-1
        if recv_idx_sizes[i+1] > 0
            for idx in recv_idx_buffers[i+1]
                push!(needed_indices, idx)
            end
        end
    end
    
    # 3. Add all indices from recvback_idx_buffers
    for i in 0:rank_sz-1
        if recvback_idx_sizes[i+1] > 0
            for idx in recvback_idx_buffers[i+1]
                push!(needed_indices, idx)
            end
        end
    end
    
    # Create compact mapping
    global_is         = collect(needed_indices)
    compact_size      = length(global_is)
    gcompact_is       = collect(1:compact_size)
    global_to_compact = Dict{Int, Int}()
    for (compact_idx, global_idx) in enumerate(global_is)
        global_to_compact[global_idx] = compact_idx
    end
    
    
    # Allocate compact sum arrays
    sum_array = zeros(Float64, compact_size*m)

    # create local to compact Dict
    local_to_compact = Dict{Int, Int}()

    for (local_idx, global_idx) in enumerate(index_a)
        compact_idx = get(global_to_compact, global_idx, 0)
        if compact_idx > 0
            local_to_compact[local_idx] = compact_idx
        end
    end
    local_is = collect(keys(local_to_compact))
    lcompact_is = collect(values(local_to_compact))

    # change recv_idx_buffers and recvback_idx_buffers to compact_idx
    for rk in 1:rank_sz
        recv_idx_buffers_rk     = recv_idx_buffers[rk]
        recvback_idx_buffers_rk = recvback_idx_buffers[rk]
        for (i,idx) in enumerate(recv_idx_buffers_rk)
            recv_idx_buffers_rk[i] = global_to_compact[idx]
        end
        for (i,idx) in enumerate(recvback_idx_buffers_rk)
            recvback_idx_buffers_rk[i] = global_to_compact[idx]
        end
    end
    # sum_array_1D = zeros(Float64, global_max_index)
    # sum_array = zeros(Float64, global_max_index, m)

    
    send_data_sizes = [send_idx_sizes[i+1] * m for i in 0:rank_sz-1]
    recv_data_sizes = [recv_idx_sizes[i+1] * m for i in 0:rank_sz-1]

    send_data_buffers = [zeros(Float64, send_idx_sizes[i+1] * m) for i in 0:rank_sz-1]
    recv_data_buffers = [zeros(Float64, recv_idx_sizes[i+1] * m) for i in 0:rank_sz-1]



    cache = AssemblerCache_v2(global_max_index, index_a, owner_a,
            recv_idx_buffers, recvback_idx_buffers,
            global_to_compact,
            global_is,
            gcompact_is,
            local_is,
            lcompact_is,
            sum_array,
            send_i,send_data_buffers,recv_data_buffers, send_data_sizes, recv_data_sizes)
    return cache
end


function assemble_mpi_v2!(a, cache::AssemblerCache_v2)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    T = eltype(a)

    is1D = ndims(a) == 1
    n = size(a, 1)
    m = is1D ? 1 : size(a, 2)

    fill!(cache.sum_array, zero(T))

    @inbounds for (i, idx) in zip(cache.local_is,cache.lcompact_is)
        owner = cache.owner_a[i]
        if owner == rank
            for j = 1:m
                cache.sum_array[(idx-1)*m+ j] += a[i, j]
            end
        end
    end


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
    # send_data_sizes = [length(cache.send_data_buffers[i+1]) for i in 0:rank_sz-1]
    # recv_data_sizes = MPI.Alltoall(MPI.UBuffer(send_data_sizes, 1), comm)
    MPI.Barrier(comm)

    # Prepare buffers for sending and receiving data
    # recv_data_buffers = [Vector{T}(undef, recv_data_sizes[i+1]) for i in 0:rank_sz-1]

    for i in 0:rank_sz-1
        fill!(cache.recv_data_buffers[i+1], zero(T))
    end


    # Communicate data
    requests = MPI.Request[]
    for i in 0:rank_sz-1
        if cache.send_data_sizes[i+1] > 0
            push!(requests, MPI.Isend(cache.send_data_buffers[i+1], i, 0, comm))
        end
        if cache.recv_data_sizes[i+1] > 0
            push!(requests, MPI.Irecv!(cache.recv_data_buffers[i+1], i, 0, comm))
        end
    end

    # Wait for all communication to complete
    MPI.Waitall!(requests)

    # Combine received data into a single vector
    # combined_recv_data = vcat(cache.recv_data_buffers...)
    # combined_recv_data = T[]
    # for i in 0:rank_sz-1
    #     if recv_data_sizes[i+1] > 0
    #         append!(combined_recv_data, cache.recv_data_buffers[i+1])
    #     end
    # end
    @inbounds for rk in 0:rank_sz-1
        if cache.recv_data_sizes[rk+1] > 0
            buffer = cache.recv_data_buffers[rk+1]
            for (i, compact_idx) in enumerate(cache.recv_idx_buffers[rk+1])
                for j = 1:m
                    cache.sum_array[(compact_idx-1)*m+j] += buffer[(i-1)*m+j]
                end
            end
        end
    end

    # send data back to original ranks
    sendback_data_buffers = cache.recv_data_buffers
    @inbounds for i in 0:rank_sz-1
        if cache.recv_data_sizes[i+1] > 0
            buf_data = sendback_data_buffers[i+1]
            for (j, compact_idx) in enumerate(cache.recv_idx_buffers[i+1])
                for k = 1:m
                    buf_data[(j-1)*m+k] = cache.sum_array[(compact_idx-1)*m+k]
                end
            end
        end
    end
    sendback_data_sizes = cache.recv_data_sizes
    recvback_data_sizes = cache.send_data_sizes
    MPI.Barrier(comm)



    # Prepare buffers for sending and receiving back data
    recvback_data_buffers = cache.send_data_buffers


    # Communicate back data
    requests_back = MPI.Request[]
    @inbounds for i in 0:rank_sz-1
        if sendback_data_sizes[i+1] > 0
            push!(requests_back, MPI.Isend(sendback_data_buffers[i+1], i, 2, comm))
        end
        if recvback_data_sizes[i+1] > 0
            push!(requests_back, MPI.Irecv!(recvback_data_buffers[i+1], i, 2, comm))
        end
    end


    # Wait for all communication to complete
    MPI.Waitall!(requests_back)


    # Combine received data into a single vector
    # combined_recv_back_data = vcat(recvback_data_buffers...)
#     @time begin 
#     combined_recv_back_data = T[]
#     for i in 0:rank_sz-1
#         if recvback_data_sizes[i+1] > 0
#             append!(combined_recv_back_data, recvback_data_buffers[i+1])
#         end
#     end
# end
    @inbounds for rk = 0:rank_sz-1
        if recvback_data_sizes[rk+1] > 0
            buffer = recvback_data_buffers[rk+1]
            for (i, compact_idx) in enumerate(cache.recvback_idx_buffers[rk+1])
                for j = 1:m
                    cache.sum_array[(compact_idx-1)*m+j] = buffer[(i-1)*m+j]
                end
            end
        end
    end

    @inbounds for (i, idx) in zip(cache.local_is,cache.lcompact_is)
        offset = (idx-1)*m
        @inbounds for j = 1:m
            a[i, j] = cache.sum_array[offset + j]
        end
    end
end

function assemble!(a, index_a, owner_a, cache)
    # Build triplets with repeated indices - sparse() will sum automatically!
    n = length(index_a)
    I = repeat(index_a, 5)  # Row indices (repeated)
    J = repeat(1:5, inner=n)  # Column indices
    V = vec(a')  # Values in column-major order

    # sparse() automatically sums repeated (i,j) entries!
    sums = sparse(I, J, V, maximum(index_a), 5, +)

    # Convert to dense and assign back
    # sums_dense = Array(sums)
    # @inbounds for (i, idx) in zip(cache.local_is,cache.lcompact_is)
    #     # if is1D
    #     #     a[i] = cache.sum_array_1D[idx]
    #     # else
    #         for j = 1:m
    #             a[i, j] = cache.sum_array[idx, j]
    #         end
    #     # end
    # end
end


Random.seed!(42)            # Set seed for reproducibility
    
    
n = 1000000
m = 5
a = rand(Float64,n,5)    
num_unique = Int(0.9 * n)      # 900,000 unique entries
num_repeated = n - num_unique   # 100,000 repeated entries

# Generate unique indices (each appears once)
index_a = collect(1:num_unique)

# Add repeated entries by copying from existing ones
append!(index_a, rand(1:num_unique, num_repeated))

# Shuffle to mix them up
shuffle!(index_a)
owner_a = zeros(Int32, n)
SD = nothing

 p_cache_1 = setup_assembler(SD, a, index_a, owner_a)
 p_cache_2 = setup_assembler_v2(SD, a, index_a, owner_a)
 @info length(keys(p_cache_1.local_to_compact))
# @btime assemble_mpi!(a, p_cache_1)
@btime assemble_mpi_v2!(a, p_cache_2)
# @btime assemble!(a, index_a, owner_a, p_cache_2)
# @btime assemble!(SD, a,  index_a, owner_a)
# MPI.Finalize()
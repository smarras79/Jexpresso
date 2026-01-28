"""
    JexpressoCoupling

Module for coupling Jexpresso with external parallel codes via MPI.
Supports one-way and two-way coupling with efficient non-blocking communication.
"""
module JexpressoCoupling

using MPI

export CouplingContext, CouplingBuffer
export initialize_coupling, finalize_coupling
export split_communicator, split_communicator_coupled
export exchange_field_data!, send_field_async, recv_field_async
export wait_all_exchanges!, synchronize_coupling
export is_coupling_root, get_coupling_info
export register_exchange_buffer!, start_buffered_exchange!, finish_buffered_exchange!
export broadcast_to_local!, scatter_to_local!, gather_from_local!
export synchronize_local, reduce_to_local_root!
export broadcast_from_code!, allreduce_across_codes!
export get_partner_root_rank


#===============================================================================
                            DATA STRUCTURES
===============================================================================#

"""
    CouplingBuffer{T}

Pre-allocated buffer for efficient repeated data exchanges.
"""
mutable struct CouplingBuffer{T<:AbstractFloat}
    send_buf::Vector{T}
    recv_buf::Vector{T}
    send_req::Union{MPI.Request, Nothing}
    recv_req::Union{MPI.Request, Nothing}
    tag::Int
    partner_rank::Int
    is_active::Bool
end

function CouplingBuffer{T}(size::Int, tag::Int, partner_rank::Int) where {T}
    return CouplingBuffer{T}(
        Vector{T}(undef, size),
        Vector{T}(undef, size),
        nothing,
        nothing,
        tag,
        partner_rank,
        false
    )
end

"""
    CouplingContext

Main context object holding all coupling-related MPI communicators and metadata.
"""
struct CouplingContext
    # Original communicator (all ranks from all codes)
    comm_world::MPI.Comm
    
    # Local communicator (ranks within this code only)
    comm_local::MPI.Comm
    
    # Inter-code communicator (for communication between codes)
    # Only valid for coupling roots; nothing for other ranks
    comm_inter::Union{MPI.Comm, Nothing}
    
    # Identification
    code_id::Int              # Which code this rank belongs to (1-indexed)
    code_name::String         # Human-readable name
    n_codes::Int              # Total number of coupled codes
    
    # Rank information
    world_rank::Int           # Rank in comm_world
    world_size::Int           # Size of comm_world
    local_rank::Int           # Rank in comm_local
    local_size::Int           # Size of comm_local
    inter_rank::Int           # Rank in comm_inter (-1 if not a member)
    inter_size::Int           # Size of comm_inter (0 if not a member)
    
    # Coupling root info
    is_root::Bool             # Is this rank the coupling root for its code?
    root_ranks::Vector{Int}   # World ranks of all coupling roots
    
    # Pre-allocated buffers for common exchanges
    buffers::Dict{Symbol, CouplingBuffer}
end

#===============================================================================
                        COMMUNICATOR SPLITTING
===============================================================================#

"""
    split_communicator(comm_initial::MPI.Comm, icolor::Int) -> Union{MPI.Comm, Nothing}

Split an MPI communicator based on color assignment.

# Arguments
- `comm_initial`: The initial MPI communicator to split
- `icolor`: Color value determining group membership
  - `icolor == 0`: Rank excluded from new communicator
  - `icolor > 0`: Rank assigned to subcommunicator with that color

# Returns
- `MPI.Comm` for participating ranks, `nothing` for excluded ranks
"""
function split_communicator(comm_initial::MPI.Comm, icolor::Int)
    rank_initial = MPI.Comm_rank(comm_initial)
    
    if icolor == 0
        jcolor = MPI.API.MPI_UNDEFINED[]
    else
        jcolor = icolor
    end
    
    ikey = rank_initial
    comm_new = MPI.Comm_split(comm_initial, jcolor, ikey)
    
    if comm_new == MPI.COMM_NULL
        return nothing
    end
    
    return comm_new
end

"""
    gather_root_ranks(comm_world::MPI.Comm, comm_local::MPI.Comm) -> Vector{Int}

Gather the world ranks of all coupling roots (local rank 0 of each code).
"""
function gather_root_ranks(comm_world::MPI.Comm, comm_local::MPI.Comm, comm_inter::Union{MPI.Comm, Nothing})
    world_rank = MPI.Comm_rank(comm_world)
    local_rank = MPI.Comm_rank(comm_local)
    
    if isnothing(comm_inter)
        # Not a coupling root - will receive result via broadcast
        return Int[]
    end
    
    # Gather world ranks of all roots
    inter_size = MPI.Comm_size(comm_inter)
    root_ranks = MPI.Allgather(world_rank, comm_inter)
    
    return root_ranks
end

#===============================================================================
                        INITIALIZATION / FINALIZATION
===============================================================================#
"""
    initialize_coupling(comm_world::MPI.Comm, code_id::Int, n_codes::Int; 
                        code_name::String="") -> CouplingContext

Initialize coupling between multiple parallel codes.
"""
function initialize_coupling(comm_world::MPI.Comm, code_id::Int, n_codes::Int;
                             code_name::String="Code_$code_id")
    
    @assert 1 <= code_id <= n_codes "code_id must be between 1 and n_codes"
    
    world_rank = MPI.Comm_rank(comm_world)
    world_size = MPI.Comm_size(comm_world)
    
    # Split by code_id to create local communicator
    comm_local = MPI.Comm_split(comm_world, code_id, world_rank)
    local_rank = MPI.Comm_rank(comm_local)
    local_size = MPI.Comm_size(comm_local)
    
    # Create inter-code communicator for coupling roots only
    inter_color = (local_rank == 0) ? 1 : 0
    comm_inter = split_communicator(comm_world, inter_color)
    
    # Determine inter-communicator rank/size
    if !isnothing(comm_inter)
        inter_rank = MPI.Comm_rank(comm_inter)
        inter_size = MPI.Comm_size(comm_inter)
    else
        inter_rank = -1
        inter_size = 0
    end
    
    is_root = (local_rank == 0)
    
    # Gather root ranks on the inter-communicator
    if is_root && !isnothing(comm_inter)
        root_ranks = MPI.Allgather(world_rank, comm_inter)
        n_roots = length(root_ranks)
    else
        root_ranks = Int[]
        n_roots = 0
    end
    
    # Broadcast n_roots to all local ranks using a buffer (not Ref)
    n_roots_buf = Int[n_roots]
    MPI.Bcast!(n_roots_buf, 0, comm_local)
    n_roots = n_roots_buf[1]
    
    # Resize and broadcast root_ranks
    if !is_root
        root_ranks = Vector{Int}(undef, n_roots)
    end
    MPI.Bcast!(root_ranks, 0, comm_local)
    
    # Verify we have the expected number of codes
    if is_root && inter_size != n_codes
        @warn "Expected $n_codes codes but found $inter_size coupling roots"
    end
    
    ctx = CouplingContext(
        comm_world,
        comm_local,
        comm_inter,
        code_id,
        code_name,
        n_codes,
        world_rank,
        world_size,
        local_rank,
        local_size,
        inter_rank,
        inter_size > 0 ? inter_size : n_codes,
        is_root,
        root_ranks,
        Dict{Symbol, CouplingBuffer}()
    )
    
    # Synchronize before returning
    MPI.Barrier(comm_world)
    
    if is_root
        @info "[$code_name] Coupling initialized: local_size=$local_size, " *
              "world_rank=$world_rank, root_ranks=$root_ranks"
    end
    
    return ctx
end


"""
    finalize_coupling(ctx::CouplingContext)

Clean up coupling resources. Call before MPI.Finalize().
"""
function finalize_coupling(ctx::CouplingContext)
    # Wait for any pending communications
    wait_all_exchanges!(ctx)
    
    # Synchronize all codes before cleanup
    MPI.Barrier(ctx.comm_world)
    
    # Free communicators (MPI.jl handles this via GC, but explicit is cleaner)
    # Note: Don't free comm_world as we didn't create it
    
    if ctx.is_root
        @info "[$(ctx.code_name)] Coupling finalized"
    end
end

#===============================================================================
                            UTILITY FUNCTIONS
===============================================================================#

"""
    is_coupling_root(ctx::CouplingContext) -> Bool

Check if this rank is the coupling root for its code.
"""
is_coupling_root(ctx::CouplingContext) = ctx.is_root

"""
    get_coupling_info(ctx::CouplingContext) -> NamedTuple

Get a summary of coupling configuration.
"""
function get_coupling_info(ctx::CouplingContext)
    return (
        code_id = ctx.code_id,
        code_name = ctx.code_name,
        n_codes = ctx.n_codes,
        local_rank = ctx.local_rank,
        local_size = ctx.local_size,
        world_rank = ctx.world_rank,
        is_root = ctx.is_root,
        root_ranks = ctx.root_ranks
    )
end

"""
    get_partner_root_rank(ctx::CouplingContext, partner_code_id::Int) -> Int

Get the world rank of the coupling root for another code.
"""
function get_partner_root_rank(ctx::CouplingContext, partner_code_id::Int)
    @assert 1 <= partner_code_id <= ctx.n_codes
    
    # root_ranks is ordered by inter_rank, which corresponds to code_id ordering
    # We need to find the root rank for the given code_id
    
    if isempty(ctx.root_ranks)
        error("Root ranks not available on this rank")
    end
    
    # Assuming root_ranks[i] corresponds to code i
    return ctx.root_ranks[partner_code_id]
end

#===============================================================================
                        SYNCHRONOUS DATA EXCHANGE
===============================================================================#

"""
    synchronize_coupling(ctx::CouplingContext)

Global barrier across all coupled codes.
"""
function synchronize_coupling(ctx::CouplingContext)
    MPI.Barrier(ctx.comm_world)
end

"""
    synchronize_local(ctx::CouplingContext)

Barrier within this code's local communicator only.
"""
function synchronize_local(ctx::CouplingContext)
    MPI.Barrier(ctx.comm_local)
end

"""
    exchange_field_data!(ctx::CouplingContext, send_data::AbstractArray{T}, 
                         recv_data::AbstractArray{T}, partner_code_id::Int;
                         tag::Int=0) where T

Synchronous exchange of field data between coupling roots.

This is a blocking operation that sends data to and receives data from 
the specified partner code. Only coupling roots participate in the actual
MPI communication.

# Arguments
- `ctx`: Coupling context
- `send_data`: Data to send (must be same size on both sides for this simple version)
- `recv_data`: Buffer to receive data into
- `partner_code_id`: Code ID of the partner to exchange with
- `tag`: MPI tag for the message

# Notes
- Only coupling roots actually communicate; other ranks return immediately
- For production use with different array sizes, use the scatter/gather variants
"""
function exchange_field_data!(ctx::CouplingContext, 
                              send_data::AbstractArray{T},
                              recv_data::AbstractArray{T}, 
                              partner_code_id::Int;
                              tag::Int=0) where T
    
    if !ctx.is_root
        return  # Only roots participate in inter-code communication
    end
    
    partner_world_rank = get_partner_root_rank(ctx, partner_code_id)
    
    # Use Sendrecv for deadlock-free exchange
    MPI.Sendrecv!(send_data, recv_data, ctx.comm_world;
                  dest=partner_world_rank, sendtag=tag,
                  source=partner_world_rank, recvtag=tag)
end

"""
    send_field_to_partner!(ctx::CouplingContext, data::AbstractArray, 
                           partner_code_id::Int; tag::Int=0)

Send field data to a partner code (blocking).
"""
function send_field_to_partner!(ctx::CouplingContext, data::AbstractArray,
                                partner_code_id::Int; tag::Int=0)
    if !ctx.is_root
        return
    end
    
    partner_world_rank = get_partner_root_rank(ctx, partner_code_id)
    MPI.Send(data, ctx.comm_world; dest=partner_world_rank, tag=tag)
end

"""
    recv_field_from_partner!(ctx::CouplingContext, data::AbstractArray,
                             partner_code_id::Int; tag::Int=0)

Receive field data from a partner code (blocking).
"""
function recv_field_from_partner!(ctx::CouplingContext, data::AbstractArray,
                                  partner_code_id::Int; tag::Int=0)
    if !ctx.is_root
        return
    end
    
    partner_world_rank = get_partner_root_rank(ctx, partner_code_id)
    MPI.Recv!(data, ctx.comm_world; source=partner_world_rank, tag=tag)
end

#===============================================================================
                        ASYNCHRONOUS DATA EXCHANGE
===============================================================================#

"""
    send_field_async(ctx::CouplingContext, data::AbstractArray,
                     partner_code_id::Int; tag::Int=0) -> Union{MPI.Request, Nothing}

Non-blocking send of field data to a partner code.

# Returns
- `MPI.Request` for coupling roots to wait on
- `nothing` for non-root ranks
"""
function send_field_async(ctx::CouplingContext, data::AbstractArray,
                          partner_code_id::Int; tag::Int=0)
    if !ctx.is_root
        return nothing
    end
    
    partner_world_rank = get_partner_root_rank(ctx, partner_code_id)
    return MPI.Isend(data, ctx.comm_world; dest=partner_world_rank, tag=tag)
end

"""
    recv_field_async(ctx::CouplingContext, data::AbstractArray,
                     partner_code_id::Int; tag::Int=0) -> Union{MPI.Request, Nothing}

Non-blocking receive of field data from a partner code.

# Returns
- `MPI.Request` for coupling roots to wait on
- `nothing` for non-root ranks
"""
function recv_field_async(ctx::CouplingContext, data::AbstractArray,
                          partner_code_id::Int; tag::Int=0)
    if !ctx.is_root
        return nothing
    end
    
    partner_world_rank = get_partner_root_rank(ctx, partner_code_id)
    return MPI.Irecv!(data, ctx.comm_world; source=partner_world_rank, tag=tag)
end

"""
    wait_request(req::Union{MPI.Request, Nothing})

Wait for a single async request to complete.
"""
function wait_request(req::Union{MPI.Request, Nothing})
    if !isnothing(req)
        MPI.Wait(req)
    end
end

"""
    wait_all_exchanges!(ctx::CouplingContext)

Wait for all pending buffer exchanges to complete.
"""
function wait_all_exchanges!(ctx::CouplingContext)
    for (name, buf) in ctx.buffers
        if buf.is_active
            wait_request(buf.send_req)
            wait_request(buf.recv_req)
            buf.send_req = nothing
            buf.recv_req = nothing
            buf.is_active = false
        end
    end
end

#===============================================================================
                    BUFFERED EXCHANGES (FOR REPEATED USE)
===============================================================================#

"""
    register_exchange_buffer!(ctx::CouplingContext, name::Symbol, size::Int,
                              partner_code_id::Int; tag::Int=0, T::Type=Float64)

Register a pre-allocated buffer for repeated exchanges.

# Arguments
- `ctx`: Coupling context (modified in place)
- `name`: Symbolic name for this buffer
- `size`: Number of elements
- `partner_code_id`: Partner code for this exchange
- `tag`: MPI tag
- `T`: Element type (default Float64)
"""
function register_exchange_buffer!(ctx::CouplingContext, name::Symbol, size::Int,
                                   partner_code_id::Int; tag::Int=0, T::Type=Float64)
    partner_rank = get_partner_root_rank(ctx, partner_code_id)
    ctx.buffers[name] = CouplingBuffer{T}(size, tag, partner_rank)
    return nothing
end

"""
    start_buffered_exchange!(ctx::CouplingContext, name::Symbol, send_data::AbstractVector)

Start an async exchange using a registered buffer.
"""
function start_buffered_exchange!(ctx::CouplingContext, name::Symbol, 
                                  send_data::AbstractVector)
    if !ctx.is_root
        return
    end
    
    buf = ctx.buffers[name]
    @assert length(send_data) == length(buf.send_buf) "Data size mismatch"
    
    # Copy data to send buffer
    copyto!(buf.send_buf, send_data)
    
    # Start async operations
    buf.send_req = MPI.Isend(buf.send_buf, ctx.comm_world; 
                             dest=buf.partner_rank, tag=buf.tag)
    buf.recv_req = MPI.Irecv!(buf.recv_buf, ctx.comm_world;
                              source=buf.partner_rank, tag=buf.tag)
    buf.is_active = true
end

"""
    finish_buffered_exchange!(ctx::CouplingContext, name::Symbol) -> Vector

Complete a buffered exchange and return received data.
"""
function finish_buffered_exchange!(ctx::CouplingContext, name::Symbol)
    if !ctx.is_root
        return eltype(ctx.buffers[name].recv_buf)[]
    end
    
    buf = ctx.buffers[name]
    
    if buf.is_active
        MPI.Wait(buf.send_req)
        MPI.Wait(buf.recv_req)
        buf.send_req = nothing
        buf.recv_req = nothing
        buf.is_active = false
    end
    
    return copy(buf.recv_buf)
end

#===============================================================================
                    COLLECTIVE OPERATIONS ACROSS CODES
===============================================================================#

"""
    broadcast_from_code!(ctx::CouplingContext, data::AbstractArray, source_code_id::Int)

Broadcast data from one code's root to all other code roots.
"""
function broadcast_from_code!(ctx::CouplingContext, data::AbstractArray, 
                              source_code_id::Int)
    if !ctx.is_root || isnothing(ctx.comm_inter)
        return
    end
    
    # In comm_inter, ranks are ordered by code_id
    source_inter_rank = source_code_id - 1  # Convert to 0-indexed
    MPI.Bcast!(data, source_inter_rank, ctx.comm_inter)
end

"""
    allreduce_across_codes!(ctx::CouplingContext, data::AbstractArray, op::MPI.Op)

Perform allreduce across all code roots.
"""
function allreduce_across_codes!(ctx::CouplingContext, data::AbstractArray, 
                                 op::MPI.Op=MPI.SUM)
    if !ctx.is_root || isnothing(ctx.comm_inter)
        return
    end
    
    MPI.Allreduce!(data, op, ctx.comm_inter)
end

#===============================================================================
                    FIELD DISTRIBUTION WITHIN CODE
===============================================================================#

"""
    scatter_to_local!(ctx::CouplingContext, global_data::AbstractArray,
                      local_data::AbstractArray; root::Int=0)

Scatter data from local root to all local ranks within this code.
"""
function scatter_to_local!(ctx::CouplingContext, global_data::AbstractArray,
                           local_data::AbstractArray; root::Int=0)
    MPI.Scatter!(global_data, local_data, ctx.comm_local; root=root)
end

"""
    gather_from_local!(ctx::CouplingContext, local_data::AbstractArray,
                       global_data::AbstractArray; root::Int=0)

Gather data from all local ranks to the local root.
"""
function gather_from_local!(ctx::CouplingContext, local_data::AbstractArray,
                            global_data::AbstractArray; root::Int=0)
    MPI.Gather!(local_data, global_data, ctx.comm_local; root=root)
end

"""
    broadcast_to_local!(ctx::CouplingContext, data::AbstractArray; root::Int=0)

Broadcast data from local root to all local ranks.
"""
function broadcast_to_local!(ctx::CouplingContext, data::AbstractArray; root::Int=0)
    MPI.Bcast!(data, root, ctx.comm_local)
end

"""
    reduce_to_local_root!(ctx::CouplingContext, send_data::AbstractArray,
                          recv_data::AbstractArray, op::MPI.Op=MPI.SUM; root::Int=0)

Reduce data from all local ranks to the local root.
"""
function reduce_to_local_root!(ctx::CouplingContext, send_data::AbstractArray,
                               recv_data::AbstractArray, op::MPI.Op=MPI.SUM; 
                               root::Int=0)
    MPI.Reduce!(send_data, recv_data, op, ctx.comm_local; root=root)
end

end # module JexpressoCoupling


#------------------------------------------------------------------------------------
# Receive the grid coordinate from Alya (structured and regular only)
#------------------------------------------------------------------------------------
function distribute_and_count!(
    rem_nx::AbstractVector{Int32},
    rem_min::AbstractVector{T},
    rem_max::AbstractVector{T},
    ndime::Int,
    nranks2::Int,
    in_my_rank::Bool,
    a::AbstractMatrix{Int32},
    wrank::Int,
    alya2world,
) where {T<:AbstractFloat}
    
    nx, ny, nz = rem_nx
    nxy  = nx * ny
    nmax = nxy * nz

    @assert ndime in (2, 3) "ndime is typically 2 or 3"
    @assert length(rem_min) ≥ ndime && length(rem_nx) ≥ ndime

    r     = mod(nmax, nranks2 - 1)
    npoin = nmax ÷ (nranks2 - 1)

    ri = zeros(Int, ndime)
    x  = zeros(T, ndime)

    rem_dx[1:ndime] = (rem_max[1:ndime] .- rem_min[1:ndime])./(rem_nx[1:ndime]-1)
    @info rem_dx
    #=
    @inbounds for ipoin in 1:nmax
        i0   = ipoin - 1
        iz   = i0 ÷ nxy
        remz = i0 - iz * nxy
        iy   = remz ÷ nx
        ix   = mod(remz - iy * nx, nx)

        if ndime ≥ 1; ri[1] = ix; end
        if ndime ≥ 2; ri[2] = iy; end
        if ndime ≥ 3; ri[3] = iz; end

        x[1:ndime] = rem_min[1:ndime] .+ T.(ri[1:ndime]) .* rem_dx[1:ndime]

        if in_my_rank
            alya_rank = if ipoin ≤ r * (npoin + 1)
                (ipoin - 1) ÷ (npoin + 1) + 1
            else
                r + (ipoin - r * (npoin + 1) - 1) ÷ npoin + 1
            end

            world_rank = alya2world[alya_rank]
            a[wrank, world_rank] += 1
        end
    end
    =#
    return a
end

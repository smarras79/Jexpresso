using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads


"""
    perform_coupling_handshake(world, nparts)

Perform the initial handshake with Alya to establish that both codes are ready.
This is a minimal synchronization step that must complete before data exchange.

# Arguments
- `world`: MPI.COMM_WORLD communicator
- `nparts`: Number of local Jexpresso ranks

# Returns
- `true` if coupled mode is active (wsize > nparts), `false` otherwise
"""
function je_perform_coupling_handshake(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    
    # Check if we are in coupled mode
    if wsize <= nparts
        return false  # Standalone mode, no coupling
    end
    
    # Exchange app identity with Alya
    local_chars = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
    MPI.Gather!(local_chars, nothing, 0, world)
    
    if wrank == nparts 
        println("[Driver] Handshake complete - Jexpresso ready.")
        flush(stdout)
    end
    
    return true  # Coupled mode active
end


"""
    receive_alya_data(world, nparts)

Receive grid metadata and coupling information from Alya.
Must be called AFTER perform_coupling_handshake() succeeds.

# Arguments
- `world`: MPI.COMM_WORLD communicator
- `nparts`: Number of local Jexpresso ranks

# Returns
- Nothing (data is stored via set_coupling_data)
"""
function je_receive_alya_data(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    
    # Verify we are in coupled mode
    if wsize <= nparts
        @warn "receive_alya_data called but not in coupled mode"
        return
    end
    
    # 1. Receive grid dimensions and bounds
    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    ndime = ndime_buf[1]
    
    rem_min = Vector{Float32}(undef, 3)
    rem_max = Vector{Float32}(undef, 3)
    rem_nx  = Vector{Int32}(undef, 3)
    
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end
    
    # 2. Synchronize rank maps across the whole world
    nranks_other = wsize - nparts
    alya2world_l = zeros(Int32, nranks_other)
    alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)
    
    a_l = zeros(Int32, wsize, wsize)
    a   = MPI.Allreduce(a_l, MPI.SUM, world)
    
    # 3. Store metadata for coupling routines
    set_coupling_data(Dict{Symbol,Any}(
        :ndime         => ndime,
        :rem_min       => rem_min,
        :rem_max       => rem_max,
        :rem_nx        => rem_nx,
        :alya2world    => alya2world,
        :couple_matrix => a,
    ))
    
    if wrank == nparts 
        println("[Driver] Data received from Alya:")
        println("         ndime=$ndime, min=$rem_min, max=$rem_max, nx=$rem_nx")
        flush(stdout)
    end
    
    return nothing
end

"""
        St_coupling

        Holds all inter-code coupling state needed during the time loop:
          - world communicator for inter-code MPI
          - local communicator for Jexpresso-internal MPI
          - rank/size info on both communicators
          - pre-allocated send/receive buffers for field exchange with Alya
          - coupling metadata received during initialization (ndime, remote grid info, etc.)
        """
mutable struct St_coupling
    comm_world::MPI.Comm
    comm_local::MPI.Comm
    wrank::Int
    wsize::Int
    lrank::Int
    lsize::Int
    nranks_alya::Int

    # Remote grid metadata from Alya
    ndime::Int32
    rem_min::Vector{Float32}
    rem_max::Vector{Float32}
    rem_nx::Vector{Int32}
    alya2world::Vector{Int32}

    # Pre-allocated exchange buffers (sized after mesh is known)
    send_buf::Vector{Float64}
    recv_buf::Vector{Float64}

    # Coupling allocation matrix
    couple::Matrix{Int64}
    bbox_locator::Any  # Store the bounding box locator
end

"""
        je_couplingSetup(inputs) -> Union{St_coupling, Nothing}

        Initialize coupling data structures using the pre-stored coupling data
        from `Jexpresso-mini-coupled.jl`.

        Returns `St_coupling` if coupling is active, `nothing` otherwise.
        """
function je_splittingSetup()
    
    world = MPI.COMM_WORLD
    wrank = MPI.Comm_rank(world)
    wsize = MPI.Comm_size(world)
    
    # Still split here so Julia ranks can talk to each other safely
    appid      = try parse(Int, get(ENV, "APPID", "2")) catch; 2 end
    local_comm = MPI.Comm_split(world, appid, wrank)
    
    set_mpi_comm(local_comm)
    set_mpi_comm_world(world) # Store this for the driver!

    return (world, wrank, wsize, appid, local_comm)
end

#------------------------------------------------------------------------------------
# Coupling exchange: send data to Alya and receive data from Alya
# Called within the time loop at every coupling step.
#
# This uses MPI_Sendrecv on COMM_WORLD between the Jexpresso local root
# and the Alya root (world rank 0).
# Non-root Jexpresso ranks participate via Bcast on the local communicator.
#------------------------------------------------------------------------------------
"""
        coupling_send_recv!(coupling::St_coupling, send_data::AbstractVector{Float64},
                                recv_data::AbstractVector{Float64};
                                alya_root::Int=0, tag::Int=100)

        Exchange field data with Alya:
          1. Jexpresso local root (lrank==0) does MPI.Sendrecv! with the Alya root on comm_world
          2. The received data is then broadcast to all Jexpresso local ranks via comm_local

        `send_data` and `recv_data` must be the same length on the root.
        On non-root ranks, `recv_data` is filled by the broadcast.
        """
function coupling_send_recv!(cpg, send_buf, recv_buf; alya_root=0, tag=1001)
    # This should ONLY be called by the leader rank (cpg.lrank == 0)
    # Use MPI.Sendrecv to match Fortran's MPI_Sendrecv
    MPI.Sendrecv!(send_buf, recv_buf, cpg.comm_world;
                  dest=alya_root, sendtag=tag,
                  source=alya_root, recvtag=tag)
end

"""
            coupling_send!(coupling::St_coupling, send_data::AbstractVector{Float64};
                           alya_root::Int=0, tag::Int=200)

        Send field data to Alya (non-blocking from Jexpresso root).
        Only the Jexpresso local root sends.
        """
function coupling_send!(coupling::St_coupling,
                        send_data::AbstractVector{Float64};
                        alya_root::Int=0, tag::Int=200)

    if coupling.lrank == 0
        MPI.Send(send_data, coupling.comm_world; dest=alya_root, tag=tag)
    end
end

"""
            coupling_recv!(coupling::St_coupling, recv_data::AbstractVector{Float64};
                           alya_root::Int=0, tag::Int=200)

        Receive field data from Alya (blocking on Jexpresso root),
        then broadcast to all local ranks.
        """
function coupling_recv!(coupling::St_coupling,
                        recv_data::AbstractVector{Float64};
                        alya_root::Int=0, tag::Int=200)

    if coupling.lrank == 0
        MPI.Recv!(recv_data, coupling.comm_world; source=alya_root, tag=tag)
    end

    MPI.Bcast!(recv_data, 0, coupling.comm_local)
end

#------------------------------------------------------------------------------------
# Receive the grid coordinate from Alya (structured and regular only)
#------------------------------------------------------------------------------------
#=function distribute_and_count!(
rem_nx,
rem_min,
rem_max,
ndime,
nranks2,
in_my_rank,
npoin_send,
npoin_recv,
wrank,
alya2world,
mesh)=#

function distribute_and_count!(coupling, mesh)
@info " DISTRIBUTE AND COUNT "
    ndime   = coupling.ndime
    nranks2 = coupling.nranks_alya

    rem_nx  = coupling.rem_nx
    rem_min = coupling.rem_min
    rem_max = coupling.rem_max
    @info rem_nx, rem_min, rem_max

    @mystop
    comm_local = get_mpi_comm()
    lrank = MPI.Comm_rank(comm_local)
    lsize = MPI.Comm_size(comm_local)
    
    # Build the bounding box locator
    locator = build_rank_bbox_locator(mesh, comm_local)
    
    #
    # is_in_my_rank must be
    # found from Jexpresso
    #
    nx, ny, nz = rem_nx
    nxy        = nx * ny
    nmax       = nxy * nz
    rem_dx     = similar(rem_max)

    @assert coupling.ndime in (2, 3) "ndime is typically 2 or 3"
    @assert length(rem_min) ≥ coupling.ndime && length(rem_nx) ≥ coupling.ndime
    
    r     = mod(nmax, nranks2 - 1)
    npoin = nmax ÷ (nranks2 - 1)

    ri = zeros(Int32,   ndime)
    x  = zeros(Float64, ndime)

    rem_dx[1:ndime] = (rem_max[1:ndime] .- rem_min[1:ndime])./(rem_nx[1:ndime] .- 1)

    #
    # Local max min on this rank:
    #
    lxmin = min(mesh.x); lxmax = max(mesh.x)
    lymin = min(mesh.y); lymax = max(mesh.y)
    lzmin = min(mesh.z); lzmax = max(mesh.z)

    my_boundingBox = St_myBoundingBox(lxmin, lxmax, lymin, lymax, lzmin, lzmax)

    @inbounds for ipoin in 1:nmax
        i0   = ipoin - 1
        iz   = i0 ÷ nxy
        remz = i0 - iz * nxy
        iy   = remz ÷ nx
        ix   = mod(remz - iy * nx, nx)

        if ndime ≥ 1; ri[1] = ix; end
        if ndime ≥ 2; ri[2] = iy; end
        if ndime ≥ 3; ri[3] = iz; end
        
        x[1:ndime] = rem_min[1:ndime] .+ Float32.(ri[1:ndime]) .* rem_dx[1:ndime]
        X = [x[1], x[2]]
        
        found, owner_ranks, is_mine = find_point_owners_bbox(locator, comm_local, X)

        for i in 1:size(X, 1)
            if found[i]
                if is_mine[i]
                    @info "Rank $locator.rank OWNS point ($( X[i,1]), $(X[i,2]))"
                    # Do interpolation here
                end
                
                # Show all owners (useful for debugging overlaps)
                if locator.myrank == 0
                    owners_str = join(owner_ranks[i], ", ")
                    @info "Point ($(X[i,1]), $(X[i,2])) owned by rank(s): $owners_str"
                end
            else
                if locator.myrank == 0
                    @warn "Point ($(X[i,1]), $(X[i,2])) is OUTSIDE all rank domains"
                end
            end
        end

        xin, yin, zin = is_on_my_rank(local::St_myBoundingBox, global::St_myBoundingBox, x, y, z; atol=default_atol(local))
             
        if is_in_my_rank
        alya_rank = if ipoin ≤ r * (npoin + 1)
        (ipoin - 1) ÷ (npoin + 1) + 1
        else
        r + (ipoin - r * (npoin + 1) - 1) ÷ npoin + 1
        end

        wrank = alya2world[alya_rank]
        npoin_send[wrank] += 1
        
        ## after all2all I need to tell Alya the local numbering of alya's points:
        #
        # send the node list 
        end
        end
        
        ### the do the all2all of npoin_send
        # ∀  mpi send list  of list
        
        @inbounds for ipoin in 1:nmax
        #ahora que tengo identificado en que rango estan los punto, ahora relleno la tabla
        
        
        end
        
        println("$a")
        return a
  
    end

    
end

struct St_myBoundingBox{T}
    xmin::T; xmax::T
    ymin::T; ymax::T
    zmin::T; zmax::T
end

# Scale-aware tolerance (important for large coordinate magnitudes)
@inline function default_atol(xlocal::St_myBoundingBox{T}) where {T<:Real}
    Lx = xlocal.xmax - xlocal.xmin
    Ly = xlocal.ymax - xlocal.ymin
    Lz = xlocal.zmax - xlocal.zmin
    scale = max(Lx, Ly, Lz, one(T))
    return sqrt(eps(T)) * scale   # ~1e-8 * scale for Float64
end

"""
        Interpolate solution at point (px, py) using Julia's mesh and basis functions.
        Returns interpolated value for the first solution component.
        """
function interpolate_at_point(px::Float64, py::Float64, mesh, u, qp)
    # Find which element contains the point
    for iel in 1:mesh.nelem
        # Get element nodes
        node_indices = mesh.conn[iel, :]
        node_coords = [(mesh.x[i], mesh.y[i]) for i in node_indices]
        
        # Check if point is in this element
        if point_in_quad(px, py, node_coords)
            # Map to reference coordinates (ξ, η) ∈ [-1, 1]²
            ξ, η = physical_to_reference(px, py, node_coords)
            
            # Evaluate basis functions at (ξ, η)
            # Assumes you have Lagrange basis functions
            basis_vals = evaluate_basis(ξ, η, mesh.ngl)
            
            # Interpolate solution
            interp_val = 0.0
            for (i_local, i_global) in enumerate(node_indices)
                interp_val += u[i_global] * basis_vals[i_local]
            end
            
            return interp_val
        end
    end
    
    # Point not found in any element
    return NaN
end

function point_in_quad(px, py, coords::Vector{Tuple{Float64,Float64}})
    # Simple bounding box check
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    return (minimum(xs) <= px <= maximum(xs) && 
        minimum(ys) <= py <= maximum(ys))
end


function map_coords_to_rank(mesh, X)
    comm_local = get_mpi_comm()
    lrank = MPI.Comm_rank(comm_local)
    lsize = MPI.Comm_size(comm_local)
    
    # Build the bounding box locator
    locator = build_rank_bbox_locator(mesh, comm_local)
    
    # Test coordinates
    #= X = [
        -5000.0   0.0;
        -3300.0   0.0;
        0.0       0.0;
        3300.0    5000.0;
        5000.0    10000.0;
    ]=#
    
    # Find owners
    found, owner_ranks, is_mine = find_point_owners_bbox(locator, comm_local, X)
    
    # Report results
    for i in 1:size(X, 1)
        if found[i]
            if is_mine[i]
                @info "Rank $locator.rank OWNS point ($( X[i,1]), $(X[i,2]))"
                # Do interpolation here
            end
            
            # Show all owners (useful for debugging overlaps)
            if locator.myrank == 0
                owners_str = join(owner_ranks[i], ", ")
                @info "Point ($(X[i,1]), $(X[i,2])) owned by rank(s): $owners_str"
            end
        else
            if locator.myrank == 0
                @warn "Point ($(X[i,1]), $(X[i,2])) is OUTSIDE all rank domains"
            end
        end
    end
    
    # Return locator for use in coupling
    return locator
end


"""
            Build a simple bounding box locator for rank ownership.
            Each rank knows its own bounding box and can test if points are inside.
            """
function build_rank_bbox_locator(mesh, comm_local)
    lrank = MPI.Comm_rank(comm_local)
    lsize = MPI.Comm_size(comm_local)
    
    # Compute local bounding box on this rank
    local_xmin = minimum(mesh.x)
    local_xmax = maximum(mesh.x)
    local_ymin = minimum(mesh.y)
    local_ymax = maximum(mesh.y)
    
    # Gather all ranks' bounding boxes
    all_xmin = zeros(Float64, lsize)
    all_xmax = zeros(Float64, lsize)
    all_ymin = zeros(Float64, lsize)
    all_ymax = zeros(Float64, lsize)
    
    # Each rank puts its bounds at its rank position
    all_xmin[lrank + 1] = local_xmin
    all_xmax[lrank + 1] = local_xmax
    all_ymin[lrank + 1] = local_ymin
    all_ymax[lrank + 1] = local_ymax
    
    # Allreduce with SUM (since zeros elsewhere)
    MPI.Allreduce!(all_xmin, MPI.SUM, comm_local)
    MPI.Allreduce!(all_xmax, MPI.SUM, comm_local)
    MPI.Allreduce!(all_ymin, MPI.SUM, comm_local)
    MPI.Allreduce!(all_ymax, MPI.SUM, comm_local)
    
    return (
        xmin = all_xmin,
        xmax = all_xmax,
        ymin = all_ymin,
        ymax = all_ymax,
        nranks = lsize,
        myrank = lrank
    )
end

#------------------------------------------------------------------------------------------
# Check which rank(s) own a batch of points based on bounding boxes.
#    Returns:
#    - found: whether point is in ANY rank's bbox
#    - owner_ranks: list of ranks that own each point (can be multiple if overlapping)
#    - is_mine: whether this rank owns each point
#
#   find_point_owners_bbox(locator, comm_local, X; margin=1e-8)
#
#    Find which rank(s) own a batch of points based on rank bounding boxes.
#
#    # Arguments
#    - `locator`: Output from `build_rank_bbox_locator` containing bbox info for all ranks
#    - `comm_local`: Julia's local MPI communicator
#    - `X::AbstractMatrix`: Nx2 matrix of coordinates [x, y] to query
#    - `margin::Float64`: Tolerance for point-on-boundary cases (default: 1e-8)
#
#    # Returns
#    - `found::Vector{Bool}`: Whether each point is in ANY rank's bbox
#    - `owner_ranks::Vector{Vector{Int32}}`: List of owner rank(s) for each point (can be multiple)
#    - `is_mine::Vector{Bool}`: Whether THIS rank owns each point
#
#    # Notes
#    - If multiple ranks claim a point (overlapping bboxes), all are returned in owner_ranks
#    - Points exactly on boundaries are included (with margin tolerance)
#    - Points outside all domains have empty owner_ranks and found[i] = false
#------------------------------------------------------------------------------------------
function find_point_owners_bbox(locator, comm_local, X::AbstractMatrix{<:Real}; margin::Float64=1e-8)
    
    # Find which rank(s) own a batch of points based on rank bounding boxes.
    #
    # Arguments:
    #   locator - Output from build_rank_bbox_locator containing bbox info for all ranks
    #   comm_local - Julia's local MPI communicator
    #   X - Nx2 matrix of coordinates [x, y] to query
    #   margin - Tolerance for point-on-boundary cases (default: 1e-8)
    #
    # Returns:
    #   found - Whether each point is in ANY rank's bbox
    #   owner_ranks - List of owner rank(s) for each point (can be multiple)
    #   is_mine - Whether THIS rank owns each point
    
    myrank = MPI.Comm_rank(comm_local)
    M = size(X, 1)
    
    # Initialize outputs
    found = falses(M)
    owner_ranks = [Int32[] for _ in 1:M]  # List of owners per point
    is_mine = falses(M)
    
    # Check each point against all ranks' bounding boxes
    for i in 1:M
        px = Float64(X[i, 1])
        py = Float64(X[i, 2])
        
        # Check against all ranks' bounding boxes
        for r in 0:(locator.nranks - 1)
            ridx = r + 1  # Convert to 1-based indexing for Julia arrays
            
            # Get this rank's bounding box
            xmin = locator.xmin[ridx] - margin
            xmax = locator.xmax[ridx] + margin
            ymin = locator.ymin[ridx] - margin
            ymax = locator.ymax[ridx] + margin
            
            # Check if point is inside this rank's bbox
            if (px >= xmin && px <= xmax && py >= ymin && py <= ymax)
                # This rank owns the point
                push!(owner_ranks[i], Int32(r))
                found[i] = true
                
                # Check if it's this rank
                if r == myrank
                    is_mine[i] = true
                end
            end
        end
    end
    
    return found, owner_ranks, is_mine
end

function je_get_alya_data!(coupling, splitting)

     if JEXPRESSO_MPI_COMM[] === nothing
        world = MPI.COMM_WORLD
        wrank = MPI.Comm_rank(world)
        wsize = MPI.Comm_size(world)

        appid = try parse(Int, get(ENV, "APPID", "2")) catch; 2 end

        local_comm = splitting.local_comm# MPI.Comm_split(world, appid, wrank)
        lrank = solitting.lrank #MPI.Comm_rank(local_comm)
        lsize = MPI.Comm_size(local_comm)

        if lsize < wsize
            #-------------------------------------------------------------
            # Coupled mode detected: we share COMM_WORLD with another code
            #-------------------------------------------------------------
            set_mpi_comm(local_comm)
            set_mpi_comm_world(world)

            if lrank == 0
                println("[Jexpresso] Coupled mode: world_size=$wsize, local_size=$lsize, appid=$appid")
                flush(stdout)
            end

            nranks_other = wsize - lsize

            # Exchange app identity (Gather names to world rank 0)
            local_chars = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
            MPI.Gather!(local_chars, nothing, 0, world)

            #-------------------------------------------------------------
            # Receive remote grid metadata from Alya via Bcast on COMM_WORLD
            # Alya broadcasts from world rank 0.
            # Use Int32/Float32 to match Fortran MPI types.
            #-------------------------------------------------------------
            ndime_buf = Vector{Int32}(undef, 1)
            MPI.Bcast!(ndime_buf, 0, world)
            ndime = ndime_buf[1]

            rem_min = Vector{Float32}(undef, 3)
            rem_max = Vector{Float32}(undef, 3)
            rem_nx  = Vector{Int32}(undef, 3)
            for idime in 1:3
                MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
                MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
                MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
            end

            alya2world_l = zeros(Int32, nranks_other)
            alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)

            a_l = zeros(Int32, wsize, wsize)
            a   = MPI.Allreduce(a_l, MPI.SUM, world)

            # Store coupling data for je_couplingSetup() to consume later
            set_coupling_data(Dict{Symbol,Any}(
                :ndime         => ndime,
                :rem_min       => rem_min,
                :rem_max       => rem_max,
                :rem_nx        => rem_nx,
                :alya2world    => alya2world,
                :couple_matrix => a,
            ))

            if lrank == 0
                println("[Jexpresso] Received from Alya: ndime=$ndime, " *
                        "rem_min=$rem_min, rem_max=$rem_max, rem_nx=$rem_nx")
                flush(stdout)
            end
        end
    end

end

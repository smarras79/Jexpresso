using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads

# Global MPI communicator - can be overridden for coupling mode
# Default is COMM_WORLD, but can be set to local communicator in coupled simulations
const JEXPRESSO_MPI_COMM = Ref{Union{MPI.Comm,Nothing}}(nothing)

# World communicator for inter-code coupling (spans all codes)
# In standalone mode this is the same as COMM_WORLD.
# In coupled mode this is MPI.COMM_WORLD while JEXPRESSO_MPI_COMM is the local split.
const JEXPRESSO_MPI_COMM_WORLD = Ref{Union{MPI.Comm,Nothing}}(nothing)

# Coupling data received during initialization (populated by Jexpresso-mini-coupled.jl)
const JEXPRESSO_COUPLING_DATA = Ref{Union{Dict{Symbol,Any},Nothing}}(nothing)

function set_mpi_comm(comm::MPI.Comm)
    JEXPRESSO_MPI_COMM[] = comm
end

function get_mpi_comm()
    return JEXPRESSO_MPI_COMM[] === nothing ? MPI.COMM_WORLD : JEXPRESSO_MPI_COMM[]
end

function set_mpi_comm_world(comm::MPI.Comm)
    JEXPRESSO_MPI_COMM_WORLD[] = comm
end

function get_mpi_comm_world()
    return JEXPRESSO_MPI_COMM_WORLD[] === nothing ? MPI.COMM_WORLD : JEXPRESSO_MPI_COMM_WORLD[]
end

function set_coupling_data(data::Dict{Symbol,Any})
    JEXPRESSO_COUPLING_DATA[] = data
end

function get_coupling_data()
    return JEXPRESSO_COUPLING_DATA[]
end

mutable struct CouplingData
    # Communication pattern
    npoin_recv::Vector{Int32}        # Points to receive from each world rank
    npoin_send::Vector{Int32}        # Points to send to each world rank
    recv_from_ranks::Vector{Int32}   # World ranks we receive from
    send_to_ranks::Vector{Int32}     # World ranks we send to

    # MPI info
    comm_world::MPI.Comm             # MPI.COMM_WORLD
    lrank::Int32                     # Local rank in Jexpresso communicator

    # Problem size
    neqs::Int                        # Number of equations per point
    ndime::Int                       # Spatial dimension (2 or 3)

    # Buffers (allocated later)
    send_bufs::Union{Nothing, Vector{Vector{Float64}}}
    recv_bufs::Union{Nothing, Vector{Vector{Float64}}}
    send_coord_bufs::Union{Nothing, Vector{Vector{Float64}}}  # coordinates sent alongside data

    # Alya coordinate information
    alya_local_coords::Union{Nothing, Matrix{Float64}}
    alya_local_ids::Union{Nothing, Vector{Int32}}
    alya_owner_ranks::Union{Nothing, Vector{Int32}}

    # Constructor
    function CouplingData(; npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
                          comm_world, lrank, neqs, ndime)
        new(npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
            comm_world, lrank, neqs, ndime, nothing, nothing, nothing)
    end
end

function setup_coupling_and_mesh(world, lsize, inputs, nranks, distribute, rank, OUTPUT_DIR, TFloat)

    # Step 2: Receive Alya's domain information
    je_receive_alya_data(world, lsize)
    coupling_data = get_coupling_data()
    local_comm = get_mpi_comm()
    lrank = MPI.Comm_rank(local_comm)
    
    # Setup SEM mesh
    sem, partitioned_model = sem_setup(inputs, nranks, distribute, rank)
    if (inputs[:backend] != CPU()) 
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs) 
    end
    
    # Extract local Alya coordinates using the binned+cropped function
    alya_local_coords, alya_local_ids, alya_owner_ranks =
        extract_local_alya_coordinates(sem.mesh, coupling_data, local_comm, world;
                                       block_size=(64,64,64), use_cropping=true)
    
    # *** NEW: Build npoin_recv DIRECTLY from alya_owner_ranks ***
    # This ensures consistency - we use the same points we extracted
    wsize = MPI.Comm_size(world)
    npoin_recv = zeros(Int32, wsize)
    
    # Count how many points we need from each Alya rank
    for owner_wrank in alya_owner_ranks
        npoin_recv[owner_wrank + 1] += 1
    end
    
    # Build recv_from_ranks list
    recv_from_ranks = Int32[]
    for i in 1:wsize
        if npoin_recv[i] > 0
            push!(recv_from_ranks, i - 1)  # 0-based world rank
        end
    end
    
    if lrank == 0
        println("[setup_coupling] npoin_recv built from alya_owner_ranks, preparing for Alltoall...")
        flush(stdout)
    end
    
    # Step 4 & 5: Alltoall to get npoin_send
    # Alya sends back the same points it received
    # so return counts = forward counts = npoin_recv
    send_counts_to_alya = Vector{Int32}(npoin_recv)
    recv_counts_from_alya = zeros(Int32, wsize)
    MPI.Barrier(world)
    MPI.Alltoall!(send_counts_to_alya, recv_counts_from_alya, 1, world)

    # Ignore recv_counts_from_alya for send counts — Alya initializes its
    # Alltoall buffer to 0 (it never initiates), so this is always zeros.
    # The return exchange uses the same points Julia sent, so counts are identical.
    npoin_send   = copy(npoin_recv)
    send_to_ranks = copy(recv_from_ranks)
    
    # Send the node list immediately after the Alltoall that sent the counts,
    # so that count and actual node list are exchanged at the same point in
    # the handshake sequence.
    ndime = coupling_data[:ndime]
    je_send_node_list(alya_local_ids, alya_owner_ranks, send_to_ranks, world)

    # Verification (should now pass!)
    verify_coupling_communication_pattern(
        npoin_recv, npoin_send,
        alya_owner_ranks, alya_local_ids,
        coupling_data, local_comm, world
    )
    
    # Print short summary
    if lrank == 0
        println("[setup_coupling] Communication pattern summary:")
        println("  npoin_recv (from Alya): ", 
                [(i-1, npoin_recv[i]) for i in 1:wsize if npoin_recv[i] > 0])
        println("  npoin_send (to Alya):   ", 
                [(i-1, npoin_send[i]) for i in 1:wsize if npoin_send[i] > 0])
        println("  recv_from_ranks: ", recv_from_ranks)
        println("  send_to_ranks:   ", send_to_ranks)
        flush(stdout)
    end
    
    if lrank == 0
        println("[setup_coupling] Building Alya point ownership map...")
        flush(stdout)
    end
    ownership = build_alya_point_ownership_map(sem.mesh, coupling_data, local_comm, world)
    print_alya_point_ownership(ownership, coupling_data; max_points_to_print=20, only_owned=true)
    verify_alya_point_distribution(ownership, coupling_data)
    
    # Initialize solution arrays
    qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    # Broadcast neqs to Alya via Allreduce(MAX).
    # Alya contributes 0; Julia contributes qp.neqs.  The collective
    # acts as the synchronisation barrier that lets Alya learn neqs
    # before it allocates its receive buffers.
    neqs_buf = Ref(Int32(qp.neqs))
    MPI.Allreduce!(neqs_buf, MPI.MAX, world)

    # Store coupling information
    ndime = coupling_data[:ndime]
    coupling = CouplingData(
        npoin_recv = npoin_recv,
        npoin_send = npoin_send,
        recv_from_ranks = recv_from_ranks,
        send_to_ranks = send_to_ranks,
        comm_world = world,
        lrank = lrank,
        neqs = qp.neqs,
        ndime = ndime
    )

    # Persist Alya coordinate data in the coupling object for interpolation
    coupling.alya_local_coords = alya_local_coords
    coupling.alya_local_ids    = alya_local_ids
    coupling.alya_owner_ranks  = alya_owner_ranks

    # Allocate coupling buffers (data + coordinates)
    coupling.send_bufs, coupling.recv_bufs, coupling.send_coord_bufs =
        allocate_coupling_buffers(npoin_recv, npoin_send, coupling.neqs, ndime)
    
    # Assert they match
 #   @assert sort(recv_from_ranks) == sort(send_to_ranks) "Coupling partners asymmetric — check alya_owner_ranks vs Alltoall result"
 #   @assert all(npoin_recv .== npoin_send) "Point counts asymmetric between send and recv"

    if lrank == 0
        println("[setup_coupling] Coupling setup complete!")
        flush(stdout)
    end
    
    return coupling, sem, partitioned_model, qp

end


function build_coupling_communication_arrays(mesh, coupling_data, local_comm, world_comm)
    
    # Extract coupling metadata from Alya
    ndime      = coupling_data[:ndime]
    rem_min    = coupling_data[:rem_min]
    rem_max    = coupling_data[:rem_max]
    rem_nx     = coupling_data[:rem_nx]
    alya2world = coupling_data[:alya2world]
    neqs       = coupling_data[:neqs]
    
    # MPI info
    wsize = MPI.Comm_size(world_comm)
    wrank = MPI.Comm_rank(world_comm)
    lsize = MPI.Comm_size(local_comm)
    lrank = MPI.Comm_rank(local_comm)
    
    # Compute Alya grid spacing
    rem_dx = zeros(Float64, 3)
    for idim in 1:ndime
        if rem_nx[idim] > 1
            rem_dx[idim] = (rem_max[idim] - rem_min[idim]) / (rem_nx[idim] - 1)
        else
            rem_dx[idim] = 0.0
        end
    end
    
    # Total number of Alya grid points
    nmax = rem_nx[1] * rem_nx[2] * rem_nx[3]
    
    # Number of Alya ranks
    nranks_alya = length(alya2world)
    
    # Distribution of points among Alya ranks
    r     = mod(nmax, nranks_alya)
    npoin = div(nmax, nranks_alya)
    
    # Initialize receive counts
    npoin_recv_local = zeros(Int32, wsize)
    
    # Local coordinate bounds for this Jexpresso rank
    xmin_local = minimum(mesh.x)
    xmax_local = maximum(mesh.x)
    ymin_local = minimum(mesh.y)
    ymax_local = maximum(mesh.y)
    zmin_local = ndime == 3 ? minimum(mesh.z) : 0.0
    zmax_local = ndime == 3 ? maximum(mesh.z) : 0.0
    
    # Tolerance for floating point comparisons
    tol = 1e-10
    
    if lrank == 0
        println("[build_coupling_arrays] Processing $nmax Alya points...")
        println("  Alya grid: nx=$(rem_nx[1]), ny=$(rem_nx[2]), nz=$(rem_nx[3])")
        println("  Alya bounds: min=$rem_min, max=$rem_max")
        println("  Grid spacing: dx=$rem_dx")
        println("  My domain bounds:")
        println("    x: [$xmin_local, $xmax_local]")
        println("    y: [$ymin_local, $ymax_local]")
        if ndime == 3
            println("    z: [$zmin_local, $zmax_local]")
        end
        flush(stdout)
    end
    
    # Counter for total points in my rank
    points_in_my_rank = 0
    
    # Loop over all Alya grid points
    for ipoin in 1:nmax
        
        # Compute 3D grid indices from flat index (0-based like Fortran)
        i0 = ipoin - 1
        
        ri3 = div(i0, rem_nx[1] * rem_nx[2])
        ri2 = div(i0 - ri3 * rem_nx[1] * rem_nx[2], rem_nx[1])
        ri1 = mod(i0 - ri3 * rem_nx[1] * rem_nx[2] - ri2 * rem_nx[1], rem_nx[1])
        
        # Compute physical coordinates
        x = zeros(Float64, 3)
        x[1] = rem_min[1] + ri1 * rem_dx[1]
        x[2] = rem_min[2] + ri2 * rem_dx[2]
        x[3] = rem_min[3] + ri3 * rem_dx[3]
        
        # Check if this point is in my rank's domain
        in_my_rank = false
        if ndime == 2
            in_my_rank = (x[1] >= xmin_local - tol && x[1] <= xmax_local + tol &&
                         x[2] >= ymin_local - tol && x[2] <= ymax_local + tol)
        else  # ndime == 3
            in_my_rank = (x[1] >= xmin_local - tol && x[1] <= xmax_local + tol &&
                         x[2] >= ymin_local - tol && x[2] <= ymax_local + tol &&
                         x[3] >= zmin_local - tol && x[3] <= zmax_local + tol)
        end
        
        if in_my_rank
            points_in_my_rank += 1
            
            # Determine which Alya rank owns this point
            alya_rank = if ipoin <= r * (npoin + 1)
                div(ipoin - 1, npoin + 1) + 1
            else
                r + div(ipoin - r * (npoin + 1) - 1, npoin) + 1
            end
            
            # Convert to world rank (0-based)
            alya_world_rank = alya2world[alya_rank]
            
            # Increment receive count from this Alya world rank
            npoin_recv_local[alya_world_rank + 1] += 1  # +1 for 1-based indexing
        end
    end
    
    if lrank == 0 || points_in_my_rank > 0
        println("[Rank $lrank] Found $points_in_my_rank Alya points in my domain")
        for i in 1:wsize
            if npoin_recv_local[i] > 0
                println("  Will receive $(npoin_recv_local[i]) points from world rank $(i-1)")
            end
        end
        flush(stdout)
    end
    
    # npoin_recv contains local counts
    npoin_recv = npoin_recv_local
    
    # npoin_send will be populated by Alltoall in the driver
    npoin_send = zeros(Int32, wsize)
    
    # Build lists of active communication partners
    recv_from_ranks = Int32[]
    
    for i in 1:wsize
        if npoin_recv[i] > 0
            push!(recv_from_ranks, i - 1)  # Store as 0-based world rank
        end
    end
    
    if lrank == 0
        total_recv = sum(npoin_recv)
        println("[build_coupling_arrays] Will receive $total_recv Alya points total")
        println("  Number of Alya ranks to receive from: $(length(recv_from_ranks))")
        flush(stdout)
    end
    
    # send_to_ranks will be populated after Alltoall
    send_to_ranks = Int32[]
    
    return npoin_recv, npoin_send, recv_from_ranks, send_to_ranks
end

# ---- Connectivity accessor using connijk (tensor-product ordered) ---------------
# Returns a function get_conn(e)::AbstractVector{Int}
# The returned nodes are in tensor-product order: i (ξ) varies fastest, then j (η),
# matching the loop convention `for j in 1:ngl, i in 1:ngl`.
function _make_conn_accessor(mesh)
    connijk = getfield(mesh, :connijk)
    nsd = mesh.nsd
    if nsd <= 2
        # connijk is [nelem × ngl × ngl × 1] for 2D
        return e -> vec(@view connijk[e, :, :, 1])
    elseif nsd == 3
        # connijk is [nelem × ngl × ngl × ngl] for 3D
        return e -> vec(@view connijk[e, :, :, :])
    else
        error("mesh.nsd must be 1, 2, or 3")
    end
end

# Number of elements in the mesh
function _num_elems(mesh)
    return mesh.nelem
end

# ---- Barycentric weights for 1D Lagrange nodes ----------------------------------
# w_j = 1 / ∏_{k≠j} (x_j - x_k)
function barycentric_weights(nodes::AbstractVector{<:Real})
    n = length(nodes)
    w = ones(Float64, n)
    for j in 1:n
        xj = nodes[j]
        denom = 1.0
        @inbounds for k in 1:n
            if k != j
                denom *= (xj - nodes[k])
            end
        end
        w[j] = 1.0 / denom
    end
    return w
end

# ---- 1D barycentric Lagrange basis and its derivative at x ----------------------
# Uses the Berrut–Trefethen formulation.
# If x matches a node (within `tol`), returns the appropriate one-hot basis and
# analytical derivatives from barycentric formula (stable handling for our use).
function bary_lagrange_and_deriv(x::Float64,
                                 nodes::AbstractVector{<:Real},
                                 w::AbstractVector{<:Real};
                                 tol=1e-12)
    n = length(nodes)
    L  = zeros(Float64, n)
    dL = zeros(Float64, n)

    # Check exact node hit
    for j in 1:n
        dx = x - nodes[j]
        if abs(dx) <= tol
            L[j] = 1.0
            # Derivative vector at node can be computed via:
            # l'_j(x_j) = -Σ_{k≠j} (1/(x_j - x_k))   and
            # l'_k(x_j) = w_k / (w_j * (x_j - x_k)) for k≠j
            s = 0.0
            @inbounds for k in 1:n
                if k != j
                    s += 1.0 / (nodes[j] - nodes[k])
                end
            end
            dL[j] = -s
            @inbounds for k in 1:n
                if k != j
                    dL[k] = w[k] / (w[j] * (nodes[j] - nodes[k]))
                end
            end
            return L, dL
        end
    end

    # General case
    S1 = 0.0           # Σ w_k / (x - x_k)
    S2 = 0.0           # Σ w_k / (x - x_k)^2
    α  = zeros(Float64, n) # α_j = w_j / (x - x_j)
    @inbounds for j in 1:n
        dx = x - nodes[j]
        aj = w[j] / dx
        α[j] = aj
        S1 += aj
        S2 += aj / dx
    end

    invS1 = 1.0 / S1
    @inbounds for j in 1:n
        L[j] = α[j] * invS1
    end

    # Derivative: L'_j = (α'_j * S1 - α_j * S1') / S1^2
    # α'_j = -w_j/(x - x_j)^2 = -(α_j)/(x - x_j)
    # S1'  = -Σ w_k/(x - x_k)^2 = -Σ α_k/(x - x_k) = -S2
    invS1sq = invS1 * invS1
    @inbounds for j in 1:n
        dx   = x - nodes[j]
        aj   = α[j]
        ajp  = -aj / dx
        dL[j] = (ajp * S1 - aj * (-S2)) * invS1sq
    end

    return L, dL
end

# ---- Uniform binning (axis-aligned) over element bounding boxes -----------------
# Returns a struct-like Dict with fields for bin search.
function _build_elem_bins(elem_bboxes::Vector{NTuple{4,Float64}}; bins_per_dim::Int=64)
    ne = length(elem_bboxes)
    xmin = +Inf; xmax = -Inf; ymin = +Inf; ymax = -Inf
    for (x0,x1,y0,y1) in elem_bboxes
        xmin = min(xmin, x0); xmax = max(xmax, x1)
        ymin = min(ymin, y0); ymax = max(ymax, y1)
    end
    nx = max(1, bins_per_dim)
    ny = max(1, bins_per_dim)
    dx = (xmax > xmin) ? (xmax - xmin) / nx : 1.0
    dy = (ymax > ymin) ? (ymax - ymin) / ny : 1.0

    bins = [Int[] for _ in 1:(nx*ny)]
    # Insert each element bbox into all bins it overlaps
    for e in 1:ne
        (x0,x1,y0,y1) = elem_bboxes[e]
        ix0 = clamp(Int(floor((x0 - xmin)/dx)), 0, nx-1)
        ix1 = clamp(Int(floor((x1 - xmin)/dx)), 0, nx-1)
        iy0 = clamp(Int(floor((y0 - ymin)/dy)), 0, ny-1)
        iy1 = clamp(Int(floor((y1 - ymin)/dy)), 0, ny-1)
        for iy in iy0:iy1, ix in ix0:ix1
            push!(bins[iy*nx + ix + 1], e)
        end
    end
    return Dict(
        :xmin => xmin, :xmax => xmax, :ymin => ymin, :ymax => ymax,
        :nx => nx, :ny => ny, :dx => dx, :dy => dy,
        :bins => bins
    )
end

# Collect candidate elements for a point using bins; if empty, fall back to all.
function _bin_candidates(bins, x, y, elem_bboxes)
    nx = bins[:nx]; ny = bins[:ny]
    dx = bins[:dx]; dy = bins[:dy]
    xmin = bins[:xmin]; ymin = bins[:ymin]
    ix = clamp(Int(floor((x - xmin)/dx)), 0, nx-1)
    iy = clamp(Int(floor((y - ymin)/dy)), 0, ny-1)
    cand = bins[:bins][iy*nx + ix + 1]
    if !isempty(cand)
        return cand
    else
        return 1:length(elem_bboxes)
    end
end

#------------------------------------------------------------------------------
#=
    extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm;
                                   block_size=(64,64,64), use_cropping=true)

Extract the subset of Alya grid coordinates that belong to this Jexpresso rank,
with optional lightweight spatial binning and exact index cropping.

# Arguments
- `mesh`: Jexpresso mesh structure (exposes `x`, `y`, and `z` if `ndime == 3`)
- `coupling_data`: Dict with Alya grid metadata:
    - `:ndime::Int`
    - `:rem_min::Vector{Float32}` length 3
    - `:rem_max::Vector{Float32}` length 3
    - `:rem_nx::Vector{Int32}`    length 3
    - `:alya2world::Vector{Int32}` (maps Alya rank → 0-based world rank)
- `local_comm`: Jexpresso local communicator
- `world_comm`: MPI.COMM_WORLD

# Keyword arguments
- `block_size`: tuple of bin sizes in Alya *index* space (default `(64,64,64)`).
  For 2D cases, the z-size is ignored.
- `use_cropping`: if `true`, compute exact index ranges that overlap local bounds
  and only scan those (highly recommended).

# Returns
- `alya_local_coords::Matrix{Float64}` [n_local_points × ndime]
- `alya_local_ids::Vector{Int32}`      [n_local_points]  (1-based Alya IDs)
- `alya_owner_ranks::Vector{Int32}`    [n_local_points]  (0-based Alya world ranks)
=#
#------------------------------------------------------------------------------
function extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm;
                                        block_size::NTuple{3,Int}=(64,64,64),
                                        use_cropping::Bool=true)

    # --- Unpack metadata ---
    ndime      = coupling_data[:ndime]
    neqs       = coupling_data[:neqs]
    rem_min_f  = coupling_data[:rem_min]   # Float32[3]
    rem_max_f  = coupling_data[:rem_max]   # Float32[3]
    rem_nx_i   = coupling_data[:rem_nx]    # Int32[3]
    alya2world = coupling_data[:alya2world]
    @assert ndime == 2 || ndime == 3 "Only ndime==2 or ndime==3 is supported"
    @assert length(alya2world) > 0 "alya2world must be non-empty"

    # Promote to Float64 for math
    rem_min = Float64.(rem_min_f)
    rem_max = Float64.(rem_max_f)
    rem_nx  = Int.(rem_nx_i)

    # Spacings
    rem_dx = zeros(Float64, 3)
    for d in 1:ndime
        rem_dx[d] = rem_nx[d] > 1 ? (rem_max[d] - rem_min[d]) / (rem_nx[d] - 1) : 0.0
    end

    # Total points and Alya rank distribution.
    # alya2world[k] (1-based) = world rank of Alya local rank k-1 (0-based).
    # Alya local rank 0 is always the driving/master rank (world rank 0);
    # it never participates in data exchange with Jexpresso.
    # Points are distributed only over the worker ranks: those whose world
    # rank is not the Alya driving rank (world rank 0).
    nmax = rem_nx[1] * rem_nx[2] * rem_nx[3]
    nranks_alya = length(alya2world)
    alya_driving_world_rank = Int32(0)   # Alya rank 0 is always the master
    # Collect 1-based Julia indices of Alya ranks that are actual workers
    alya_worker_indices = [k for k in 1:nranks_alya
                           if alya2world[k] != alya_driving_world_rank]
    nworkers_alya = length(alya_worker_indices)
    r_w   = mod(nmax, nworkers_alya)
    np_w  = div(nmax, nworkers_alya)

    # Local Jexpresso bounds
    xmin_local = minimum(mesh.x); xmax_local = maximum(mesh.x)
    ymin_local = minimum(mesh.y); ymax_local = maximum(mesh.y)
    zmin_local = ndime == 3 ? minimum(mesh.z) : 0.0
    zmax_local = ndime == 3 ? maximum(mesh.z) : 0.0
    tol = 1e-10

    lrank = MPI.Comm_rank(local_comm)
    wrank = MPI.Comm_rank(world_comm)

    # --- Helpers ---
    # Map structured indices (0-based) to physical coords
    @inline function idx_to_xyz(i1::Int, i2::Int, i3::Int)
        x = rem_min[1] + i1 * rem_dx[1]
        y = rem_min[2] + i2 * rem_dx[2]
        z = rem_min[3] + i3 * rem_dx[3]
        return x, y, z
    end

    # Map structured indices (0-based) to global 1-based Alya ID
    @inline function idx_to_ipoin(i1::Int, i2::Int, i3::Int)
        # i0 = i1 + nx1*(i2 + nx2*i3)
        i0 = i1 + rem_nx[1] * (i2 + rem_nx[2] * i3)
        return i0 + 1
    end

    # Check if a single point is inside local bounds
    @inline function inside_local(x::Float64, y::Float64, z::Float64)
        if ndime == 2
            return (x >= xmin_local - tol && x <= xmax_local + tol &&
                    y >= ymin_local - tol && y <= ymax_local + tol)
        else
            return (x >= xmin_local - tol && x <= xmax_local + tol &&
                    y >= ymin_local - tol && y <= ymax_local + tol &&
                    z >= zmin_local - tol && z <= zmax_local + tol)
        end
    end

    # Map ipoin (1-based) to the 1-based Julia index into alya2world of the
    # worker that owns it.  The driving rank (world rank 0) is excluded;
    # only alya_worker_indices are considered.
    @inline function owner_alya_rank(ipoin::Int)
        iworker = if ipoin <= r_w * (np_w + 1)
            div(ipoin - 1, np_w + 1)       # 0-based index into alya_worker_indices
        else
            r_w + div(ipoin - r_w * (np_w + 1) - 1, np_w)
        end
        return alya_worker_indices[iworker + 1]  # 1-based Julia index into alya2world
    end

    # --- Exact cropping in index space (optional but very effective) ---
    # Find i_lo/i_hi in each dim that may intersect local bounds.
    @inline function crop_1d(minA, maxA, dx, nx, minB, maxB)
        if nx <= 1
            return 0, 0
        else
            # Map physical → index with tolerance, clamp to [0, nx-1]
            ilo = Int(clamp(floor((minB - minA) / dx + 1e-12), 0, nx - 1))
            ihi = Int(clamp( ceil((maxB - minA) / dx - 1e-12), 0, nx - 1))
            # Expand by 1 index to be safe with tolerance if needed
            ilo = max(0, ilo - 1)
            ihi = min(nx - 1, ihi + 1)
            return ilo, ihi
        end
    end

    # Default full ranges
    i1_lo, i1_hi = 0, rem_nx[1] - 1
    i2_lo, i2_hi = 0, rem_nx[2] - 1
    i3_lo, i3_hi = 0, rem_nx[3] - 1

    if use_cropping
        if rem_nx[1] > 1
            i1_lo, i1_hi = crop_1d(rem_min[1], rem_max[1], rem_dx[1], rem_nx[1], xmin_local - tol, xmax_local + tol)
        end
        if rem_nx[2] > 1
            i2_lo, i2_hi = crop_1d(rem_min[2], rem_max[2], rem_dx[2], rem_nx[2], ymin_local - tol, ymax_local + tol)
        end
        if ndime == 3 && rem_nx[3] > 1
            i3_lo, i3_hi = crop_1d(rem_min[3], rem_max[3], rem_dx[3], rem_nx[3], zmin_local - tol, zmax_local + tol)
        else
            i3_lo = 0; i3_hi = 0
        end
        # Early-out if no overlap
        if i1_lo > i1_hi || i2_lo > i2_hi || i3_lo > i3_hi
            if lrank == 0
                println("[extract_local_alya_coordinates] (lrank=$lrank, wrank=$wrank) no Alya points overlap local domain.")
                flush(stdout)
            end
            return zeros(Float64, 0, ndime), Int32[], Int32[]
        end
    else
        if ndime == 2
            i3_lo = 0; i3_hi = 0
        end
    end

    # --- Binning setup ---
    Bx, By, Bz = block_size
    Bx = max(1, Bx); By = max(1, By); Bz = max(1, (ndime == 3 ? Bz : 1))

    # Pre-size hints (upper bound for this rank)
    est_count = (i1_hi - i1_lo + 1) * (i2_hi - i2_lo + 1) * (i3_hi - i3_lo + 1)
    # For 2D, i3 range is 1, so it's fine.

    X = Float64[]; sizehint!(X, est_count)
    Y = Float64[]; sizehint!(Y, est_count)
    Z = ndime == 3 ? (sizehint!(Float64[], est_count)) : Float64[]
    ids    = Int32[]; sizehint!(ids, est_count)
    owners = Int32[]; sizehint!(owners, est_count)

    # --- Iterate over bins with AABB pruning ---
    # Utility to compute bin bbox in physical space
    @inline function bin_bbox(i1s::Int, i1e::Int, i2s::Int, i2e::Int, i3s::Int, i3e::Int)
        # Coordinates are monotonic (dx >= 0 assumed)
        x_min = rem_min[1] + i1s * rem_dx[1]
        x_max = rem_min[1] + i1e * rem_dx[1]
        y_min = rem_min[2] + i2s * rem_dx[2]
        y_max = rem_min[2] + i2e * rem_dx[2]
        z_min = rem_min[3] + i3s * rem_dx[3]
        z_max = rem_min[3] + i3e * rem_dx[3]
        return x_min, x_max, y_min, y_max, z_min, z_max
    end

    @inline function bbox_outside_local(xmin::Float64, xmax::Float64,
                                        ymin::Float64, ymax::Float64,
                                        zmin::Float64, zmax::Float64)
        if ndime == 2
            return (xmax < xmin_local - tol || xmin > xmax_local + tol ||
                    ymax < ymin_local - tol || ymin > ymax_local + tol)
        else
            return (xmax < xmin_local - tol || xmin > xmax_local + tol ||
                    ymax < ymin_local - tol || ymin > ymax_local + tol ||
                    zmax < zmin_local - tol || zmin > zmax_local + tol)
        end
    end

    @inline function bbox_inside_local(xmin::Float64, xmax::Float64,
                                       ymin::Float64, ymax::Float64,
                                       zmin::Float64, zmax::Float64)
        if ndime == 2
            return (xmin >= xmin_local - tol && xmax <= xmax_local + tol &&
                    ymin >= ymin_local - tol && ymax <= ymax_local + tol)
        else
            return (xmin >= xmin_local - tol && xmax <= xmax_local + tol &&
                    ymin >= ymin_local - tol && ymax <= ymax_local + tol &&
                    zmin >= zmin_local - tol && zmax <= zmax_local + tol)
        end
    end

    # Iterate over cropped index ranges in bins
    for i3s in i3_lo: Bz : i3_hi
        i3e = min(i3s + Bz - 1, i3_hi)
        for i2s in i2_lo: By : i2_hi
            i2e = min(i2s + By - 1, i2_hi)
            for i1s in i1_lo: Bx : i1_hi
                i1e = min(i1s + Bx - 1, i1_hi)

                # Bin bbox in physical space
                bxmin, bxmax, bymin, bymax, bzmin, bzmax = bin_bbox(i1s, i1e, i2s, i2e, i3s, i3e)

                # Reject bin quickly if entirely outside
                if bbox_outside_local(bxmin, bxmax, bymin, bymax, bzmin, bzmax)
                    continue
                end

                # If bin is entirely inside local domain, add all bin points fast
                if bbox_inside_local(bxmin, bxmax, bymin, bymax, bzmin, bzmax)
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x, y, z = idx_to_xyz(i1, i2, i3)
                        ipoin = idx_to_ipoin(i1, i2, i3)
                        alya_rank = owner_alya_rank(ipoin)
                        alya_world_rank = alya2world[alya_rank]  # 0-based

                        push!(X, x); push!(Y, y)
                        if ndime == 3
                            push!(Z, z)
                        end
                        push!(ids,    Int32(ipoin))
                        push!(owners, Int32(alya_world_rank))
                    end
                else
                    # Partial overlap: check points individually
                    @inbounds for i3 in i3s:i3e, i2 in i2s:i2e, i1 in i1s:i1e
                        x, y, z = idx_to_xyz(i1, i2, i3)
                        if inside_local(x, y, z)
                            ipoin = idx_to_ipoin(i1, i2, i3)
                            alya_rank = owner_alya_rank(ipoin)
                            alya_world_rank = alya2world[alya_rank]  # 0-based

                            push!(X, x); push!(Y, y)
                            if ndime == 3
                                push!(Z, z)
                            end
                            push!(ids,    Int32(ipoin))
                            push!(owners, Int32(alya_world_rank))
                        end
                    end
                end
            end
        end
    end

    # Assemble matrix
    n_local = length(ids)
    alya_local_coords = zeros(Float64, n_local, ndime)
    if ndime == 2
        @inbounds for i in 1:n_local
            alya_local_coords[i,1] = X[i]
            alya_local_coords[i,2] = Y[i]
        end
    else
        @inbounds for i in 1:n_local
            alya_local_coords[i,1] = X[i]
            alya_local_coords[i,2] = Y[i]
            alya_local_coords[i,3] = Z[i]
        end
    end

    # Sort by (owner_rank, global_id) so that data sent to each Alya rank
    # arrives in ascending flat-index order, matching the Fortran Gatherv
    # displacement layout.  Without this sort the bin-iteration order of
    # extract_local_alya_coordinates may differ from flat-index order for
    # grids spanning more than one spatial bin.
    if n_local > 1
        perm = sortperm(collect(zip(owners, ids)))
        alya_local_coords = alya_local_coords[perm, :]
        ids    = ids[perm]
        owners = owners[perm]
    end

    if lrank == 0 || n_local > 0
        println("[extract_local_alya_coordinates] (lrank=$lrank, wrank=$wrank) ",
                "collected $n_local Alya points using cropping=$(use_cropping) ",
                "and block_size=$(block_size).")
        flush(stdout)
    end

    return alya_local_coords, ids, owners
end

function allocate_coupling_buffers(npoin_recv, npoin_send, neqs, ndime)
    wsize = length(npoin_recv)

    send_bufs = Vector{Vector{Float64}}(undef, wsize)
    recv_bufs = Vector{Vector{Float64}}(undef, wsize)
    send_coord_bufs = Vector{Vector{Float64}}(undef, wsize)

    for i in 1:wsize
        # Allocate send buffer (neqs values per point)
        nsend = npoin_send[i] * neqs
        send_bufs[i] = zeros(Float64, nsend)

        # Allocate receive buffer
        nrecv = npoin_recv[i] * neqs
        recv_bufs[i] = zeros(Float64, nrecv)

        # Allocate coordinate send buffer (ndime coords per point)
        send_coord_bufs[i] = zeros(Float64, npoin_send[i] * ndime)
    end

    return send_bufs, recv_bufs, send_coord_bufs
end

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

#----------------------------------------------------------------------------------------
#    je_send_node_list(mesh, send_to_ranks, world)
#
# Send the global node ID list (from connijk) to each Alya communication partner.
# Called once during coupling setup, immediately after the Alltoall that exchanges
# node counts. Alya already knows the count from the Alltoall result, so only the
# sorted unique global node IDs are sent here (one message per partner).
#----------------------------------------------------------------------------------------
function je_send_node_list(alya_local_ids::Vector{Int64},
                           alya_owner_ranks::Vector{Int32},
                           send_to_ranks::Vector{Int32},
                           world::MPI.Comm)
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    wrank = MPI.Comm_rank(world)

    send_requests = MPI.Request[]
    for dest_rank in send_to_ranks
        # Send only the Alya grid-point IDs whose owner is dest_rank
        mask   = alya_owner_ranks .== dest_rank
        gid_buf = Int64.(alya_local_ids[mask])
        push!(send_requests, MPI.Isend(gid_buf, dest_rank, 0, world))
        println("[je_send_node_list] lrank=$lrank (wrank=$wrank): sending $(length(gid_buf)) Alya point IDs to world rank $dest_rank")
    end
    isempty(send_requests) || MPI.Waitall(send_requests)
    flush(stdout)

    return nothing
end

function je_receive_alya_data(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    
    if wsize <= nparts
        @warn "je_receive_alya_data called but not in coupled mode"
        return
    end
    
    # 1. ndime  (Fortran STEP 2, first broadcast)
    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    ndime = Int(ndime_buf[1])

    # 2. rem_min, rem_max, rem_nx  (Fortran STEP 2)
    # Broadcast the full 3-element arrays in one call each to avoid
    # SubArray aliasing issues that occur when using per-element views.
    rem_min = Vector{Float64}(undef, 3)
    rem_max = Vector{Float64}(undef, 3)
    rem_nx  = Vector{Int32}(undef, 3)
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end

    # 3. neqs  (Fortran: after rem_nx loop)
    neqs_buf = Vector{Int32}(undef, 1)
    # MPI.Bcast!(neqs_buf, 0, world)
    neqs = Int(neqs_buf[1])

    # 4. nsteps  (Fortran: after neqs)
    nsteps_buf = Vector{Int32}(undef, 1)
    # MPI.Bcast!(nsteps_buf, 0, world)
    nsteps = Int(nsteps_buf[1])

    # 5. Alya->world rank map  (Fortran: MPI_AllReduce)
    nranks_alya  = wsize - nparts
    alya2world_l = zeros(Int32, nranks_alya)

    flush(stdout)
    
    alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)
    
    set_coupling_data(Dict{Symbol,Any}(
        :ndime   => ndime,
        :neqs    => neqs,
        :nsteps  => nsteps,
        :rem_min => rem_min,
        :rem_max => rem_max,
        :rem_nx  => rem_nx,
        :alya2world => alya2world,
    ))
    
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    if lrank == 0
        println("[je_receive_alya_data] ndime=$ndime, neqs=$neqs, nsteps=$nsteps")
        println("  min=$rem_min, max=$rem_max, nx=$rem_nx")
        flush(stdout)
    end
    
    return nothing
end

#------------------------------------------------------------------------------------
# Functions for Step 6: Data exchange during time loop
# These will be used when implementing velocity interpolation and exchange
#------------------------------------------------------------------------------------

function interpolate_at_point(px::Float64, py::Float64, mesh, u, qp)
    # TODO: Implement proper interpolation using basis functions
    # For now, this is a placeholder
    return NaN
end

function point_in_quad(px, py, coords::Vector{Tuple{Float64,Float64}})
    # Simple bounding box check
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    return (minimum(xs) <= px <= maximum(xs) && 
            minimum(ys) <= py <= maximum(ys))
end


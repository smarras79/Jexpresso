using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads

"""
    CouplingData

Structure to store coupling information between Alya and Jexpresso.
"""
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
    
    # Buffers (allocated later)
    send_bufs::Union{Nothing, Vector{Vector{Float64}}}
    recv_bufs::Union{Nothing, Vector{Vector{Float64}}}
    
    # Constructor
    function CouplingData(; npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
                          comm_world, lrank, neqs)
        new(npoin_recv, npoin_send, recv_from_ranks, send_to_ranks,
            comm_world, lrank, neqs, nothing, nothing)
    end
end


"""
    setup_coupling_and_mesh(world, lsize, inputs, nranks, distribute, rank, OUTPUT_DIR, TFloat)

Complete coupling initialization including:
- Receiving Alya's domain information
- Setting up SEM mesh
- Building coupling communication pattern
- Exchanging pattern with Alya via MPI_Alltoall
- Initializing solution arrays
- Allocating coupling buffers

# Arguments
- `world`: MPI.COMM_WORLD communicator
- `lsize`: Number of local Jexpresso ranks
- `inputs`: Input parameters dictionary
- `nranks`: Number of ranks for mesh setup
- `distribute`: Distribution function for partitioning
- `rank`: Current rank
- `OUTPUT_DIR`: Output directory path
- `TFloat`: Floating point type (Float32 or Float64)

# Returns
- `coupling`: CouplingData structure with communication pattern and buffers
- `sem`: SEM mesh and solver structures
- `partitioned_model`: Partitioned mesh model
- `qp`: Initialized solution arrays
"""
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
    
    # Step 3: Build coupling communication pattern
    npoin_recv, npoin_send, recv_from_ranks, send_to = 
        build_coupling_communication_arrays(sem.mesh, coupling_data, 
                                            local_comm, world)
    
    if lrank == 0
        println("[setup_coupling] Coupling arrays built, preparing for Alltoall...")
        flush(stdout)
    end
    
    # Step 4 & 5: Communicate the pattern back to Alya via Alltoall
    send_counts_to_alya = Vector{Int32}(npoin_recv)
    recv_counts_from_alya = zeros(Int32, MPI.Comm_size(world))
    
    if lrank == 0
        println("[setup_coupling] Julia: Before barrier before Alltoall")
        flush(stdout)
    end
    
    # Barrier to synchronize with Alya
    MPI.Barrier(world)
    
    if lrank == 0
        println("[setup_coupling] Julia: About to call MPI.Alltoall!")
        println("  Sending counts: ", send_counts_to_alya)
        flush(stdout)
    end
    
    # Alltoall: Tell Alya how many points I need from each of its ranks
    MPI.Alltoall!(send_counts_to_alya, recv_counts_from_alya, 1, world)
    
    if lrank == 0
        println("[setup_coupling] Julia: MPI.Alltoall completed!")
        flush(stdout)
    end
    
    # Update npoin_send from received data
    npoin_send = recv_counts_from_alya
    
    # Update send_to_ranks based on actual send pattern
    send_to_ranks = Int32[]
    for i in 1:length(npoin_send)
        if npoin_send[i] > 0
            push!(send_to_ranks, i - 1)  # 0-based world rank
        end
    end
    
    # VERIFICATION: Print detailed coupling pattern
    print_coupling_pattern(lrank, npoin_recv, npoin_send, recv_from_ranks, 
                          send_to_ranks, coupling_data, world)
    
    # Initialize solution arrays
    qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)
    
    # Store coupling information
    coupling = CouplingData(
        npoin_recv = npoin_recv,
        npoin_send = npoin_send,
        recv_from_ranks = recv_from_ranks,
        send_to_ranks = send_to_ranks,
        comm_world = world,
        lrank = lrank,
        neqs = qp.neqs
    )
    
    # Allocate coupling buffers
    coupling.send_bufs, coupling.recv_bufs = 
        allocate_coupling_buffers(npoin_recv, npoin_send, coupling.neqs)
    
    if lrank == 0
        println("[setup_coupling] Coupling setup complete!")
        flush(stdout)
    end
    
    return coupling, sem, partitioned_model, qp
end


"""
    print_coupling_pattern(lrank, npoin_recv, npoin_send, recv_from_ranks, 
                          send_to_ranks, coupling_data, world)

Print detailed verification of the coupling communication pattern.
"""
function print_coupling_pattern(lrank, npoin_recv, npoin_send, recv_from_ranks, 
                               send_to_ranks, coupling_data, world)
    
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    nranks_alya = length(coupling_data[:alya2world])
    
    println("=========================================")
    println("Julia local rank $lrank, world rank $wrank")
    println("=========================================")
    println("Total points to RECV from Alya: $(sum(npoin_recv))")
    println("Total points to SEND to Alya: $(sum(npoin_send))")
    println("-----------------------------------------")
    
    # Print world rank breakdown
    println("World rank breakdown (size=$wsize):")
    for i in 0:nranks_alya-1
        println("  World rank $(coupling_data[:alya2world][i+1]) = Alya rank $i (ALYA)")
    end
    for i in 0:(wsize - nranks_alya - 1)
        println("  World rank $(nranks_alya + i) = Julia rank $i (JULIA)")
    end
    println("-----------------------------------------")
    
    # Print detailed receive pattern
    if sum(npoin_recv) > 0
        println("RECV pattern (npoin_recv array) - receive FROM these ranks:")
        for i in 1:wsize
            if npoin_recv[i] > 0
                world_rank = i - 1
                if world_rank < nranks_alya
                    println("  <- World rank $world_rank (Alya rank $world_rank): $(npoin_recv[i]) points")
                else
                    julia_rank = world_rank - nranks_alya
                    println("  <- World rank $world_rank (Julia rank $julia_rank): $(npoin_recv[i]) points")
                end
            end
        end
    else
        println("RECV pattern: No points to receive")
    end
    
    # Print detailed send pattern
    if sum(npoin_send) > 0
        println("SEND pattern (npoin_send array) - send TO these ranks:")
        for i in 1:wsize
            if npoin_send[i] > 0
                world_rank = i - 1
                if world_rank < nranks_alya
                    println("  -> World rank $world_rank (Alya rank $world_rank): $(npoin_send[i]) points")
                else
                    julia_rank = world_rank - nranks_alya
                    println("  -> World rank $world_rank (Julia rank $julia_rank): $(npoin_send[i]) points")
                end
            end
        end
    else
        println("SEND pattern: No points to send")
    end
    
    println("-----------------------------------------")
    println("Active communication partners:")
    println("  recv_from_ranks: $recv_from_ranks")
    println("  send_to_ranks: $send_to_ranks")
    println("=========================================")
    flush(stdout)
end


"""
    build_coupling_communication_arrays(mesh, coupling_data, local_comm, world_comm)

Build send/receive arrays for coupling communication between Alya and Jexpresso.
Determines which Alya grid points fall into which Jexpresso rank domains.

# Arguments
- `mesh`: Jexpresso mesh structure with local coordinates
- `coupling_data`: Dict containing remote grid metadata from Alya
- `local_comm`: MPI communicator for Jexpresso ranks only
- `world_comm`: MPI.COMM_WORLD (shared with Alya)

# Returns
- `npoin_recv`: Array[wsize] - number of Alya points in my domain from each world rank
- `npoin_send`: Array[wsize] - number of points to send to each world rank (filled by Alltoall)
- `recv_from_ranks`: Vector of world ranks we receive from (non-zero entries)
- `send_to_ranks`: Vector of world ranks we send to (non-zero entries)
"""
function build_coupling_communication_arrays(mesh, coupling_data, local_comm, world_comm)
    
    # Extract coupling metadata from Alya
    ndime      = coupling_data[:ndime]
    rem_min    = coupling_data[:rem_min]
    rem_max    = coupling_data[:rem_max]
    rem_nx     = coupling_data[:rem_nx]
    alya2world = coupling_data[:alya2world]
    
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

"""
    allocate_coupling_buffers(npoin_recv, npoin_send, neqs)

Allocate send and receive buffers based on communication pattern.

# Arguments
- `npoin_recv`: Array of receive counts per rank
- `npoin_send`: Array of send counts per rank  
- `neqs`: Number of equations (variables per point)

# Returns
- `send_bufs`: Vector of buffers for sending to each rank
- `recv_bufs`: Vector of buffers for receiving from each rank
"""
function allocate_coupling_buffers(npoin_recv, npoin_send, neqs)
    wsize = length(npoin_recv)
    
    send_bufs = Vector{Vector{Float64}}(undef, wsize)
    recv_bufs = Vector{Vector{Float64}}(undef, wsize)
    
    for i in 1:wsize
        # Allocate send buffer
        nsend = npoin_send[i] * neqs
        send_bufs[i] = zeros(Float64, nsend)
        
        # Allocate receive buffer
        nrecv = npoin_recv[i] * neqs
        recv_bufs[i] = zeros(Float64, nrecv)
    end
    
    return send_bufs, recv_bufs
end

"""
    je_perform_coupling_handshake(world, nparts)

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
    je_receive_alya_data(world, nparts)

Receive grid metadata from Alya via Bcast.
"""
function je_receive_alya_data(world, nparts)
    wsize = MPI.Comm_size(world)
    wrank = MPI.Comm_rank(world)
    
    # Verify we are in coupled mode
    if wsize <= nparts
        @warn "je_receive_alya_data called but not in coupled mode"
        return
    end
    
    # Receive grid dimensions
    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    ndime = ndime_buf[1]
    
    # Receive grid bounds and resolution
    rem_min = Vector{Float32}(undef, 3)
    rem_max = Vector{Float32}(undef, 3)
    rem_nx  = Vector{Int32}(undef, 3)
    
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end
    
    # Build Alya to world rank mapping
    nranks_other = wsize - nparts
    alya2world_l = zeros(Int32, nranks_other)
    alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)
    
    # Store coupling data
    set_coupling_data(Dict{Symbol,Any}(
        :ndime         => ndime,
        :rem_min       => rem_min,
        :rem_max       => rem_max,
        :rem_nx        => rem_nx,
        :alya2world    => alya2world,
    ))
    
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    
    if lrank == 0
        println("[je_receive_alya_data] Received from Alya:")
        println("  ndime=$ndime, min=$rem_min, max=$rem_max, nx=$rem_nx")
        flush(stdout)
    end
    
    return nothing
end

#------------------------------------------------------------------------------------
# Functions for Step 6: Data exchange during time loop
# These will be used when implementing velocity interpolation and exchange
#------------------------------------------------------------------------------------

"""
    interpolate_at_point(px, py, mesh, u, qp)

Interpolate solution at point (px, py) using Julia's mesh and basis functions.
Returns interpolated value for the first solution component.

NOTE: This is a placeholder for Step 6 implementation.
"""
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


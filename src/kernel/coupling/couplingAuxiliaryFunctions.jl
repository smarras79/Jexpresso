#------------------------------------------------------------------------------------
#    Alya_Point_Ownership
#
#Structure to store which Jexpresso rank owns each Alya grid point.
#------------------------------------------------------------------------------------
struct Alya_Point_Ownership
    alya_coords::Matrix{Float64}       # Alya point coordinates [npoints × ndime]
    owner_jrank::Vector{Int32}         # Jexpresso local rank that owns each point (-1 if not owned)
    owner_wrank::Vector{Int32}         # World rank that owns each point (-1 if not owned)
    alya_point_id::Vector{Int32}       # Original Alya point ID (1-based)
    ndime::Int32                       # Spatial dimension
end

#------------------------------------------------------------------------------------
#
#    build_alya_point_ownership_map(mesh, coupling_data, local_comm, world_comm)
#
#Build a complete map of which Jexpresso ranks own which Alya grid points.
#This is useful for debugging and verification.
#
## Arguments
#- `mesh`: Jexpresso mesh structure with local coordinates
#- `coupling_data`: Dict containing remote grid metadata from Alya
#- `local_comm`: MPI communicator for Jexpresso ranks only
#- `world_comm`: MPI.COMM_WORLD (shared with Alya)
#
## Returns
#- `ownership`: Alya_Point_Ownership structure with complete mapping
#------------------------------------------------------------------------------------
function build_alya_point_ownership_map(mesh, coupling_data, local_comm, world_comm)
    
    # Extract coupling metadata from Alya
    ndime      = coupling_data[:ndime]
    neqs       = coupling_data[:neqs]
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
    
    # Local coordinate bounds for this Jexpresso rank
    xmin_local = minimum(mesh.x)
    xmax_local = maximum(mesh.x)
    ymin_local = minimum(mesh.y)
    ymax_local = maximum(mesh.y)
    zmin_local = ndime == 3 ? minimum(mesh.z) : 0.0
    zmax_local = ndime == 3 ? maximum(mesh.z) : 0.0
    
    # Tolerance for floating point comparisons
    tol = 1e-10
    
    # Arrays to store ownership information (locally)
    local_coords      = zeros(Float64, nmax, ndime)
    local_owner_jrank = fill(Int32(-1), nmax)  # -1 means not owned by any rank
    local_owner_wrank = fill(Int32(-1), nmax)
    
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
        
        # Store coordinates
        for idim in 1:ndime
            local_coords[ipoin, idim] = x[idim]
        end
        
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
            # This point is in my domain
            local_owner_jrank[ipoin] = lrank
            local_owner_wrank[ipoin] = wrank
        end
    end
    
    # Gather ownership information from all Jexpresso ranks
    # Use Allreduce with MAX to get the owner (since only one rank claims each point)
    global_owner_jrank = MPI.Allreduce(local_owner_jrank, MPI.MAX, local_comm)
    global_owner_wrank = MPI.Allreduce(local_owner_wrank, MPI.MAX, local_comm)
    
    # Create point IDs (1-based)
    alya_point_ids = collect(Int32(1):Int32(nmax))
    
    # Create ownership structure
    ownership = Alya_Point_Ownership(local_coords,
                                     global_owner_jrank,
                                     global_owner_wrank,
                                     alya_point_ids,
                                     Int32(ndime)
                                     )
    
    return ownership
end


#------------------------------------------------------------------------------------
#    print_alya_point_ownership(ownership, coupling_data; 
#                               max_points_to_print=20,
#                               only_owned=false)
#
#Print which Jexpresso rank owns each Alya grid point.
#
## Arguments
#- `ownership`: Alya_Point_Ownership structure
#- `coupling_data`: Dict containing coupling metadata
#- `max_points_to_print`: Maximum number of points to print (default: 20, use -1 for all)
#- `only_owned`: If true, only print points that are owned by some Jexpresso rank
#
#------------------------------------------------------------------------------------
function print_alya_point_ownership(ownership, coupling_data; 
                                   max_points_to_print=20,
                                   only_owned=false)
    
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    
    # Only rank 0 prints
    if lrank != 0
        return
    end
    
    npoints = length(ownership.alya_point_id)
    ndime = ownership.ndime
    nranks_alya = length(coupling_data[:alya2world])
    
    println("\n" * "="^70)
    println("Alya Point Ownership Map")
    println("="^70)
    println("Total Alya points: $npoints")
    println("Spatial dimension: $ndime")
    println("-"^70)
    
    # Count owned vs unowned
    n_owned = count(x -> x >= 0, ownership.owner_jrank)
    n_unowned = npoints - n_owned
    
    println("Points owned by Jexpresso: $n_owned")
    println("Points not owned: $n_unowned")
    println("-"^70)
    
    # Determine how many points to print
    n_to_print = max_points_to_print
    if max_points_to_print < 0
        n_to_print = npoints
    end
    
    # Print header
    if ndime == 2
        println("Point ID | Coordinates (x, y)          | Owner")
    else
        println("Point ID | Coordinates (x, y, z)                | Owner")
    end
    println("-"^70)
    
    # Print point ownership
    printed_count = 0
    for i in 1:npoints
        # Skip if only showing owned points and this point is not owned
        if only_owned && ownership.owner_jrank[i] < 0
            continue
        end
        
        # Stop if we've printed enough
        if printed_count >= n_to_print
            if only_owned
                remaining = count(x -> x >= 0, ownership.owner_jrank[i+1:end])
            else
                remaining = npoints - i
            end
            if remaining > 0
                println("... ($(remaining) more points not shown)")
            end
            break
        end
        
        # Format coordinates
        if ndime == 2
            coord_str = @sprintf("(%10.2f, %10.2f)", 
                                ownership.alya_coords[i, 1], 
                                ownership.alya_coords[i, 2])
        else
            coord_str = @sprintf("(%10.2f, %10.2f, %10.2f)", 
                                ownership.alya_coords[i, 1], 
                                ownership.alya_coords[i, 2],
                                ownership.alya_coords[i, 3])
        end
        
        # Format owner information
        jrank = ownership.owner_jrank[i]
        wrank = ownership.owner_wrank[i]
        
        if jrank >= 0
            owner_str = @sprintf("Jexpresso rank %d (world rank %d)", jrank, wrank)
        else
            owner_str = "NOT OWNED"
        end
        
        # Print the line
        println(@sprintf("%8d | %s | %s", 
                        ownership.alya_point_id[i], 
                        coord_str, 
                        owner_str))
        
        printed_count += 1
    end
    
    println("="^70)
    println()
    flush(stdout)
end

#------------------------------------------------------------------------------------
#    verify_alya_point_distribution(ownership, coupling_data)
#
#  Print statistics about how Alya points are distributed across Jexpresso ranks.
#------------------------------------------------------------------------------------
function verify_alya_point_distribution(ownership, coupling_data)
    
    lcomm = get_mpi_comm()
    lrank = MPI.Comm_rank(lcomm)
    lsize = MPI.Comm_size(lcomm)
    
    # Only rank 0 prints
    if lrank != 0
        return
    end
    
    npoints = length(ownership.alya_point_id)
    nranks_alya = length(coupling_data[:alya2world])
    
    println("\n" * "="^70)
    println("Alya Point Distribution Statistics")
    println("="^70)
    
    # Count points per Jexpresso rank
    points_per_rank = zeros(Int, lsize)
    for i in 1:npoints
        jrank = ownership.owner_jrank[i]
        if jrank >= 0
            points_per_rank[jrank + 1] += 1  # +1 for 1-based indexing
        end
    end
    
    println("Points per Jexpresso rank:")
    for r in 0:lsize-1
        println(@sprintf("  Jexpresso rank %d: %d points (%.1f%%)", 
                        r, 
                        points_per_rank[r+1],
                        100.0 * points_per_rank[r+1] / npoints))
    end
    
    println("-"^70)
    
    # Check for unowned points
    n_unowned = count(x -> x < 0, ownership.owner_jrank)
    if n_unowned > 0
        println("WARNING: $n_unowned Alya points are not owned by any Jexpresso rank!")
        println("This may indicate a domain mismatch or tolerance issue.")
    else
        println("✓ All Alya points are owned by Jexpresso ranks")
    end
    
    println("="^70)
    println()
    flush(stdout)
end

function verify_extracted_vs_counts!(alya_owner_ranks::Vector{Int32}, npoin_recv::Vector{Int32})
    # Tally by 0-based world rank
    by_world = Dict{Int32,Int}()
    for r in alya_owner_ranks
        by_world[r] = get(by_world, r, 0) + 1
    end
    # Compare to npoin_recv (1-based with (world_rank+1) indexing)
    for (wr, cnt) in by_world
        expected = npoin_recv[Int(wr) + 1]
        @assert cnt == expected "Mismatch for world_rank=$wr: extracted=$cnt vs expected=$expected"
    end
    return nothing
end


#------------------------------------------------------------------------------------
#=
    coupling_exchange_data!(cpg::CouplingData)

Perform MPI data exchange between Jexpresso and Alya using non-blocking communication.
Exchanges data stored in cpg.send_bufs with partners and receives into cpg.recv_bufs.

Uses the communication pattern established during setup:
- Send to ranks in send_to_ranks with data from send_bufs
- Receive from ranks in recv_from_ranks into recv_bufs
=#
#------------------------------------------------------------------------------------
"""
    coupling_exchange_data!(cpg::CouplingData)

Exchange data AND coordinates with Alya using non-blocking communication.
- Data   send to rank R: tag = TAG_DATA  + R
- Coord  send to rank R: tag = TAG_COORD + R
- Data receive from R:   tag = TAG_DATA  + MY_RANK
"""
function coupling_exchange_data!(cpg::CouplingData)
    wsize = length(cpg.npoin_send)

    # Get my world rank
    lcomm = get_mpi_comm()
    wrank = MPI.Comm_rank(cpg.comm_world)

    # Allocate arrays for MPI requests
    send_requests = MPI.Request[]
    recv_requests = MPI.Request[]

    # Tags (must match Alya Fortran code)
    TAG_DATA  = 2000
    TAG_COORD = 3000

    # Post all non-blocking receives (data from Alya → Julia)
    for src_rank in cpg.recv_from_ranks
        if cpg.npoin_recv[src_rank + 1] > 0
            buf = cpg.recv_bufs[src_rank + 1]
            tag = TAG_DATA + wrank
            req = MPI.Irecv!(buf, src_rank, tag, cpg.comm_world)
            push!(recv_requests, req)
        end
    end

    # Post all non-blocking sends (data + coordinates, Julia → Alya)
    for dest_rank in cpg.send_to_ranks
        if cpg.npoin_send[dest_rank + 1] > 0
            # Send interpolated variable data
            buf = cpg.send_bufs[dest_rank + 1]
            tag = TAG_DATA + dest_rank
            req = MPI.Isend(buf, dest_rank, tag, cpg.comm_world)
            push!(send_requests, req)

            # Send coordinates so Alya can map values to grid positions
            cbuf = cpg.send_coord_bufs[dest_rank + 1]
            ctag = TAG_COORD + dest_rank
            creq = MPI.Isend(cbuf, dest_rank, ctag, cpg.comm_world)
            push!(send_requests, creq)
        end
    end

    # Wait for all receives to complete
    if !isempty(recv_requests)
        MPI.Waitall(recv_requests)
    end

    # Wait for all sends to complete
    if !isempty(send_requests)
        MPI.Waitall(send_requests)
    end

    return nothing
end

function coupling_exchange_data_old!(cpg::CouplingData)
    wsize = length(cpg.npoin_send)
    
    # Allocate arrays for MPI requests
    send_requests = MPI.Request[]
    recv_requests = MPI.Request[]
    
    # Post all non-blocking receives first (good practice)
    for (idx, src_rank) in enumerate(cpg.recv_from_ranks)
        if cpg.npoin_recv[src_rank + 1] > 0
            buf = cpg.recv_bufs[src_rank + 1]
            tag = 2000 + src_rank  # Unique tag per sender
            req = MPI.Irecv!(buf, src_rank, tag, cpg.comm_world)
            push!(recv_requests, req)
        end
    end
    
    # Post all non-blocking sends
    for (idx, dest_rank) in enumerate(cpg.send_to_ranks)
        if cpg.npoin_send[dest_rank + 1] > 0
            buf = cpg.send_bufs[dest_rank + 1]
            tag = 2000 + cpg.lrank  # Unique tag per sender (use local rank converted to world rank)
            # Need to convert local rank to world rank
            lcomm = get_mpi_comm()
            lrank = MPI.Comm_rank(lcomm)
            # Find my world rank (Jexpresso ranks start after Alya ranks)
            nranks_alya = length(cpg.npoin_send) - MPI.Comm_size(lcomm)
            my_world_rank = nranks_alya + lrank
            tag = 2000 + my_world_rank
            req = MPI.Isend(buf, dest_rank, tag, cpg.comm_world)
            push!(send_requests, req)
        end
    end
    
    # Wait for all receives to complete
    if !isempty(recv_requests)
        MPI.Waitall(recv_requests)
    end
    
    # Wait for all sends to complete
    if !isempty(send_requests)
        MPI.Waitall(send_requests)
    end
    
    return nothing
end


#------------------------------------------------------------------------------------
#=

    pack_interpolated_data!(cpg::CouplingData, interp_values::Matrix{Float64}, 
                           alya_owner_ranks::Vector{Int32})

Pack interpolated values into send buffers according to which Alya rank owns each point.

# Arguments
- `cpg`: CouplingData structure with send_bufs allocated
- `interp_values`: [n_local_points × neqs] interpolated values at Alya coordinates
- `alya_owner_ranks`: [n_local_points] world rank that owns each Alya point (0-based)
=#
#------------------------------------------------------------------------------------
function pack_interpolated_data!(cpg::CouplingData, interp_values::Matrix{Float64},
                                 alya_owner_ranks::Vector{Int32},
                                 alya_local_coords::Matrix{Float64})
    n_local = size(interp_values, 1)
    neqs  = cpg.neqs
    ndime = cpg.ndime
    wsize = length(cpg.npoin_send)

    # Zero out send buffers (data + coordinates)
    for i in 1:wsize
        fill!(cpg.send_bufs[i], 0.0)
        fill!(cpg.send_coord_bufs[i], 0.0)
    end

    # Track how many values we've packed for each destination rank
    send_offsets = zeros(Int, wsize)
    coord_offsets = zeros(Int, wsize)

    # Pack data AND coordinates for each Alya point
    @inbounds for i in 1:n_local
        owner_rank = alya_owner_ranks[i]  # 0-based world rank
        buf_idx = owner_rank + 1  # 1-based indexing

        if cpg.npoin_send[buf_idx] > 0
            # Pack variable values
            offset = send_offsets[buf_idx]
            buf = cpg.send_bufs[buf_idx]
            for q in 1:neqs
                buf[offset + q] = interp_values[i, q]
            end
            send_offsets[buf_idx] += neqs

            # Pack coordinates so Alya can map each value to its grid position
            coffset = coord_offsets[buf_idx]
            cbuf = cpg.send_coord_bufs[buf_idx]
            for d in 1:ndime
                cbuf[coffset + d] = alya_local_coords[i, d]
            end
            coord_offsets[buf_idx] += ndime
        end
    end

    return nothing
end


#------------------------------------------------------------------------------------
#= 
    unpack_received_data!(cpg::CouplingData, u::AbstractVector, mesh, 
                         alya_local_coords::Matrix{Float64}, alya_local_ids::Vector{Int32})

Unpack received data from Alya and apply to Jexpresso solution.
This is a placeholder - implement based on your physics coupling requirements.

# Arguments
- `cpg`: CouplingData with filled recv_bufs
- `u`: Jexpresso solution vector to modify
- `mesh`: Jexpresso mesh
- `alya_local_coords`: Local Alya coordinates [n_local × ndime]
- `alya_local_ids`: Global Alya IDs for local points (1-based)
=#
#------------------------------------------------------------------------------------
function unpack_received_data!(cpg::CouplingData, u::AbstractVector, mesh,
                              alya_local_coords::Matrix{Float64}, 
                              alya_local_ids::Vector{Int32})
    # This function depends on your specific coupling physics
    # Example: apply received Alya data as forcing or boundary conditions
    
    neqs = cpg.neqs
    wsize = length(cpg.npoin_recv)
    
    # Reconstruct full received data organized by Alya point ID
    n_local = size(alya_local_coords, 1)
    received_values = zeros(Float64, n_local, neqs)
    
    # Unpack from receive buffers
    recv_offsets = zeros(Int, wsize)
    
    for src_rank in cpg.recv_from_ranks
        buf_idx = src_rank + 1
        if cpg.npoin_recv[buf_idx] > 0
            buf = cpg.recv_bufs[buf_idx]
            n_from_src = cpg.npoin_recv[buf_idx]
            
            # Find which local points came from this rank
            for i in 1:n_local
                # Check if this point should come from this rank
                # (This is simplified - you may need a more sophisticated mapping)
                if recv_offsets[buf_idx] < length(buf)
                    offset = recv_offsets[buf_idx]
                    for q in 1:neqs
                        received_values[i, q] = buf[offset + q]
                    end
                    recv_offsets[buf_idx] += neqs
                end
            end
        end
    end
    
    # TODO: Apply received data to Jexpresso solution
    # This depends on your coupling approach (forcing, BC, etc.)
    # Example:
    # for i in 1:n_local
    #     px, py = alya_local_coords[i, 1], alya_local_coords[i, 2]
    #     # Find nearest Jexpresso mesh point
    #     distances = sqrt.((mesh.x .- px).^2 .+ (mesh.y .- py).^2)
    #     nearest_idx = argmin(distances)
    #     # Apply coupling (example: weak forcing)
    #     u[nearest_idx] += 0.01 * received_values[i, 1]
    # end
    
    return nothing
end

#------------------------------------------------------------------------------------
#=
     view_state_matrix(u::AbstractVector, npoin::Int, neqs::Int)

Reshape the flat solution vector into a matrix view [npoin × neqs].

# Arguments
- `u`: Flat solution vector of length npoin*neqs
- `npoin`: Number of spatial points
- `neqs`: Number of equations per point

# Returns
- Matrix view of u with shape [npoin × neqs]
=#
#------------------------------------------------------------------------------------
function view_state_matrix(u::AbstractVector, npoin::Int, neqs::Int)
    @assert length(u) == npoin * neqs "Solution vector size mismatch"
    return reshape(view(u, :), neqs, npoin)'
end

#------------------------------------------------------------------------------------
#    interpolate_solution_to_alya_coords(alya_coords::Matrix{Float64}, mesh, 
#                                       u_mat::AbstractMatrix, basis, neqs;
#                                       use_bins::Bool=true, bins_per_dim::Int=64)
#
#Interpolate SEM solution to Alya grid coordinates.
#
## Arguments
#- `alya_coords`: [n_points × ndime] coordinates where to interpolate
#- `mesh`: Jexpresso SEM mesh
#- `u_mat`: [npoin × neqs] solution matrix
#- `basis`: SEM basis functions
#- `neqs`: Number of equations
#- `use_bins`: Whether to use spatial binning for faster element search
#- `bins_per_dim`: Number of bins per dimension if use_bins=true
#
## Returns
#- `u_interp`: [n_points × neqs] interpolated values
#
#------------------------------------------------------------------------------------
function interpolate_solution_to_alya_coords(alya_coords::Matrix{Float64}, mesh, 
                                             u_mat::AbstractMatrix, basis, ξ, neqs, inputs;
                                             use_bins::Bool=true, bins_per_dim::Int=64)

    
    n_points = size(alya_coords, 1)
    u_interp = zeros(Float64, n_points, neqs)
    
    # Extract mesh info
    get_conn = _make_conn_accessor(mesh)
    nelem = _num_elems(mesh)
    ngl = mesh.ngl

    # Get LGL nodes from basis (where solution is defined)
    ξ_nodes = Vector{Float64}(ξ)
    ω = barycentric_weights(ξ_nodes)

    # Build element bounding boxes
    elem_bboxes = Vector{NTuple{4,Float64}}(undef, nelem)
    @inbounds for e in 1:nelem
        nodes = get_conn(e)
        xs = mesh.x[nodes]
        ys = mesh.y[nodes]
        elem_bboxes[e] = (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
    end

    bins = use_bins ? _build_elem_bins(elem_bboxes; bins_per_dim=bins_per_dim) : nothing
    
    # Interpolate each point
    @inbounds for ipt in 1:n_points
        px, py = alya_coords[ipt, 1], alya_coords[ipt, 2]
        
        candidates = bins !== nothing ? _bin_candidates(bins, px, py, elem_bboxes) : (1:nelem)
        
        found = false
        for e in candidates
            nodes = get_conn(e)
            x_elem = mesh.x[nodes]
            y_elem = mesh.y[nodes]
            
            # Quick bbox check
            if px < minimum(x_elem) - 1e-10 || px > maximum(x_elem) + 1e-10 ||
               py < minimum(y_elem) - 1e-10 || py > maximum(y_elem) + 1e-10
                continue
            end
            
            # Map to reference coordinates
            ξ_ref, η_ref, converged = physical_to_reference(
                px, py, x_elem, y_elem, ξ_nodes, ω, ngl
            )
            
            if !converged || abs(ξ_ref) > 1.0 + 1e-10 || abs(η_ref) > 1.0 + 1e-10
                continue
            end
            
            # Evaluate basis at (ξ_ref, η_ref)
            ψξ = evaluate_lagrange_1d(ξ_ref, ξ_nodes, ω)
            ψη = evaluate_lagrange_1d(η_ref, ξ_nodes, ω)
            
            # Interpolate solution
            for q in 1:neqs
                val = 0.0
                idx = 1
                for j in 1:ngl, i in 1:ngl
                    val += ψξ[i] * ψη[j] * u_mat[nodes[idx], q]
                    idx += 1
                end
                u_interp[ipt, q] = val
            end
            
            found = true
            break
        end
        
        if !found
            # Fallback: nearest neighbor
            distances = (mesh.x .- px).^2 .+ (mesh.y .- py).^2
            nearest = argmin(distances)
            for q in 1:neqs
                u_interp[ipt, q] = u_mat[nearest, q]
            end
        end
    end

    return u_interp
end

# Evaluate 1D Lagrange basis at arbitrary point using barycentric formula
function evaluate_lagrange_1d(ξ::Float64, ξ_nodes::Vector{Float64}, ω::Vector{Float64})
    n = length(ξ_nodes)
    ψ = zeros(Float64, n)
    
    # Check for exact node match
    @inbounds for i in 1:n
        if abs(ξ - ξ_nodes[i]) < 1e-14
            ψ[i] = 1.0
            return ψ
        end
    end
    
    # Barycentric formula
    sum_val = 0.0
    @inbounds for i in 1:n
        ψ[i] = ω[i] / (ξ - ξ_nodes[i])
        sum_val += ψ[i]
    end
    
    inv_sum = 1.0 / sum_val
    @inbounds for i in 1:n
        ψ[i] *= inv_sum
    end
    
    return ψ
end

# Map physical to reference coordinates using Newton iteration
function physical_to_reference(px::Float64, py::Float64,
                               x_elem::AbstractVector, y_elem::AbstractVector,
                               ξ_nodes::Vector{Float64}, ω::Vector{Float64}, ngl::Int)
    max_iter = 20
    tol = 1e-12
    ξ, η = 0.0, 0.0
    
    for iter in 1:max_iter
        # Evaluate basis and derivatives
        ψξ = evaluate_lagrange_1d(ξ, ξ_nodes, ω)
        ψη = evaluate_lagrange_1d(η, ξ_nodes, ω)
        dψξ = evaluate_lagrange_1d_derivative(ξ, ξ_nodes, ω)
        dψη = evaluate_lagrange_1d_derivative(η, ξ_nodes, ω)
        
        # Compute current position and Jacobian
        x_curr = y_curr = dxdξ = dxdη = dydξ = dydη = 0.0
        
        idx = 1
        @inbounds for j in 1:ngl, i in 1:ngl
            ψ_val = ψξ[i] * ψη[j]
            dξ_val = dψξ[i] * ψη[j]
            dη_val = ψξ[i] * dψη[j]
            
            x_curr += ψ_val * x_elem[idx]
            y_curr += ψ_val * y_elem[idx]
            dxdξ += dξ_val * x_elem[idx]
            dxdη += dη_val * x_elem[idx]
            dydξ += dξ_val * y_elem[idx]
            dydη += dη_val * y_elem[idx]
            idx += 1
        end
        
        # Residual
        rx = px - x_curr
        ry = py - y_curr
        
        if sqrt(rx*rx + ry*ry) < tol
            return ξ, η, true
        end
        
        # Newton update
        det_J = dxdξ * dydη - dxdη * dydξ
        abs(det_J) < 1e-15 && return ξ, η, false
        
        inv_det = 1.0 / det_J
        ξ += inv_det * ( dydη * rx - dxdη * ry)
        η += inv_det * (-dydξ * rx + dxdξ * ry)
        
        (abs(ξ) > 10.0 || abs(η) > 10.0) && return ξ, η, false
    end
    
    return ξ, η, false
end

# Evaluate 1D Lagrange derivative using barycentric formula
function evaluate_lagrange_1d_derivative(ξ::Float64, ξ_nodes::Vector{Float64}, 
                                        ω::Vector{Float64})
    n = length(ξ_nodes)
    dψ = zeros(Float64, n)
    
    # Check for exact node match
    @inbounds for j in 1:n
        if abs(ξ - ξ_nodes[j]) < 1e-14
            # Derivative at node j
            for k in 1:n
                k == j && continue
                dψ[k] = ω[k] / (ω[j] * (ξ_nodes[j] - ξ_nodes[k]))
            end
            for k in 1:n
                k == j && continue
                dψ[j] -= 1.0 / (ξ_nodes[j] - ξ_nodes[k])
            end
            return dψ
        end
    end
    
    # General case
    α = zeros(Float64, n)
    sum_α = 0.0
    @inbounds for i in 1:n
        α[i] = ω[i] / (ξ - ξ_nodes[i])
        sum_α += α[i]
    end
    
    sum_dα = 0.0
    @inbounds for i in 1:n
        sum_dα -= α[i] / (ξ - ξ_nodes[i])
    end
    
    inv_sum2 = 1.0 / (sum_α * sum_α)
    @inbounds for i in 1:n
        dα_i = -α[i] / (ξ - ξ_nodes[i])
        dψ[i] = (dα_i * sum_α - α[i] * sum_dα) * inv_sum2
    end
    
    return dψ
end

#------------------------------------------------------------------------------------
#    perform_coupling_exchange!(integrator, cpg::CouplingData, basis)
#
#Complete coupling exchange: interpolate Jexpresso solution to Alya grid,
#exchange with Alya, and apply received data.
#
#This is the main function to be called from the coupling callback.
#------------------------------------------------------------------------------------
function je_perform_coupling_exchange(u, u_mat, t, cpg::CouplingData, mesh, basis, inputs, ξ, neqs)
    
    # 1. Prepare solution view
    npoin = mesh.npoin
    neqs  = neqs
    
    qout = zeros(Float64, npoin, neqs)
    u2uaux!(u_mat, u, neqs, npoin)
    call_user_uout(qout, u_mat, u_mat, 0, inputs[:SOL_VARS_TYPE], npoin, neqs, neqs)
    
    # 2. Interpolate to local Alya coordinates
    u_interp_local = interpolate_solution_to_alya_coords(
        cpg.alya_local_coords, mesh, qout, basis, 
        ξ, neqs, inputs;
        use_bins=true, bins_per_dim=64
    )

    # 3. Pack interpolated data AND coordinates into send buffers
    pack_interpolated_data!(cpg, u_interp_local, cpg.alya_owner_ranks, cpg.alya_local_coords)

    # 4. Exchange data with Alya
    coupling_exchange_data!(cpg)

    # 5. Unpack and apply received data from Alya
    unpack_received_data!(cpg, u, mesh,
                          cpg.alya_local_coords, cpg.alya_local_ids)

    return nothing
end

#-------------------------------------------------------------------------------------
#    verify_coupling_communication_pattern(npoin_recv, npoin_send, 
#                                         alya_owner_ranks, alya_local_ids,
#                                        coupling_data, local_comm, world_comm)
#
#Verify that the coupling communication pattern is correct by checking:
#1. npoin_recv matches the distribution of alya_owner_ranks
#2. Global point conservation (no duplicates/missing points)
#3. Communication symmetry after Alltoall
#4. Detailed breakdown of send/receive patterns
#-------------------------------------------------------------------------------------
function verify_coupling_communication_pattern(npoin_recv, npoin_send,
                                               alya_owner_ranks, alya_local_ids,
                                               coupling_data, local_comm, world_comm)
    
    lrank = MPI.Comm_rank(local_comm)
    wrank = MPI.Comm_rank(world_comm)
    wsize = MPI.Comm_size(world_comm)
    lsize = MPI.Comm_size(local_comm)
    
    nranks_alya = length(coupling_data[:alya2world])
    
    println("="^80)
    println("[VERIFY] Julia lrank=$lrank, wrank=$wrank")
    println("="^80)
    
    # -------------------------------------------------------------------------
    # CHECK 1: Verify npoin_recv matches alya_owner_ranks distribution
    # -------------------------------------------------------------------------
    println("\n[CHECK 1] Verifying npoin_recv against alya_owner_ranks...")
    
    # Count how many points this rank owns for each Alya world rank
    local_counts = zeros(Int32, wsize)
    for owner_wrank in alya_owner_ranks
        local_counts[owner_wrank + 1] += 1
    end
    
    # Compare with npoin_recv
    mismatch = false
    for i in 1:wsize
        if local_counts[i] != npoin_recv[i]
            println("  ❌ MISMATCH at world rank $(i-1):")
            println("     alya_owner_ranks count: $(local_counts[i])")
            println("     npoin_recv value: $(npoin_recv[i])")
            mismatch = true
        end
    end
    
    if !mismatch
        println("  ✓ npoin_recv matches alya_owner_ranks distribution")
    else
        @warn "npoin_recv does NOT match alya_owner_ranks!"
    end
    
    # -------------------------------------------------------------------------
    # CHECK 2: Verify total points to receive
    # -------------------------------------------------------------------------
    println("\n[CHECK 2] Verifying total points to receive...")
    
    total_local_points = length(alya_owner_ranks)
    total_recv = sum(npoin_recv)
    
    println("  Local Alya points extracted: $total_local_points")
    println("  Sum of npoin_recv: $total_recv")
    
    if total_local_points == total_recv
        println("  ✓ Total counts match")
    else
        @warn "Total counts do NOT match!"
    end
    
    # -------------------------------------------------------------------------
    # CHECK 3: Global point conservation (USE LOCAL_COMM ONLY)
    # -------------------------------------------------------------------------
    println("\n[CHECK 3] Checking global point conservation...")
    
    # Total Alya grid points
    rem_nx = coupling_data[:rem_nx]
    total_alya_points = rem_nx[1] * rem_nx[2] * rem_nx[3]
    
    # Sum across all Julia ranks ONLY (use local_comm, not world_comm)
    global_owned = MPI.Allreduce(total_local_points, MPI.SUM, local_comm)
    
    println("  Total Alya grid points: $total_alya_points")
    println("  Total points owned by all Julia ranks: $global_owned")
    
    if global_owned == total_alya_points
        println("  ✓ All Alya points are accounted for (no duplicates/missing)")
    else
        @warn "Point conservation failed! Expected $total_alya_points, got $global_owned"
    end
    
    # -------------------------------------------------------------------------
    # CHECK 4: Detailed send/receive breakdown by rank
    # -------------------------------------------------------------------------
    println("\n[CHECK 4] Detailed communication pattern...")
    
    println("  RECEIVE (npoin_recv) - I will receive FROM:")
    total = 0
    for i in 1:wsize
        if npoin_recv[i] > 0
            wrank_from = i - 1
            if wrank_from < nranks_alya
                println("    $(npoin_recv[i]) points FROM world rank $wrank_from (Alya rank $wrank_from)")
            else
                jrank = wrank_from - nranks_alya
                println("    $(npoin_recv[i]) points FROM world rank $wrank_from (Julia rank $jrank)")
            end
            total += npoin_recv[i]
        end
    end
    println("  Total to RECEIVE: $total")
    
    println("\n  SEND (npoin_send) - I will send TO:")
    total = 0
    for i in 1:wsize
        if npoin_send[i] > 0
            wrank_to = i - 1
            if wrank_to < nranks_alya
                println("    $(npoin_send[i]) points TO world rank $wrank_to (Alya rank $wrank_to)")
            else
                jrank = wrank_to - nranks_alya
                println("    $(npoin_send[i]) points TO world rank $wrank_to (Julia rank $jrank)")
            end
            total += npoin_send[i]
        end
    end
    println("  Total to SEND: $total")
    
    # -------------------------------------------------------------------------
    # CHECK 5: Communication symmetry
    # -------------------------------------------------------------------------
    println("\n[CHECK 5] Checking communication symmetry...")
    
    recv_partners = [i-1 for i in 1:wsize if npoin_recv[i] > 0]
    send_partners = [i-1 for i in 1:wsize if npoin_send[i] > 0]
    
    println("  Ranks I receive from: $recv_partners")
    println("  Ranks I send to: $send_partners")
    
    # -------------------------------------------------------------------------
    # CHECK 6: Verify specific Alya point IDs
    # -------------------------------------------------------------------------
    println("\n[CHECK 6] Sample Alya point ownership...")
    
     n_sample = min(10, length(alya_local_ids))
    #n_sample = length(alya_local_ids)
    if n_sample > 0
        println("  First $n_sample local Alya points:")
        for i in 1:n_sample
            gid = alya_local_ids[i]
            owner_wrank = alya_owner_ranks[i]
            if owner_wrank < nranks_alya
                println("    Point ID $gid → owned by world rank $owner_wrank (Alya rank $owner_wrank)")
            else
                jrank = owner_wrank - nranks_alya
                println("    Point ID $gid → owned by world rank $owner_wrank (Julia rank $jrank)")
            end
        end
    end
    
    # -------------------------------------------------------------------------
    # CHECK 7: Julia-only communication matrix (USE LOCAL_COMM ONLY)
    # -------------------------------------------------------------------------
    println("\n[CHECK 7] Julia-only communication pattern check...")
    
    # Build communication statistics for Julia ranks only
    julia_send_to_alya = zeros(Int32, nranks_alya)
    julia_recv_from_alya = zeros(Int32, nranks_alya)
    
    for i in 1:nranks_alya
        julia_send_to_alya[i] = npoin_send[i]
        julia_recv_from_alya[i] = npoin_recv[i]
    end
    
    # Sum across Julia ranks using local_comm
    global_send_to_alya = MPI.Allreduce(julia_send_to_alya, MPI.SUM, local_comm)
    global_recv_from_alya = MPI.Allreduce(julia_recv_from_alya, MPI.SUM, local_comm)
    
    if lrank == 0
        println("\n  Global Julia → Alya communication:")
        for i in 1:nranks_alya
            if global_send_to_alya[i] > 0 || global_recv_from_alya[i] > 0
                println("    Alya rank $(i-1): Julia sends $(global_send_to_alya[i]), receives $(global_recv_from_alya[i])")
            end
        end
    end
    
    println("\n" * "="^80)
    println("[VERIFY] Verification complete for Julia lrank=$lrank")
    println("="^80)
    flush(stdout)
    
    return true
end

"""
    Alya_Point_Ownership

Structure to store which Jexpresso rank owns each Alya grid point.
"""
struct Alya_Point_Ownership
    alya_coords::Matrix{Float64}       # Alya point coordinates [npoints × ndime]
    owner_jrank::Vector{Int32}         # Jexpresso local rank that owns each point (-1 if not owned)
    owner_wrank::Vector{Int32}         # World rank that owns each point (-1 if not owned)
    alya_point_id::Vector{Int32}       # Original Alya point ID (1-based)
    ndime::Int32                       # Spatial dimension
end


"""
    build_alya_point_ownership_map(mesh, coupling_data, local_comm, world_comm)

Build a complete map of which Jexpresso ranks own which Alya grid points.
This is useful for debugging and verification.

# Arguments
- `mesh`: Jexpresso mesh structure with local coordinates
- `coupling_data`: Dict containing remote grid metadata from Alya
- `local_comm`: MPI communicator for Jexpresso ranks only
- `world_comm`: MPI.COMM_WORLD (shared with Alya)

# Returns
- `ownership`: Alya_Point_Ownership structure with complete mapping
"""
function build_alya_point_ownership_map(mesh, coupling_data, local_comm, world_comm)
    
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
    local_coords = zeros(Float64, nmax, ndime)
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
    ownership = Alya_Point_Ownership(
        local_coords,
        global_owner_jrank,
        global_owner_wrank,
        alya_point_ids,
        Int32(ndime)
    )
    
    return ownership
end


"""
    print_alya_point_ownership(ownership, coupling_data; 
                               max_points_to_print=20,
                               only_owned=false)

Print which Jexpresso rank owns each Alya grid point.

# Arguments
- `ownership`: Alya_Point_Ownership structure
- `coupling_data`: Dict containing coupling metadata
- `max_points_to_print`: Maximum number of points to print (default: 20, use -1 for all)
- `only_owned`: If true, only print points that are owned by some Jexpresso rank

# Example output:
```
Alya Point Ownership Map:
========================
Point    1: (  -5000.0,     0.0) -> Jexpresso rank 1 (world rank 3)
Point    2: (  -3888.9,     0.0) -> Jexpresso rank 1 (world rank 3)
Point   11: (  -5000.0,  1111.1) -> Jexpresso rank 1 (world rank 3)
...
```
"""
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


"""
    verify_alya_point_distribution(ownership, coupling_data)

Print statistics about how Alya points are distributed across Jexpresso ranks.
"""
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

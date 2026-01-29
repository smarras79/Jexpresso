struct AngularElementGhostInfo
    # Element identification
    spatial_elem_id::Int
    owner_rank::Int
    
    # Angular mesh structure (still need full angular mesh for interpolation)
    nelem_ang::Int
    nop_ang::Vector{Int}
    npoin_ang::Int
    
    # Angular coordinates (full angular mesh)
    coords_ang::Matrix{Float64}
    connijk_ang::Array{Int, 3}
    
    # Angular metrics (full angular mesh - needed for interpolation weights)
    Je_ang::Array{Float64, 3}
    dξdx_ang::Array{Float64, 3}
    dξdy_ang::Array{Float64, 3}
    dηdx_ang::Array{Float64, 3}
    dηdy_ang::Array{Float64, 3}
    
    # Refinement information
    ref_level::Vector{Int}
    
    # NEW: Only store interface nodes, not entire element
    interface_spatial_nodes::Vector{Int}  # Local spatial node IDs on interface
    interface_spatial_global_ids::Vector{Int}  # Global spatial node IDs
    
    # Map: (local_spatial_idx, e_ext, iθ, jθ) -> global spatial-angular ID
    interface_node_gids::Dict{Tuple{Int,Int,Int,Int}, Int}
end

struct NonConformingGhostLayer
    ghost_elements::Dict{Tuple{Int,Int}, AngularElementGhostInfo}
    owned_to_ghost_map::Dict{Int, Vector{AngularElementGhostInfo}}
    send_angular_data_to::Dict{Int, Vector{Int}}
    recv_angular_data_from::Dict{Int, Vector{Int}}
    interface_hanging_nodes::Set{Int}
    parent_search_cache::Dict{Int, AngularElementGhostInfo}
    n_ghost_nodes::Int
    ghost_node_offset::Int
end

function build_nonconforming_ghost_layer_corrected(
    mesh, connijk_spa, ip2gip, ip2gip_spa, gip2owner_spa,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
    extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
    extra_meshes_ref_level,
    n_spa, neighbors
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)
    
    nelem = mesh.nelem  # Only locally owned elements
    ngl = mesh.ngl
    
    @info "[Rank $rank] Building non-conforming ghost layer (corrected)..."
    
    # =========================================================================
    # PHASE 1: Identify partition boundary nodes
    # =========================================================================
    # These are nodes that are shared with other processors
    # (nodes on faces/edges/vertices of the partition boundary)
    
    partition_boundary_nodes = Set{Int}()  # Local node IDs
    node_to_owner = Dict{Int, Int}()        # Local node ID -> owning rank
    
    for ip = 1:mesh.npoin
        owner = mesh.gip2owner[ip]
        node_to_owner[ip] = owner
        
        # A node is on partition boundary if another processor owns it
        # OR if we own it but it appears on other processors
        if owner != rank
            # We have this node but don't own it -> definitely on boundary
            push!(partition_boundary_nodes, ip)
        end
    end
    
    # Also find nodes we own that appear on other processors
    # This requires communication
    owned_boundary_nodes = find_owned_nodes_on_boundary(mesh, ip2gip, rank, comm)
    union!(partition_boundary_nodes, owned_boundary_nodes)
    
    @info "[Rank $rank] Found $(length(partition_boundary_nodes)) nodes on partition boundary"

    # =========================================================================
    # PHASE 2: Identify locally owned elements that touch partition boundary
    # =========================================================================
    # These are the elements that might have ghost neighbors
    
    boundary_elements = Set{Int}()
    elem_boundary_nodes = Dict{Int, Set{Int}}()  # elem -> boundary nodes it contains
    
    for iel = 1:nelem  # All elements in mesh are locally owned
        boundary_nodes_in_elem = Set{Int}()
        
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            
            if ip in partition_boundary_nodes
                push!(boundary_nodes_in_elem, ip)
            end
        end
        
        if !isempty(boundary_nodes_in_elem)
            push!(boundary_elements, iel)
            elem_boundary_nodes[iel] = boundary_nodes_in_elem
        end
    end
    
    @info "[Rank $rank] $(length(boundary_elements)) locally owned elements touch partition boundary"
    
    # =========================================================================
    # PHASE 3: Exchange ghost neighbor information
    # =========================================================================
    
    ghost_neighbor_map = exchange_ghost_neighbor_info(
        boundary_elements, elem_boundary_nodes, ip2gip, rank, comm
    )
    
    @info "[Rank $rank] Found $(sum(length(v) for v in values(ghost_neighbor_map))) ghost neighbor relationships"
    
    # =========================================================================
    # PHASE 3.5: Exchange interface node requests (NEW)
    # =========================================================================
    
    interface_requests = exchange_interface_node_requests(
        ghost_neighbor_map, elem_boundary_nodes, ip2gip, rank, comm
    )
    
    interface_to_send = exchange_interface_requests(interface_requests, rank, comm)
    
    @info "[Rank $rank] Received interface requests from $(length(interface_to_send)) ranks"
    
    # Log interface sizes for verification
    for (dest_rank, elem_dict) in interface_to_send
        for (elem_id, node_set) in elem_dict
            @info "[Rank $rank] Will send $(length(node_set)) interface nodes for element $elem_id to rank $dest_rank"
            # Should be 1, 5, or 25 for 4th order elements
        end
    end
    
    # =========================================================================
    # PHASE 4: Request angular data for ghost elements
    # =========================================================================
    
    needed_ghost_elements = Set{Tuple{Int, Int}}()
    
    for neighbors in values(ghost_neighbor_map)
        for (ghost_iel, ghost_owner) in neighbors
            push!(needed_ghost_elements, (ghost_iel, ghost_owner))
        end
    end
    
    requests_by_rank = Dict{Int, Vector{Int}}()
    
    for (ghost_iel, ghost_owner) in needed_ghost_elements
        if !haskey(requests_by_rank, ghost_owner)
            requests_by_rank[ghost_owner] = Int[]
        end
        push!(requests_by_rank[ghost_owner], ghost_iel)
    end
    
    @info "[Rank $rank] Requesting angular data from $(length(requests_by_rank)) ranks"
    
    # =========================================================================
    # PHASE 5: Exchange element requests
    # =========================================================================
    
    send_elem_data_to = exchange_element_requests(requests_by_rank, rank, comm)
    
    # =========================================================================
    # PHASE 6: Pack and send angular element data (UPDATED)
    # =========================================================================
    
    send_buffers = Dict{Int, Vector{UInt8}}()
    
    for (dest_rank, elem_ids) in send_elem_data_to
        # Get interface nodes for this destination
        interface_nodes = if haskey(interface_to_send, dest_rank)
            interface_to_send[dest_rank]
        else
            @warn "[Rank $rank] No interface info for rank $dest_rank, sending empty"
            Dict{Int, Set{Int}}()
        end
        
        buffer = pack_angular_element_data_with_interface(
            elem_ids, ngl,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_extra_Je, extra_meshes_extra_dξdx, 
            extra_meshes_extra_dξdy, extra_meshes_extra_dηdx, 
            extra_meshes_extra_dηdy, extra_meshes_ref_level,
            connijk_spa, ip2gip_spa,
            mesh, ip2gip,
            interface_nodes
        )
        send_buffers[dest_rank] = buffer
    end

    @info "[Rank $rank] has packed data for sending"
    # =========================================================================
    # PHASE 7: Receive angular element data
    # =========================================================================
    @info "[Rank $rank] is exchanging buffers"
    recv_buffers = exchange_buffers(send_buffers, requests_by_rank, comm)
    @info "[Rank $rank] finished exchanging buffers"
    # =========================================================================
    # PHASE 8: Unpack received data into ghost element structures
    # =========================================================================
    @info "[Rank $rank] is unpacking data after receiving"
    ghost_elements = Dict{Tuple{Int,Int}, AngularElementGhostInfo}()  # (ghost_iel, owner) -> info
    ghost_node_counter = 0
    
    for (src_rank, buffer) in recv_buffers
        elem_infos = unpack_angular_element_data_optimized(buffer, src_rank)
        
        for elem_info in elem_infos
            key = (elem_info.spatial_elem_id, elem_info.owner_rank)
            ghost_elements[key] = elem_info
            ghost_node_counter += count_ghost_nodes_in_element(elem_info, ngl)
        end
    end
    
    @info "[Rank $rank] Received $(length(ghost_elements)) ghost elements with ~$ghost_node_counter ghost nodes"

    # =========================================================================
    # PHASE 9: Map owned elements to their ghost neighbors
    # =========================================================================
    
    owned_to_ghost_map = build_owned_to_ghost_map(
        ghost_neighbor_map, ghost_elements
    )

    # =========================================================================
    # PHASE 10: Identify hanging nodes on interfaces
    # =========================================================================
    
    interface_hanging_nodes, parent_search_cache = identify_interface_hanging_nodes_final(
        owned_to_ghost_map, ghost_elements,
        mesh, connijk_spa,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops, extra_meshes_extra_nelems,
        ip2gip_spa, ngl, rank
    )
    
    @info "[Rank $rank] Found $(length(interface_hanging_nodes)) hanging nodes on processor interfaces"
    
    ghost_layer = NonConformingGhostLayer(
        ghost_elements, 
        owned_to_ghost_map,
        send_elem_data_to,
        requests_by_rank,
        interface_hanging_nodes,
        parent_search_cache,
        ghost_node_counter,
        n_spa
    )
    
    return ghost_layer
end

# =========================================================================
# Helper: Find owned nodes that appear on other processors
# =========================================================================

function find_owned_nodes_on_boundary(mesh, ip2gip, rank, comm)
    """
    Find nodes that:
    1. This processor owns (gip2owner[ip] == rank)
    2. Also appear on other processors (are in their mesh.connijk)
    
    Strategy:
    - Gather all owned nodes from this processor
    - Gather all non-owned nodes from all other processors
    - Intersection gives us owned nodes on partition boundary
    """
    
    nproc = MPI.Comm_size(comm)
    
    # Collect global IDs of nodes we own
    owned_nodes_global = Set{Int}()
    # Collect global IDs of nodes we have but don't own
    not_owned_nodes_global = Set{Int}()
    for ip = 1:mesh.npoin
        gip = ip2gip[ip]
        owner = mesh.gip2owner[ip]
        
        if owner == rank
            push!(owned_nodes_global, gip)
        else
            push!(not_owned_nodes_global, gip)
        end
    end

    # Exchange lists
    local_owned_list = collect(owned_nodes_global)
    local_not_owned_list = collect(not_owned_nodes_global)
    
    n_owned = Int32(length(local_owned_list))
    n_not_owned = Int32(length(local_not_owned_list))
    
    # Gather counts
    owned_counts = MPI.Allgather([n_owned], comm)
    not_owned_counts = MPI.Allgather([n_not_owned], comm)
    
    # Gather actual data
    if isempty(local_owned_list)
        local_owned_list = Int[]
    end
    if isempty(local_not_owned_list)
        local_not_owned_list = Int[]
    end
    
    total_owned = sum(owned_counts)
    total_not_owned = sum(not_owned_counts)
    
    all_owned = if total_owned > 0
        MPI.Allgatherv(local_owned_list, owned_counts, comm)
    else
        Int[]
    end
    
    all_not_owned = if total_not_owned > 0
        MPI.Allgatherv(local_not_owned_list, not_owned_counts, comm)
    else
        Int[]
    end
    
    # Build sets of not-owned nodes for each processor
    not_owned_by_proc = Dict{Int, Set{Int}}()
    offset = 0
    for r = 0:(nproc-1)
        count = not_owned_counts[r+1]
        not_owned_by_proc[r] = Set(all_not_owned[offset+1:offset+count])
        offset += count
    end
    
    # Find intersection: nodes we own that appear as not-owned on other processors
    boundary_owned_global = Set{Int}()
    
    for other_rank = 0:(nproc-1)
        if other_rank == rank
            continue
        end
        
        # Intersection of our owned nodes with other rank's not-owned nodes
        shared = intersect(owned_nodes_global, not_owned_by_proc[other_rank])
        union!(boundary_owned_global, shared)
    end
    
    @info "[Rank $rank] Found $(length(boundary_owned_global)) owned nodes on partition boundary"
    
    # Convert to local IDs
    boundary_owned_local = Set{Int}()
    for ip = 1:mesh.npoin
        gip = ip2gip[ip]
        if gip in boundary_owned_global
            push!(boundary_owned_local, ip)
        end
    end
    
    return boundary_owned_local
end

# =========================================================================
# Helper: Exchange ghost neighbor information
# =========================================================================

function exchange_ghost_neighbor_info(
    boundary_elements, elem_boundary_nodes, ip2gip, rank, comm
)
    """
    For each boundary element, find which elements on other processors
    share its boundary nodes.
    
    Returns: Dict{Int, Vector{Tuple{Int,Int}}}
             owned_elem_id -> [(ghost_elem_id, ghost_owner_rank), ...]
    """
    
    nproc = MPI.Comm_size(comm)
    
    # Package boundary node information to send
    # For each boundary element, send its global boundary node IDs
    local_boundary_info = Dict{Int, Vector{Int}}()  # elem_id -> [global node IDs]
    
    for iel in boundary_elements
        global_nodes = [ip2gip[ip] for ip in elem_boundary_nodes[iel]]
        local_boundary_info[iel] = global_nodes
    end
    
    # Serialize and exchange
    local_buffer = IOBuffer()
    serialize(local_buffer, local_boundary_info)
    local_data = take!(local_buffer)
    
    n_local = Int32(length(local_data))
    counts = MPI.Allgather([n_local], comm)
    total_count = sum(counts)
    
    all_boundary_info = Dict{Int, Dict{Int, Vector{Int}}}()  # rank -> boundary_info
    
    if total_count > 0
        all_data = MPI.Allgatherv(local_data, counts, comm)
        
        offset = 1
        for src_rank = 0:(nproc-1)
            count = counts[src_rank+1]
            if count > 0
                chunk = all_data[offset:offset+count-1]
                all_boundary_info[src_rank] = deserialize(IOBuffer(chunk))
            else
                all_boundary_info[src_rank] = Dict{Int, Vector{Int}}()
            end
            offset += count
        end
    else
        for src_rank = 0:(nproc-1)
            all_boundary_info[src_rank] = Dict{Int, Vector{Int}}()
        end
    end
    
    # Now check which remote elements share nodes with our boundary elements
    ghost_neighbor_map = Dict{Int, Vector{Tuple{Int,Int}}}()
    
    for iel_owned in boundary_elements
        owned_global_nodes = Set([ip2gip[ip] for ip in elem_boundary_nodes[iel_owned]])
        ghost_neighbors = Tuple{Int,Int}[]
        
        # Check elements from other processors
        for (remote_rank, remote_boundary_info) in all_boundary_info
            if remote_rank == rank
                continue  # Skip self
            end
            
            for (remote_iel, remote_global_nodes) in remote_boundary_info
                # Check if they share nodes
                shared = intersect(owned_global_nodes, Set(remote_global_nodes))
                
                if !isempty(shared)
                    # These elements are neighbors across partition boundary
                    push!(ghost_neighbors, (remote_iel, remote_rank))
                end
            end
        end
        
        if !isempty(ghost_neighbors)
            ghost_neighbor_map[iel_owned] = ghost_neighbors
        end
    end
    
    return ghost_neighbor_map
end

# =========================================================================
# Helper: Build owned-to-ghost mapping
# =========================================================================

function build_owned_to_ghost_map(
    ghost_neighbor_map, ghost_elements
)
    """
    Create a more convenient mapping structure
    """
    
    owned_to_ghost = Dict{Int, Vector{AngularElementGhostInfo}}()
    
    for (iel_owned, ghost_list) in ghost_neighbor_map
        ghost_infos = AngularElementGhostInfo[]
        
        for (ghost_iel, ghost_owner) in ghost_list
            key = (ghost_iel, ghost_owner)
            if haskey(ghost_elements, key)
                push!(ghost_infos, ghost_elements[key])
            else
                @warn "Ghost element ($ghost_iel, $ghost_owner) not found in received data"
            end
        end
        
        if !isempty(ghost_infos)
            owned_to_ghost[iel_owned] = ghost_infos
        end
    end
    
    return owned_to_ghost
end

# =========================================================================
# Helper: Identify interface hanging nodes (final version)
# =========================================================================

function identify_interface_hanging_nodes_final(
    owned_to_ghost_map, ghost_elements,
    mesh, connijk_spa,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    ip2gip_spa, ngl, rank
)
    """
    For each owned element with ghost neighbors:
    1. Compare angular mesh with each ghost neighbor
    2. Identify child-parent relationships
    3. Mark hanging nodes on the interface
    """
    
    interface_hanging = Set{Int}()
    parent_cache = Dict{Int, AngularElementGhostInfo}()
    
    for (iel_owned, ghost_infos) in owned_to_ghost_map
        
        # Check each angular element in the owned element
        for e_ext_owned = 1:extra_meshes_extra_nelems[iel_owned]
            θmin_owned, θmax_owned, ϕmin_owned, ϕmax_owned = get_element_bounds_fast(
                iel_owned, e_ext_owned, extra_meshes_coords,
                extra_meshes_connijk, extra_meshes_extra_nops
            )
            
            # Check against each ghost neighbor
            for ghost_info in ghost_infos
                for e_ext_ghost = 1:ghost_info.nelem_ang
                    θmin_ghost, θmax_ghost, ϕmin_ghost, ϕmax_ghost = 
                        get_ghost_element_bounds(ghost_info, e_ext_ghost)
                    
                    # Check if owned angular element is a child of ghost angular element
                    if is_child_element(θmin_owned, θmax_owned, ϕmin_owned, ϕmax_owned,
                                       θmin_ghost, θmax_ghost, ϕmin_ghost, ϕmax_ghost)
                        
                        # Mark hanging nodes in owned element
                        mark_hanging_nodes_in_element_advanced!(
                            interface_hanging, parent_cache,
                            iel_owned, e_ext_owned,
                            ghost_info, e_ext_ghost,
                            connijk_spa, mesh.connijk, ngl,
                            extra_meshes_coords, extra_meshes_connijk,
                            extra_meshes_extra_nops
                        )
                    end
                end
            end
        end
    end
    
    return interface_hanging, parent_cache
end

function mark_hanging_nodes_in_element_advanced!(
    interface_hanging, parent_cache,
    iel_owned, e_ext_owned,
    ghost_info, e_ext_ghost,
    connijk_spa, spatial_connijk, ngl,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops
)
    """
    Mark all non-vertex nodes in the child (owned) element as hanging
    using coordinate-based checking to determine parent vertices
    """
    
    nop_owned = extra_meshes_extra_nops[iel_owned][e_ext_owned]
    
    # Mark all nodes in this angular element that are not parent vertices
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        for jθ = 1:(nop_owned+1), iθ = 1:(nop_owned+1)
            ip_spa = connijk_spa[iel_owned][i, j, k, e_ext_owned, iθ, jθ]
            
            # Use advanced coordinate-based checking
            is_hanging = !is_parent_vertex_advanced(
                iθ, jθ, nop_owned,
                iel_owned, e_ext_owned,
                ghost_info, e_ext_ghost,
                extra_meshes_coords, extra_meshes_connijk
            )
            
            if is_hanging
                push!(interface_hanging, ip_spa)
                parent_cache[ip_spa] = ghost_info
            end
        end
    end
end

# =========================================================================
# Helper: Pack angular element data for sending
# =========================================================================

function pack_angular_element_data_with_interface(
    elem_ids, ngl,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je, extra_meshes_extra_dξdx, 
    extra_meshes_extra_dξdy, extra_meshes_extra_dηdx, 
    extra_meshes_extra_dηdy, extra_meshes_ref_level,
    connijk_spa, ip2gip_spa,
    mesh, ip2gip,
    interface_nodes_to_send  # Dict{elem_id -> Set(global_spatial_ids)}
)
    """
    Pack only the interface nodes specified by the requester
    """
    
    buffer = IOBuffer()
    
    # Write number of elements
    write(buffer, Int32(length(elem_ids)))
    
    for iel in elem_ids
        # Element ID
        write(buffer, Int32(iel))
        write(buffer, Int32(ngl))
        
        # Angular mesh structure (full - needed for interpolation)
        write(buffer, Int32(extra_meshes_extra_nelems[iel]))
        write(buffer, Int32(length(extra_meshes_extra_nops[iel])))
        
        for nop in extra_meshes_extra_nops[iel]
            write(buffer, Int32(nop))
        end
        
        # Angular coordinates (full)
        npoin_ang = size(extra_meshes_coords[iel], 2)
        write(buffer, Int32(npoin_ang))
        write(buffer, extra_meshes_coords[iel][:])
        
        # Angular connectivity (full)
        nelem_ang = extra_meshes_extra_nelems[iel]
        for e_ext = 1:nelem_ang
            nop = extra_meshes_extra_nops[iel][e_ext]
            for j = 1:(nop+1)
                for i = 1:(nop+1)
                    write(buffer, Int32(extra_meshes_connijk[iel][e_ext, i, j]))
                end
            end
        end
        
        # Angular metrics (full)
        for e_ext = 1:nelem_ang
            nop = extra_meshes_extra_nops[iel][e_ext]
            for j = 1:(nop+1)
                for i = 1:(nop+1)
                    write(buffer, Float64(extra_meshes_extra_Je[iel][e_ext, i, j]))
                    write(buffer, Float64(extra_meshes_extra_dξdx[iel][e_ext, i, j]))
                    write(buffer, Float64(extra_meshes_extra_dξdy[iel][e_ext, i, j]))
                    write(buffer, Float64(extra_meshes_extra_dηdx[iel][e_ext, i, j]))
                    write(buffer, Float64(extra_meshes_extra_dηdy[iel][e_ext, i, j]))
                end
            end
        end
        
        # Refinement levels
        for level in extra_meshes_ref_level[iel]
            write(buffer, Int32(level))
        end
        
        # ===================================================================
        # Find interface nodes based on provided list
        # ===================================================================
        
        requested_global_ids = interface_nodes_to_send[iel]
        interface_spatial_indices = Tuple{Int,Int,Int}[]
        
        # Find which (i,j,k) in this element correspond to requested global IDs
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip_local = mesh.connijk[iel, i, j, k]
            gip = ip2gip[ip_local]
            
            if gip in requested_global_ids
                push!(interface_spatial_indices, (i, j, k))
            end
        end
        
        # Sort for consistent ordering
        sort!(interface_spatial_indices)
        
        n_interface = length(interface_spatial_indices)
        write(buffer, Int32(n_interface))
        
        @debug "Packing element $iel: $n_interface interface nodes"
        
        # Write interface node data
        for (i, j, k) in interface_spatial_indices
            # Write spatial indices
            write(buffer, Int32(i))
            write(buffer, Int32(j))
            write(buffer, Int32(k))
            
            # Write global spatial ID
            ip_spatial = mesh.connijk[iel, i, j, k]
            gip_spatial = ip2gip[ip_spatial]
            write(buffer, Int64(gip_spatial))
            
            # Write spatial-angular global IDs for this spatial node
            for e_ext = 1:nelem_ang
                nop = extra_meshes_extra_nops[iel][e_ext]
                for jθ = 1:(nop+1)
                    for iθ = 1:(nop+1)
                        ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        gid = ip2gip_spa[ip_spa]
                        write(buffer, Int64(gid))
                    end
                end
            end
        end
    end
    
    return take!(buffer)
end

# =========================================================================
# Updated Unpack Function - Only Unpack Interface Nodes
# =========================================================================

function unpack_angular_element_data_optimized(buffer, owner_rank)
    """
    Unpack ghost element data with only interface nodes
    """
    
    io = IOBuffer(buffer)
    n_elems = read(io, Int32)
    
    elem_infos = AngularElementGhostInfo[]
    
    for elem_idx = 1:n_elems
        elem_id = read(io, Int32)
        ngl = read(io, Int32)
        
        # Read angular mesh structure
        nelem_ang = read(io, Int32)
        n_nops = read(io, Int32)
        nop_ang = [read(io, Int32) for _ = 1:n_nops]
        
        npoin_ang = read(io, Int32)
        coords_data = [read(io, Float64) for _ = 1:(2*npoin_ang)]
        coords_ang = reshape(coords_data, 2, npoin_ang)
        
        # Read angular connectivity
        max_nop = maximum(nop_ang)
        connijk_ang = zeros(Int, nelem_ang, max_nop+1, max_nop+1)
        
        for e_ext = 1:nelem_ang
            nop = nop_ang[e_ext]
            for j = 1:(nop+1)
                for i = 1:(nop+1)
                    connijk_ang[e_ext, i, j] = read(io, Int32)
                end
            end
        end
        
        # Read angular metrics
        Je_ang = zeros(Float64, nelem_ang, max_nop+1, max_nop+1)
        dξdx_ang = zeros(Float64, nelem_ang, max_nop+1, max_nop+1)
        dξdy_ang = zeros(Float64, nelem_ang, max_nop+1, max_nop+1)
        dηdx_ang = zeros(Float64, nelem_ang, max_nop+1, max_nop+1)
        dηdy_ang = zeros(Float64, nelem_ang, max_nop+1, max_nop+1)
        
        for e_ext = 1:nelem_ang
            nop = nop_ang[e_ext]
            for j = 1:(nop+1)
                for i = 1:(nop+1)
                    Je_ang[e_ext, i, j] = read(io, Float64)
                    dξdx_ang[e_ext, i, j] = read(io, Float64)
                    dξdy_ang[e_ext, i, j] = read(io, Float64)
                    dηdx_ang[e_ext, i, j] = read(io, Float64)
                    dηdy_ang[e_ext, i, j] = read(io, Float64)
                end
            end
        end
        
        # Read refinement levels
        ref_level = [read(io, Int32) for _ = 1:nelem_ang]
        
        # ===================================================================
        # NEW: Read only interface nodes
        # ===================================================================
        
        n_interface = read(io, Int32)
        
        interface_spatial_nodes = Int[]
        interface_spatial_global_ids = Int[]
        interface_node_gids = Dict{Tuple{Int,Int,Int,Int}, Int}()
        
        for interface_idx = 1:n_interface
            # Read spatial indices
            i = read(io, Int32)
            j = read(io, Int32)
            k = read(io, Int32)
            
            # Read global spatial ID
            gip_spatial = read(io, Int64)
            
            push!(interface_spatial_nodes, interface_idx)  # Store sequential index
            push!(interface_spatial_global_ids, gip_spatial)
            
            # Read spatial-angular global IDs
            for e_ext = 1:nelem_ang
                nop = nop_ang[e_ext]
                for jθ = 1:(nop+1)
                    for iθ = 1:(nop+1)
                        gid = read(io, Int64)
                        # Key: (interface_spatial_idx, e_ext, iθ, jθ)
                        interface_node_gids[(interface_idx, e_ext, iθ, jθ)] = gid
                    end
                end
            end
        end
        
        elem_info = AngularElementGhostInfo(
            elem_id,
            owner_rank,
            nelem_ang,
            nop_ang,
            npoin_ang,
            coords_ang,
            connijk_ang,
            Je_ang,
            dξdx_ang,
            dξdy_ang,
            dηdx_ang,
            dηdy_ang,
            ref_level,
            interface_spatial_nodes,
            interface_spatial_global_ids,
            interface_node_gids
        )
        
        push!(elem_infos, elem_info)
    end
    
    return elem_infos
end

function get_face_spatial_indices(face_id, ngl)
    """Get (i,j,k) indices for nodes on a specific face"""
    if face_id == 1
        return [(1, j, k) for j = 1:ngl, k = 1:ngl]
    elseif face_id == 2
        return [(ngl, j, k) for j = 1:ngl, k = 1:ngl]
    elseif face_id == 3
        return [(i, 1, k) for i = 1:ngl, k = 1:ngl]
    elseif face_id == 4
        return [(i, ngl, k) for i = 1:ngl, k = 1:ngl]
    elseif face_id == 5
        return [(i, j, 1) for i = 1:ngl, j = 1:ngl]
    else  # face_id == 6
        return [(i, j, ngl) for i = 1:ngl, j = 1:ngl]
    end
end

# =========================================================================
# MPI Communication helpers
# =========================================================================

function exchange_element_requests(requests_by_rank, rank, comm)
    """Exchange which elements each rank needs"""
    # Similar to previous MPI exchange patterns
    # Returns: Dict(rank -> [elem_ids that rank needs from me])
    
    # Serialize requests
    local_buffer = IOBuffer()
    serialize(local_buffer, requests_by_rank)
    local_data = take!(local_buffer)
    
    n_local = Int32(length(local_data))
    counts = MPI.Allgather([n_local], comm)
    
    total_count = sum(counts)
    if total_count > 0
        all_data = MPI.Allgatherv(local_data, counts, comm)
        
        # Deserialize all requests
        all_requests = Dict{Int, Dict{Int, Vector{Int}}}()
        offset = 1
        for src_rank = 0:(MPI.Comm_size(comm)-1)
            count = counts[src_rank+1]
            if count > 0
                chunk = all_data[offset:offset+count-1]
                all_requests[src_rank] = deserialize(IOBuffer(chunk))
            else
                all_requests[src_rank] = Dict{Int, Vector{Int}}()
            end
            offset += count
        end
    else
        all_requests = Dict{Int, Dict{Int, Vector{Int}}}()
    end
    
    # Build send list
    send_to = Dict{Int, Vector{Int}}()
    for (requesting_rank, their_requests) in all_requests
        if requesting_rank == rank
            continue
        end
        
        if haskey(their_requests, rank)
            send_to[requesting_rank] = their_requests[rank]
        end
    end
    
    return send_to
end

function exchange_buffers(send_buffers, requests_by_rank, comm)
    """
    Exchange variable-length data buffers between processors
    
    Protocol:
    1. Exchange buffer sizes (so receivers know how much data to expect)
    2. Exchange actual buffer data
    
    Uses non-blocking communication to avoid deadlocks
    """
    
    nproc = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    
    recv_buffers = Dict{Int, Vector{UInt8}}()
    
    # =========================================================================
    # PHASE 1: Exchange buffer sizes
    # =========================================================================
    
    # Prepare send size data
    send_sizes = Dict{Int, Int32}()
    for (dest_rank, buffer) in send_buffers
        send_sizes[dest_rank] = Int32(length(buffer))
    end
    
    # Post non-blocking receives for sizes
    size_requests = MPI.Request[]
    recv_sizes = Dict{Int, Ref{Int32}}()
    
    for (src_rank, _) in requests_by_rank
        size_buf = Ref{Int32}(0)
        recv_sizes[src_rank] = size_buf
        
        req = MPI.Irecv!(size_buf, comm, source=src_rank, tag=src_rank*1000)
        push!(size_requests, req)
    end
    
    # Post non-blocking sends for sizes
    for (dest_rank, size) in send_sizes
        MPI.Isend(Ref(size), comm, dest=dest_rank, tag=rank*1000)
    end
    
    # Wait for all size exchanges to complete
    MPI.Waitall(size_requests)
    
    @debug "[Rank $rank] Phase 1 complete: exchanged sizes with $(length(recv_sizes)) ranks"
    
    # =========================================================================
    # PHASE 2: Allocate receive buffers based on received sizes
    # =========================================================================
    
    for (src_rank, size_ref) in recv_sizes
        size = size_ref[]
        if size > 0
            recv_buffers[src_rank] = Vector{UInt8}(undef, size)
        else
            recv_buffers[src_rank] = Vector{UInt8}()
        end
    end
    
    @debug "[Rank $rank] Phase 2 complete: allocated receive buffers"
    
    # =========================================================================
    # PHASE 3: Exchange actual buffer data
    # =========================================================================
    
    data_requests = MPI.Request[]
    
    # Post non-blocking receives for data
    for (src_rank, buffer) in recv_buffers
        if !isempty(buffer)
            req = MPI.Irecv!(buffer, comm, source=src_rank, tag=src_rank*2000)
            push!(data_requests, req)
        end
    end
    
    # Post non-blocking sends for data
    for (dest_rank, buffer) in send_buffers
        if !isempty(buffer)
            MPI.Isend(buffer, comm, dest=dest_rank, tag=rank*2000)
        end
    end
    
    # Wait for all data exchanges to complete
    MPI.Waitall(data_requests)
    
    @debug "[Rank $rank] Phase 3 complete: exchanged data buffers"
    
    return recv_buffers
end

# =========================================================================
# Helper: Get angular element bounds from ghost info
# =========================================================================

function get_ghost_element_bounds(ghost_info::AngularElementGhostInfo, e_ext)
    """
    Extract angular element bounds from ghost element info
    Returns: (θmin, θmax, ϕmin, ϕmax)
    """
    
    nop = ghost_info.nop_ang[e_ext]
    
    # Get corner coordinates
    ip_11 = ghost_info.connijk_ang[e_ext, 1, 1]
    ip_nn = ghost_info.connijk_ang[e_ext, nop+1, nop+1]
    
    θ_11 = ghost_info.coords_ang[1, ip_11]
    ϕ_11 = ghost_info.coords_ang[2, ip_11]
    
    θ_nn = ghost_info.coords_ang[1, ip_nn]
    ϕ_nn = ghost_info.coords_ang[2, ip_nn]
    
    # Handle periodic boundary in ϕ
    if ϕ_nn == 0.0 && ϕ_11 > π
        ϕ_nn = 2π
    end
    
    θmin, θmax = minmax(θ_11, θ_nn)
    ϕmin, ϕmax = minmax(ϕ_11, ϕ_nn)
    
    return θmin, θmax, ϕmin, ϕmax
end

# =========================================================================
# Updated Count Function
# =========================================================================

function count_ghost_nodes_in_element(elem_info::AngularElementGhostInfo, ngl)
    """
    Count only interface nodes, not entire element
    """
    
    n_interface_spatial = length(elem_info.interface_spatial_nodes)
    
    # Each interface spatial node has all angular nodes
    total_angular_nodes_per_spatial = sum((nop+1)^2 for nop in elem_info.nop_ang)
    
    return n_interface_spatial * total_angular_nodes_per_spatial
end

# =========================================================================
# Helper: Check if angular node is a parent vertex
# =========================================================================

function is_parent_vertex_advanced(
    iθ_child, jθ_child, nop_child,
    iel_child, e_ext_child,
    ghost_info::AngularElementGhostInfo, e_ext_parent,
    coords_ang_child, connijk_ang_child
)
    """
    Advanced version that checks actual coordinates to determine
    if a child node aligns with a parent vertex
    
    This handles arbitrary refinement patterns correctly by comparing
    actual (θ, φ) coordinates
    """
    
    nop_parent = ghost_info.nop_ang[e_ext_parent]
    
    # Get child node coordinates
    ip_child = connijk_ang_child[iel_child][e_ext_child, iθ_child, jθ_child]
    θ_child = coords_ang_child[iel_child][1, ip_child]
    ϕ_child = coords_ang_child[iel_child][2, ip_child]
    
    # Tolerance for coordinate comparison
    tol = 1e-12
    
    # Check against all parent vertices
    for jθ_parent = 1:(nop_parent+1), iθ_parent = 1:(nop_parent+1)
        ip_parent = ghost_info.connijk_ang[e_ext_parent, iθ_parent, jθ_parent]
        θ_parent = ghost_info.coords_ang[1, ip_parent]
        ϕ_parent = ghost_info.coords_ang[2, ip_parent]
        
        # Direct coordinate match
        if abs(θ_child - θ_parent) < tol && abs(ϕ_child - ϕ_parent) < tol
            return true
        end
        
        # Handle periodic boundary in ϕ (at 0/2π)
        if abs(θ_child - θ_parent) < tol
            # Check if coordinates match with 2π offset
            if abs(ϕ_child - (ϕ_parent + 2π)) < tol || 
               abs(ϕ_child - (ϕ_parent - 2π)) < tol
                return true
            end
            
            # Also check for ϕ=0 vs ϕ=2π equivalence
            if (abs(ϕ_child) < tol && abs(ϕ_parent - 2π) < tol) ||
               (abs(ϕ_child - 2π) < tol && abs(ϕ_parent) < tol)
                return true
            end
        end
    end
    
    return false
end

# Helper to get ghost element by tuple key
function get_ghost_element(ghost_layer::NonConformingGhostLayer, elem_id::Int, owner_rank::Int)
    return ghost_layer.ghost_elements[(elem_id, owner_rank)]
end

# Helper to check if ghost element exists
function has_ghost_element(ghost_layer::NonConformingGhostLayer, elem_id::Int, owner_rank::Int)
    return haskey(ghost_layer.ghost_elements, (elem_id, owner_rank))
end

# Iterator over all ghost elements
function iterate_ghost_elements(ghost_layer::NonConformingGhostLayer)
    return values(ghost_layer.ghost_elements)
end

# =========================================================================
# Helper: Identify interface spatial nodes
# =========================================================================

function identify_interface_spatial_nodes(
    iel_owned, iel_ghost,
    mesh, ip2gip, ngl,
    shared_nodes_set::Set{Int}
)
    """
    Find which spatial nodes in the ghost element are on the interface
    with the owned element.
    
    Returns:
    - interface_nodes_local: Local spatial indices in ghost element
    - interface_nodes_global: Global spatial node IDs
    - spatial_index_map: Maps ghost element local index -> interface index
    """
    
    interface_nodes_local = Int[]
    interface_nodes_global = Int[]
    spatial_index_map = Dict{Tuple{Int,Int,Int}, Int}()  # (i,j,k) -> interface_idx
    
    interface_idx = 1
    
    # Check all nodes in ghost element
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip_ghost = mesh.connijk[iel_ghost, i, j, k]
        gip_ghost = ip2gip[ip_ghost]
        
        # Is this node on the interface (shared with owned element)?
        if gip_ghost in shared_nodes_set
            push!(interface_nodes_local, ip_ghost)
            push!(interface_nodes_global, gip_ghost)
            spatial_index_map[(i, j, k)] = interface_idx
            interface_idx += 1
        end
    end
    
    return interface_nodes_local, interface_nodes_global, spatial_index_map
end

# =========================================================================
# Helper: Get ghost node global ID (updated)
# =========================================================================

function get_ghost_interface_node_gid(
    elem_info::AngularElementGhostInfo,
    interface_spatial_idx::Int,
    e_ext::Int,
    iθ::Int,
    jθ::Int
)
    """
    Get global ID for a ghost interface node
    """
    key = (interface_spatial_idx, e_ext, iθ, jθ)
    
    if !haskey(elem_info.interface_node_gids, key)
        error("Interface node not found: interface_idx=$interface_spatial_idx, e_ext=$e_ext, iθ=$iθ, jθ=$jθ")
    end
    
    return elem_info.interface_node_gids[key]
end

# =========================================================================
# Helper: Find interface spatial index from global ID
# =========================================================================

function find_interface_spatial_index(
    elem_info::AngularElementGhostInfo,
    gip_spatial::Int
)
    """
    Find which interface index corresponds to a global spatial node ID
    Returns: interface index or nothing
    """
    
    for (idx, gid) in enumerate(elem_info.interface_spatial_global_ids)
        if gid == gip_spatial
            return idx
        end
    end
    
    return nothing
end

# =========================================================================
# Helper: Exchange interface node information 
# =========================================================================

function exchange_interface_node_requests(
    ghost_neighbor_map, elem_boundary_nodes, ip2gip, rank, comm
)
    """
    For each ghost element we need, tell the owner which spatial nodes
    we actually need (the interface nodes).
    
    Returns: Dict{Int, Dict{Int, Set{Int}}}
             dest_rank -> elem_id -> Set(global_spatial_node_ids)
    """
    
    # Build requests: for each ghost element, which nodes do we need?
    interface_requests = Dict{Int, Dict{Int, Set{Int}}}()
    
    for (iel_owned, ghost_list) in ghost_neighbor_map
        # Get global IDs of boundary nodes in owned element
        owned_boundary_global = Set(ip2gip[ip] for ip in elem_boundary_nodes[iel_owned])
        
        for (ghost_iel, ghost_owner) in ghost_list
            # We need the nodes from ghost_iel that overlap with iel_owned's boundary
            if !haskey(interface_requests, ghost_owner)
                interface_requests[ghost_owner] = Dict{Int, Set{Int}}()
            end
            
            if !haskey(interface_requests[ghost_owner], ghost_iel)
                interface_requests[ghost_owner][ghost_iel] = Set{Int}()
            end
            
            # The interface nodes are the boundary nodes of the owned element
            # (these will be shared with the ghost element)
            union!(interface_requests[ghost_owner][ghost_iel], owned_boundary_global)
        end
    end
    
    return interface_requests
end

function exchange_interface_requests(interface_requests, rank, comm)
    """
    Exchange interface node requests to determine what to send
    
    Returns: Dict{Int, Dict{Int, Set{Int}}}
             src_rank -> elem_id -> Set(global_spatial_node_ids to send)
    """
    
    nproc = MPI.Comm_size(comm)
    
    # Serialize requests
    local_buffer = IOBuffer()
    serialize(local_buffer, interface_requests)
    local_data = take!(local_buffer)
    
    n_local = Int32(length(local_data))
    counts = MPI.Allgather([n_local], comm)
    total_count = sum(counts)
    
    all_requests = Dict{Int, Dict{Int, Dict{Int, Set{Int}}}}()
    
    if total_count > 0
        all_data = MPI.Allgatherv(local_data, counts, comm)
        
        offset = 1
        for src_rank = 0:(nproc-1)
            count = counts[src_rank+1]
            if count > 0
                chunk = all_data[offset:offset+count-1]
                all_requests[src_rank] = deserialize(IOBuffer(chunk))
            else
                all_requests[src_rank] = Dict{Int, Dict{Int, Set{Int}}}()
            end
            offset += count
        end
    else
        for src_rank = 0:(nproc-1)
            all_requests[src_rank] = Dict{Int, Dict{Int, Set{Int}}}()
        end
    end
    
    # Extract what other processors need from me
    interface_to_send = Dict{Int, Dict{Int, Set{Int}}}()
    
    for (requesting_rank, their_requests) in all_requests
        if requesting_rank == rank
            continue
        end
        
        if haskey(their_requests, rank)
            interface_to_send[requesting_rank] = their_requests[rank]
        end
    end
    
    return interface_to_send
end
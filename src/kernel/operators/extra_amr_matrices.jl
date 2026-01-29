# =========================================================================
# Build Extended Local Numbering with Ghost Nodes
# =========================================================================

function build_extended_local_numbering(
    n_spa, ghost_layer::NonConformingGhostLayer,
    ip2gip_spa, rank
)
    """
    Build mapping from global IDs to extended local indices (including ghosts)
    
    Local indices: 1:n_spa (owned nodes)
    Ghost indices: (n_spa+1):(n_spa+n_ghost) (ghost nodes)
    
    Returns:
    - gid_to_extended_local: Dict{Int, Int} mapping global ID -> extended local index
    - extended_local_to_gid: Vector{Int} inverse mapping
    - n_total: Total number of nodes (owned + ghost)
    """
    
    @info "[Rank $rank] Building extended local numbering..."
    
    gid_to_extended_local = Dict{Int, Int}()
    extended_local_to_gid = Int[]
    
    # Reserve space
    sizehint!(gid_to_extended_local, n_spa + ghost_layer.n_ghost_nodes)
    sizehint!(extended_local_to_gid, n_spa + ghost_layer.n_ghost_nodes)
    
    # =========================================================================
    # Phase 1: Map owned nodes (1:n_spa)
    # =========================================================================
    
    for ip_local = 1:n_spa
        gid = ip2gip_spa[ip_local]
        gid_to_extended_local[gid] = ip_local
        push!(extended_local_to_gid, gid)
    end
    
    @info "[Rank $rank] Mapped $n_spa owned nodes"
    
    # =========================================================================
    # Phase 2: Map ghost nodes (n_spa+1:n_total)
    # =========================================================================
    
    ghost_idx = n_spa + 1
    ghost_gids_added = Set{Int}()  # Track to avoid duplicates
    
    for ((ghost_iel, ghost_owner), ghost_info) in ghost_layer.ghost_elements
        # Iterate through all interface nodes in this ghost element
        for (key, gid) in ghost_info.interface_node_gids
            # key is (interface_spatial_idx, e_ext, iθ, jθ)
            
            # Only add if we haven't seen this global ID yet
            if !haskey(gid_to_extended_local, gid) && !(gid in ghost_gids_added)
                gid_to_extended_local[gid] = ghost_idx
                push!(extended_local_to_gid, gid)
                push!(ghost_gids_added, gid)
                ghost_idx += 1
            end
        end
    end
    
    n_ghost_mapped = ghost_idx - n_spa - 1
    n_total = ghost_idx - 1
    
    @info "[Rank $rank] Mapped $n_ghost_mapped ghost nodes"
    @info "[Rank $rank] Total extended system size: $n_total"
    
    # Verify count
    if n_ghost_mapped != ghost_layer.n_ghost_nodes
        @warn "[Rank $rank] Ghost node count mismatch: expected $(ghost_layer.n_ghost_nodes), got $n_ghost_mapped"
    end
    
    return gid_to_extended_local, extended_local_to_gid, n_total
end

# =========================================================================
# Updated: Build Restriction Matrix with Extended Mapping
# =========================================================================

function build_restriction_matrix_with_interface_hanging(
    connijk_spa, nc_non_global_nodes, n_spa, n_total,
    ghost_layer::NonConformingGhostLayer,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je,
    mesh, ngl, nelem,
    neighbors,  # Add neighbors
    ip2gip_spa, gid_to_extended_local,
    rank
)
    """
    Build restriction matrix with neighbors infrastructure
    """
    
    n_free = n_spa - length(nc_non_global_nodes)
    hanging_node_set = Set(nc_non_global_nodes)
    
    @info "[Rank $rank] Building restriction matrix:"
    @info "  Free DOFs: $n_free"
    @info "  Hanging DOFs: $(length(nc_non_global_nodes))"
    @info "  Interface hanging: $(length(ghost_layer.interface_hanging_nodes))"
    @info "  Total columns (including ghosts): $n_total"
    
    # Triplet storage
    I_vec = Int[]
    J_vec = Int[]
    V_vec = Float64[]
    
    estimated_entries = n_free + length(nc_non_global_nodes) * 20
    sizehint!(I_vec, estimated_entries)
    sizehint!(J_vec, estimated_entries)
    sizehint!(V_vec, estimated_entries)
    
    # =========================================================================
    # Phase 1: Identity for free nodes
    # =========================================================================
    
    for ip_g = 1:n_free
        push!(I_vec, ip_g)
        push!(J_vec, ip_g)
        push!(V_vec, 1.0)
    end
    
    @info "[Rank $rank] Added identity for $n_free free nodes"
    
    # =========================================================================
    # Phase 2: Constraints for hanging nodes
    # =========================================================================
    
    # Build interpolation cache
    interpolation_cache = Dict{NTuple{4,Int}, Tuple{Matrix{Float64}, Matrix{Float64}}}()
    
    # Separate interior and interface hanging nodes
    interior_hanging = setdiff(hanging_node_set, ghost_layer.interface_hanging_nodes)
    interface_hanging = ghost_layer.interface_hanging_nodes
    
    @info "[Rank $rank] Interior hanging: $(length(interior_hanging)), Interface hanging: $(length(interface_hanging))"
    
    # Process interior hanging nodes (parents are local)
    n_interior_constraints = 0
    for ip_hanging in interior_hanging
        if ip_hanging > n_free
            # Find parent element and build constraint
            constraint_entries = build_interior_hanging_constraint(
                ip_hanging, connijk_spa, extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_extra_Je,
                mesh, ngl, nelem, neighbors, interpolation_cache  # Pass neighbors
            )
            
            for (parent_idx, weight) in constraint_entries
                if parent_idx <= n_free  # Parent must be a free node
                    push!(I_vec, parent_idx)
                    push!(J_vec, ip_hanging)
                    push!(V_vec, weight)
                    n_interior_constraints += 1
                end
            end
        end
    end
    
    @info "[Rank $rank] Added $n_interior_constraints interior hanging node constraints"
    
    # Process interface hanging nodes (parents may be ghost)
    n_interface_constraints = 0
    n_interface_to_ghost = 0
    n_interface_to_free = 0
    
    for ip_hanging in interface_hanging
        if ip_hanging > n_free
            # Get parent ghost element info
            if !haskey(ghost_layer.parent_search_cache, ip_hanging)
                @warn "[Rank $rank] Hanging node $ip_hanging has no parent cache entry"
                continue
            end
            
            ghost_info = ghost_layer.parent_search_cache[ip_hanging]
            
            # Build constraint using ghost element data with extended mapping
            constraint_entries = build_interface_hanging_constraint(
                ip_hanging, ghost_info, connijk_spa,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                extra_meshes_extra_Je,
                mesh, ngl, nelem,
                ip2gip_spa, gid_to_extended_local,
                ghost_offset, n_free,
                interpolation_cache,
                rank
            )
            
            for (parent_extended_idx, weight) in constraint_entries
                if parent_extended_idx <= n_free
                    # Parent is a free node
                    push!(I_vec, parent_extended_idx)
                    push!(J_vec, ip_hanging)
                    push!(V_vec, weight)
                    n_interface_to_free += 1
                elseif parent_extended_idx <= n_total
                    # Parent is a ghost node
                    push!(I_vec, n_free)  # Placeholder - will need special handling
                    push!(J_vec, parent_extended_idx)
                    push!(V_vec, weight)
                    n_interface_to_ghost += 1
                else
                    @warn "[Rank $rank] Parent index $parent_extended_idx out of range [1, $n_total]"
                end
                
                n_interface_constraints += 1
            end
        end
    end
    
    @info "[Rank $rank] Added $n_interface_constraints interface hanging node constraints"
    @info "[Rank $rank]   -> $n_interface_to_free to free nodes"
    @info "[Rank $rank]   -> $n_interface_to_ghost to ghost nodes"
    
    # Build sparse matrix
    nc_mat = sparse(I_vec, J_vec, V_vec, n_free, n_total)
    
    # Normalize columns
    @info "[Rank $rank] Normalizing constraint matrix columns..."
    for j = 1:n_total
        col_sum = 0.0
        for idx in nzrange(nc_mat, j)
            col_sum += nonzeros(nc_mat)[idx]
        end
        
        if abs(col_sum - 1.0) > 1e-13 && abs(col_sum) > 1e-14
            for idx in nzrange(nc_mat, j)
                nonzeros(nc_mat)[idx] /= col_sum
            end
        end
    end
    
    @info "[Rank $rank] Restriction matrix complete: nnz=$(nnz(nc_mat))"
    
    return nc_mat
end

# =========================================================================
# Build Constraint for Interior Hanging Node
# =========================================================================

function build_interior_hanging_constraint(
    ip_hanging, connijk_spa, 
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je,
    mesh, ngl, nelem,
    neighbors,  # Add neighbors
    interpolation_cache
)
    """
    Build interpolation constraint for a hanging node whose parent element
    is locally owned (not a ghost).
    
    Updated to use neighbors infrastructure.
    """
    
    # Find which element and indices this hanging node belongs to
    iel_child, i_child, j_child, k_child, e_ext_child, iθ_child, jθ_child = 
        find_node_location(ip_hanging, connijk_spa, ngl, nelem)
    
    if iel_child === nothing
        @warn "Could not find location for hanging node $ip_hanging"
        return Tuple{Int, Float64}[]
    end
    
    # Find parent element using neighbors infrastructure
    iel_parent, e_ext_parent = find_parent_element_local(
        iel_child, e_ext_child,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops, extra_meshes_extra_nelems,
        neighbors, nelem
    )
    
    if iel_parent === nothing
        @warn "Could not find parent element for hanging node $ip_hanging in element $iel_child, angular element $e_ext_child"
        return Tuple{Int, Float64}[]
    end
    
    # Build or retrieve interpolation matrices
    cache_key = (iel_child, e_ext_child, iel_parent, e_ext_parent)
    
    if !haskey(interpolation_cache, cache_key)
        Lθ, Lϕ = build_interpolation_matrices_local!(
            iel_child, e_ext_child, iel_parent, e_ext_parent,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops
        )
        interpolation_cache[cache_key] = (Lθ, Lϕ)
    end
    
    Lθ, Lϕ = interpolation_cache[cache_key]
    
    # Compute interpolation weights with Jacobian scaling
    constraint_entries = Tuple{Int, Float64}[]
    
    Je_child = extra_meshes_extra_Je[iel_child][e_ext_child, iθ_child, jθ_child]
    nop_parent = extra_meshes_extra_nops[iel_parent][e_ext_parent]
    
    # Determine spatial location in parent element
    # Key point: if parent is in a different spatial element, we need to find
    # the corresponding spatial location
    i_parent, j_parent, k_parent = find_corresponding_spatial_location(
        iel_child, i_child, j_child, k_child,
        iel_parent, mesh, ngl
    )
    
    for jθ_parent = 1:(nop_parent+1)
        for iθ_parent = 1:(nop_parent+1)
            Lθ_val = Lθ[iθ_child, iθ_parent]
            Lϕ_val = Lϕ[jθ_child, jθ_parent]
            
            if abs(Lθ_val) < 1e-14 || abs(Lϕ_val) < 1e-14
                continue
            end
            
            # Get parent node index
            ip_parent = connijk_spa[iel_parent][i_parent, j_parent, k_parent, 
                                                 e_ext_parent, iθ_parent, jθ_parent]
            
            Je_parent = extra_meshes_extra_Je[iel_parent][e_ext_parent, 
                                                           iθ_parent, jθ_parent]
            
            # Weight includes interpolation and Jacobian ratio
            weight = Lθ_val * Lϕ_val * (Je_child / Je_parent)
            
            push!(constraint_entries, (ip_parent, weight))
        end
    end
    
    return constraint_entries
end

# =========================================================================
# Updated: Build Interface Hanging Constraint with Mapping
# =========================================================================

function build_interface_hanging_constraint(
    ip_hanging, ghost_info::AngularElementGhostInfo,
    connijk_spa,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je,
    mesh, ngl, nelem,
    ip2gip_spa, gid_to_extended_local,  # Added mappings
    ghost_offset, n_free,
    interpolation_cache,
    rank
)
    """
    Build interpolation constraint for a hanging node whose parent element
    is a ghost element on another processor.
    
    Now uses gid_to_extended_local to convert global IDs to extended local indices.
    
    Returns: Vector{Tuple{Int, Float64}} - list of (extended_local_index, weight) pairs
    """
    
    # Find which element and indices this hanging node belongs to
    iel_child, i_child, j_child, k_child, e_ext_child, iθ_child, jθ_child = 
        find_node_location(ip_hanging, connijk_spa, ngl, nelem)
    
    if iel_child === nothing
        @warn "[Rank $rank] Could not find location for interface hanging node $ip_hanging"
        return Tuple{Int, Float64}[]
    end
    
    # Find which angular element in ghost is the parent
    e_ext_parent = find_parent_angular_element_in_ghost(
        iel_child, e_ext_child, ghost_info,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops
    )
    
    if e_ext_parent === nothing
        @warn "[Rank $rank] Could not find parent angular element in ghost for hanging node $ip_hanging"
        return Tuple{Int, Float64}[]
    end
    
    # Build interpolation matrices (child to ghost parent)
    Lθ, Lϕ = build_ghost_interpolation_matrices(
        iel_child, e_ext_child, ghost_info, e_ext_parent,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops
    )
    
    # Get child node Jacobian
    Je_child = extra_meshes_extra_Je[iel_child][e_ext_child, iθ_child, jθ_child]
    
    # Build constraint entries
    constraint_entries = Tuple{Int, Float64}[]
    
    nop_parent = ghost_info.nop_ang[e_ext_parent]
    
    # Find which interface spatial index corresponds to (i_child, j_child, k_child)
    ip_spatial = mesh.connijk[iel_child, i_child, j_child, k_child]
    gip_spatial = mesh.ip2gip[ip_spatial]
    
    # Find this in ghost's interface nodes
    interface_spatial_idx = find_interface_spatial_index(ghost_info, gip_spatial)
    
    if interface_spatial_idx === nothing
        @warn "[Rank $rank] Spatial node (gip=$gip_spatial) not found in ghost interface for hanging node $ip_hanging"
        return Tuple{Int, Float64}[]
    end
    
    for jθ_parent = 1:(nop_parent+1)
        for iθ_parent = 1:(nop_parent+1)
            Lθ_val = Lθ[iθ_child, iθ_parent]
            Lϕ_val = Lϕ[jθ_child, jθ_parent]
            
            if abs(Lθ_val) < 1e-14 || abs(Lϕ_val) < 1e-14
                continue
            end
            
            # Get ghost parent node global ID
            gid_parent = get_ghost_interface_node_gid(
                ghost_info, interface_spatial_idx, e_ext_parent, 
                iθ_parent, jθ_parent
            )
            
            # Convert global ID to extended local index using mapping
            if !haskey(gid_to_extended_local, gid_parent)
                @warn "[Rank $rank] Ghost parent node gid=$gid_parent not in extended mapping"
                continue
            end
            
            extended_local_idx = gid_to_extended_local[gid_parent]
            
            # Get parent Jacobian
            Je_parent = ghost_info.Je_ang[e_ext_parent, iθ_parent, jθ_parent]
            
            # Compute weight with Jacobian scaling
            weight = Lθ_val * Lϕ_val * (Je_child / Je_parent)
            
            push!(constraint_entries, (extended_local_idx, weight))
        end
    end
    
    @debug "[Rank $rank] Built constraint for hanging node $ip_hanging: $(length(constraint_entries)) parent nodes"
    
    return constraint_entries
end

# =========================================================================
# Helper: Find node location in connectivity
# =========================================================================

function find_node_location(ip_spa, connijk_spa, ngl, nelem)
    """
    Find which element and indices a spatial-angular node belongs to
    
    Returns: (iel, i, j, k, e_ext, iθ, jθ) or (nothing, ...) if not found
    """
    
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            for e_ext = 1:size(connijk_spa[iel], 4)
                nop = size(connijk_spa[iel], 5) - 1
                for jθ = 1:(nop+1), iθ = 1:(nop+1)
                    if connijk_spa[iel][i, j, k, e_ext, iθ, jθ] == ip_spa
                        return iel, i, j, k, e_ext, iθ, jθ
                    end
                end
            end
        end
    end
    
    return nothing, 0, 0, 0, 0, 0, 0
end

# =========================================================================
# Helper: Find parent element using neighbors infrastructure
# =========================================================================

function find_parent_element_local(
    iel_child, e_ext_child,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    neighbors, nelem
)
    """
    Find the parent angular element using the neighbors infrastructure.
    
    Key insight: The hanging node's spatial element (iel_child) has the same
    angular mesh at all spatial points. We need to find which neighboring
    spatial element has a coarser angular element that contains e_ext_child.
    
    Returns: (iel_parent, e_ext_parent) or (nothing, nothing)
    """
    
    # Get child element bounds
    θmin_child, θmax_child, ϕmin_child, ϕmax_child = get_element_bounds_fast(
        iel_child, e_ext_child, extra_meshes_coords,
        extra_meshes_connijk, extra_meshes_extra_nops
    )
    
    # Check neighbors (up to 26 in 3D)
    for ineighbor = 1:26
        iel_neighbor = neighbors[iel_child, ineighbor, 1]
        is_nonconforming = neighbors[iel_child, ineighbor, 2]
        
        if iel_neighbor == 0
            # No neighbor in this direction
            continue
        end
        
        if is_nonconforming == 0
            # Conforming neighbor - no parent here
            continue
        end
        
        # This is a non-conforming neighbor - check its angular elements
        for e_ext_neighbor = 1:extra_meshes_extra_nelems[iel_neighbor]
            θmin_neighbor, θmax_neighbor, ϕmin_neighbor, ϕmax_neighbor = 
                get_element_bounds_fast(
                    iel_neighbor, e_ext_neighbor, extra_meshes_coords,
                    extra_meshes_connijk, extra_meshes_extra_nops
                )
            
            # Check if child element is contained in neighbor element
            if is_child_element(θmin_child, θmax_child, ϕmin_child, ϕmax_child,
                               θmin_neighbor, θmax_neighbor, ϕmin_neighbor, ϕmax_neighbor)
                # Found the parent!
                return iel_neighbor, e_ext_neighbor
            end
        end
    end
    
    # Parent not found
    return nothing, nothing
end

# =========================================================================
# Helper: Find parent angular element in ghost
# =========================================================================

function find_parent_angular_element_in_ghost(
    iel_child, e_ext_child, ghost_info::AngularElementGhostInfo,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops
)
    """
    Find which angular element in the ghost is the parent of the child element
    """
    
    # Get child element bounds
    θmin_child, θmax_child, ϕmin_child, ϕmax_child = get_element_bounds_fast(
        iel_child, e_ext_child, extra_meshes_coords,
        extra_meshes_connijk, extra_meshes_extra_nops
    )
    
    # Check each ghost angular element
    for e_ext_ghost = 1:ghost_info.nelem_ang
        θmin_ghost, θmax_ghost, ϕmin_ghost, ϕmax_ghost = 
            get_ghost_element_bounds(ghost_info, e_ext_ghost)
        
        if is_child_element(θmin_child, θmax_child, ϕmin_child, ϕmax_child,
                           θmin_ghost, θmax_ghost, ϕmin_ghost, ϕmax_ghost)
            return e_ext_ghost
        end
    end
    
    return nothing
end

# =========================================================================
# Helper: Build standard interpolation matrices (local to local)
# =========================================================================

function build_interpolation_matrices_local!(
    iel_child, e_ext_child, iel_parent, e_ext_parent,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops
)
    """
    Build interpolation matrices from child to parent (both local)
    """
    
    nop_child = extra_meshes_extra_nops[iel_child][e_ext_child]
    nop_parent = extra_meshes_extra_nops[iel_parent][e_ext_parent]
    
    # Get child grid points
    θ_child = zeros(Float64, nop_child+1)
    ϕ_child = zeros(Float64, nop_child+1)
    
    for i = 1:(nop_child+1)
        ip = extra_meshes_connijk[iel_child][e_ext_child, i, i]
        θ_child[i] = extra_meshes_coords[iel_child][1, ip]
        ϕ_child[i] = extra_meshes_coords[iel_child][2, ip]
    end
    
    # Get parent grid points
    θ_parent = zeros(Float64, nop_parent+1)
    ϕ_parent = zeros(Float64, nop_parent+1)
    
    for i = 1:(nop_parent+1)
        ip = extra_meshes_connijk[iel_parent][e_ext_parent, i, i]
        θ_parent[i] = extra_meshes_coords[iel_parent][1, ip]
        ϕ_parent[i] = extra_meshes_coords[iel_parent][2, ip]
    end
    
    # Build interpolation matrices
    Lθ = zeros(Float64, nop_child+1, nop_parent+1)
    Lϕ = zeros(Float64, nop_child+1, nop_parent+1)
    
    ωθ = zeros(Float64, nop_parent+1)
    ωϕ = zeros(Float64, nop_parent+1)
    
    BarycentricWeights!(θ_parent, ωθ)
    BarycentricWeights!(ϕ_parent, ωϕ)
    
    PolynomialInterpolationMatrix!(θ_parent, ωθ, θ_child, Lθ)
    PolynomialInterpolationMatrix!(ϕ_parent, ωϕ, ϕ_child, Lϕ)
    
    # Normalize rows
    for i = 1:(nop_child+1)
        row_sum_θ = sum(Lθ[i, :])
        if abs(row_sum_θ) > 1e-14
            Lθ[i, :] ./= row_sum_θ
        end
        
        row_sum_ϕ = sum(Lϕ[i, :])
        if abs(row_sum_ϕ) > 1e-14
            Lϕ[i, :] ./= row_sum_ϕ
        end
    end
    
    return Lθ, Lϕ
end

# =========================================================================
# Helper: Build ghost interpolation matrices
# =========================================================================

function build_ghost_interpolation_matrices(
    iel_child, e_ext_child, ghost_info::AngularElementGhostInfo, e_ext_parent,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops
)
    """
    Build interpolation matrices from child (local) to parent (ghost)
    """
    
    nop_child = extra_meshes_extra_nops[iel_child][e_ext_child]
    nop_parent = ghost_info.nop_ang[e_ext_parent]
    
    # Get child grid points
    θ_child = zeros(Float64, nop_child+1)
    ϕ_child = zeros(Float64, nop_child+1)
    
    for i = 1:(nop_child+1)
        ip = extra_meshes_connijk[iel_child][e_ext_child, i, i]
        θ_child[i] = extra_meshes_coords[iel_child][1, ip]
        ϕ_child[i] = extra_meshes_coords[iel_child][2, ip]
    end
    
    # Get parent (ghost) grid points
    θ_parent = zeros(Float64, nop_parent+1)
    ϕ_parent = zeros(Float64, nop_parent+1)
    
    for i = 1:(nop_parent+1)
        ip = ghost_info.connijk_ang[e_ext_parent, i, i]
        θ_parent[i] = ghost_info.coords_ang[1, ip]
        ϕ_parent[i] = ghost_info.coords_ang[2, ip]
    end
    
    # Build interpolation matrices
    Lθ = zeros(Float64, nop_child+1, nop_parent+1)
    Lϕ = zeros(Float64, nop_child+1, nop_parent+1)
    
    ωθ = zeros(Float64, nop_parent+1)
    ωϕ = zeros(Float64, nop_parent+1)
    
    BarycentricWeights!(θ_parent, ωθ)
    BarycentricWeights!(ϕ_parent, ωϕ)
    
    PolynomialInterpolationMatrix!(θ_parent, ωθ, θ_child, Lθ)
    PolynomialInterpolationMatrix!(ϕ_parent, ωϕ, ϕ_child, Lϕ)
    
    # Normalize rows
    for i = 1:(nop_child+1)
        row_sum_θ = sum(Lθ[i, :])
        if abs(row_sum_θ) > 1e-14
            Lθ[i, :] ./= row_sum_θ
        end
        
        row_sum_ϕ = sum(Lϕ[i, :])
        if abs(row_sum_ϕ) > 1e-14
            Lϕ[i, :] ./= row_sum_ϕ
        end
    end
    
    return Lθ, Lϕ
end

# =========================================================================
# Helper: Find corresponding spatial location in parent element
# =========================================================================

function find_corresponding_spatial_location(
    iel_child, i_child, j_child, k_child,
    iel_parent, mesh, ngl
)
    """
    Find the corresponding spatial location in the parent element.
    
    If iel_parent == iel_child: return same (i, j, k)
    If iel_parent != iel_child: find the shared spatial node
    """
    
    if iel_parent == iel_child
        # Same spatial element - same location
        return i_child, j_child, k_child
    end
    
    # Different spatial elements - find shared spatial node
    ip_child = mesh.connijk[iel_child, i_child, j_child, k_child]
    gip_child = mesh.ip2gip[ip_child]
    
    # Find this node in parent element
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip_parent = mesh.connijk[iel_parent, i, j, k]
        gip_parent = mesh.ip2gip[ip_parent]
        
        if gip_parent == gip_child
            return i, j, k
        end
    end
    
    # Should not reach here if elements actually share this node
    @warn "Could not find corresponding spatial node in parent element"
    return i_child, j_child, k_child  # Fallback
end
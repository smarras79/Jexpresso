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
# Build Restriction Matrices - CORRECTED
# =========================================================================

function build_restriction_matrices_local_and_ghost(
    connijk_spa, nc_non_global_nodes, n_spa, 
    ghost_layer::NonConformingGhostLayer,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    extra_meshes_extra_Je,
    mesh, ngl, nelem,
    neighbors,
    ip2gip_spa, gip2owner_extra, gid_to_extended_local, extended_local_to_gid,
    rank
)
    """
    Build restriction matrices with identity for ALL hanging nodes initially.
    
    nc_mat: (n_spa × n_spa)
    
    Rows:
    1. Free nodes (1:n_free): Identity
    2. Interior hanging nodes (n_free+1:n_spa, not interface): Constraint equation
    3. Interface hanging nodes: Identity (will be used to send data, then removed)
    
    After all operations, we manually remove ALL hanging node rows/columns.
    """
    
    #These are used to store information about pure ghost parents
    parent_ghost_idx = 1
    gid_to_extended_parents = Dict{Int, Int}()
    extended_parents_to_gid = Int[]
    extended_parents_x = Float64[]
    extended_parents_y = Float64[]
    extended_parents_z = Float64[]
    extended_parents_θ = Float64[]
    extended_parents_ϕ = Float64[]
    extended_parents_ip = Int[]


    #store parent data for later use
    local_parent_indices = Set{Int}()
    local_non_owned_parents = Set{Int}()
    nonowned_parent_gids = Set{Int}()

    # Hanging node set definition
    n_free = n_spa - length(nc_non_global_nodes)
    hanging_node_set = Set(nc_non_global_nodes)
    
    @info "[Rank $rank] Building restriction matrices:"
    @info "  Free DOFs: $n_free"
    @info "  Total hanging: $(length(nc_non_global_nodes))"
    @info "  Interior hanging: $(length(setdiff(hanging_node_set, ghost_layer.interface_hanging_nodes)))"
    @info "  Interface hanging: $(length(ghost_layer.interface_hanging_nodes))"
    
    # Triplet storage for nc_mat
    I_local = Int[]
    J_local = Int[]
    V_local = Float64[]

    I_rhs = Int[]
    J_rhs = Int[]
    V_rhs = Float64[]
    
    # Ghost constraint data
    ghost_constraint_data = Dict{Int, Vector{Tuple{Int, Float64}}}()
    ghost_constraint_data_rhs = Dict{Int, Vector{Tuple{Int, Float64}}}()
    # Track ALL hanging nodes for final removal
    # CRITICAL: Include both interior AND interface-only hanging nodes
    all_hanging_nodes = Set{Int}(nc_non_global_nodes)
    union!(all_hanging_nodes, ghost_layer.interface_hanging_nodes)
    
    estimated_entries = n_spa + length(nc_non_global_nodes) * 20
    sizehint!(I_local, estimated_entries)
    sizehint!(J_local, estimated_entries)
    sizehint!(V_local, estimated_entries)

    sizehint!(I_rhs, estimated_entries)
    sizehint!(J_rhs, estimated_entries)
    sizehint!(V_rhs, estimated_entries)
    
    # =========================================================================
    # Phase 1: Free nodes - Identity
    # =========================================================================
    
    for ip_free = 1:n_free
        if !(ip_free in all_hanging_nodes )
            push!(I_local, ip_free)
            push!(J_local, ip_free)
            push!(V_local, 1.0)
            push!(I_rhs, ip_free)
            push!(J_rhs, ip_free)
            push!(V_rhs, 1.0)
        end
    end
    
    @info "[Rank $rank] Added identity for $n_free free nodes"
    
    # =========================================================================
    # Phase 2: Interior hanging nodes - Constraint equations
    # =========================================================================
    
    interpolation_cache = Dict{NTuple{4,Int}, Tuple{Matrix{Float64}, Matrix{Float64}}}()
    interior_hanging = setdiff(hanging_node_set, ghost_layer.interface_hanging_nodes)
    
    @info "[Rank $rank] Processing interior hanging: $(length(interior_hanging))"
    n_processed = 0
    for ip_hanging in interior_hanging
        
        constraint_entries = build_interior_hanging_constraint(
            ip_hanging, connijk_spa, extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_extra_Je,
            mesh, ngl, nelem, neighbors, interpolation_cache
        )
        
        # Add constraint equation to nc_mat (replaces identity)
        for (parent_local_idx, weight) in constraint_entries
            if parent_local_idx <= n_spa
                #push!(I_local, ip_hanging)
                #push!(J_local, parent_local_idx)
                push!(I_local, parent_local_idx)
                push!(J_local, ip_hanging)
                push!(V_local, weight)
                push!(I_rhs, parent_local_idx)
                push!(J_rhs, ip_hanging)
                push!(V_rhs, weight)
            end
        end
        n_processed += 1
    end
    
    @info "[Rank $rank] Interior hanging constraint equations added to nc_mat"
    @info "max of nc_mat with interior", rank, maximum(sparse(I_local, J_local, V_local, n_spa, n_spa))
    # =========================================================================
    # Phase 3: Interface hanging nodes - Identity + Store constraints
    # =========================================================================
    
    interface_hanging = ghost_layer.interface_hanging_nodes
    
    @info "[Rank $rank] Processing interface hanging: $(length(interface_hanging))"
    
    for ip_hanging in interface_hanging
        # Add IDENTITY to nc_mat (will use this to extract data before removal)
        #=if !(ip_hanging in interior_hanging) && !(ip_hanging <= n_free)
            push!(I_local, ip_hanging)
            push!(J_local, ip_hanging)
            push!(V_local, 1.0)
            push!(I_rhs, ip_hanging)
            push!(J_rhs, ip_hanging)
            push!(V_rhs, 1.0)
        end=#
        
        # Get and store constraint information
        if !haskey(ghost_layer.parent_search_cache, ip_hanging)
            @warn "[Rank $rank] Hanging node $ip_hanging has no parent cache"
            continue
        end
        
        ghost_info = ghost_layer.parent_search_cache[ip_hanging]
        
        constraint_entries = build_interface_hanging_constraint(
            ip_hanging, ghost_info, connijk_spa,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_extra_Je,
            mesh, ngl, nelem,
            ip2gip_spa, gid_to_extended_local,
            n_spa, n_free,
            interpolation_cache,
            rank
        )
        
        # Separate local and ghost parents
        local_parents = Tuple{Int, Float64}[]
        ghost_parents = Tuple{Int, Float64}[]

        for (parent_extended_idx, weight, θ, ϕ, ip) in constraint_entries
            if parent_extended_idx <= n_spa
                push!(local_parents, (parent_extended_idx, weight))
            else
                ghost_gid = extended_local_to_gid[parent_extended_idx]
                push!(ghost_parents, (ghost_gid, weight))
                if !(haskey(gid_to_extended_parents,ghost_gid))
                    gid_to_extended_parents[ghost_gid] = n_spa + parent_ghost_idx
                    push!(extended_parents_to_gid,ghost_gid)
                    push!(extended_parents_x, mesh.x[ip])
                    push!(extended_parents_y, mesh.y[ip])
                    push!(extended_parents_z, mesh.z[ip])
                    push!(extended_parents_θ, θ)
                    push!(extended_parents_ϕ, ϕ)
                    push!(extended_parents_ip, ip)
                    parent_ghost_idx += 1
                    push!(I_local, n_spa + parent_ghost_idx-1)
                    push!(J_local, n_spa + parent_ghost_idx-1)
                    push!(V_local, 1.0)           
                end
            end
        end
        
        # CRITICAL: Separate owned vs non-owned parents to avoid double-counting
        # - Parents OWNED by this rank go to nc_mat for local constraint application
        # - Parents NOT owned by this rank go to ghost_constraint_data for remote communication

        owned_local_parents = Tuple{Int, Float64}[]
        nonowned_local_parents = Tuple{Int, Float64}[]

        for (parent_local_idx, weight) in local_parents
            if parent_local_idx <= n_spa
                if gip2owner_extra[parent_local_idx] == rank
                    push!(owned_local_parents, (parent_local_idx, weight))
                else
                    push!(nonowned_local_parents, (parent_local_idx, weight))
                end
            end
        end

        # Add OWNED local parents to nc_mat and nc_mat_rhs
        for (parent_local_idx, weight) in owned_local_parents
            push!(I_local, parent_local_idx)
            push!(J_local, ip_hanging)
            push!(V_local, weight)
            push!(I_rhs, parent_local_idx)
            push!(J_rhs, ip_hanging)
            push!(V_rhs, weight)
            if !(parent_local_idx in local_parent_indices)
                push!(local_parent_indices, parent_local_idx)
            end
        end

        # Add NON-OWNED local parents to nc_mat only
        for (parent_local_idx, weight) in nonowned_local_parents
            push!(I_local, parent_local_idx)
            push!(J_local, ip_hanging)
            push!(V_local, weight)
            if !(parent_local_idx in local_parent_indices)
                push!(local_parent_indices, parent_local_idx)
                push!(local_non_owned_parents, parent_local_idx)
                push!(nonowned_parent_gids, ip2gip_spa[parent_local_idx])
            end
        end


        # Store non-owned parent constraints for inter-processor effects exchange for rhs
        # Store owned and non-owned for matrix exchange.
        # This includes:
        # - Local parents NOT owned by this rank (needed by other ranks)
        # - Ghost parents (from remote ranks)
        all_parent_constraints = [(gid, w) for (gid, w) in vcat(
            [(ip2gip_spa[p_idx], w) for (p_idx, w) in local_parents],
            ghost_parents
        )]

        all_parent_constraints_rhs = [(gid, w) for (gid, w) in vcat(
            [(ip2gip_spa[p_idx], w) for (p_idx, w) in nonowned_local_parents],
            ghost_parents
        )]
        if !isempty(all_parent_constraints)
            ghost_constraint_data[ip_hanging] = all_parent_constraints
        end

        if !isempty(all_parent_constraints_rhs)
            ghost_constraint_data_rhs[ip_hanging] = all_parent_constraints_rhs
        end

    end
    
    @info "[Rank $rank] Interface hanging: $(length(ghost_constraint_data)) nodes with constraints"
    #@info "max of nc_mat before construction", rank, maximum(sparse(I_local, J_local, V_local, n_spa, n_spa))
    # =========================================================================
    # Build nc_mat: (n_spa × n_spa)
    # =========================================================================
    @info maximum(I_local), maximum(J_local), n_free, n_spa
    nc_mat = sparse(I_local, J_local, V_local, n_spa+parent_ghost_idx-1, n_spa+parent_ghost_idx-1)
    nc_mat_rhs = sparse(I_rhs, J_rhs, V_rhs, n_spa + parent_ghost_idx-1, n_spa + parent_ghost_idx-1)
    @info "maximum of nc_mat", rank, maximum(nc_mat)
    # Normalize constraint rows (skip free nodes and interface hanging identities)
    for i = 1:n_spa
        if i in interior_hanging
            col_sum = 0.0
            for idx in nzrange(nc_mat, i)
                col_sum += nonzeros(nc_mat)[idx]
            end
        
            if abs(col_sum - 1.0) > 1e-13 && abs(col_sum) > 1e-14
                @info "Column $j sum = $col_sum"
                for idx in nzrange(nc_mat, i)
                    nonzeros(nc_mat)[idx] /= col_sum
                end
            end
        end
    end

    for i = 1:n_spa
        if i in interior_hanging
            col_sum = 0.0
            for idx in nzrange(nc_mat_rhs, i)
                col_sum += nonzeros(nc_mat_rhs)[idx]
            end
        
            if abs(col_sum - 1.0) > 1e-13 && abs(col_sum) > 1e-14
                @info "Column $j sum = $col_sum"
                for idx in nzrange(nc_mat_rhs, i)
                    nonzeros(nc_mat_rhs)[idx] /= col_sum
                end
            end
        end
    end

    @info "maximum of nc_mat after column normalization", rank, maximum(nc_mat)
    P = nc_mat'
    P_vec = nc_mat_rhs'
    @info "[Rank $rank] nc_mat: $(size(nc_mat)), nnz=$(nnz(nc_mat))"
    #Separate matrix to used for rhs restriction from one used for matrix restriction
    return nc_mat, P, nc_mat_rhs, P_vec, ghost_constraint_data, ghost_constraint_data_rhs,
        all_hanging_nodes, gid_to_extended_parents, extended_parents_to_gid,
        extended_parents_x, extended_parents_y, extended_parents_z, extended_parents_θ, extended_parents_ϕ, extended_parents_ip,
        local_parent_indices, local_non_owned_parents, nonowned_parent_gids
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
        iel_child, i_child, j_child, k_child, e_ext_child, mesh.connijk,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops, extra_meshes_extra_nelems,
        neighbors, nelem
    )
    
    if iel_parent === nothing
        @warn "Could not find parent element for hanging node $ip_hanging in element $iel_child, angular element $e_ext_child"
        return Tuple{Int, Float64}[]
    end
    
    #=if (iel_child == 14 && i_child == 5 && j_child == 3 && k_child == 1)
        @info iel_child, i_child, j_child, k_child, e_ext_child, iθ_child, jθ_child, iel_parent, e_ext_parent
    end=#
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
    if (iel_child == 14 && i_child == 5 && j_child == 3 && k_child == 1 && e_ext_child == 3 
        && iθ_child == 3 && jθ_child == 3)
        @info "e_ext=3 hanging parent" iel_parent, e_ext_parent, i_parent, j_parent, k_parent
    end
    if (iel_child == 14 && i_child == 5 && j_child == 3 && k_child == 2 && e_ext_child == 3 
        && iθ_child == 3 && jθ_child == 3)
        @info "e_ext=3 hanging parent k=2" iel_parent, e_ext_parent, i_parent, j_parent, k_parent
    end

    for jθ_parent = 1:(nop_parent+1)
        for iθ_parent = 1:(nop_parent+1)
            Lθ_val = Lθ[iθ_child, iθ_parent]
            Lϕ_val = Lϕ[jθ_child, jθ_parent]
            
            if abs(Lθ_val) < 1e-14 || abs(Lϕ_val) < 1e-14
                continue
            end

            # Skip if this is identity (vertex node)
            if abs(Lθ_val - 1.0) < 1e-14 && abs(Lϕ_val - 1.0) < 1e-14
                continue
            end
            
            # Get parent node index
            ip_parent = connijk_spa[iel_parent][i_parent, j_parent, k_parent, 
                                                 e_ext_parent, iθ_parent, jθ_parent]
            
            Je_parent = extra_meshes_extra_Je[iel_parent][e_ext_parent, 
                                                           iθ_parent, jθ_parent]
            
            # Weight includes interpolation and Jacobian ratio
            
            weight = Lθ_val * Lϕ_val #* (Je_child / Je_parent)
            
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
    Lθ, Lϕ, θ_parent, ϕ_parent = build_ghost_interpolation_matrices(
        iel_child, e_ext_child, ghost_info, e_ext_parent,
        extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops
    )
    
    # Get child node Jacobian
    Je_child = extra_meshes_extra_Je[iel_child][e_ext_child, iθ_child, jθ_child]
    
    # Build constraint entries
    constraint_entries = Tuple{Int, Float64, Float64, Float64, Int}[]
    
    nop_parent = ghost_info.nop_ang[e_ext_parent]
    
    # Find which interface spatial index corresponds to (i_child, j_child, k_child)
    ip_spatial = mesh.connijk[iel_child, i_child, j_child, k_child]
    gip_spatial = mesh.ip2gip[ip_spatial]
    
    # Verify spatial node is in ghost's interface
    if !(gip_spatial in ghost_info.interface_spatial_global_ids)
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

            # Get ghost parent node global ID using gip_spatial directly
            gid_parent = get_ghost_interface_node_gid(
                ghost_info, gip_spatial, e_ext_parent,
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
            weight = Lθ_val * Lϕ_val #* (Je_child / Je_parent)
            
            push!(constraint_entries, (extended_local_idx, weight ,θ_parent[iθ_parent], ϕ_parent[jθ_parent], ip_spatial))
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
    iel_child, i_child, j_child, k_child, e_ext_child, connijk,
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
    ip_child = connijk[iel_child, i_child, j_child, k_child]
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

        if !(ip_child in connijk[iel_neighbor, :, :, :])
            # non-conforming neighbor but the spatial node is not on this particular neighbor
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
            if is_child_element(θmin_neighbor, θmax_neighbor, ϕmin_neighbor, ϕmax_neighbor,
                               θmin_child, θmax_child, ϕmin_child, ϕmax_child)
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
        
        if is_child_element(θmin_ghost, θmax_ghost, ϕmin_ghost, ϕmax_ghost,
                           θmin_child, θmax_child, ϕmin_child, ϕmax_child)
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
        if i == nop_child+1 && ϕ_child[i] < 1e-10 && ϕ_child[i-1] > π
            ϕ_child[i] = 2π
        end
    end
    
    # Get parent grid points
    θ_parent = zeros(Float64, nop_parent+1)
    ϕ_parent = zeros(Float64, nop_parent+1)
    
    for i = 1:(nop_parent+1)
        ip = extra_meshes_connijk[iel_parent][e_ext_parent, i, i]
        θ_parent[i] = extra_meshes_coords[iel_parent][1, ip]
        ϕ_parent[i] = extra_meshes_coords[iel_parent][2, ip]
        if i == nop_parent+1 && ϕ_parent[i] < 1e-10 && ϕ_parent[i-1] > π
            ϕ_parent[i] = 2π
        end
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
        if i == nop_child+1 && ϕ_child[i] < 1e-10 && ϕ_child[i-1] > π
            ϕ_child[i] = 2π
        end
    end
    
    # Get parent (ghost) grid points
    θ_parent = zeros(Float64, nop_parent+1)
    ϕ_parent = zeros(Float64, nop_parent+1)
    
    for i = 1:(nop_parent+1)
        ip = ghost_info.connijk_ang[e_ext_parent, i, i]
        θ_parent[i] = ghost_info.coords_ang[1, ip]
        ϕ_parent[i] = ghost_info.coords_ang[2, ip]
        if i == nop_parent+1 && ϕ_parent[i] < 1e-10 && ϕ_parent[i-1] > π
            ϕ_parent[i] = 2π
        end
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
    
    return Lθ, Lϕ, θ_parent, ϕ_parent
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
        @warn "parent and child are the same somehow"
        return i_child, j_child, k_child
    end
    
    # Different spatial elements - find shared spatial node
    ip_child = mesh.connijk[iel_child, i_child, j_child, k_child]
    
    # Find this node in parent element
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip_parent = mesh.connijk[iel_parent, i, j, k]
        if ip_parent == ip_child
            return i, j, k
        end
    end
    
    # Should not reach here if elements actually share this node
    @warn "Could not find corresponding spatial node in parent element"
    return i_child, j_child, k_child  # Fallback
end

# =========================================================================
# Apply LEFT Ghost Effects
# =========================================================================

function apply_left_ghost_effects(
    A_local, ghost_effects_left, 
    ip2gip_spa, n_free, n_spa, rank
)
    """
    Apply ghost constraint effects from LEFT multiplication (restriction).
    
    For each hanging node that depends on ghost parents:
    - Add the weighted contribution from ghost parent rows
    """
    
    I_vec, J_vec, V_vec = findnz(A_local)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)
    
    # Ghost effects represent: for hanging node i, add sum(weight * row_ghost_parent)
    for (gid_hanging, row_contributions) in ghost_effects_left
        # Find local index of hanging node
        ip_hanging_local = findfirst(x -> x == gid_hanging, ip2gip_spa)
        
        if ip_hanging_local === nothing || ip_hanging_local > n_free
            continue  # Not owned or not a free node index
        end
        
        # Add contributions to this row
        for (gid_col, value) in row_contributions
            ip_col_local = findfirst(x -> x == gid_col, ip2gip_spa)
            
            if ip_col_local !== nothing && ip_col_local <= n_spa
                push!(I_new, ip_hanging_local)
                push!(J_new, ip_col_local)
                push!(V_new, value)
            end
        end
    end
    
    # Rebuild matrix (sums duplicates)
    A_with_ghost = sparse(I_new, J_new, V_new, n_free, n_spa)
    
    return A_with_ghost
end

# =========================================================================
# Exchange Ghost Constraints for RIGHT multiplication (Transpose)
# =========================================================================

function exchange_ghost_constraints_transpose(
    ghost_constraints, A_current,
    ip2gip_spa, gip2owner_spa,
    n_free, n_spa, rank, comm
)
    """
    Exchange ghost constraints for RIGHT prolongation (P^T operation).
    
    This is the transpose operation: we need columns instead of rows.
    For each hanging node that depends on ghost parents:
    - Get the column of A_current for the hanging node
    - Distribute it to ghost parents (weighted)
    """
    
    # Organize what to send
    send_to_rank = Dict{Int, Vector{Tuple{Int, Int, Float64}}}()
    
    for (ip_hanging_local, ghost_parents) in ghost_constraints
        gid_hanging = ip2gip_spa[ip_hanging_local]
        
        # Get column of A_current for this hanging node
        # A_current is (n_free × n_spa), so column ip_hanging_local
        col_entries = Tuple{Int, Float64}[]
        
        # Extract column from CSC matrix
        if ip_hanging_local <= n_spa
            for idx in nzrange(A_current, ip_hanging_local)
                row = rowvals(A_current)[idx]
                val = nonzeros(A_current)[idx]
                gid_row = ip2gip_spa[row]
                push!(col_entries, (gid_row, val))
            end
        end
        
        # For each ghost parent, send weighted column
        for (gid_ghost_parent, weight) in ghost_parents
            owner_ghost = gip2owner_spa[gid_ghost_parent]
            
            if !haskey(send_to_rank, owner_ghost)
                send_to_rank[owner_ghost] = Tuple{Int, Int, Float64}[]
            end
            
            # Send: (gid_ghost_parent_col, gid_row, weighted_value)
            for (gid_row, value) in col_entries
                push!(send_to_rank[owner_ghost], (gid_ghost_parent, gid_row, weight * value))
            end
        end
    end
    
    # Exchange
    send_buffers = Dict{Int, Vector{UInt8}}()
    for (dest_rank, data) in send_to_rank
        buffer = IOBuffer()
        serialize(buffer, data)
        send_buffers[dest_rank] = take!(buffer)
    end
    
    recv_buffers = exchange_buffers(send_buffers,
                                    Dict(k => [] for k in keys(send_to_rank)),
                                    comm)
    
    # Organize received contributions
    ghost_effects_right = Dict{Int, Vector{Tuple{Int, Float64}}}()
    # gid_col -> [(gid_row, value)]
    
    for (src_rank, buffer) in recv_buffers
        data = deserialize(IOBuffer(buffer))
        for (gid_col, gid_row, value) in data
            if !haskey(ghost_effects_right, gid_col)
                ghost_effects_right[gid_col] = Tuple{Int, Float64}[]
            end
            push!(ghost_effects_right[gid_col], (gid_row, value))
        end
    end
    
    return ghost_effects_right
end

# =========================================================================
# Apply RIGHT Ghost Effects
# =========================================================================

function apply_right_ghost_effects(
    A_local, ghost_effects_right,
    ip2gip_spa, n_free, rank
)
    """
    Apply ghost constraint effects from RIGHT multiplication (prolongation).
    
    For each ghost parent that has hanging children:
    - Add the weighted contributions to the corresponding column
    """
    
    I_vec, J_vec, V_vec = findnz(A_local)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)
    
    # Ghost effects represent: for ghost parent column, add contributions
    for (gid_ghost_col, row_contributions) in ghost_effects_right
        # Find local index of ghost parent (should be a free node here)
        ip_col_local = findfirst(x -> x == gid_ghost_col, ip2gip_spa)
        
        if ip_col_local === nothing || ip_col_local > n_free
            continue
        end
        
        # Add contributions to this column
        for (gid_row, value) in row_contributions
            ip_row_local = findfirst(x -> x == gid_row, ip2gip_spa)
            
            if ip_row_local !== nothing && ip_row_local <= n_free
                push!(I_new, ip_row_local)
                push!(J_new, ip_col_local)
                push!(V_new, value)
            end
        end
    end
    
    # Rebuild matrix
    A_final = sparse(I_new, J_new, V_new, n_free, n_free)
    
    return A_final
end

function exchange_ghost_constraints(
    ghost_constraints, M_inv_LHS,
    ip2gip_spa, gip2owner_spa,
    n_spa, rank, comm
)
    """
    Exchange ghost constraint information between processors.
    
    For each hanging node that depends on ghost parents:
    1. Compute the effect: weight * (row of M_inv_LHS corresponding to ghost parent)
    2. Send to the processor that owns the hanging node
    
    Returns: Dict{Int, SparseVector} - gid_hanging -> constraint_row
    """
    
    nproc = MPI.Comm_size(comm)
    
    # Organize what to send to each processor
    send_to_rank = Dict{Int, Vector{Tuple{Int, Int, Float64}}}()  # rank -> [(gid_hanging, gid_col, value)]
    
    for (ip_hanging_local, ghost_parents) in ghost_constraints
        # Get global ID of hanging node
        gid_hanging = ip2gip_spa[ip_hanging_local]
        
        # For each ghost parent this hanging node depends on
        for (gid_ghost_parent, weight) in ghost_parents
            # The ghost parent is owned by some other processor
            # We need to get the row of M_inv_LHS corresponding to this ghost
            # But we don't have it locally!
            
            # Instead, we send a request to the owner of the ghost parent
            # asking them to compute: weight * (M_inv_LHS row for ghost_parent)
            # and send it to the owner of the hanging node
            
            owner_ghost = gip2owner_spa[gid_ghost_parent]
            
            if !haskey(send_to_rank, owner_ghost)
                send_to_rank[owner_ghost] = Tuple{Int, Int, Float64}[]
            end
            
            # Request: (hanging_gid, ghost_parent_gid, weight)
            push!(send_to_rank[owner_ghost], (gid_hanging, gid_ghost_parent, weight))
        end
    end
    
    @info "[Rank $rank] Sending ghost constraint requests to $(length(send_to_rank)) ranks"
    
    # Exchange requests
    received_requests = exchange_ghost_requests(send_to_rank, rank, comm)
    
    @info "[Rank $rank] Received requests from $(length(received_requests)) ranks"
    
    # Process requests and compute responses
    responses = Dict{Int, Vector{Tuple{Int, Vector{Tuple{Int, Float64}}}}}()  
    # rank -> [(gid_hanging, [(gid_col, value)])]
    
    for (src_rank, requests) in received_requests
        rank_responses = Tuple{Int, Vector{Tuple{Int, Float64}}}[]
        
        for (gid_hanging, gid_ghost_parent, weight) in requests
            # Find local index of ghost_parent
            ip_ghost_local = findfirst(x -> x == gid_ghost_parent, ip2gip_spa)
            
            if ip_ghost_local === nothing
                @warn "[Rank $rank] Requested ghost parent gid=$gid_ghost_parent not found locally"
                continue
            end
            
            # Get row of M_inv_LHS for this ghost parent
            row_entries = Tuple{Int, Float64}[]
            rows = rowvals(M_inv_LHS)
            vals = nonzeros(M_inv_LHS)
            
            for idx in nzrange(M_inv_LHS, ip_ghost_local)
                col = rows[idx]
                val = vals[idx]
                gid_col = ip2gip_spa[col]
                
                # Scale by weight
                push!(row_entries, (gid_col, weight * val))
            end
            
            push!(rank_responses, (gid_hanging, row_entries))
        end
        
        if !isempty(rank_responses)
            responses[src_rank] = rank_responses
        end
    end
    
    # Exchange responses
    ghost_effects = exchange_ghost_responses(responses, rank, comm)
    
    return ghost_effects
end

# =========================================================================
# Helper: Exchange Ghost Requests
# =========================================================================

function exchange_ghost_requests(send_to_rank, rank, comm)
    """Exchange ghost constraint requests"""
    
    # Serialize and exchange
    send_buffers = Dict{Int, Vector{UInt8}}()
    
    for (dest_rank, requests) in send_to_rank
        buffer = IOBuffer()
        serialize(buffer, requests)
        send_buffers[dest_rank] = take!(buffer)
    end
    
    # Use the exchange_buffers function we already have
    recv_buffers = exchange_buffers(send_buffers, 
                                    Dict(k => [] for k in keys(send_to_rank)),
                                    comm)
    
    # Deserialize
    received_requests = Dict{Int, Vector{Tuple{Int, Int, Float64}}}()
    for (src_rank, buffer) in recv_buffers
        received_requests[src_rank] = deserialize(IOBuffer(buffer))
    end
    
    return received_requests
end

# =========================================================================
# Helper: Exchange Ghost Responses
# =========================================================================

function exchange_ghost_responses(responses, rank, comm)
    """Exchange computed ghost effects back"""
    
    # Serialize and exchange
    send_buffers = Dict{Int, Vector{UInt8}}()
    
    for (dest_rank, response_data) in responses
        buffer = IOBuffer()
        serialize(buffer, response_data)
        send_buffers[dest_rank] = take!(buffer)
    end
    
    # Build requests_by_rank for ALL ranks (including those not in send_buffers)
    # This prevents deadlock when only some ranks have data to send
    nproc = MPI.Comm_size(comm)
    requests_by_rank = Dict{Int, Vector}()
    for src_rank = 0:(nproc-1)
        requests_by_rank[src_rank] = []
    end

    recv_buffers = exchange_buffers(send_buffers, requests_by_rank, comm)
    
    # Deserialize and organize by hanging node gid
    ghost_effects = Dict{Int, Vector{Tuple{Int, Float64}}}()  # gid_hanging -> [(gid_col, value)]
    
    for (src_rank, buffer) in recv_buffers
        response_list = deserialize(IOBuffer(buffer))
        
        for (gid_hanging, row_entries) in response_list
            if !haskey(ghost_effects, gid_hanging)
                ghost_effects[gid_hanging] = Tuple{Int, Float64}[]
            end
            append!(ghost_effects[gid_hanging], row_entries)
        end
    end
    
    return ghost_effects
end

# =========================================================================
# Compute Row Effects BEFORE Restriction
# =========================================================================

function compute_hanging_row_effects_before_restriction(
    ghost_constraint_data, M_inv_LHS, ip2gip_spa, gip2owner_spa, rank,
    gid_to_extended_parents, extended_parents_to_gid, gip_to_local
)
    n_spa = size(M_inv_LHS, 1)
    I_vec, J_vec, V_vec = findnz(M_inv_LHS)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)
    n_ghost_parents = length(extended_parents_to_gid)

    # extended_parents_to_gid is an array of GIDs — set of values for O(1) lookup
    extended_gid_set = Set(extended_parents_to_gid)

    effects_by_owner = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()

    # Build row lookup from findnz — O(nnz) once, avoids dense [row,col] indexing
    row_index = Dict{Int, Vector{Tuple{Int,Float64}}}()
    for k = 1:length(I_vec)
        list = get!(row_index, I_vec[k], Tuple{Int,Float64}[])
        push!(list, (J_vec[k], V_vec[k]))
    end

    for (ip_hanging, parent_constraints) in ghost_constraint_data

        hanging_row = Tuple{Int,Float64}[]
        for (col, val) in get(row_index, ip_hanging, Tuple{Int,Float64}[])
            col > n_spa && continue
            push!(hanging_row, (ip2gip_spa[col], val))
        end

        for (parent_gid, weight) in parent_constraints
            owner = gip2owner_spa[parent_gid]
            owner == rank && continue

            if parent_gid in extended_gid_set
                ip_parent = gid_to_extended_parents[parent_gid]
                for (gid_col, value) in hanging_row
                    ip_col = get(gip_to_local, gid_col, nothing)
                    ip_col === nothing && continue
                    push!(I_new, ip_parent)
                    push!(J_new, ip_col)
                    push!(V_new, weight * value)
                end
            else
                dest = get!(effects_by_owner, owner, Tuple{Int,Int,Float64}[])
                for (gid_col, value) in hanging_row
                    push!(dest, (parent_gid, gid_col, weight * value))
                end
            end
        end
    end

    M_inv_LHS = sparse(I_new, J_new, V_new,
                       n_spa + n_ghost_parents, n_spa + n_ghost_parents)
    return effects_by_owner, M_inv_LHS
end

# =========================================================================
# Compute Column Effects BEFORE Prolongation
# =========================================================================

function compute_hanging_col_effects_before_prolongation(
    ghost_constraint_data, A_matrix, ip2gip_spa, gip2owner_spa,
    n_spa, n_spa_g, all_hanging_nodes, rank,
    gid_to_extended_parents, extended_parents_to_gid, gip_to_local
)
    n_spa_g   = size(A_matrix, 1)
    I_vec, J_vec, V_vec = findnz(A_matrix)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)

    extended_gid_set  = Set(extended_parents_to_gid)
    hanging_set       = all_hanging_nodes isa Set ? all_hanging_nodes :
                        Set(all_hanging_nodes)

    # Build extended reverse lookup
    gip_to_extended = copy(gip_to_local)
    for (gid, local_idx) in gid_to_extended_parents
        gip_to_extended[gid] = local_idx
    end

    effects_by_owner = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()

    for (ip_hanging, parent_constraints) in ghost_constraint_data
        ip_hanging > size(A_matrix, 2) && continue

        # Extract column ip_hanging using nzrange — correct for CSC
        hanging_col = Tuple{Int, Float64}[]
        for idx in nzrange(A_matrix, ip_hanging)
            row = rowvals(A_matrix)[idx]
            val = nonzeros(A_matrix)[idx]
            row in hanging_set && continue

            gid_row = if row <= n_spa
                ip2gip_spa[row]
            else
                extended_parents_to_gid[row - n_spa]
            end
            push!(hanging_col, (gid_row, val))
        end

        for (parent_gid, weight) in parent_constraints
            owner = gip2owner_spa[parent_gid]
            owner == rank && continue

            if parent_gid in extended_gid_set
                ip_parent = gid_to_extended_parents[parent_gid]
                for (gid_row, value) in hanging_col
                    ip_row = get(gip_to_extended, gid_row, nothing)
                    ip_row === nothing && continue
                    push!(I_new, ip_row)
                    push!(J_new, ip_parent)
                    push!(V_new, weight * value)
                end
            else
                dest = get!(effects_by_owner, owner, Tuple{Int,Int,Float64}[])
                for (gid_row, value) in hanging_col
                    push!(dest, (parent_gid, gid_row, weight * value))
                end
            end
        end
    end

    A_matrix = sparse(I_new, J_new, V_new, n_spa_g, n_spa_g)
    return effects_by_owner, A_matrix
end

# =========================================================================
# Exchange Hanging Effects
# =========================================================================

function exchange_hanging_effects(effects_to_send, rank, comm)
    nproc = MPI.Comm_size(comm)

    # Pack as flat Int/Float64 arrays: [parent_gid, col_gid, value, ...]
    send_int   = Dict{Int, Vector{Int}}()
    send_float = Dict{Int, Vector{Float64}}()
    for dest = 0:nproc-1
        effects = get(effects_to_send, dest, Tuple{Int,Int,Float64}[])
        ints   = Vector{Int}(undef,     2 * length(effects))
        floats = Vector{Float64}(undef, 1 * length(effects))
        for (k, (pg, cg, v)) in enumerate(effects)
            ints[2k-1]  = pg
            ints[2k]    = cg
            floats[k]   = v
        end
        send_int[dest]   = ints
        send_float[dest] = floats
    end

    # Exchange counts first
    send_counts = Int32[length(get(effects_to_send, d,
                    Tuple{Int,Int,Float64}[])) for d = 0:nproc-1]
    recv_counts = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)

    # Exchange int data
    send_int_flat   = vcat([send_int[d]   for d = 0:nproc-1]...)
    recv_int_flat   = Vector{Int}(undef, 2 * sum(recv_counts))
    MPI.Alltoallv!(MPI.VBuffer(send_int_flat,   Int32.(2 .* send_counts)),
               MPI.VBuffer(recv_int_flat,   Int32.(2 .* recv_counts)), comm)

    # Exchange float data  
    send_float_flat = vcat([send_float[d] for d = 0:nproc-1]...)
    recv_float_flat = Vector{Float64}(undef, sum(recv_counts))
    MPI.Alltoallv!(MPI.VBuffer(send_float_flat, Int32.(send_counts)),
               MPI.VBuffer(recv_float_flat, Int32.(recv_counts)), comm)

    # Unpack
    received = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()
    offset_i = 0; offset_f = 0
    for src = 0:nproc-1
        n = recv_counts[src+1]
        effects = Vector{Tuple{Int,Int,Float64}}(undef, n)
        for k = 1:n
            pg = recv_int_flat[offset_i + 2k-1]
            cg = recv_int_flat[offset_i + 2k]
            v  = recv_float_flat[offset_f + k]
            effects[k] = (pg, cg, v)
        end
        received[src] = effects
        offset_i += 2n
        offset_f += n
    end
    return received
end

function exchange_hanging_effects_vector(effects_to_send, rank, comm)
    nproc = MPI.Comm_size(comm)

    send_int   = Dict{Int, Vector{Int}}()
    send_float = Dict{Int, Vector{Float64}}()
    for dest = 0:nproc-1
        effects = get(effects_to_send, dest, Tuple{Int,Float64}[])
        ints   = Vector{Int}(undef,     length(effects))
        floats = Vector{Float64}(undef, length(effects))
        for (k, (pg, v)) in enumerate(effects)
            ints[k]   = pg
            floats[k] = v
        end
        send_int[dest]   = ints
        send_float[dest] = floats
    end
    send_counts = Int32[length(get(effects_to_send, d,
                    Tuple{Int,Int,Float64}[])) for d = 0:nproc-1]
    recv_counts = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)

    # Exchange int data
    send_int_flat   = vcat([send_int[d]   for d = 0:nproc-1]...)
    recv_int_flat   = Vector{Int}(undef, sum(recv_counts))
    MPI.Alltoallv!(MPI.VBuffer(send_int_flat,   Int32.(send_counts)),
               MPI.VBuffer(recv_int_flat,   Int32.(recv_counts)), comm)

    # Exchange float data  
    send_float_flat = vcat([send_float[d] for d = 0:nproc-1]...)
    recv_float_flat = Vector{Float64}(undef, sum(recv_counts))
    MPI.Alltoallv!(MPI.VBuffer(send_float_flat, Int32.(send_counts)),
               MPI.VBuffer(recv_float_flat, Int32.(recv_counts)), comm)

    received = Dict{Int, Vector{Tuple{Int,Float64}}}()
    offset = 0
    for src = 0:nproc-1
        n = recv_counts[src+1]
        effects = Vector{Tuple{Int,Float64}}(undef, n)
        for k = 1:n
            effects[k] = (recv_int_flat[offset+k], recv_float_flat[offset+k])
        end
        received[src] = effects
        offset += n
    end
    return received
end

## =========================================================================
# Add Row Effects (Updated for n_spa dimension)
# =========================================================================

function add_hanging_row_effects(A_local, received_effects, ip2gip_spa,
                                  n_spa, rank, gip_to_local)
    I_vec, J_vec, V_vec = findnz(A_local)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)

    for (src_rank, effects) in received_effects
        for (parent_gid, gid_col, value) in effects
            ip_parent = get(gip_to_local, parent_gid, nothing)
            ip_col    = get(gip_to_local, gid_col,    nothing)
            ip_parent === nothing && continue
            ip_col    === nothing && continue
            ip_parent > n_spa    && continue
            ip_col    > n_spa    && continue
            push!(I_new, ip_parent)
            push!(J_new, ip_col)
            push!(V_new, value)
        end
    end

    return sparse(I_new, J_new, V_new, n_spa, n_spa)
end

# =========================================================================
# Add Column Effects (Updated for n_spa dimension)
# =========================================================================

function add_hanging_col_effects(A_local, received_effects, ip2gip_spa,
                                  n_spa, rank, gip_to_local)
    I_vec, J_vec, V_vec = findnz(A_local)
    I_new = copy(I_vec)
    J_new = copy(J_vec)
    V_new = copy(V_vec)

    for (src_rank, effects) in received_effects
        for (parent_gid, gid_row, value) in effects
            ip_parent = get(gip_to_local, parent_gid, nothing)
            ip_row    = get(gip_to_local, gid_row,    nothing)
            ip_parent === nothing && continue
            ip_row    === nothing && continue
            ip_parent > n_spa    && continue
            ip_row    > n_spa    && continue
            push!(I_new, ip_row)
            push!(J_new, ip_parent)
            push!(V_new, value)
        end
    end

    return sparse(I_new, J_new, V_new, n_spa, n_spa)
end

function extract_free_submatrix_remove_all_hanging(
    A_full, all_hanging_nodes, n_free, n_spa_g, rank
)
    # Ensure O(1) lookup
    hanging_set = all_hanging_nodes isa Set ? all_hanging_nodes :
                  Set(all_hanging_nodes)
    
    I_full, J_full, V_full = findnz(A_full)
    n = length(I_full)

    I_free = sizehint!(Int[],     n)
    J_free = sizehint!(Int[],     n)
    V_free = sizehint!(Float64[], n)

    for idx = 1:n
        i = I_full[idx]; j = J_full[idx]
        (i in hanging_set || j in hanging_set) && continue
        push!(I_free, i)
        push!(J_free, j)
        push!(V_free, V_full[idx])
    end

    return sparse(I_free, J_free, V_free, n_spa_g, n_spa_g)
end

# =========================================================================
# Parallel RHS Restriction (Vector Version)
# =========================================================================

function compute_hanging_rhs_effects_before_restriction(
    ghost_constraint_data, RHS_vector, ip2gip_spa, gip2owner_spa, rank
)
    """
    Compute RHS effects for interface hanging nodes BEFORE applying nc_mat.

    For each interface hanging node:
    - Get its RHS value (identity preserved it)
    - Apply constraint weights
    - Send to ghost parent owners

    This is the vector version of compute_hanging_row_effects_before_restriction.
    """
    
    effects_by_owner = Dict{Int, Vector{Tuple{Int, Float64}}}()

    for (ip_hanging, parent_constraints) in ghost_constraint_data
        gid_hanging = ip2gip_spa[ip_hanging]
        # Get RHS value for hanging node
        rhs_value = RHS_vector[ip_hanging]

        # Apply constraint: for each parent, send weighted RHS value
        for (parent_gid, weight) in parent_constraints
            owner = gip2owner_spa[parent_gid]
            
            if !haskey(effects_by_owner, owner)
                effects_by_owner[owner] = Tuple{Int, Float64}[]
            end

            # Send: (parent_gid, weight * rhs_value)
            push!(effects_by_owner[owner], (parent_gid, weight * rhs_value))
        end
    end

    return effects_by_owner
end

function add_hanging_rhs_effects(RHS_local, received_effects, ip2gip_spa,
                                  n_spa, rank, gip_to_local)
    RHS_updated = copy(RHS_local)
    for (src_rank, effects) in received_effects
        for (parent_gid, value) in effects
            ip_parent = get(gip_to_local, parent_gid, nothing)
            ip_parent === nothing && continue
            ip_parent > n_spa    && continue
            RHS_updated[ip_parent] += value
        end
    end
    return RHS_updated
end

function extract_free_rhs_subvector(
    RHS_full, all_hanging_nodes, n_free, n_spa, n_spa_g, rank
)
    """
    Extract the (n_free) subvector by removing ALL hanging node entries.

    Keep only: entries 1:n_free (free nodes only)
    """

    RHS_free = zeros(eltype(RHS_full), n_spa_g)

    
    for i = 1:n_spa
        if !(i in all_hanging_nodes) && i <= n_free
            RHS_free[i] = RHS_full[i]
        elseif i <= n_free
            RHS_free[i] = 0
        end

    end

    @info "[Rank $rank] Extracted RHS free subvector: $n_free entries"

    return RHS_free
end

# =========================================================================
# Parallel Solution Prolongation (Vector Version)
# =========================================================================

function build_reverse_ghost_constraint_map(
    ghost_constraint_data, ip2gip_spa, gip2owner_spa, rank, gip_to_local
)
    reverse_map = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()

    for (ip_hanging, parent_constraints) in ghost_constraint_data
        gid_hanging   = ip2gip_spa[ip_hanging]
        owner_hanging = gip2owner_spa[gid_hanging]

        for (parent_gid, weight) in parent_constraints
            gip2owner_spa[parent_gid] == rank || continue
            ip_parent = get(gip_to_local, parent_gid, nothing)
            ip_parent === nothing && continue
            dest = get!(reverse_map, ip_parent, Tuple{Int,Int,Float64}[])
            push!(dest, (gid_hanging, owner_hanging, weight))
        end
    end
    return reverse_map
end

function compute_solution_prolongation_contributions(
    reverse_ghost_map, solution_free, ip2gip_spa, n_free, rank
)
    """
    For free nodes that are ghost parents of hanging nodes on other processors,
    compute weighted contributions to send.

    For each free parent node:
    - solution_value = solution_free[parent]
    - For each hanging child: send weight * solution_value to child's owner

    Returns: Dict{Int, Vector{Tuple{Int, Float64}}}
             owner_rank -> [(hanging_gid, contribution)]
    """

    contributions_by_owner = Dict{Int, Vector{Tuple{Int, Float64}}}()

    for (ip_parent, hanging_children) in reverse_ghost_map
        # Get solution value at this parent node
        if ip_parent <= n_free
            parent_value = solution_free[ip_parent]
            
            # Send weighted contributions to each hanging child
            for (gid_hanging, owner_hanging, weight) in hanging_children
                if !haskey(contributions_by_owner, owner_hanging)
                    contributions_by_owner[owner_hanging] = Tuple{Int, Float64}[]
                end

                contribution = weight * parent_value
                push!(contributions_by_owner[owner_hanging], (gid_hanging, contribution))
            end
        end
    end

    return contributions_by_owner
end

function add_solution_prolongation_contributions(
    solution_extended, received_contributions, ip2gip_spa,
    n_spa, rank, gip_to_local
)
    solution_updated = copy(solution_extended)
    for (src_rank, contributions) in received_contributions
        for (gid_hanging, value) in contributions
            ip_hanging = get(gip_to_local, gid_hanging, nothing)
            ip_hanging === nothing && continue
            ip_hanging > n_spa    && continue
            solution_updated[ip_hanging] += value
        end
    end
    return solution_updated
end

# =========================================================================
# All-Reduce Parent-Parent Matrix Entries
# =========================================================================

function allreduce_parent_parent_entries(
    MLHS, nc_mat, gip2owner_extra, n_spa, rank, comm, ip2gip_spa,
    local_parent_indices, local_non_owned_parents, nonowned_parent_gids,
    gip_to_local
)
    nproc = MPI.Comm_size(comm)

    # Exchange non-owned parent GIDs as flat Int array
    local_gids  = collect(nonowned_parent_gids)
    n_local     = Int32(length(local_gids))
    gid_counts  = MPI.Allgather([n_local], comm)
    all_gids    = MPI.Allgatherv(local_gids, gid_counts, comm)
    all_parent_gids = Set(all_gids)

    # Extract relevant (row_gid, col_gid, val) triples
    rows_sp = rowvals(MLHS)
    vals_sp = nonzeros(MLHS)
    local_I = Int[]
    local_J = Int[]
    local_V = Float64[]

    for col_idx = 1:n_spa
        col_gid = ip2gip_spa[col_idx]
        col_gid in all_parent_gids || continue
        for idx in nzrange(MLHS, col_idx)
            row_idx = rows_sp[idx]
            row_gid = ip2gip_spa[row_idx]
            row_gid in all_parent_gids || continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end

    # Pack as flat arrays: [row_gid, col_gid, ...] and [val, ...]
    n_local_entries = Int32(length(local_I))
    entry_counts    = MPI.Allgather([n_local_entries], comm)

    send_ij  = Vector{Int}(undef,     2 * n_local_entries)
    send_v   = Vector{Float64}(undef, n_local_entries)
    for k = 1:n_local_entries
        send_ij[2k-1] = local_I[k]
        send_ij[2k]   = local_J[k]
        send_v[k]     = local_V[k]
    end

    all_ij = MPI.Allgatherv(send_ij, Int32.(entry_counts .* 2), comm)
    all_v  = MPI.Allgatherv(send_v,  entry_counts,              comm)

    # Accumulate
    global_entries = Dict{Tuple{Int,Int}, Float64}()
    for k = 1:length(all_v)
        rg = all_ij[2k-1]; cg = all_ij[2k]; v = all_v[k]
        key = (rg, cg)
        global_entries[key] = get(global_entries, key, 0.0) + v
    end

    # Update MLHS in-place
    MLHS_updated = copy(MLHS)
    for ((row_gid, col_gid), reduced_val) in global_entries
        row_idx = get(gip_to_local, row_gid, nothing)
        col_idx = get(gip_to_local, col_gid, nothing)
        row_idx === nothing && continue
        col_idx === nothing && continue
        MLHS_updated[row_idx, col_idx] = reduced_val
    end

    return MLHS_updated
end
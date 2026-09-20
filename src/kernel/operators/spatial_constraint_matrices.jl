# ============================================================================
# Spatial Constraint Matrix Construction
# ============================================================================
# Build Lagrange interpolation constraint equations for spatial hanging nodes

"""
    build_spatial_constraint_matrices(
        mesh, spatial_amr_cache, extra_meshes_coords, extra_meshes_connijk,
        extra_meshes_extra_nops, extra_meshes_extra_nelems,
        ngl, rank
    ) -> SpatialAMRCache

Build Lagrange interpolation constraint matrices for spatial hanging nodes.

CRITICAL: Each spatial element couples with ALL associated angular nodes.

For each spatial non-conforming facet with a child element:
1. Get child and parent spatial element IDs
2. Extract spatial element node coordinates from mesh
3. Extract ALL angular mesh data for both elements (nelem_ang elements per spatial elem)
4. Build 1D Lagrange interpolation matrices for spatial dimensions
5. For each hanging spatial node, create constraint for ALL angular DOFs at that location
6. Store full constraint including spatial + angular dimensions

Populates cache.parent_weights with constraints for all spatial-angular nodes:
Dict{spatial_hanging_node_id, [(parent_spatial_node_id, weight), ...]}

For each spatial hanging node location, ALL angular nodes at that location inherit
the spatial constraints (with identity on angular dimension).

Returns updated cache with spatially-constrained angular nodes.
"""
function build_spatial_constraint_matrices(
    mesh::St_mesh, spatial_amr_cache::SpatialAMRCache,
    extra_meshes_coords, extra_meshes_connijk,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    ngl::Int, rank::Int
)
    @rankinfo rank "Building spatial constraint matrices with angular coupling..."

    num_ncf = spatial_amr_cache.num_spatial_hanging_facets

    if num_ncf == 0
        @rankinfo rank "No spatial hanging facets found, skipping constraint matrix construction"
        return spatial_amr_cache
    end

    nelem = length(extra_meshes_coords)
    @rankinfo rank "Processing $(num_ncf) spatial non-conforming facets..."
    @rankinfo rank "Spatial elements: $nelem, each with angular mesh"
    @rankinfo rank "BUILDING: Full spatial-angular constraint system"

    # Cache for interpolation matrices to avoid recomputation
    interp_cache = Dict{Tuple{Int,Int,Int}, Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}()

    # Track statistics
    num_spatial_hanging = 0
    num_spatial_angular_hanging = 0

    # Process each spatial non-conforming facet
    for ncf_idx = 1:num_ncf
        # Get child and parent element IDs
        child_elem_id = Int(mesh.cip[ncf_idx])
        parent_elem_id = Int(mesh.pip[ncf_idx])

        # Skip if this facet has no hanging nodes on child side
        # (can happen with ghost parent facets)
        if child_elem_id == 0 || parent_elem_id == 0
            continue
        end

        # Get facet information
        local_facet_id = Int(mesh.lfid[ncf_idx])

        try
            # Get node indices on this facet
            IPc_nodes = mesh.IPc_list[:, ncf_idx]  # Child face nodes
            IPp_nodes = mesh.IPp_list[:, ncf_idx]  # Parent face nodes

            # Build interpolation matrices from actual facet node coordinates
            L1, L2= build_spatial_interpolation_matrices_from_facet_nodes(
                child_elem_id, parent_elem_id, local_facet_id,
                IPc_nodes, IPp_nodes, mesh, ngl
            )

            # Build constraint equations for each hanging node on child facet
            # INCLUDING all associated angular nodes
            num_spatial_added = build_hanging_node_constraints_with_angular(
                child_elem_id, parent_elem_id, local_facet_id,
                IPc_nodes, IPp_nodes, L1, L2,
                spatial_amr_cache, mesh, ngl,
                extra_meshes_extra_nelems, extra_meshes_extra_nops, rank
            )

            num_spatial_hanging += num_spatial_added
            # Each spatial node couples with all angular nodes
            nelem_ang_child = extra_meshes_extra_nelems[child_elem_id]
            nop_child = extra_meshes_extra_nops[child_elem_id][1]  # Assume uniform in first element
            num_spatial_angular_hanging += num_spatial_added * nelem_ang_child * (nop_child+1) * (nop_child+1)

        catch err
            @rankinfo rank "Error processing facet $ncf_idx (child=$child_elem_id, parent=$parent_elem_id): $err"
            rethrow(err)
        end
    end

    @rankinfo rank "✓ Spatial constraint matrices built"
    @rankinfo rank "  - $(num_spatial_hanging) spatial hanging nodes"
    @rankinfo rank "  - $(num_spatial_angular_hanging) spatial-angular hanging DOFs (including all angular elements)"
    return spatial_amr_cache
end

# ============================================================================
# Build Lagrange Interpolation Matrices for One Spatial Facet
# ============================================================================

"""
    build_spatial_interpolation_matrices_from_facet_nodes(
        child_elem_id, parent_elem_id, local_facet_id,
        IPc_nodes, IPp_nodes, mesh, ngl
    ) -> (L1, L2)

Build 1D Lagrange interpolation matrices for the TWO varying directions on a facet.

CRITICAL: Extract coordinates from IPc_nodes and IPp_nodes arrays (already structured as ngl×ngl grids).

Extract 1D lines along each varying direction:
- Line in direction 1: take nodes with fixed j index (j=1), vary i from 1:ngl
- Line in direction 2: take nodes with fixed i index (i=1), vary j from 1:ngl

Each L matrix is ngl×ngl mapping one 1D line to corresponding 1D line.

Tensor product: weight[child(i,j), parent(p1,p2)] = L1[i,p1] * L2[j,p2]

This ensures correct high-order interpolation while maintaining partition of unity.
"""
function build_spatial_interpolation_matrices_from_facet_nodes(
    child_elem_id::Int, parent_elem_id::Int, local_facet_id::Int,
    IPc_nodes::Vector, IPp_nodes::Vector,
    mesh::St_mesh, ngl::Int
)
    # IPc_nodes and IPp_nodes are ngl² arrays representing ngl×ngl grids on the facet
    # Grid position (i,j) maps to linear index: (i-1)*ngl + j

    # Extract 1D lines along direction 1 (i-direction, j=1)
    # and direction 2 (j-direction, i=1)

    parent_line_1 = Vector{Int}()  # i-direction: j=1 (fixed)
    parent_line_2 = Vector{Int}()  # j-direction: i=1 (fixed)
    child_line_1 = Vector{Int}()
    child_line_2 = Vector{Int}()

    for i = 1:ngl
        # Line 1: (i, 1) for i in 1:ngl
        idx = (i - 1) * ngl + 1
        if idx <= length(IPp_nodes)
            push!(parent_line_1, Int(IPp_nodes[idx]))
            push!(child_line_1, Int(IPc_nodes[idx]))
        end
    end

    for j = 1:ngl
        # Line 2: (1, j) for j in 1:ngl
        idx = 0 * ngl + j  # i=1
        if idx <= length(IPp_nodes)
            push!(parent_line_2, Int(IPp_nodes[idx]))
            push!(child_line_2, Int(IPc_nodes[idx]))
        end
    end

    if isempty(parent_line_1) || isempty(parent_line_2) ||
       isempty(child_line_1) || isempty(child_line_2)
        error("Failed to extract 1D node lines from facet (IPp_nodes length=$(length(IPp_nodes)), IPc_nodes length=$(length(IPc_nodes)))")
    end

    # Now extract the appropriate coordinate from each direction based on facet type
    if local_facet_id in [1, 2]
        # z-faces: x and y vary, z constant
        # Line 1 varies in x, Line 2 varies in y
        parent_x = [mesh.x[n] for n in parent_line_2 if n > 0]
        child_x = [mesh.x[n] for n in child_line_2 if n > 0]
        parent_z = [mesh.z[n] for n in parent_line_1 if n > 0]
        child_z = [mesh.z[n] for n in child_line_1 if n > 0]
        
        L1 = build_1d_lagrange_matrix(parent_x, child_x)  # x-direction
        L2 = build_1d_lagrange_matrix(parent_z, child_z)  # y-direction
    elseif local_facet_id in [3, 4]
        # y-faces: x and z vary, y constant
        # Line 1 varies in x, Line 2 varies in z
        parent_x = [mesh.x[n] for n in parent_line_2 if n > 0]
        child_x = [mesh.x[n] for n in child_line_2 if n > 0]
        parent_y = [mesh.y[n] for n in parent_line_1 if n > 0]
        child_y = [mesh.y[n] for n in child_line_1 if n > 0]
        
        L1 = build_1d_lagrange_matrix(parent_x, child_x)  # x-direction
        L2 = build_1d_lagrange_matrix(parent_y, child_y)  # z-direction
    elseif local_facet_id in [5, 6]
        # x-faces: y and z vary, x constant
        # Line 1 varies in y, Line 2 varies in z
        parent_y = [mesh.y[n] for n in parent_line_2 if n > 0]
        child_y = [mesh.y[n] for n in child_line_2 if n > 0]
        parent_z = [mesh.z[n] for n in parent_line_1 if n > 0]
        child_z = [mesh.z[n] for n in child_line_1 if n > 0]
        
        L1 = build_1d_lagrange_matrix(parent_y, child_y)  # y-direction
        L2 = build_1d_lagrange_matrix(parent_z, child_z)  # z-direction
    else
        error("Invalid local_facet_id: $local_facet_id")
    end

    return L1, L2
end

# ============================================================================
# Coord-Cache Variant of Interpolation Matrix Builder (for cross-rank NCFs)
# ============================================================================

"""
    build_spatial_interpolation_matrices_with_coord_cache(
        local_facet_id, IPc_nodes, parent_global_ips, mesh, coord_cache, ngl
    ) -> (L1, L2)

Like `build_spatial_interpolation_matrices_from_facet_nodes` but uses a
`coord_cache::Dict{Int,NTuple{3,Float64}}` (keyed by GLOBAL IP) for parent
node coordinates instead of `mesh.x/y/z[local_ip]`.

Child node coordinates still come from the local mesh (child is always local).
"""
function build_spatial_interpolation_matrices_with_coord_cache(
    local_facet_id::Int,
    IPc_nodes::Vector{<:Integer},
    parent_global_ips::Vector{<:Integer},
    mesh::St_mesh,
    coord_cache::Dict{Int, NTuple{3, Float64}},
    ngl::Int
)
    # Extract 1D lines: line1 = i-direction (j=1), line2 = j-direction (i=1)
    parent_line_1 = Int[]
    parent_line_2 = Int[]
    child_line_1  = Int[]
    child_line_2  = Int[]

    for i = 1:ngl
        idx = (i - 1) * ngl + 1      # (i, j=1)
        idx <= length(parent_global_ips) || continue
        push!(parent_line_1, Int(parent_global_ips[idx]))
        push!(child_line_1,  Int(IPc_nodes[idx]))
    end
    for j = 1:ngl
        idx = j                       # (i=1, j)
        idx <= length(parent_global_ips) || continue
        push!(parent_line_2, Int(parent_global_ips[idx]))
        push!(child_line_2,  Int(IPc_nodes[idx]))
    end

    if isempty(parent_line_1) || isempty(parent_line_2)
        error("Failed to extract 1D node lines from parent_global_ips (length=$(length(parent_global_ips)))")
    end

    get_coord(gip::Int, dim::Int) = coord_cache[gip][dim]
    get_child(lip::Int, dim::Int) = dim == 1 ? mesh.x[lip] : (dim == 2 ? mesh.y[lip] : mesh.z[lip])

    if local_facet_id in [1, 2]       # z-faces: x and y vary
        parent_x = [get_coord(g, 1) for g in parent_line_2 if g > 0]
        child_x  = [get_child(c, 1) for c in child_line_2  if c > 0]
        parent_z = [get_coord(g, 3) for g in parent_line_1 if g > 0]
        child_z  = [get_child(c, 3) for c in child_line_1  if c > 0]
        L1 = build_1d_lagrange_matrix(parent_x, child_x)
        L2 = build_1d_lagrange_matrix(parent_z, child_z)
    elseif local_facet_id in [3, 4]   # y-faces: x and z vary
        parent_x = [get_coord(g, 1) for g in parent_line_2 if g > 0]
        child_x  = [get_child(c, 1) for c in child_line_2  if c > 0]
        parent_y = [get_coord(g, 2) for g in parent_line_1 if g > 0]
        child_y  = [get_child(c, 2) for c in child_line_1  if c > 0]
        L1 = build_1d_lagrange_matrix(parent_x, child_x)
        L2 = build_1d_lagrange_matrix(parent_y, child_y)
    elseif local_facet_id in [5, 6]   # x-faces: y and z vary
        parent_y = [get_coord(g, 2) for g in parent_line_2 if g > 0]
        child_y  = [get_child(c, 2) for c in child_line_2  if c > 0]
        parent_z = [get_coord(g, 3) for g in parent_line_1 if g > 0]
        child_z  = [get_child(c, 3) for c in child_line_1  if c > 0]
        L1 = build_1d_lagrange_matrix(parent_y, child_y)
        L2 = build_1d_lagrange_matrix(parent_z, child_z)
    else
        error("Invalid local_facet_id: $local_facet_id")
    end

    return L1, L2
end

# ============================================================================
# Cross-Rank Hanging Node Constraint Builder
# ============================================================================

"""
    build_cross_rank_hanging_node_constraints(
        IPc_nodes, parent_global_ips, L1, L2,
        cache, mesh, ngl,
        extra_meshes_extra_nelems, extra_meshes_extra_nops, rank
    ) -> Int

Build constraint equations for hanging nodes whose parents are on another rank.

Stores entries in `cache.cross_rank_parent_weights` using GLOBAL spatial IPs
for the parent nodes (unlike `parent_weights` which uses local IPs).

Returns the number of new hanging nodes added.
"""
function build_cross_rank_hanging_node_constraints(
    IPc_nodes::Vector{<:Integer},
    parent_global_ips::Vector{<:Integer},
    L1::Matrix{Float64}, L2::Matrix{Float64},
    cache::SpatialAMRCache, mesh::St_mesh, ngl::Int,
    extra_meshes_extra_nelems, extra_meshes_extra_nops, rank::Int
)
    num_added = 0

    for node_idx = 1:length(IPc_nodes)
        child_node_id = Int(IPc_nodes[node_idx])
        child_node_id <= 0 && continue

        i = div(node_idx - 1, ngl) + 1
        j = mod(node_idx - 1, ngl) + 1
        (i < 1 || i > size(L1, 1) || j < 1 || j > size(L2, 1)) && continue

        w_1 = L1[j, :]
        w_2 = L2[i, :]

        parent_weights_raw = Tuple{Int, Float64}[]
        for p1 = 1:ngl
            for p2 = 1:ngl
                parent_face_idx = (p2 - 1) * ngl + p1
                parent_face_idx <= length(parent_global_ips) || continue
                parent_gip = Int(parent_global_ips[parent_face_idx])
                parent_gip <= 0 && continue
                weight = w_1[p1] * w_2[p2]
                abs(weight) > 1e-14 || continue
                push!(parent_weights_raw, (parent_gip, weight))
            end
        end

        isempty(parent_weights_raw) && continue

        # Skip if already constrained by a previous NCF face.
        if haskey(cache.cross_rank_parent_weights, child_node_id) ||
           haskey(cache.parent_weights, child_node_id)           ||
           haskey(cache.coincident_nodes, child_node_id)         ||
           haskey(cache.cross_rank_coincident_nodes, child_node_id)
            continue
        end

        # Coincident cross-rank: single parent with weight ≈ 1.
        # Track separately — zero child's row/RHS, copy parent solution after solve.
        if length(parent_weights_raw) == 1 && abs(parent_weights_raw[1][2] - 1.0) < 1e-14
            parent_gip = parent_weights_raw[1][1]
            cache.cross_rank_coincident_nodes[child_node_id] = parent_gip
            num_added += 1
            continue
        end

        # Store: child_node_id → [(parent_GLOBAL_ip, weight), ...]
        cache.cross_rank_parent_weights[child_node_id] = parent_weights_raw
        num_added += 1
    end

    return num_added
end

# ============================================================================
# Build 1D Lagrange Interpolation Matrix
# ============================================================================

"""
    build_1d_lagrange_matrix(x_parent::Vector, x_child::Vector) -> Matrix

Build 1D Lagrange interpolation matrix.

Given parent and child node coordinate arrays in one dimension:
- Compute barycentric weights for parent nodes
- Build Lagrange interpolation matrix from parent to child
- Normalize rows to satisfy partition of unity

Returns (num_child × num_parent) matrix where:
- L[i,j] = weight of parent node j in constraint for child node i
"""
function build_1d_lagrange_matrix(x_parent::Vector{Float64}, x_child::Vector{Float64})
    n_parent = length(x_parent)
    n_child = length(x_child)

    # Initialize matrix
    L = zeros(Float64, n_child, n_parent)
    # Compute barycentric weights for parent nodes
    ω = zeros(Float64, n_parent)
    BarycentricWeights!(x_parent, ω)
    
    # Build interpolation matrix: columns are parent nodes, rows are child nodes
    PolynomialInterpolationMatrix!(x_parent, ω, x_child, L)

    # CRITICAL: Do NOT normalize rows!
    # Lanrange interpolation already satisfies partition of unity (row sums = 1).
    # Normalizing breaks column sums (each parent's total contribution).
    # Trust the Lagrange formula.

    return L
end

# ============================================================================
# Build Hanging Node Constraints
# ============================================================================

"""
    build_hanging_node_constraints_with_angular(
        child_elem_id, parent_elem_id, local_facet_id,
        IPc_nodes, IPp_nodes, Lx, Ly, Lz,
        cache, mesh, ngl,
        extra_meshes_extra_nelems, extra_meshes_extra_nops, rank
    ) -> Int

Build constraint equations for spatial hanging nodes WITH angular coupling.

CRITICAL: Each spatial hanging node constraint applies to ALL associated angular nodes.

For each spatial hanging node at location (i,j,k):
- Gets spatial constraint weights: u_spatial_hanging = Σ w_s * u_spatial_parent
- For each angular element e in 1:nelem_ang:
  - For each angular node (e, iθ, iφ):
    - Same spatial constraint applies to all angular DOFs at that location

Returns: Number of spatial hanging nodes processed
"""
function build_hanging_node_constraints_with_angular(
    child_elem_id::Int, parent_elem_id::Int, local_facet_id::Int,
    IPc_nodes::Vector, IPp_nodes::Vector,
    L1::Matrix{Float64}, L2::Matrix{Float64},
    cache::SpatialAMRCache, mesh::St_mesh, ngl::Int,
    extra_meshes_extra_nelems, extra_meshes_extra_nops, rank::Int
)
    num_spatial_hanging = 0

    # Process ALL ngl² nodes on the facet (not just unique ones)
    for node_idx = 1:length(IPc_nodes)
        child_node_id = Int(IPc_nodes[node_idx])
        parent_node_id = Int(IPp_nodes[node_idx])

        # Skip invalid nodes
        if child_node_id <= 0
            continue
        end

        # Convert linear index to (i,j) coordinates on the ngl×ngl facet
        # Node ordering in mesh: column-major
        i = div(node_idx - 1, ngl) + 1  # Row (0-indexed -> 1-indexed)
        j = mod(node_idx - 1, ngl) + 1  # Column (0-indexed -> 1-indexed)

        # Check bounds
        if i < 1 || i > size(L1, 1) || j < 1 || j > size(L2, 1)
            continue
        end

        # Get spatial interpolation weights from correct dimensions
        # L1[i, :] gives interpolation weight for child position i from all parent positions
        # L2[j, :] gives interpolation weight for child position j from all parent positions
        w_1 = L1[j, :]  # Weights for dimension 1
        w_2 = L2[i, :]  # Weights for dimension 2

        # Build spatial constraint: tensor product of 1D weights
        # For each parent node position (p1, p2) on the parent face:
        #   weight = w_1[p1] * w_2[p2]
        parent_weights_raw = Tuple{Int, Float64}[]
        for p1 = 1:ngl
            for p2 = 1:ngl
                parent_face_idx = (p2 - 1) * ngl + p1
                if parent_face_idx <= length(IPp_nodes)
                    parent_node_id = Int(IPp_nodes[parent_face_idx])
                    if parent_node_id > 0
                        weight = w_1[p1] * w_2[p2]
                        if abs(weight) > 1e-14
                            push!(parent_weights_raw, (parent_node_id, weight))
                        end
                    end
                end
            end
        end

        # If no weights found, use uniform interpolation
        if isempty(parent_weights_raw)
            unique_parents = unique(filter(x -> x > 0, IPp_nodes))
            uniform_weight = 1.0 / max(1, length(unique_parents))
            for p in unique_parents
                push!(parent_weights_raw, (p, uniform_weight))
            end
        end

        # CRITICAL: Do NOT normalize weights!
        # Lagrange interpolation weights already satisfy partition of unity.
        # Normalizing here breaks the column sum property (each parent's total contribution = 1).
        # Just use the raw weights from tensor product of 1D Lagrange matrices.

        if isempty(parent_weights_raw)
            # Fallback only if absolutely no weights computed
            unique_parents = unique(filter(x -> x > 0, IPp_nodes))
            uniform_weight = 1.0 / max(1, length(unique_parents))
            spatial_weights = [(p, uniform_weight) for p in unique_parents]
        else
            # Use weights directly from Lagrange tensor product
            spatial_weights = parent_weights_raw
        end

        # Classify: coincident (single parent, weight=1) vs truly-interpolated hanging node.
        # Coincident nodes sit at exactly the same location as a parent node.
        # They must NOT be assembled as free nodes (causes double-counting with the parent)
        # and their contribution must NOT be accumulated into the parent (would double-count
        # the parent's own assembled value). Instead, track them separately: zero their
        # assembled row (identity equation, zero RHS) and copy parent solution after solve.
        if length(spatial_weights) == 1 && abs(spatial_weights[1][2] - 1.0) < 1e-14
            parent_ip = spatial_weights[1][1]
            if !haskey(cache.coincident_nodes, child_node_id)
                cache.coincident_nodes[child_node_id] = parent_ip
                num_spatial_hanging += 1
            end
        elseif !haskey(cache.parent_weights, child_node_id)
            cache.parent_weights[child_node_id] = spatial_weights
            num_spatial_hanging += 1
        end
    end

    return num_spatial_hanging
end

# ============================================================================
# Verification of Constraint Matrices
# ============================================================================

"""
    verify_spatial_constraints(cache::SpatialAMRCache, rank::Int, n_spatial_nodes::Int,
                              extra_meshes_extra_nelems=nothing, extra_meshes_extra_nops=nothing)

Verify that spatial constraint matrices satisfy partition of unity by checking column sums of restriction matrix.

Builds the full restriction matrix M (n_spatial_nodes × n_spatial_nodes) with:
- Identity rows for free nodes
- Interpolation weight rows for hanging nodes

Checks: Each column sums to 1.0 (partition of unity - each parent conserves its contribution)

Parameters:
- cache: SpatialAMRCache with parent_weights
- rank: MPI rank for logging
- n_spatial_nodes: Total spatial nodes (free + hanging)
- extra_meshes_extra_nelems: Number of angular elements per spatial element (optional)
- extra_meshes_extra_nops: Angular polynomial orders per element (optional)
"""
function verify_spatial_constraints(cache::SpatialAMRCache, rank::Int, n_spatial_nodes::Int,
                                    extra_meshes_extra_nelems=nothing,
                                    extra_meshes_extra_nops=nothing)
    if length(cache.parent_weights) == 0
        @rankinfo rank "WARNING: No constraint weights found in cache"
        return
    end

    num_hanging = length(cache.parent_weights)

    # Build full n_spatial_nodes × n_spatial_nodes matrix with:
    # - Identity rows for free nodes
    # - Interpolation weight rows for hanging nodes
    M_full = sparse(1.0I(n_spatial_nodes))  # Start with identity

    # For each hanging node, replace its row with interpolation weights
    for (hanging_node_id, weights) in cache.parent_weights
        # Clear the identity row for this hanging node
        for j in 1:n_spatial_nodes
            M_full[hanging_node_id, j] = 0.0
        end
        # Fill in the parent weights
        for (parent_id, weight) in weights
            if parent_id > 0 && parent_id <= n_spatial_nodes
                M_full[parent_id, hanging_node_id] = weight
            end
        end
    end

    # Check column sums (partition of unity)
    num_violated_col_unity = 0
    for col in 1:n_spatial_nodes
        col_sum = sum(M_full[:, col])
        if abs(col_sum - 1.0) > 1e-10
            num_violated_col_unity += 1
            if num_violated_col_unity <= 5  # Print first 5 violations
                @rankinfo rank "WARNING: Column $col sum = $col_sum (should be 1.0)"
            end
        end
    end

    @rankinfo rank "Spatial Constraint Verification Summary:"
    @rankinfo rank "  - Spatial hanging nodes: $num_hanging"
    @rankinfo rank "  - Restriction matrix size: $n_spatial_nodes × $n_spatial_nodes"
    @rankinfo rank "  - Columns with sums ≠ 1.0: $num_violated_col_unity"

    # If angular mesh provided, verify spatial-angular coupling awareness
    if extra_meshes_extra_nelems !== nothing && extra_meshes_extra_nops !== nothing
        total_angular_nelems = sum(extra_meshes_extra_nelems)
        avg_nop = Int(round(mean([nops[1] for nops in extra_meshes_extra_nops])))

        expected_spatial_angular = num_hanging * total_angular_nelems * (avg_nop + 1) * (avg_nop + 1)
        @rankinfo rank "  - Expected spatial-angular DOFs (with coupling): ~$expected_spatial_angular"
        @rankinfo rank "  - NOTE: RHS/Matrix assembly will expand these constraints to all angular nodes"
    end

    if num_violated_col_unity > 0
        error("[Rank $rank] Spatial constraint verification FAILED: " *
              "$num_violated_col_unity column(s) have sum ≠ 1.0 (partition of unity violated)")
    end

    @rankinfo rank "✓ Spatial constraints verified: partition of unity satisfied"
end

# ============================================================================
# Build Sparse Constraint Matrices for Assembly
# ============================================================================

"""
    build_spatial_restriction_and_prolongation(
        spatial_amr_cache::SpatialAMRCache,
        n_spatial_nodes::Int,
        n_angular_nodes::Int
    ) -> (SparseMatrix, SparseMatrix)

Build sparse restriction and prolongation matrices from parent_weights with angular coupling.

CRITICAL: Spatial constraints couple with ALL angular DOFs at each location (Kronecker product structure).

Restriction matrix R_spatial structure:
- Size: (n_free*n_angular) × (n_total*n_angular) where n_total = n_spatial, n_free = n_free_spatial
- Layout: [I_free | W] where I_free is identity on free DOFs, W is weight block
- Free DOF rows: identity (pass through to free DOFs)
- Hanging DOF columns: contain interpolation weights from parents
- For each spatial hanging node s with parents {p1, p2, ...} and weights {w1, w2, ...}:
  - For each angular node a ∈ [1, n_angular]:
    - For each parent p:
      - R[(p-1)*n_ang + a, n_free_ang + (s-1)*n_ang + a] = weight_p

Kronecker product form: R = R_spatial ⊗ I_angular (with partial rank structure)

Prolongation matrix P_spatial:
- P_spatial = R_spatial' (transpose of restriction)
- Size: (n_total*n_angular) × (n_free*n_angular)

Returns: (R_spatial, P_spatial)
- R_spatial is (n_free*n_angular) × (n_total*n_angular)
- P_spatial is (n_total*n_angular) × (n_free*n_angular)
- Matrix product: A_free = R * A * P

For no spatial hanging nodes, returns (I, I) of size (n_spatial*n_angular).
"""
function build_spatial_restriction_and_prolongation(
    spatial_amr_cache::SpatialAMRCache,
    n_spa_new::Int,
    spatial_hanging_nodes_all_angular::Set{Int},
    point_dict_spa::Dict{NTuple{5,Float64}, Int},
    mesh,
    xyz_ang_map::Dict{NTuple{3,Float64}, Vector{Tuple{Float64,Float64,Int}}}
)::Tuple{SparseMatrixCSC, SparseMatrixCSC}

    if isempty(spatial_amr_cache.parent_weights)
        return sparse(I, n_spa_new, n_spa_new),
               sparse(I, n_spa_new, n_spa_new)
    end

    rows_M = Int[]
    cols_M = Int[]
    vals_M = Float64[]

    # Phase 1: identity rows for free DOFs
    for dof_idx = 1:n_spa_new
        if !(dof_idx in spatial_hanging_nodes_all_angular)
            push!(rows_M, dof_idx)
            push!(cols_M, dof_idx)
            push!(vals_M, 1.0)
        end
    end

    # Phase 2: constraint weight rows for interpolated hanging DOFs.
    # Use xyz_ang_map (keyed by rounded (x,y,z)) to iterate over ALL angular DOFs
    # present at the hanging location — including newly-refined ones from angular
    # AMR that are absent from the base extra_mesh.
    for (hanging_node_id, parent_list) in spatial_amr_cache.parent_weights
        x_child = mesh.x[hanging_node_id]
        y_child = mesh.y[hanging_node_id]
        z_child = mesh.z[hanging_node_id]
        xyz_key_h = (round(x_child, digits=12), round(y_child, digits=12), round(z_child, digits=12))
        ang_dofs_h = get(xyz_ang_map, xyz_key_h, Tuple{Float64,Float64,Int}[])
        for (θ, ϕ, col_hanging) in ang_dofs_h
            col_hanging > 0 || continue
            for (parent_node_id, weight) in parent_list
                x_p = mesh.x[parent_node_id]; y_p = mesh.y[parent_node_id]; z_p = mesh.z[parent_node_id]
                key_p = (round(x_p, digits=12), round(y_p, digits=12),
                         round(z_p, digits=12), round(θ, digits=12), round(ϕ, digits=12))
                row_parent = get(point_dict_spa, key_p, 0)
                row_parent == 0 && continue
                push!(rows_M, row_parent)
                push!(cols_M, col_hanging)
                push!(vals_M, weight)
            end
        end
    end

    # Phase 3: cross-rank parents that are locally present.
    gip_to_local_r = Dict{Int,Int}(Int(mesh.ip2gip[ip]) => ip for ip = 1:mesh.npoin)
    for (hanging_node_id, cross_weights) in spatial_amr_cache.cross_rank_parent_weights
        x_child = mesh.x[hanging_node_id]
        y_child = mesh.y[hanging_node_id]
        z_child = mesh.z[hanging_node_id]
        xyz_key_h = (round(x_child, digits=12), round(y_child, digits=12), round(z_child, digits=12))
        ang_dofs_h = get(xyz_ang_map, xyz_key_h, Tuple{Float64,Float64,Int}[])
        for (θ, ϕ, col_hanging) in ang_dofs_h
            col_hanging > 0 || continue
            for (parent_gip, weight) in cross_weights
                local_ip = get(gip_to_local_r, parent_gip, 0)
                local_ip == 0 && continue
                x_p = mesh.x[local_ip]; y_p = mesh.y[local_ip]; z_p = mesh.z[local_ip]
                key_p = (round(x_p, digits=12), round(y_p, digits=12),
                         round(z_p, digits=12), round(θ, digits=12), round(ϕ, digits=12))
                row_parent = get(point_dict_spa, key_p, 0)
                row_parent == 0 && continue
                push!(rows_M, row_parent)
                push!(cols_M, col_hanging)
                push!(vals_M, weight)
            end
        end
    end

    R_spatial = sparse(rows_M, cols_M, vals_M, n_spa_new, n_spa_new)
    P_spatial = R_spatial'
    return R_spatial, P_spatial
end

# =============================================================================
# MPI exchange: combined DOF GIDs for cross-rank spatial parents
# =============================================================================

"""
    _exchange_spatial_parent_combined_gids(...)

Allgatherv-based exchange that mirrors the uniform angular path's `global_sp_ang_to_gid`
construction (build_rad_3d.jl:1060-1103).

Round 1: every rank broadcasts the set of cross-rank parent spatial GIPs it needs so
that ALL ranks know the full needed set.

Round 2: every rank contributes `(parent_spatial_gip, angular_index_k, compact_gid)`
triples for any needed GIP present in its local spatial mesh, then ALL ranks receive
all triples and build the result locally.

This avoids the point-to-point owner-lookup approach, which silently drops parents
whose owner is not recorded in `mesh.pgip_ghost` (causing wrong results in parallel).

Returns `Dict{Int,Vector{Int}}`: parent_global_spatial_ip → [combined_gid_1,...,combined_gid_n_ang_base].
"""
function _exchange_spatial_parent_combined_gids(
    spatial_amr_cache::SpatialAMRCache,
    ip2gip_spa::Vector{Int},
    xyz_ang_map::Dict{NTuple{3,Float64}, Vector{Tuple{Float64,Float64,Int}}},
    mesh::St_mesh,
    rank::Int, comm::MPI.Comm
)
    # Round 1: Allgather all needed cross-rank parent GIPs across ranks.
    local_needed = Int[]
    seen_needed  = Set{Int}()
    for (_, weights) in spatial_amr_cache.cross_rank_parent_weights
        for (parent_gip, _) in weights
            parent_gip in seen_needed && continue
            push!(seen_needed, parent_gip)
            push!(local_needed, parent_gip)
        end
    end
    nn         = Int32(length(local_needed))
    nn_all     = MPI.Allgather([nn], comm)
    all_needed = MPI.Allgatherv(local_needed, nn_all, comm)
    all_needed_set = Set(all_needed)

    # Round 2: for every needed GIP locally visible, contribute quadruples
    # (parent_gip, θ, ϕ, compact_combined_gid).  Using xyz_ang_map covers ALL
    # adapted angular DOFs at the parent location, not just the base-mesh ones.
    local_quads = Float64[]
    seen_spa_ang = Set{Tuple{Int,Float64,Float64}}()
    for lip = 1:mesh.npoin
        sp_gip = Int(mesh.ip2gip[lip])
        sp_gip in all_needed_set || continue
        xyz_key = (round(mesh.x[lip], digits=12),
                   round(mesh.y[lip], digits=12),
                   round(mesh.z[lip], digits=12))
        ang_dofs = get(xyz_ang_map, xyz_key, Tuple{Float64,Float64,Int}[])
        for (θ_k, ϕ_k, ip_comb) in ang_dofs
            key_t = (sp_gip, θ_k, ϕ_k)
            key_t in seen_spa_ang && continue
            push!(seen_spa_ang, key_t)
            push!(local_quads, Float64(sp_gip), θ_k, ϕ_k, Float64(ip2gip_spa[ip_comb]))
        end
    end
    nq      = Int32(length(local_quads) ÷ 4)
    nq_all  = MPI.Allgather([nq], comm)
    all_q   = MPI.Allgatherv(local_quads, Int32.(nq_all .* 4), comm)

    # Build result: (parent_gip, (θ, ϕ)) → compact_combined_gid
    result = Dict{Tuple{Int, NTuple{2,Float64}}, Int}()
    for k = 1:length(all_q) ÷ 4
        gip = Int(round(all_q[4k-3]))
        θ_k = all_q[4k-2]
        ϕ_k = all_q[4k-1]
        gid = Int(round(all_q[4k]))
        result[(gip, (θ_k, ϕ_k))] = gid
    end
    return result
end

# =============================================================================
# Combine angular and spatial restriction operators into one nc_mat
# =============================================================================

"""
    combine_spatial_angular_restrictions(
        nc_mat, nc_mat_rhs,
        ghost_constraint_data, ghost_constraint_data_rhs,
        gid_to_extended_parents, extended_parents_to_gid,
        extended_parents_x, extended_parents_y, extended_parents_z,
        extended_parents_θ, extended_parents_ϕ, extended_parents_ip,
        all_hanging_nodes,
        R_spatial_ext, R_spatial_rhs_ext,
        ghost_constraint_data_spa, ghost_constraint_data_spa_rhs,
        extended_parents_to_gid_spa,
        spatial_hanging_nodes,
        n_spa, rank
    )

Combine an already-built angular restriction matrix `nc_mat` with already-built
spatial constraint outputs (from `build_spatial_constraints_for_combined_path`).

The spatial outputs use extended-parent indices starting at `n_spa + 1`.  This
function re-indexes them to start at `n_spa + n_ang_ext + 1` (after the existing
angular extended parents) before assembling the combined matrix.

Mathematical identity (disjoint hanging-node sets by `can_refine_angular` mask):

    nc_mat_combined[i, j] = nc_mat[i, j]
                           + R_spatial_weight_only[i, j]     (hanging rows, off-diag)
                           - δ_{ij} for i in spatial_hanging  (remove angular identity rows)

Mutates in-place: `ghost_constraint_data`, `ghost_constraint_data_rhs`,
`gid_to_extended_parents`, `extended_parents_to_gid`, and the six
`extended_parents_*` coordinate/index vectors.
"""
function combine_spatial_angular_restrictions(
    nc_mat::SparseMatrixCSC,
    nc_mat_rhs::SparseMatrixCSC,
    ghost_constraint_data::Dict{Int,Vector{Tuple{Int,Float64}}},
    ghost_constraint_data_rhs::Dict{Int,Vector{Tuple{Int,Float64}}},
    gid_to_extended_parents::Dict{Int,Int},
    extended_parents_to_gid::Vector{Int},
    extended_parents_x::Vector{Float64},
    extended_parents_y::Vector{Float64},
    extended_parents_z::Vector{Float64},
    extended_parents_θ::Vector{Float64},
    extended_parents_ϕ::Vector{Float64},
    extended_parents_ip::Vector{Int},
    all_hanging_nodes::Set{Int},
    R_spatial_ext::SparseMatrixCSC,
    R_spatial_rhs_ext::SparseMatrixCSC,
    ghost_constraint_data_spa::Dict{Int,Vector{Tuple{Int,Float64}}},
    ghost_constraint_data_spa_rhs::Dict{Int,Vector{Tuple{Int,Float64}}},
    extended_parents_to_gid_spa::Vector{Int},
    spatial_hanging_nodes::Set{Int},
    n_spa::Int,
    rank::Int
)
    n_ang_ext = length(extended_parents_to_gid)   # angular extended parents already present
    n_spa_ext = length(extended_parents_to_gid_spa)
    n_new     = n_spa + n_ang_ext + n_spa_ext

    # ── 1. Merge ghost constraint dicts (spatial entries into angular ones) ───
    for (k, v) in ghost_constraint_data_spa
        existing = get!(ghost_constraint_data, k, Tuple{Int,Float64}[])
        for entry in v
            any(g == entry[1] for (g, _) in existing) || push!(existing, entry)
        end
    end
    for (k, v) in ghost_constraint_data_spa_rhs
        existing = get!(ghost_constraint_data_rhs, k, Tuple{Int,Float64}[])
        for entry in v
            any(g == entry[1] for (g, _) in existing) || push!(existing, entry)
        end
    end

    # ── 2. Merge extended parents: spatial ones go after angular, with offset ─
    # Spatial extended parents from build_spatial_constraints_for_combined_path
    # have indices n_spa+1..n_spa+n_spa_ext.  Re-index them to
    # n_spa+n_ang_ext+1..n_new so they don't collide with angular ext parents.
    for (k_idx, gid) in enumerate(extended_parents_to_gid_spa)
        new_ext_idx = n_spa + n_ang_ext + k_idx
        gid_to_extended_parents[gid] = new_ext_idx
        push!(extended_parents_to_gid, gid)
        # Spatial ext parents carry no angular coordinate data; use ip=0 so the
        # ghost-BC block (guarded by !isempty(extended_parents_ip)) can safely
        # iterate over them without attempting to apply BC logic.
        push!(extended_parents_ip, 0)
        push!(extended_parents_x,  0.0)
        push!(extended_parents_y,  0.0)
        push!(extended_parents_z,  0.0)
        push!(extended_parents_θ,  0.0)
        push!(extended_parents_ϕ,  0.0)
    end

    # ── 3. Extract spatial triplets (filter free-DOF identities; shift ext indices) ─
    # R_spatial_ext contains:
    #   identity rows for free DOFs (i == j, i ≤ n_spa, i ∉ spatial_hanging_nodes)
    #   weight rows   for hanging DOFs (off-diagonal, or ext-parent diagonal)
    #   identity rows for spatial extended parents (i == j, i > n_spa)
    # The combining formula already has identity from nc_mat for free DOFs, so we
    # only add weight rows and spatial extended-parent identity rows.
    Is_r, Js_r, Vs_r = findnz(R_spatial_ext)
    I_spa = Int[]; J_spa = Int[]; V_spa = Float64[]
    for k in eachindex(Is_r)
        i, j, v = Is_r[k], Js_r[k], Vs_r[k]
        # Skip free-DOF identity (not in hanging set and within base DOF range)
        (i == j && i <= n_spa && !(i in spatial_hanging_nodes)) && continue
        push!(I_spa, i > n_spa ? i + n_ang_ext : i)
        push!(J_spa, j > n_spa ? j + n_ang_ext : j)
        push!(V_spa, v)
    end

    Is_rr, Js_rr, Vs_rr = findnz(R_spatial_rhs_ext)
    I_spa_rhs = Int[]; J_spa_rhs = Int[]; V_spa_rhs = Float64[]
    for k in eachindex(Is_rr)
        i, j, v = Is_rr[k], Js_rr[k], Vs_rr[k]
        (i == j && i <= n_spa && !(i in spatial_hanging_nodes)) && continue
        push!(I_spa_rhs, i > n_spa ? i + n_ang_ext : i)
        push!(J_spa_rhs, j > n_spa ? j + n_ang_ext : j)
        push!(V_spa_rhs, v)
    end

    # ── 4. Build combined matrices ────────────────────────────────────────────
    # nc_mat has identity rows for DOFs free in angular (including spatial-hanging).
    # Subtracting I_{spatial_hanging} removes those identity entries so the combined
    # row carries only the spatial weights.
    I_nc,     J_nc,     V_nc     = findnz(nc_mat)
    I_nc_rhs, J_nc_rhs, V_nc_rhs = findnz(nc_mat_rhs)

    I_sub = collect(spatial_hanging_nodes)
    J_sub = collect(spatial_hanging_nodes)
    V_sub = fill(-1.0, length(I_sub))

    nc_mat_combined = sparse(
        vcat(I_nc,     I_spa,     I_sub),
        vcat(J_nc,     J_spa,     J_sub),
        vcat(Float64.(V_nc),     V_spa,     V_sub), n_new, n_new)
    nc_mat_rhs_combined = sparse(
        vcat(I_nc_rhs, I_spa_rhs, I_sub),
        vcat(J_nc_rhs, J_spa_rhs, J_sub),
        vcat(Float64.(V_nc_rhs), V_spa_rhs, V_sub), n_new, n_new)

    P_combined     = sparse(nc_mat_combined')
    P_vec_combined = sparse(nc_mat_rhs_combined')

    all_hanging_combined = union(all_hanging_nodes, spatial_hanging_nodes)

    @rankinfo rank "Combined angular+spatial restriction: " *
                   "$(length(spatial_hanging_nodes)) spatial hanging DOFs, " *
                   "$n_ang_ext angular ext + $n_spa_ext spatial ext parents → n_new=$n_new"

    return nc_mat_combined, P_combined, nc_mat_rhs_combined, P_vec_combined,
           all_hanging_combined, spatial_hanging_nodes
end

# =============================================================================
# Spatial constraints for the combined adaptive path (spatial-only case)
# Mirrors the proven uniform-angular-mesh path (build_rad_3d.jl lines 952-1429)
# instead of using combine_spatial_angular_restrictions, which has parallel bugs.
# =============================================================================

"""
    build_spatial_constraints_for_combined_path(...)

Build spatial AMR constraint structures for the combined adaptive path.

Works for both the spatial-only case (`global_max_ref == 0`) and the combined
angular+spatial case (`global_max_ref > 0`).  Uses the same logic as the proven
uniform-angular-mesh path but drives node lookup through `xyz_ang_map` built from
`point_dict_combined_adapted`, so adapted angular DOFs are covered correctly.

Returns:
  (R_spatial_ext, P_spatial_ext, R_spatial_rhs_ext,
   ghost_constraint_data_spa, ghost_constraint_data_spa_rhs,
   gid_to_extended_parents_spa, extended_parents_to_gid_spa,
   gip_to_local_spa_ext, spatial_hanging_nodes_all_angular, n_ext_spa)
"""
function build_spatial_constraints_for_combined_path(
    spatial_amr_cache::SpatialAMRCache,
    mesh,
    ip2gip_spa::Vector{Int},
    gip2owner_extra::Vector{Int},
    gip2owner_spa,
    point_dict_combined_adapted::Dict{NTuple{5,Float64}, Int},
    n_spa::Int,
    rank::Int,
    comm::MPI.Comm
)

    # ── Build (x,y,z) → [(θ,ϕ,combined_DOF)] from point_dict_combined_adapted ──
    # This reverse map covers ALL adapted angular DOFs, including refined ones from
    # angular AMR that are absent from extra_mesh_base.extra_coords.  Both the
    # spatial-only path (global_max_ref==0, no angular adaptation) and the combined
    # path (global_max_ref>0) are handled correctly: for the spatial-only case the
    # map simply reproduces what extra_mesh_base would give.
    xyz_ang_map = Dict{NTuple{3,Float64}, Vector{Tuple{Float64,Float64,Int}}}()
    for (key, ip_comb) in point_dict_combined_adapted
        x, y, z, θ, ϕ = key
        xyz_key = (x, y, z)  # keys are already rounded at point_dict construction
        push!(get!(xyz_ang_map, xyz_key, Tuple{Float64,Float64,Int}[]), (θ, ϕ, ip_comb))
    end

    # ── spatial_hanging_nodes_all_angular ──────────────────────────────────────
    nc_non_global_nodes_spa = Int[]
    nc_set_spa = Set{Int}()
    for hanging_ip in keys(spatial_amr_cache.parent_weights)
        xyz_key = (round(mesh.x[hanging_ip], digits=12),
                   round(mesh.y[hanging_ip], digits=12),
                   round(mesh.z[hanging_ip], digits=12))
        for (_, _, ip_spa_h) in get(xyz_ang_map, xyz_key, Tuple{Float64,Float64,Int}[])
            ip_spa_h in nc_set_spa && continue
            push!(nc_non_global_nodes_spa, ip_spa_h)
            push!(nc_set_spa, ip_spa_h)
        end
    end
    for hanging_ip in keys(spatial_amr_cache.cross_rank_parent_weights)
        xyz_key = (round(mesh.x[hanging_ip], digits=12),
                   round(mesh.y[hanging_ip], digits=12),
                   round(mesh.z[hanging_ip], digits=12))
        for (_, _, ip_spa_h) in get(xyz_ang_map, xyz_key, Tuple{Float64,Float64,Int}[])
            ip_spa_h in nc_set_spa && continue
            push!(nc_non_global_nodes_spa, ip_spa_h)
            push!(nc_set_spa, ip_spa_h)
        end
    end
    spatial_hanging_nodes_all_angular = Set(nc_non_global_nodes_spa)
    @rankinfo rank "$(length(spatial_hanging_nodes_all_angular)) spatial-angular hanging DOFs"

    # ── R_spatial: Phases 1, 2, 3 ─────────────────────────────────────────────
    R_spatial, _ = build_spatial_restriction_and_prolongation(
        spatial_amr_cache, n_spa, spatial_hanging_nodes_all_angular,
        point_dict_combined_adapted, mesh, xyz_ang_map
    )

    # ── global (spatial_GIP, (θ,ϕ)) → compact_GID ────────────────────────────
    # Exchange cross-rank parent combined GIDs using xyz_ang_map so adapted
    # angular DOFs are covered.  Returns Dict{Tuple{Int,NTuple{2,Float64}},Int}.
    global_sp_ang_to_gid = _exchange_spatial_parent_combined_gids(
        spatial_amr_cache, ip2gip_spa, xyz_ang_map, mesh, rank, comm
    )

    # ── gip_to_local_spa ──────────────────────────────────────────────────────
    gip_to_local_spa = Dict{Int,Int}(ip2gip_spa[ip] => ip for ip = 1:n_spa)

    # ── _rhs_handled_dofs: combined DOFs whose RHS is covered by local parents ─
    _rhs_handled_dofs = Set{Int}()
    for (hsp, _) in spatial_amr_cache.parent_weights
        xyz_key = (round(mesh.x[hsp], digits=12),
                   round(mesh.y[hsp], digits=12),
                   round(mesh.z[hsp], digits=12))
        for (_, _, d) in get(xyz_ang_map, xyz_key, Tuple{Float64,Float64,Int}[])
            d == 0 && continue
            push!(_rhs_handled_dofs, d)
        end
    end

    # ── Ghost constraint dicts + Source 1 R_spatial augmentation ─────────────
    ghost_constraint_data_spa     = Dict{Int, Vector{Tuple{Int, Float64}}}()
    ghost_constraint_data_spa_rhs = Dict{Int, Vector{Tuple{Int, Float64}}}()
    gip_to_local_r = Dict{Int,Int}(Int(ip2gip_spa[ip]) => ip for ip = 1:n_spa)

    # Source 1: cross-rank parents — iterate over all adapted angular DOFs at the
    # hanging location via xyz_ang_map; match parent by (θ, ϕ) key.
    for (hanging_local_spa, cross_weights) in spatial_amr_cache.cross_rank_parent_weights
        xyz_key_h = (round(mesh.x[hanging_local_spa], digits=12),
                     round(mesh.y[hanging_local_spa], digits=12),
                     round(mesh.z[hanging_local_spa], digits=12))
        for (θ_k, ϕ_k, hanging_dof) in get(xyz_ang_map, xyz_key_h, Tuple{Float64,Float64,Int}[])
            (hanging_dof == 0 || hanging_dof > n_spa) && continue
            ang_key = (θ_k, ϕ_k)
            parent_constraints = Tuple{Int, Float64}[]
            for (parent_global_spa, weight) in cross_weights
                parent_gid = get(global_sp_ang_to_gid, (parent_global_spa, ang_key), 0)
                parent_gid == 0 && continue
                push!(parent_constraints, (parent_gid, weight))
                local_ip = get(gip_to_local_r, parent_gid, 0)
                local_ip == 0 && continue
                R_spatial[local_ip, hanging_dof] = weight
            end
            isempty(parent_constraints) && continue
            ghost_constraint_data_spa[hanging_dof] = copy(parent_constraints)
            if !(hanging_dof in _rhs_handled_dofs)
                ghost_constraint_data_spa_rhs[hanging_dof] = copy(parent_constraints)
            end
        end
    end
    P_spatial = sparse(R_spatial')

    # Source 2: local parents split by ownership
    I_rhs_sp = Int[]; J_rhs_sp = Int[]; V_rhs_sp = Float64[]
    for dof = 1:n_spa
        if !(dof in spatial_hanging_nodes_all_angular)
            push!(I_rhs_sp, dof); push!(J_rhs_sp, dof); push!(V_rhs_sp, 1.0)
        end
    end
    for (hanging_local_spa, local_weights) in spatial_amr_cache.parent_weights
        xyz_key_h = (round(mesh.x[hanging_local_spa], digits=12),
                     round(mesh.y[hanging_local_spa], digits=12),
                     round(mesh.z[hanging_local_spa], digits=12))
        for (θ_k, ϕ_k, hanging_dof) in get(xyz_ang_map, xyz_key_h, Tuple{Float64,Float64,Int}[])
            (hanging_dof == 0 || hanging_dof > n_spa) && continue
            for (parent_local_spa, weight) in local_weights
                key_p = (round(mesh.x[parent_local_spa], digits=12),
                         round(mesh.y[parent_local_spa], digits=12),
                         round(mesh.z[parent_local_spa], digits=12),
                         round(θ_k, digits=12), round(ϕ_k, digits=12))
                parent_dof = get(point_dict_combined_adapted, key_p, 0)
                (parent_dof == 0 || parent_dof > n_spa) && continue
                parent_gid = ip2gip_spa[parent_dof]
                mat_entry = get!(ghost_constraint_data_spa, hanging_dof, Tuple{Int,Float64}[])
                if !any(g == parent_gid for (g, _) in mat_entry)
                    push!(mat_entry, (parent_gid, weight))
                end
                if gip2owner_extra[parent_dof] == rank
                    push!(I_rhs_sp, parent_dof); push!(J_rhs_sp, hanging_dof); push!(V_rhs_sp, weight)
                else
                    rhs_entry = get!(ghost_constraint_data_spa_rhs, hanging_dof, Tuple{Int,Float64}[])
                    if !any(g == parent_gid for (g, _) in rhs_entry)
                        push!(rhs_entry, (parent_gid, weight))
                    end
                end
            end
        end
    end
    R_spatial_rhs = sparse(I_rhs_sp, J_rhs_sp, V_rhs_sp, n_spa, n_spa)

    # ── Extended parent system (mirrors uniform path lines 1394-1429) ─────────
    # Uses gip2owner_spa (Vector by compact GID, already covers all global DOFs)
    # instead of gip2owner_spa_gid (Dict) — no ghost augmentation loop needed.
    gid_to_extended_parents_spa = Dict{Int,Int}()
    extended_parents_to_gid_spa = Int[]
    n_ghost_ext_spa = 0
    for (_, parent_list) in ghost_constraint_data_spa
        for (parent_gid, _) in parent_list
            owner = gip2owner_spa[parent_gid]
            if owner != rank &&
               !haskey(gid_to_extended_parents_spa, parent_gid) &&
               !haskey(gip_to_local_spa, parent_gid)
                n_ghost_ext_spa += 1
                gid_to_extended_parents_spa[parent_gid] = n_spa + n_ghost_ext_spa
                push!(extended_parents_to_gid_spa, parent_gid)
            end
        end
    end
    n_ext_spa = n_spa + n_ghost_ext_spa
    @rankinfo rank "$n_ghost_ext_spa extended parents → n_ext_spa=$n_ext_spa"

    gip_to_local_spa_ext = copy(gip_to_local_spa)
    for (pgid, ext_idx) in gid_to_extended_parents_spa
        gip_to_local_spa_ext[pgid] = ext_idx
    end

    if n_ghost_ext_spa > 0
        I_re, J_re, V_re = findnz(R_spatial)
        for k = 1:n_ghost_ext_spa
            ext_idx = n_spa + k
            push!(I_re, ext_idx); push!(J_re, ext_idx); push!(V_re, 1.0)
        end
        R_spatial_ext = sparse(I_re, J_re, V_re, n_ext_spa, n_ext_spa)
        P_spatial_ext = sparse(R_spatial_ext')
        I_rr, J_rr, V_rr = findnz(R_spatial_rhs)
        for k = 1:n_ghost_ext_spa
            ext_idx = n_spa + k
            push!(I_rr, ext_idx); push!(J_rr, ext_idx); push!(V_rr, 1.0)
        end
        R_spatial_rhs_ext = sparse(I_rr, J_rr, V_rr, n_ext_spa, n_ext_spa)
    else
        R_spatial_ext     = R_spatial
        P_spatial_ext     = P_spatial
        R_spatial_rhs_ext = R_spatial_rhs
    end

    return (R_spatial_ext, P_spatial_ext, R_spatial_rhs_ext,
            ghost_constraint_data_spa, ghost_constraint_data_spa_rhs,
            gid_to_extended_parents_spa, extended_parents_to_gid_spa,
            gip_to_local_spa_ext, spatial_hanging_nodes_all_angular, n_ext_spa)
end

export build_spatial_constraint_matrices, verify_spatial_constraints,
       build_spatial_restriction_and_prolongation,
       combine_spatial_angular_restrictions,
       build_spatial_constraints_for_combined_path

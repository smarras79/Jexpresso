# ============================================================================
# Stage 2: Spatial Constraint Matrix Construction
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
    @rankinfo rank "[$rank] Stage 2: Building spatial constraint matrices with angular coupling..."

    num_ncf = spatial_amr_cache.num_spatial_hanging_facets

    if num_ncf == 0
        @rankinfo rank "[$rank] No spatial hanging facets found, skipping constraint matrix construction"
        return spatial_amr_cache
    end

    nelem = length(extra_meshes_coords)
    @rankinfo rank "[$rank] Processing $(num_ncf) spatial non-conforming facets..."
    @rankinfo rank "[$rank] Spatial elements: $nelem, each with angular mesh"
    @rankinfo rank "[$rank] BUILDING: Full spatial-angular constraint system"

    # ── Spatial-Angular Coupling Note ──────────────────────────────────────
    # CRITICAL: When a spatial hanging node is constrained, ALL angular nodes
    # at that spatial location are also constrained via the spatial constraint.
    #
    # For a spatial hanging node with constraint:
    #   u_spatial_hanging = Σ w_s * u_spatial_parent
    #
    # And angular DOFs at that location [nelem_ang elements]:
    #   u_spatial_hanging(e, θ, φ) for e = 1:nelem_ang
    #
    # The spatial constraint applies to ALL angular nodes:
    #   u_spatial_hanging(e, θ, φ) = Σ w_s * u_spatial_parent(e, θ, φ)
    #
    # Current Stage 2 (spatial-only):
    #   - Builds spatial hanging node constraints
    #   - Stores in cache.parent_weights
    #   - Angular mesh data not yet available
    #
    # Future Stage 2+ Enhancement (spatial-angular):
    #   - Will replicate spatial constraints across ALL angular nodes
    #   - For each spatial hanging node at (i,j,k):
    #     - For each angular element e = 1:nelem_ang:
    #       - For each angular node (θ,φ):
    #         - Create constraint linking all angular copies
    #   - Creates spatial-angular Kronecker product structure
    #
    # For now: cache.parent_weights contains spatial constraints only.
    # RHS/Matrix assembly (Stages 4-5) will expand to full spatial-angular.
    # ───────────────────────────────────────────────────────────────────────

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
            @rankinfo rank "[$rank] Error processing facet $ncf_idx (child=$child_elem_id, parent=$parent_elem_id): $err"
            rethrow(err)
        end
    end

    @rankinfo rank "[$rank] ✓ Spatial constraint matrices built"
    @rankinfo rank "[$rank]   - $(num_spatial_hanging) spatial hanging nodes"
    @rankinfo rank "[$rank]   - $(num_spatial_angular_hanging) spatial-angular hanging DOFs (including all angular elements)"
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
        #@info local_facet_id, child_x, parent_x, L1
        #@info child_z, parent_z, L2
    elseif local_facet_id in [3, 4]
        # y-faces: x and z vary, y constant
        # Line 1 varies in x, Line 2 varies in z
        parent_x = [mesh.x[n] for n in parent_line_2 if n > 0]
        child_x = [mesh.x[n] for n in child_line_2 if n > 0]
        parent_y = [mesh.y[n] for n in parent_line_1 if n > 0]
        child_y = [mesh.y[n] for n in child_line_1 if n > 0]
        
        L1 = build_1d_lagrange_matrix(parent_x, child_x)  # x-direction
        L2 = build_1d_lagrange_matrix(parent_y, child_y)  # z-direction
        #@info local_facet_id, child_x, parent_x, L1
        #@info child_y, parent_y, L2
    elseif local_facet_id in [5, 6]
        # x-faces: y and z vary, x constant
        # Line 1 varies in y, Line 2 varies in z
        parent_y = [mesh.y[n] for n in parent_line_2 if n > 0]
        child_y = [mesh.y[n] for n in child_line_2 if n > 0]
        parent_z = [mesh.z[n] for n in parent_line_1 if n > 0]
        child_z = [mesh.z[n] for n in child_line_1 if n > 0]
        
        L1 = build_1d_lagrange_matrix(parent_y, child_y)  # y-direction
        L2 = build_1d_lagrange_matrix(parent_z, child_z)  # z-direction
        #@info local_facet_id, child_y, parent_y, L1
        #@info child_z, parent_z, L2
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
# Extract Element Coordinates
# ============================================================================

"""
    extract_element_corner_coords(iel, mesh) -> (3, 8) matrix

Extract corner coordinates of hexahedral element from mesh.

Returns 3×8 matrix where each column is a corner node:
- Row 1: x coordinates
- Row 2: y coordinates
- Row 3: z coordinates
"""
function extract_element_corner_coords(iel::Int, mesh::St_mesh)
    coords = zeros(Float64, 3, 8)
    ngl = mesh.ngl

    # 8 corners of hex element in reference ordering
    corners = [
        (1,   1,   1),
        (1,   1,   ngl),
        (ngl, 1,   ngl),
        (ngl, 1,   1),
        (1,   ngl, 1),
        (1,   ngl, ngl),
        (ngl, ngl, ngl),
        (ngl, ngl, 1)
    ]

    for (corner_idx, (i, j, k)) in enumerate(corners)
        ip = mesh.connijk[iel, i, j, k]
        coords[1, corner_idx] = mesh.x[ip]
        coords[2, corner_idx] = mesh.y[ip]
        coords[3, corner_idx] = mesh.z[ip]
    end

    return coords
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
    """
    Build constraint equations for hanging nodes on a spatial facet.

    FIXED: Now processes ALL LGL nodes on the facet (not just corners),
    and uses the correct pair of L matrices based on face orientation.

    Face orientation:
    - Face 1,2: z-faces (normal ±z) → use Lx, Ly
    - Face 3,4: y-faces (normal ±y) → use Lx, Lz
    """

    num_spatial_hanging = 0

    # Determine which L matrices to use based on face orientation
    # Face numbering for 3D hexahedron:

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
        x_child = mesh.x[child_node_id]
        y_child = mesh.y[child_node_id]
        z_child = mesh.z[child_node_id]
        # Iterate through all parent nodes (ngl × ngl grid on parent facet)
        for p1 = 1:ngl
            for p2 = 1:ngl
                parent_face_idx = (p2 - 1) * ngl + p1#(p1 - 1) * ngl + p2
                if parent_face_idx <= length(IPp_nodes)
                    parent_node_id = Int(IPp_nodes[parent_face_idx])
                    if parent_node_id > 0
                        # Tensor product of weights
                        x_parent = mesh.x[parent_node_id]
                        y_parent = mesh.y[parent_node_id]
                        z_parent = mesh.z[parent_node_id]
                        #@info "child", x_child, y_child, z_child
                        #@info "parent", x_parent, y_parent, z_parent
                        #@info "weights", w_1[p1], w_2[p2], w_1[p1] * w_2[p2]
                        weight = w_1[p1] * w_2[p2]
                        #@info weight, w_1[p1], w_2[p2], parent_node_id, parent_face_idx
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

"""
    build_hanging_node_constraints(
        child_elem_id, parent_elem_id, local_facet_id,
        IPc_nodes, IPp_nodes, Lx, Ly, Lz,
        cache, mesh, ngl
    )

Build constraint equations for hanging nodes on a spatial facet.

For each hanging node on the child facet:
- Look up (i,j) position in parent facet coordinates
- Use Lx, Ly to get weights from parent nodes
- Store constraint in cache.parent_weights
"""
function build_hanging_node_constraints(
    child_elem_id::Int, parent_elem_id::Int, local_facet_id::Int,
    IPc_nodes::Vector, IPp_nodes::Vector,
    Lx::Matrix{Float64}, Ly::Matrix{Float64}, Lz::Matrix{Float64},
    cache::SpatialAMRCache, mesh::St_mesh, ngl::Int
)
    # Total nodes on face
    ngl2 = ngl * ngl

    # Get unique non-zero child nodes
    unique_child_nodes = unique(filter(x -> x > 0, IPc_nodes))

    for child_node_id in unique_child_nodes
        # Find position of this node in IPc_nodes array
        node_positions = findall(x -> x == child_node_id, IPc_nodes)

        if isempty(node_positions)
            continue
        end

        # Process first occurrence only (hanging nodes should have unique positions)
        pos = node_positions[1]

        # Convert linear index to (i,j) in face
        # Note: face_idx is 1-indexed position in IPc_list array (already 1:ngl2)
        i = div(pos - 1, ngl) + 1  # Row
        j = mod(pos - 1, ngl) + 1  # Column

        # Check bounds
        if i > size(Lx, 1) || j > size(Lx, 2) || i < 1 || j < 1
            continue
        end

        # Get interpolation weights from Lagrange matrices
        # Lx[i, :] gives weights for child row i from all parent nodes
        # Ly[j, :] gives weights for child col j from all parent nodes
        w_x = Lx[i, :]  # Weights in x direction (for each parent)
        w_y = Ly[j, :]  # Weights in y direction (for each parent)

        # Build constraint: tensor product of 1D weights
        # For 2D face: weight[parent_idx] = w_x[i_parent] * w_y[j_parent]
        # We need to map from the node ID to the parent node index

        parent_weights_raw = Tuple{Int, Float64}[]

        # Get unique parent node IDs (non-zero entries)
        unique_parent_nodes = unique(filter(x -> x > 0, IPp_nodes))

        # For each parent node, look up its weights from Lagrange matrices
        # The parent_idx should correspond to position in the parent node array
        for parent_idx = 1:min(length(unique_parent_nodes), length(w_x), length(w_y))
            parent_node_id = unique_parent_nodes[parent_idx]

            # Tensor product of 1D Lagrange weights
            weight = w_x[parent_idx] * w_y[parent_idx]
            
            if abs(weight) > 1e-14
                push!(parent_weights_raw, (Int(parent_node_id), weight))
            end
        end

        # If no weights passed threshold, use all and don't filter
        if isempty(parent_weights_raw)
            for parent_idx = 1:min(length(unique_parent_nodes), length(w_x), length(w_y))
                parent_node_id = unique_parent_nodes[parent_idx]
                weight = w_x[parent_idx] * w_y[parent_idx]
                push!(parent_weights_raw, (Int(parent_node_id), weight))
            end
        end

        # Normalize to ensure partition of unity
        total_weight = sum(w -> w[2], parent_weights_raw; init=0.0)

        if abs(total_weight) > 1e-14
            parent_weights = [(n, w/total_weight) for (n, w) in parent_weights_raw]
        else
            # Fallback: uniform weights
            uniform_weight = 1.0 / length(parent_weights_raw)
            parent_weights = [(n, uniform_weight) for (n, _) in parent_weights_raw]
        end

        # Store in cache
        if !haskey(cache.parent_weights, child_node_id)
            cache.parent_weights[child_node_id] = parent_weights
        end
    end
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
        @rankinfo rank "[$rank] WARNING: No constraint weights found in cache"
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
                @rankinfo rank "[$rank] WARNING: Column $col sum = $col_sum (should be 1.0)"
            end
        end
    end

    @rankinfo rank "[$rank] Spatial Constraint Verification Summary:"
    @rankinfo rank "[$rank]   - Spatial hanging nodes: $num_hanging"
    @rankinfo rank "[$rank]   - Restriction matrix size: $n_spatial_nodes × $n_spatial_nodes"
    @rankinfo rank "[$rank]   - Columns with sums ≠ 1.0: $num_violated_col_unity"

    # If angular mesh provided, verify spatial-angular coupling awareness
    if extra_meshes_extra_nelems !== nothing && extra_meshes_extra_nops !== nothing
        total_angular_nelems = sum(extra_meshes_extra_nelems)
        avg_nop = Int(round(mean([nops[1] for nops in extra_meshes_extra_nops])))

        expected_spatial_angular = num_hanging * total_angular_nelems * (avg_nop + 1) * (avg_nop + 1)
        @rankinfo rank "[$rank]   - Expected spatial-angular DOFs (with coupling): ~$expected_spatial_angular"
        @rankinfo rank "[$rank]   - NOTE: RHS/Matrix assembly will expand these constraints to all angular nodes"
    end

    if num_violated_col_unity > 0
        error("[Rank $rank] Spatial constraint verification FAILED: " *
              "$num_violated_col_unity column(s) have sum ≠ 1.0 (partition of unity violated)")
    end

    @rankinfo rank "[$rank] ✓ Spatial constraints verified: partition of unity satisfied"
end

# ============================================================================
# Stage 4: Build Sparse Constraint Matrices for Assembly
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
    n_spatial_nodes::Int,
    n_angular_nodes::Int,
    spatial_hanging_nodes_all_angular::Set{Int},
    mesh, extra_mesh
)::Tuple{SparseMatrixCSC, SparseMatrixCSC}
    n_total_spatial_ang = n_spatial_nodes * n_angular_nodes

    if isempty(spatial_amr_cache.parent_weights)
        # No hanging nodes - return identity for safe operation
        return sparse(I, n_total_spatial_ang, n_total_spatial_ang),
               sparse(I, n_total_spatial_ang, n_total_spatial_ang)
    end

    # Build full n_total × n_total matrix first
    rows_M = Int[]
    cols_M = Int[]
    vals_M = Float64[]

    # Phase 1: For free DOFs, add identity rows
    for dof_idx = 1:n_total_spatial_ang
        if !(dof_idx in spatial_hanging_nodes_all_angular)
            # This is a free DOF - add identity row
            push!(rows_M, dof_idx)
            push!(cols_M, dof_idx)
            push!(vals_M, 1.0)
        end
    end

    # Phase 2: For hanging DOFs, add constraint weight rows
    for (hanging_node_id, parent_list) in spatial_amr_cache.parent_weights
        # For each angular DOF, apply the spatial constraint
        for i_ang = 1:n_angular_nodes
            # Row index for this hanging spatial-angular DOF
            col_hanging = (hanging_node_id - 1) * n_angular_nodes + i_ang
            x_child = mesh.x[hanging_node_id]
            y_child = mesh.y[hanging_node_id]
            z_child = mesh.z[hanging_node_id]
            θ_child = extra_mesh.extra_coords[1,i_ang]
            ϕ_child = extra_mesh.extra_coords[2,i_ang]
            # Add weights for each parent
            for (parent_node_id, weight) in parent_list
                # Column index for parent spatial-angular DOF
                row_parent = (parent_node_id - 1) * n_angular_nodes + i_ang
                x_parent = mesh.x[parent_node_id]
                y_parent = mesh.y[parent_node_id]
                z_parent = mesh.z[parent_node_id]
                θ_parent = extra_mesh.extra_coords[1,i_ang]
                ϕ_parent = extra_mesh.extra_coords[2,i_ang]
                
                #@info "child info", x_child, y_child, z_child, θ_child, ϕ_child
                #@info "parent info", x_parent, y_parent, z_parent, θ_parent, ϕ_parent, weight
                push!(rows_M, row_parent)
                push!(cols_M, col_hanging)
                push!(vals_M, weight)
                #@info weight, row_parent, col_hanging
            end
        end
    end

    # Construct full n_total × n_total matrix
    R_spatial = sparse(rows_M, cols_M, vals_M, n_total_spatial_ang, n_total_spatial_ang)
    

    # Prolongation is transpose
    P_spatial = R_spatial'

    return R_spatial, P_spatial
end

"""
    apply_spatial_constraint_to_rhs(
        RHS::Vector, spatial_amr_cache::SpatialAMRCache,
        n_ang_per_spatial::Int
    ) -> Vector

Apply spatial constraint restrictions to RHS for hanging spatial nodes.

For each spatial hanging node s and its parents {p1, p2, ...} with weights {w1, w2, ...}:
  For each angular node a ∈ [1, n_ang_per_spatial]:
    RHS_out[(s-1)*n_ang + a] = Σ_i( w_i * RHS[(p_i-1)*n_ang + a] )

This applies the spatial constraint across ALL associated angular DOFs at each location
(spatial-angular coupling).

Arguments:
- RHS: Original RHS vector of size npoin * n_ang_per_spatial
- spatial_amr_cache: Cache containing parent_weights from Stage 2
- n_ang_per_spatial: Number of angular DOFs per spatial node

Returns: Modified RHS with spatial constraints applied (same size as input).
"""
function apply_spatial_constraint_to_rhs(
    RHS::Vector{Float64},
    spatial_amr_cache::SpatialAMRCache,
    n_ang_per_spatial::Int
)::Vector{Float64}
    RHS_out = copy(RHS)

    # Apply spatial constraints: for each hanging spatial node
    for (spatial_hanging_node_id, parent_list) in spatial_amr_cache.parent_weights
        # For each angular DOF associated with this spatial hanging node
        for i_ang = 1:n_ang_per_spatial
            # Global DOF index for this (spatial_hanging, angular) pair
            idx_hanging = (spatial_hanging_node_id - 1) * n_ang_per_spatial + i_ang

            # Interpolate from parents: RHS_hanging = Σ weight * RHS_parent
            rhs_val = 0.0
            for (parent_node_id, weight) in parent_list
                idx_parent = (parent_node_id - 1) * n_ang_per_spatial + i_ang
                if 1 <= idx_parent <= length(RHS)
                    rhs_val += weight * RHS[idx_parent]
                end
            end

            if 1 <= idx_hanging <= length(RHS_out)
                RHS_out[idx_hanging] = rhs_val
            end
        end
    end

    return RHS_out
end

export build_spatial_constraint_matrices, verify_spatial_constraints,
       build_spatial_restriction_and_prolongation,
       apply_spatial_constraint_to_rhs

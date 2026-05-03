# ============================================================================
# Stage 1: Spatial AMR Cache Structure - Task 1.2
# ============================================================================
# Cache for spatial AMR data needed by Stage 2+ (constraint matrix construction)
# This structure holds pre-computed parent element information and interpolation data

"""
    SpatialAMRCache

Cache structure for spatial AMR hanging node information to be used during
constraint matrix construction and RHS/matrix assembly.

Fields:
- `parent_map::Dict{Int, Int}`: Mapping from hanging node global index to parent element
  (Will be populated in Stage 1b, used in Stage 2 for constraint equations)

- `parent_weights::Dict{Int, Vector{Tuple{Int, Float64}}}`: Lagrange interpolation weights
  for each hanging node. Maps hanging_node_id → [(parent_node_id, weight), ...]
  (Populated in Stage 2 during spatial constraint matrix construction)

- `ghost_spatial_parents::Dict{Int, SpatialElementGhostInfo}`: Information about
  spatial parent elements on ghost processors (for MPI communication, Stage 3)

- `spatial_interp_cache::Dict{Tuple{Int,Int,Int}, Matrix{Float64}}`: Cache of
  Lagrange interpolation matrices. Key: (child_ref_level, parent_ref_level, ngl)
  (Populated in Stage 2 to avoid recomputation)

- `element_refinement_levels::Vector{Int}`: Copy of mesh.ad_lvl for quick lookup

- `num_spatial_hanging_facets::Int`: Count of spatial non-conforming facets
  (Copy of mesh.num_ncf)
"""
Base.@kwdef struct SpatialAMRCache
    # Parent element mapping for hanging nodes
    parent_map::Dict{Int, Int} = Dict{Int, Int}()

    # Interpolation weights for LOCAL (same-rank) constraints.
    # hanging_local_spatial_ip → [(parent_local_spatial_ip, weight), ...]
    parent_weights::Dict{Int, Vector{Tuple{Int, Float64}}} =
        Dict{Int, Vector{Tuple{Int, Float64}}}()

    # Interpolation weights for CROSS-RANK constraints.
    # Parent lives on another rank — keyed by GLOBAL spatial IP (not local).
    # hanging_local_spatial_ip → [(parent_GLOBAL_spatial_ip, weight), ...]
    # Used in Stage 5 MPI pattern to exchange row/col effects to the parent rank.
    cross_rank_parent_weights::Dict{Int, Vector{Tuple{Int, Float64}}} =
        Dict{Int, Vector{Tuple{Int, Float64}}}()

    # Ghost spatial parent info (for inter-processor communication)
    ghost_spatial_parents::Dict{Int, Any} = Dict{Int, Any}()

    # Cached Lagrange matrices to avoid recomputation during RHS/matrix assembly
    spatial_interp_cache::Dict{Tuple{Int,Int,Int}, Matrix{Float64}} =
        Dict{Tuple{Int,Int,Int}, Matrix{Float64}}()

    # Quick lookup for refinement levels
    element_refinement_levels::Vector{Int} = Int[]

    # Count of spatial non-conforming facets
    num_spatial_hanging_facets::Int = 0

    # Coincident nodes: child face nodes that sit exactly on a parent node (weight=1,
    # single parent). They are NOT constrained via interpolation — instead their
    # assembled row is zeroed (identity equation, zero RHS) and after the solve the
    # parent's solution value is copied into them.
    # Maps local spatial child IP → local spatial parent IP (same-rank case).
    coincident_nodes::Dict{Int, Int} = Dict{Int, Int}()

    # Cross-rank coincident nodes: child is local, parent is on another rank.
    # Maps local spatial child IP → GLOBAL spatial parent IP.
    cross_rank_coincident_nodes::Dict{Int, Int} = Dict{Int, Int}()
end

"""
    SpatialElementGhostInfo

Information about a spatial element on a ghost processor.
Used in Stage 3 for ghost layer extension.

Fields:
- `spatial_elem_id::Int`: Local element ID on ghost processor
- `owner_rank::Int`: Rank that owns this element
- `parent_element_id::Int`: ID of parent element (if this is a hanging element)
- `parent_spatial_elem_id::Int`: Parent element ID on owner's rank
- `spatial_jacobian::Float64`: Jacobian of spatial element (or SVector for 3D)
- `spatial_coords::Matrix{Float64}`: Node coordinates on ghost processor
- `spatial_hanging_nodes::Dict{Int, Vector{Tuple{Int, Float64}}}`: Hanging node constraints
"""
Base.@kwdef mutable struct SpatialElementGhostInfo
    spatial_elem_id::Int = 0
    owner_rank::Int = 0
    parent_element_id::Int = 0
    parent_spatial_elem_id::Int = 0
    spatial_jacobian::Float64 = 0.0
    spatial_coords::Matrix{Float64} = zeros(3, 8)  # 8 nodes per 3D hex element
    spatial_hanging_nodes::Dict{Int, Vector{Tuple{Int, Float64}}} =
        Dict{Int, Vector{Tuple{Int, Float64}}}()
end

export SpatialAMRCache, SpatialElementGhostInfo

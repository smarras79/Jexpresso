# ============================================================================
# Stage 1: Spatial Mesh Infrastructure Integration - Query Functions
# ============================================================================
# These functions provide access to spatial AMR structures for Stage 1+ implementation
# They query the mesh data structures populated by p4est/Gridap for spatial hanging nodes

"""
    get_spatial_hanging_nodes_for_element(mesh::St_mesh, iel::Int) -> Vector{Tuple{Int, Int}}

Return list of spatial non-conforming facets involving element `iel`.

For each non-conforming facet where `iel` appears as child or parent:
- Returns tuple `(ncf_idx, other_elem_id)` where `ncf_idx` is the index in mesh.cip/pip arrays
- `other_elem_id` is the parent if `iel` is child, or child if `iel` is parent

Returns empty vector if element has no spatial hanging node neighbors.
Those facets can be queried further using get_spatial_parent_element() and
get_spatial_parent_element_child_facet_indices().
"""
function get_spatial_hanging_nodes_for_element(mesh::St_mesh, iel::Int)
    hanging_facets = Tuple{Int, Int}[]

    # Iterate through all spatial non-conforming facets
    for idx in 1:mesh.num_ncf
        cid = Int(mesh.cip[idx])  # Child element ID
        pid = Int(mesh.pip[idx])  # Parent element ID

        # Check if this facet involves our element
        if cid == iel
            # Our element is the child: we have hanging nodes on this facet
            push!(hanging_facets, (idx, pid))
        elseif pid == iel
            # Our element is the parent: neighbor is refined
            push!(hanging_facets, (idx, cid))
        end
    end

    return hanging_facets
end

"""
    get_spatial_parent_element(mesh::St_mesh, ncf_idx::Int) -> Int

Return parent spatial element ID for non-conforming facet at index `ncf_idx`.

The parent element is the coarser element. Its DOFs will appear in
the constraint equation for hanging nodes on the child element side of the facet.
"""
function get_spatial_parent_element(mesh::St_mesh, ncf_idx::Int)
    if ncf_idx < 1 || ncf_idx > mesh.num_ncf
        error("Invalid non-conforming facet index: $ncf_idx (valid range: 1:$(mesh.num_ncf))")
    end
    return Int(mesh.pip[ncf_idx])
end

"""
    get_spatial_child_element(mesh::St_mesh, ncf_idx::Int) -> Int

Return child spatial element ID for non-conforming facet at index `ncf_idx`.

The child element is the finer element. Its DOFs will be constrained to
the parent element's DOFs via Lagrange interpolation on this facet.
"""
function get_spatial_child_element(mesh::St_mesh, ncf_idx::Int)
    if ncf_idx < 1 || ncf_idx > mesh.num_ncf
        error("Invalid non-conforming facet index: $ncf_idx (valid range: 1:$(mesh.num_ncf))")
    end
    return Int(mesh.cip[ncf_idx])
end

"""
    get_spatial_child_facet_indices(mesh::St_mesh, ncf_idx::Int) -> (Int, Int, Int)

Return information about the child facet for non-conforming facet `ncf_idx`.

Returns tuple `(local_facet_id, refinement_half_1, refinement_half_2)`:
- `local_facet_id`: Which face of the child element (1-6 in 3D)
- `refinement_half_1`, `refinement_half_2`: Info about refinement pattern on that face

This info is needed to build the spatial constraint matrix in Stage 2.
"""
function get_spatial_child_facet_indices(mesh::St_mesh, ncf_idx::Int)
    if ncf_idx < 1 || ncf_idx > mesh.num_ncf
        error("Invalid non-conforming facet index: $ncf_idx (valid range: 1:$(mesh.num_ncf))")
    end
    return (Int(mesh.lfid[ncf_idx]), Int(mesh.half1[ncf_idx]), Int(mesh.half2[ncf_idx]))
end

"""
    get_spatial_parent_child_node_mapping(mesh::St_mesh, ncf_idx::Int) -> (Vector{Int}, Vector{Int})

Return node mappings between parent and child faces for non-conforming facet `ncf_idx`.

Returns tuple `(IPc_nodes, IPp_nodes)` where:
- `IPc_nodes`: Global node indices on child element's side of the facet
- `IPp_nodes`: Global node indices on parent element's side of the facet

These node lists are pre-computed and stored in mesh.IPc_list and mesh.IPp_list.
They will be used in Stage 2 to build the actual Lagrange interpolation constraint matrix.
"""
function get_spatial_parent_child_node_mapping(mesh::St_mesh, ncf_idx::Int)
    if ncf_idx < 1 || ncf_idx > mesh.num_ncf
        error("Invalid non-conforming facet index: $ncf_idx (valid range: 1:$(mesh.num_ncf))")
    end
    IPc = mesh.IPc_list[:, ncf_idx]
    IPp = mesh.IPp_list[:, ncf_idx]
    return (Vector{Int}(IPc), Vector{Int}(IPp))
end

"""
    get_spatial_hanging_nodes_degrees(mesh::St_mesh) -> Vector{Int}

Return refinement level (polynomial degree) for each spatial element.

mesh.ad_lvl[iel] gives the refinement level. Elements with higher levels
are refined and may have hanging nodes on their boundaries.
"""
function get_spatial_hanging_nodes_degrees(mesh::St_mesh)
    return Vector{Int}(mesh.ad_lvl)
end

export get_spatial_hanging_nodes_for_element
export get_spatial_parent_element
export get_spatial_child_element
export get_spatial_child_facet_indices
export get_spatial_parent_child_node_mapping
export get_spatial_hanging_nodes_degrees

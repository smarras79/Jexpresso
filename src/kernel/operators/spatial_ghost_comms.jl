# ============================================================================
# Stage 3: Ghost Layer Extension for Spatial Elements
# ============================================================================
# MPI communication of spatial constraint data across rank boundaries.
#
# The mesh already performs ghost exchange for NCF (non-conforming facet) nodes:
#   - mesh.pgip_ghost: global IPs of parent interface nodes for each
#     parent-ghost NCF (child local, parent ghost). These were received from
#     the parent rank during mesh setup via get_ghost_ips + send_and_receive.
#   - mesh.IPc_list_pg: child-side local interface node IPs for each parent-ghost NCF.
#   - mesh.IPp_list_cg: parent-side local interface node IPs for each child-ghost NCF.
#
# For purely spatial AMR, the Lagrange interpolation weights needed to express
# hanging child nodes in terms of parent nodes can be built entirely from:
#   1. Local child interface nodes (always available).
#   2. Ghost parent interface nodes (global IPs in mesh.pgip_ghost, convertible
#      to local IPs via mesh.gip2ip).
#
# After building these constraints, the GMRES AllReduce correctly combines
# contributions from all ranks for both the LHS (R*A*P) and RHS (R*b) because:
#   - Ghost parent DOFs appear in ip2gip_spa with valid global IDs.
#   - The GMRES global solution x covers all global DOFs, so P_spatial*solution
#     correctly fills hanging node values from parent values (including ghost parents).
#
# No additional MPI communication is required beyond what the mesh setup provides.

"""
    exchange_spatial_ghosts(
        mesh, spatial_amr_cache, extra_meshes_extra_nops, extra_meshes_extra_nelems, rank, comm
    ) -> SpatialAMRCache

Build spatial constraint matrices for cross-rank non-conforming facets.

Uses pre-exchanged ghost node data from the mesh setup (mesh.pgip_ghost,
mesh.IPc_list_pg, etc.) to build Lagrange interpolation constraints for
parent-ghost NCFs (child element local, parent element on another rank).

Returns updated cache with all spatial hanging node constraints populated,
including cross-rank ones.
"""
function exchange_spatial_ghosts(
    mesh::St_mesh, spatial_amr_cache::SpatialAMRCache,
    extra_meshes_extra_nops, extra_meshes_extra_nelems,
    rank::Int, comm::MPI.Comm
)
    @info rank "[$rank] Stage 3: Building spatial ghost constraints..."

    nprocs = MPI.Comm_size(comm)

    if nprocs == 1
        @rankinfo rank "[$rank] Single rank, no ghost exchange needed"
        return spatial_amr_cache
    end

    num_ncf_pg = mesh.num_ncf_pg
    @rankinfo rank "[$rank] Processing $num_ncf_pg parent-ghost NCFs (child local, parent ghost)..."

    # ALL ranks must call build_spatial_constraints_parent_ghost regardless of
    # num_ncf_pg. Ranks without child NCFs may still own parent nodes that other
    # ranks need coordinates from — skipping here would deadlock MPI.Alltoall
    # inside fetch_parent_face_coordinates.
    spatial_amr_cache = build_spatial_constraints_parent_ghost(
        mesh, spatial_amr_cache,
        extra_meshes_extra_nelems, extra_meshes_extra_nops,
        mesh.ngl, rank, comm
    )

    n_cross_rank = length(spatial_amr_cache.cross_rank_parent_weights)
    @rankinfo rank "[$rank] ✓ Cross-rank constraints built for $n_cross_rank hanging nodes"

    return spatial_amr_cache
end

# ============================================================================
# MPI: Fetch Parent Face Node Coordinates from Owning Ranks
# ============================================================================

"""
    fetch_parent_face_coordinates(mesh, rank, comm) -> Dict{Int, NTuple{3,Float64}}

Two-round point-to-point coordinate exchange for cross-rank NCF parent nodes.

Round 1: each rank sends the global IPs it needs (from mesh.pgip_ghost) to the
         owning ranks (from mesh.pgip_owner).
Round 2: each owner looks up its local coordinates and sends (gip, x, y, z) back.

Returns a Dict mapping parent_global_ip → (x, y, z).
Local nodes that ARE already in the local mesh are also added to the cache.
"""
function fetch_parent_face_coordinates(
    mesh::St_mesh, rank::Int, comm::MPI.Comm
)
    coord_cache = Dict{Int, NTuple{3,Float64}}()

    # ── Add locally available parent nodes (gip2ip gives a valid index) ──────
    for gip in mesh.pgip_ghost
        gip <= 0 && continue
        gip > length(mesh.gip2ip) && continue
        lip = Int(mesh.gip2ip[gip])
        lip >= 1 && lip <= mesh.npoin || continue
        coord_cache[gip] = (mesh.x[lip], mesh.y[lip], mesh.z[lip])
    end

    nprocs = MPI.Comm_size(comm)
    if nprocs == 1
        return coord_cache
    end

    # ── Filter: GIPs not already resolved locally ────────────────────────────
    needed_gips   = Int[]
    needed_owners = Int[]
    for k = 1:length(mesh.pgip_ghost)
        gip   = Int(mesh.pgip_ghost[k])
        owner = Int(mesh.pgip_owner[k])
        gip   <= 0 && continue
        owner == rank && continue           # we own it (already in cache above)
        haskey(coord_cache, gip) && continue
        push!(needed_gips,   gip)
        push!(needed_owners, owner)
    end

    # NOTE: do NOT return early here even if needed_gips is empty.
    # Ranks with no child NCFs may still own parent nodes that other ranks are
    # requesting. All ranks must participate in every send_and_receive call
    # (each uses MPI.Alltoall internally) to avoid deadlock.

    # ── Round 1: send GIP requests to owners ─────────────────────────────────
    recv_gip_requests, request_senders =
        send_and_receive(needed_gips, needed_owners, comm)

    # ── Build coordinate responses for received requests ─────────────────────
    n_recv = length(recv_gip_requests)
    resp_gips = Vector{Int}(undef, n_recv)
    resp_x    = Vector{Float64}(undef, n_recv)
    resp_y    = Vector{Float64}(undef, n_recv)
    resp_z    = Vector{Float64}(undef, n_recv)

    for k = 1:n_recv
        gip = Int(recv_gip_requests[k])
        resp_gips[k] = gip
        lip = (gip >= 1 && gip <= length(mesh.gip2ip)) ? Int(mesh.gip2ip[gip]) : 0
        if lip >= 1 && lip <= mesh.npoin
            resp_x[k] = mesh.x[lip]
            resp_y[k] = mesh.y[lip]
            resp_z[k] = mesh.z[lip]
        else
            # GIP requested from us but not in our local mesh — should not happen
            # if pgip_owner is correct. Respond with zeros as fallback.
            resp_x[k] = 0.0; resp_y[k] = 0.0; resp_z[k] = 0.0
        end
    end

    # ── Round 2: send coordinates back to requesters ─────────────────────────
    recv_gips_back, _ = send_and_receive(resp_gips, request_senders, comm)
    recv_x,         _ = send_and_receive(resp_x,    request_senders, comm)
    recv_y,         _ = send_and_receive(resp_y,    request_senders, comm)
    recv_z,         _ = send_and_receive(resp_z,    request_senders, comm)

    for k = 1:length(recv_gips_back)
        gip = Int(recv_gips_back[k])
        gip > 0 || continue
        coord_cache[gip] = (recv_x[k], recv_y[k], recv_z[k])
    end

    return coord_cache
end

# ============================================================================
# Build Constraints for Parent-Ghost NCFs
# ============================================================================

"""
    build_spatial_constraints_parent_ghost(
        mesh, spatial_amr_cache,
        extra_meshes_extra_nelems, extra_meshes_extra_nops,
        ngl, rank, comm
    ) -> SpatialAMRCache

Build Lagrange interpolation constraints for cross-rank parent-ghost NCFs.

For each NCF where the child element is LOCAL and the parent is on ANOTHER rank:
  - Fetches parent face node coordinates via MPI (two-round exchange).
  - Builds Lagrange interpolation weights from child/parent face coordinates.
  - Stores constraints in spatial_amr_cache.cross_rank_parent_weights using
    GLOBAL spatial IPs for the parent (no valid local IP exists on this rank).

The Stage 5 MPI pattern later communicates these effects to the parent rank.
"""
function build_spatial_constraints_parent_ghost(
    mesh::St_mesh, spatial_amr_cache::SpatialAMRCache,
    extra_meshes_extra_nelems, extra_meshes_extra_nops,
    ngl::Int, rank::Int, comm::MPI.Comm
)
    num_ncf_pg = mesh.num_ncf_pg
    ngl2 = ngl * ngl  # nodes per face (3D)

    num_constraints_added = 0
    num_skipped = 0

    # ── Pre-fetch all parent face node coordinates via MPI ────────────────────
    coord_cache = fetch_parent_face_coordinates(mesh, rank, comm)
    @info rank "[$rank] coord_cache populated with $(length(coord_cache)) parent node coordinates"

    for pg_idx = 1:num_ncf_pg
        local_facet_id = Int(mesh.lfid_pg[pg_idx])

        # ── Child interface nodes (local) ─────────────────────────────────────
        IPc_nodes = Vector{Int}(mesh.IPc_list_pg[:, pg_idx])

        if all(x -> x <= 0, IPc_nodes)
            num_skipped += 1
            continue
        end

        # ── Parent interface node GLOBAL IPs ─────────────────────────────────
        offset_start = (pg_idx - 1) * ngl2 + 1
        offset_end   = pg_idx * ngl2

        if offset_end > length(mesh.pgip_ghost)
            num_skipped += 1
            continue
        end

        parent_global_ips = Vector{Int}(mesh.pgip_ghost[offset_start:offset_end])

        if any(x -> x <= 0, parent_global_ips)
            num_skipped += 1
            continue
        end

        # Check that all parent coordinates are available in the cache
        if !all(gip -> haskey(coord_cache, gip), parent_global_ips)
            missing = count(gip -> !haskey(coord_cache, gip), parent_global_ips)
            @info rank "[$rank] pg_idx=$pg_idx: $missing parent GIPs missing from coord_cache, skipping"
            num_skipped += 1
            continue
        end

        # ── Build Lagrange interpolation matrices using coordinate cache ───────
        try
            L1, L2 = build_spatial_interpolation_matrices_with_coord_cache(
                local_facet_id, IPc_nodes, parent_global_ips, mesh, coord_cache, ngl
            )

            # ── Store cross-rank constraints (global parent IPs) ──────────────
            num_added = build_cross_rank_hanging_node_constraints(
                IPc_nodes, parent_global_ips, L1, L2,
                spatial_amr_cache, mesh, ngl,
                extra_meshes_extra_nelems, extra_meshes_extra_nops, rank
            )
            num_constraints_added += num_added

        catch err
            @info rank "[$rank] pg_idx=$pg_idx: error building cross-rank constraints: $err"
            num_skipped += 1
            continue
        end
    end

    @info rank "[$rank] Parent-ghost NCFs: $num_constraints_added cross-rank constraints added, $num_skipped skipped"
    return spatial_amr_cache
end

# ============================================================================
# Build Spatial Ghost Constraint Data (for MPI effects exchange)
# ============================================================================

"""
    build_spatial_ghost_constraint_data(
        spatial_amr_cache, mesh, ip2gip_spa, gip2owner_extra, n_ang, rank
    ) -> (ghost_constraint_data, reverse_ghost_map)

Identify cross-rank spatial hanging nodes and build data structures needed for:
  1. ghost_constraint_data: local hanging DOF → [(global parent DOF, weight)]
     For hanging DOFs whose parents are on other ranks.
     Used to detect cross-rank structure (informational).
  2. reverse_ghost_map: local owned parent DOF → [(global hanging DOF, owner, weight)]
     Used for solution prolongation: after GMRES solve, owned parents send their
     solution values to the hanging node owners.

NOTE: In practice, for the non-adaptive spatial case, the GMRES AllReduce
already handles the LHS and the global solution already contains correct values
for ghost parent DOFs. So these structures are mainly for diagnostics and
for potential future explicit RHS communication.
"""
function build_spatial_ghost_constraint_data(
    spatial_amr_cache::SpatialAMRCache,
    ip2gip_spa::Vector{Int},
    gip2owner_extra::Vector{Int},
    n_ang::Int,
    rank::Int
)
    ghost_constraint_data = Dict{Int, Vector{Tuple{Int, Float64}}}()
    reverse_ghost_map = Dict{Int, Vector{Tuple{Int, Int, Float64}}}()

    n_total_dofs = length(ip2gip_spa)

    for (hanging_spatial_ip, weights) in spatial_amr_cache.parent_weights
        has_cross_rank_parent = false
        cross_rank_parents = Tuple{Int, Float64}[]

        for (parent_spatial_ip, weight) in weights
            # Parent DOF owner (check first angular DOF)
            parent_dof = (parent_spatial_ip - 1) * n_ang + 1
            if parent_dof < 1 || parent_dof > n_total_dofs
                continue
            end
            parent_owner = gip2owner_extra[parent_dof]

            if parent_owner != rank
                has_cross_rank_parent = true
                push!(cross_rank_parents, (parent_spatial_ip, weight))
            end
        end

        if has_cross_rank_parent
            # For each angular DOF of this hanging spatial node
            for i_ang = 1:n_ang
                dof_h = (hanging_spatial_ip - 1) * n_ang + i_ang
                if dof_h < 1 || dof_h > n_total_dofs
                    continue
                end

                # Build full parent constraint list (global DOF, weight)
                all_parent_dofs = Tuple{Int, Float64}[]
                for (parent_spatial_ip, weight) in weights
                    parent_dof = (parent_spatial_ip - 1) * n_ang + i_ang
                    if parent_dof >= 1 && parent_dof <= n_total_dofs
                        parent_gid = ip2gip_spa[parent_dof]
                        push!(all_parent_dofs, (parent_gid, weight))
                    end
                end

                if !isempty(all_parent_dofs)
                    ghost_constraint_data[dof_h] = all_parent_dofs
                end
            end
        end
    end

    # Build reverse_ghost_map: for owned parents, what hanging nodes depend on them
    # (on other ranks). Used for explicit prolongation communication if needed.
    for (dof_h, parent_constraints) in ghost_constraint_data
        dof_h <= n_total_dofs || continue
        gid_h = ip2gip_spa[dof_h]
        owner_h = gip2owner_extra[dof_h]

        for (parent_gid, weight) in parent_constraints
            # Find local IP for this parent
            # We need to check if this rank OWNS this parent
            parent_dof_local = nothing
            for ip = 1:n_total_dofs
                if ip2gip_spa[ip] == parent_gid
                    parent_dof_local = ip
                    break
                end
            end

            parent_dof_local === nothing && continue
            gip2owner_extra[parent_dof_local] == rank || continue

            # This rank owns the parent → add to reverse map
            dest = get!(reverse_ghost_map, parent_dof_local, Tuple{Int, Int, Float64}[])
            push!(dest, (gid_h, owner_h, weight))
        end
    end

    n_interface = length(unique(dof_h ÷ n_ang + 1 for dof_h in keys(ghost_constraint_data)))
    @rankinfo rank "[$rank] Spatial ghost constraint data: $n_interface cross-rank hanging spatial nodes"
    @rankinfo rank "[$rank] Reverse ghost map: $(length(reverse_ghost_map)) owned parent DOFs have cross-rank children"

    return ghost_constraint_data, reverse_ghost_map
end

# ============================================================================
# Verification
# ============================================================================

"""
    verify_spatial_ghost_exchange(cache, rank) -> Bool

Verify that spatial ghost exchange completed successfully.
"""
function verify_spatial_ghost_exchange(cache::SpatialAMRCache, rank::Int)
    @rankinfo rank "[$rank] Verifying spatial ghost exchange..."

    n_local      = length(cache.parent_weights)
    n_cross_rank = length(cache.cross_rank_parent_weights)
    @rankinfo rank "[$rank] Local hanging nodes: $n_local, cross-rank hanging nodes: $n_cross_rank"

    # Basic integrity check for local constraints
    for (hanging_ip, weights) in cache.parent_weights
        if isempty(weights)
            error("[Rank $rank] Hanging node $hanging_ip has no parent weights")
        end
        w_sum = sum(w for (_, w) in weights)
        if abs(w_sum - 1.0) > 1e-8
            @warn "[$rank] Hanging node $hanging_ip weight sum = $w_sum (expected 1.0)"
        end
    end

    # Basic integrity check for cross-rank constraints
    for (hanging_ip, weights) in cache.cross_rank_parent_weights
        if isempty(weights)
            error("[Rank $rank] Cross-rank hanging node $hanging_ip has no parent weights")
        end
        w_sum = sum(w for (_, w) in weights)
        if abs(w_sum - 1.0) > 1e-8
            @warn "[$rank] Cross-rank hanging node $hanging_ip weight sum = $w_sum (expected 1.0)"
        end
    end

    @rankinfo rank "[$rank] ✓ Ghost exchange verification passed"
    return true
end

export exchange_spatial_ghosts, verify_spatial_ghost_exchange,
       build_spatial_ghost_constraint_data, fetch_parent_face_coordinates

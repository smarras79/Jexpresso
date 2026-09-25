function find_separated_rank_cut_elements(mesh; min_hops::Int=3)


    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
 

    # ── Same-rank adjacency (already validated) ─────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

     # ── Rank-cut identification (already validated via two-sided handshake) ──
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)


    # Every element on THIS rank that has at least one rank-cut face.
    rank_cut_elements = unique(iel for (iel,fid) in keys(neighbor_rank))
    println("Rank $rank: $(length(rank_cut_elements)) elements with a rank-cut face")



    if length(rank_cut_elements) < 2
        println("⚠ fewer than 2 rank-cut elements on this rank — cannot select a pair")
        return nothing
    end

    # BFS graph distance (via same-rank neighbor_elem edges only) from a
    # given starting element to every other element — capped, since we
    # only need to confirm distance >= min_hops, not exact shortest paths
    # everywhere.
    function bfs_distances(start_iel, max_depth)
        dist = Dict{Int,Int}(start_iel => 0)
        frontier = [start_iel]
        depth = 0
        while !isempty(frontier) && depth < max_depth
            depth += 1
            next_frontier = Int[]
            for iel in frontier
                for fid = 1:6
                    jel = neighbor_elem[iel][fid]
                    jel == zero(TInt) && continue
                    jel = Int(jel)
                    haskey(dist, jel) && continue
                    dist[jel] = depth
                    push!(next_frontier, jel)
                end
            end
            frontier = next_frontier
        end
        return dist
    end

    # Try each rank-cut element as a candidate "first" target, check if any
    # OTHER rank-cut element is at least min_hops away.
    for iel1 in rank_cut_elements
        dist = bfs_distances(iel1, min_hops + 2)   # search slightly beyond min_hops
        for iel2 in rank_cut_elements
            iel2 == iel1 && continue
            d = get(dist, iel2, typemax(Int))   # not found within search depth => "far enough"
            if d >= min_hops
                println("Rank $rank: found pair (iel1=$iel1, iel2=$iel2), graph distance >= $min_hops " *
                        "(actual: $(d == typemax(Int) ? "beyond search depth" : d))")
                return (iel1, iel2)
            end
        end
    end

    println("⚠ No pair of rank-cut elements found with graph distance >= $min_hops on this rank")
    return nothing
end

"""
    resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts, comm)

Generalizes resolve_rank_cut_neighbors to handle SPATIALLY NON-CONFORMING
rank-cut faces: a coarse face on one rank may correspond to 4 fine
children faces on another rank (or vice versa), not just a 1-to-1 match.

Returns:
  neighbor_rank      :: Dict{(iel,fid), rank}           — SAME as before (target rank)
  neighbor_kind      :: Dict{(iel,fid), Symbol}         — :conforming, :i_am_coarse, :i_am_fine
  neighbor_gel_group :: Dict{(iel,fid), Vector{TInt}}   — remote element(s): length 1 (conforming
                                                            or i_am_fine), or length 4 (i_am_coarse)
"""
function resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts::Vector{Tuple{TInt,TInt}}, comm;
                                         tol::Float64 = 1e-6)
    rank    = MPI.Comm_rank(comm)
    ngl     = mesh.ngl

    # ── Phase 1: pack (4 corner GIDs sorted [for a QUICK exact-match probe],
    # 4 corner COORDINATES [for containment], rank, iel, fid, gel) ──────────
    local_keys   = NTuple{4,TInt}[]
    local_coords = Vector{NTuple{4,NTuple{3,Float64}}}()
    local_iel = TInt[]; local_fid = TInt[]; local_gel = TInt[]

    for (iel, fid) in suspected_cuts
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
        key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        coords = ((mesh.x[n1],mesh.y[n1],mesh.z[n1]), (mesh.x[n2],mesh.y[n2],mesh.z[n2]),
                  (mesh.x[n3],mesh.y[n3],mesh.z[n3]), (mesh.x[n4],mesh.y[n4],mesh.z[n4]))
        push!(local_keys, key); push!(local_coords, coords)
        push!(local_iel, iel); push!(local_fid, fid); push!(local_gel, mesh.el2gel[iel])
    end

    # ── Phase 2: Allgatherv EVERYTHING (small — ~100 faces/rank) ─────────────
    n_local = TInt(length(local_keys))
    all_n = MPI.Allgather(n_local, comm)

    # Flatten: 4 keys + 12 coords (4 corners x 3) + rank,iel,fid,gel = 20 values
    flat_local = Float64[]
    for i = eachindex(local_keys)
        k = local_keys[i]; c = local_coords[i]
        append!(flat_local, Float64[k[1],k[2],k[3],k[4]])
        for corner in c
            append!(flat_local, Float64[corner[1],corner[2],corner[3]])
        end
        append!(flat_local, Float64[rank, local_iel[i], local_fid[i], local_gel[i]])
    end
    counts = all_n .* TInt(20)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    n_total = length(all_flat) ÷ 20
    entries = Vector{NamedTuple}(undef, n_total)
    for m = 1:n_total
        b = (m-1)*20
        key = (TInt(all_flat[b+1]),TInt(all_flat[b+2]),TInt(all_flat[b+3]),TInt(all_flat[b+4]))
        corners = ((all_flat[b+5],all_flat[b+6],all_flat[b+7]), (all_flat[b+8],all_flat[b+9],all_flat[b+10]),
                   (all_flat[b+11],all_flat[b+12],all_flat[b+13]), (all_flat[b+14],all_flat[b+15],all_flat[b+16]))
        entries[m] = (key=key, corners=corners, rank=TInt(all_flat[b+17]), iel=TInt(all_flat[b+18]),
                      fid=TInt(all_flat[b+19]), gel=TInt(all_flat[b+20]))
    end

    # ── Phase 3: for each of MY faces, find matches among ALL entries ───────
    neighbor_rank = Dict{Tuple{TInt,TInt}, TInt}()
    neighbor_kind = Dict{Tuple{TInt,TInt}, Symbol}()
    neighbor_gel_group = Dict{Tuple{TInt,TInt}, Vector{TInt}}()
    n_unresolved = 0

    function bbox(corners)
        xs = [c[1] for c in corners]; ys = [c[2] for c in corners]; zs = [c[3] for c in corners]
        return (minimum(xs),maximum(xs), minimum(ys),maximum(ys), minimum(zs),maximum(zs))
    end

    for i = eachindex(local_keys)
        my_key = local_keys[i]; my_iel = local_iel[i]; my_fid = local_fid[i]
        my_bbox = bbox(local_coords[i])

        # Exact match probe (conforming case — same as before)
        exact_matches = [e for e in entries if e.key == my_key && e.rank != rank]
        if length(exact_matches) == 1
            neighbor_rank[(my_iel,my_fid)] = exact_matches[1].rank
            neighbor_kind[(my_iel,my_fid)] = :conforming
            neighbor_gel_group[(my_iel,my_fid)] = [exact_matches[1].gel]
            continue
        end

        # Containment probes: find OTHER-rank entries whose bbox is
        # STRICTLY inside mine (I am coarse), or whose bbox STRICTLY
        # contains mine (I am fine) — same principle as
        # build_angular_ncf_correspondence, applied to (x,y,z).
        candidates_inside = NamedTuple[]   # I am coarse
        candidates_containing = NamedTuple[]  # I am fine
        for e in entries
            e.rank == rank && continue
            ebbox = bbox(e.corners)
            inside = (ebbox[1] >= my_bbox[1]-tol && ebbox[2] <= my_bbox[2]+tol &&
                      ebbox[3] >= my_bbox[3]-tol && ebbox[4] <= my_bbox[4]+tol &&
                      ebbox[5] >= my_bbox[5]-tol && ebbox[6] <= my_bbox[6]+tol)
            containing = (my_bbox[1] >= ebbox[1]-tol && my_bbox[2] <= ebbox[2]+tol &&
                          my_bbox[3] >= ebbox[3]-tol && my_bbox[4] <= ebbox[4]+tol &&
                          my_bbox[5] >= ebbox[5]-tol && my_bbox[6] <= ebbox[6]+tol)
            inside && push!(candidates_inside, e)
            containing && push!(candidates_containing, e)
        end

        if length(candidates_inside) == 4
            neighbor_rank[(my_iel,my_fid)] = candidates_inside[1].rank
            neighbor_kind[(my_iel,my_fid)] = :i_am_coarse
            neighbor_gel_group[(my_iel,my_fid)] = [c.gel for c in candidates_inside]
        elseif length(candidates_containing) == 1
            neighbor_rank[(my_iel,my_fid)] = candidates_containing[1].rank
            neighbor_kind[(my_iel,my_fid)] = :i_am_fine
            neighbor_gel_group[(my_iel,my_fid)] = [candidates_containing[1].gel]
        else
            n_unresolved += 1
            @warn "resolve_rank_cut_neighbors_ncf: face (iel=$my_iel,fid=$my_fid) matched " *
                  "$(length(candidates_inside)) 'inside' and $(length(candidates_containing)) " *
                  "'containing' candidates — expected exactly 4 or 1. Unresolved."
        end
    end

    return neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved
end

function test_resolve_rank_cut_neighbors_ncf(inputs, nparts, distribute, sem, partitioned_model;
                                              spatial_target_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    println("Rank $rank: applying spatial refinement at iel=$spatial_target_iel " *
            "(only on rank 0, per test design)...")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    # ── Same-rank adjacency + rank-cut identification (already validated) ───
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        push!(ncf_as_child_keys, (child_iel, child_fid))
        push!(ncf_as_parent_keys, (parent_iel, parent_fid))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
    	end
    end
    
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    # ── The actual new function ───────────────────────────────────────────
    neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved =
        resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts, comm)

    println("Rank $rank: $n_unresolved unresolved faces  " * (n_unresolved==0 ? "✓" : "✗"))
    println("Rank $rank: kind breakdown:")
    for kind in (:conforming, :i_am_coarse, :i_am_fine)
        n = count(==(kind), values(neighbor_kind))
        println("  $kind: $n")
    end

    # Print the specific NCF entries for inspection
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming && continue
        println("  (iel=$iel, fid=$fid): kind=$kind, target_rank=$(neighbor_rank[(iel,fid)]), " *
                "remote_gels=$(neighbor_gel_group[(iel,fid)])")
    end

    return mesh, neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved
end

"""
    exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)

Exchanges each rank-cut face's FULL angular-mesh description (extra_coords,
extra_connijk, extra_nop, extra_nelem) across ranks — enough to run the
SAME per-node containment logic build_angular_dependency_graph uses
locally, against a remote element's angular mesh.

Reuses the SELF-DESCRIBING matching pattern already validated in
build_rank_cut_exchange: every exchanged float is tagged with
(sender_gel, position_within_payload), and the receiver matches incoming
data back to its own faces by that tag — NOT by positional trust, which
was proven unsafe for exactly this kind of one-shot MPI exchange earlier
in this effort.
"""
#=function exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, comm)
    rank = MPI.Comm_rank(comm)

    send_keys    = collect(keys(neighbor_rank))
    send_targets = [Int(neighbor_rank[k]) for k in send_keys]

    # ── Build each face's flattened payload + remember MY OWN gel (so the
    # remote side can tag its reply, and so I can match incoming data back
    # to MY OWN expected remote gel — the SAME disambiguation principle as
    # build_rank_cut_exchange's sender-gel tagging). ────────────────────────
    payloads   = Vector{Vector{Float64}}(undef, length(send_keys))
    my_gels    = Vector{TInt}(undef, length(send_keys))
    for (idx,(iel,fid)) in enumerate(send_keys)
        em = mesh.extra_mesh[iel]
        my_gels[idx] = mesh.el2gel[iel]
        payload = Float64[Float64(em.extra_nelem), Float64(em.extra_npoin)]
        for p = 1:em.extra_npoin
            append!(payload, (em.extra_coords[1,p], em.extra_coords[2,p]))
        end
        for e = 1:em.extra_nelem
            n = em.extra_nop[e]+1
            push!(payload, Float64(n))
            for a=1:n, b=1:n
                push!(payload, Float64(em.extra_connijk[e,a,b]))
            end
        end
        payloads[idx] = payload
    end

    # ── Flatten: one entry per (face, position-within-payload), tagged
    # with MY gel and MY position — self-describing, exactly like
    # build_rank_cut_exchange's (global_id, sender_gel) tuple. ──────────────
    flat_targets = Int[]
    flat_vals    = Float64[]
    flat_gel_tag = TInt[]     # MY gel, repeated for every float in this face's payload
    flat_pos_tag = Int[]      # position within THIS face's payload
    flat_send_keys = Tuple{Int,Int,Int}[]   # (face_idx, position) -> for MY OWN bookkeeping

    for (idx, payload) in enumerate(payloads)
        for (pos, v) in enumerate(payload)
            push!(flat_targets, send_targets[idx])
            push!(flat_vals, v)
            push!(flat_gel_tag, my_gels[idx])
            push!(flat_pos_tag, pos)
            push!(flat_send_keys, (idx, pos, 0))
        end
    end

    # ── Three parallel caches: value, gel-tag, position-tag — SAME shape
    # as the value+gel two-cache pattern already validated for spatial NCF. ──
    val_cache = SendReceiveCache(comm, flat_vals, flat_targets)
    send_and_receive!(val_cache, flat_vals)
    gel_cache = SendReceiveCache(comm, Float64.(flat_gel_tag), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_gel_tag))
    pos_cache = SendReceiveCache(comm, Float64.(flat_pos_tag), flat_targets)
    send_and_receive!(pos_cache, Float64.(flat_pos_tag))

    # ── Reconstruct: group received (value, sender_gel, position) triples
    # by sender_gel, sort by position, rebuild each remote face's payload
    # in order — self-describing, no positional trust in the RECEIVE order
    # across different senders/faces. ────────────────────────────────────────
    n_recv = length(val_cache.combined_recv_data)
    by_gel = Dict{TInt, Vector{Tuple{Int,Float64}}}()  # remote_gel -> [(position, value), ...]
    for m = 1:n_recv
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        rpos = Int(round(Int, pos_cache.combined_recv_data[m]))
        rval = val_cache.combined_recv_data[m]
        push!(get!(by_gel, rgel, Tuple{Int,Float64}[]), (rpos, rval))
    end

    remote_angular_mesh = Dict{TInt, NamedTuple}()  # remote_gel -> (extra_nelem, extra_npoin, extra_coords, extra_connijk, extra_nop)
    for (rgel, pairs) in by_gel
        sorted_pairs = sort(pairs, by = x -> x[1])
        vals = [p[2] for p in sorted_pairs]

        nelem = Int(round(Int, vals[1]))
        npoin = Int(round(Int, vals[2]))
        coords = zeros(2, npoin)
        offset = 2
        for p = 1:npoin
            coords[1,p] = vals[offset+1]; coords[2,p] = vals[offset+2]
            offset += 2
        end
        nop = zeros(Int, nelem)
        connijk = Vector{Array{Int,2}}(undef, nelem)
        for e = 1:nelem
            n = Int(round(Int, vals[offset+1])); offset += 1
            nop[e] = n - 1
            grid = zeros(Int, n, n)
            for a=1:n, b=1:n
                offset += 1
                grid[a,b] = Int(round(Int, vals[offset]))
            end
            connijk[e] = grid
        end

        remote_angular_mesh[rgel] = (extra_nelem=nelem, extra_npoin=npoin, extra_coords=coords,
                                       extra_connijk=connijk, extra_nop=nop)
    end

    return remote_angular_mesh
end=#

function exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, comm)
    rank = MPI.Comm_rank(comm)

    send_keys    = collect(keys(neighbor_rank))
    send_targets = [Int(neighbor_rank[k]) for k in send_keys]

    payloads   = Vector{Vector{Float64}}(undef, length(send_keys))
    my_gels    = Vector{TInt}(undef, length(send_keys))
    for (idx,(iel,fid)) in enumerate(send_keys)
        em = mesh.extra_mesh[iel]
        my_gels[idx] = mesh.el2gel[iel]
        payload = Float64[Float64(em.extra_nelem), Float64(em.extra_npoin)]
        for p = 1:em.extra_npoin
            append!(payload, (em.extra_coords[1,p], em.extra_coords[2,p]))
        end
        for e = 1:em.extra_nelem
            n = em.extra_nop[e]+1
            push!(payload, Float64(n))
            for a=1:n, b=1:n
                push!(payload, Float64(em.extra_connijk[e,a,b]))
            end
        end
        payloads[idx] = payload
    end

    flat_targets = Int[]
    flat_vals    = Float64[]
    flat_gel_tag = TInt[]
    flat_pos_tag = Int[]

    for (idx, payload) in enumerate(payloads)
        for (pos, v) in enumerate(payload)
            push!(flat_targets, send_targets[idx])
            push!(flat_vals, v)
            push!(flat_gel_tag, my_gels[idx])
            push!(flat_pos_tag, pos)
        end
    end

    val_cache = SendReceiveCache(comm, flat_vals, flat_targets)
    send_and_receive!(val_cache, flat_vals)
    gel_cache = SendReceiveCache(comm, Float64.(flat_gel_tag), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_gel_tag))
    pos_cache = SendReceiveCache(comm, Float64.(flat_pos_tag), flat_targets)
    send_and_receive!(pos_cache, Float64.(flat_pos_tag))

    # ── Reconstruct with DEDUPLICATION: a remote element bordering my rank
    # via MULTIPLE faces sends its full angular-mesh payload once PER
    # FACE, each restarting its own position counter at 1 — meaning the
    # SAME (rgel,pos) pair can legitimately arrive more than once, with
    # IDENTICAL values. Deduplicate by (rgel,pos) BEFORE reconstructing,
    # so each position is used exactly once regardless of how many
    # faces triggered the send. ──────────────────────────────────────────────
    n_recv = length(val_cache.combined_recv_data)
    by_gel_pos = Dict{Tuple{TInt,Int}, Float64}()
    n_dup_detected = 0
    for m = 1:n_recv
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        rpos = Int(round(Int, pos_cache.combined_recv_data[m]))
        rval = val_cache.combined_recv_data[m]
        key = (rgel, rpos)
        if haskey(by_gel_pos, key)
            isapprox(by_gel_pos[key], rval; atol=1e-8) ||
                @warn "Rank $rank: duplicate (gel=$rgel,pos=$rpos) with DIFFERING values: " *
                      "$(by_gel_pos[key]) vs $rval"
            n_dup_detected += 1
        else
            by_gel_pos[key] = rval
        end
    end
    n_dup_detected > 0 && println("Rank $rank: deduplicated $n_dup_detected repeated (gel,pos) entries " *
                                    "(expected when a remote element borders this rank via multiple faces)")

    by_gel = Dict{TInt, Vector{Tuple{Int,Float64}}}()
    for ((rgel,rpos), rval) in by_gel_pos
        push!(get!(by_gel, rgel, Tuple{Int,Float64}[]), (rpos, rval))
    end

    remote_angular_mesh = Dict{TInt, NamedTuple}()
    for (rgel, pairs) in by_gel
        sorted_pairs = sort(pairs, by = x -> x[1])
        vals = [p[2] for p in sorted_pairs]

        nelem = Int(round(Int, vals[1]))
        npoin = Int(round(Int, vals[2]))
        coords = zeros(2, npoin)
        offset = 2
        for p = 1:npoin
            coords[1,p] = vals[offset+1]; coords[2,p] = vals[offset+2]
            offset += 2
        end
        nop = zeros(Int, nelem)
        connijk = Vector{Array{Int,2}}(undef, nelem)
        for e = 1:nelem
            n = Int(round(Int, vals[offset+1])); offset += 1
            nop[e] = n - 1
            grid = zeros(Int, n, n)
            for a=1:n, b=1:n
                offset += 1
                grid[a,b] = Int(round(Int, vals[offset]))
            end
            connijk[e] = grid
        end

        remote_angular_mesh[rgel] = (extra_nelem=nelem, extra_npoin=npoin, extra_coords=coords,
                                       extra_connijk=connijk, extra_nop=nop)
    end

    return remote_angular_mesh
end

function test_exchange_angular_mesh_full(inputs, nparts, distribute, sem, partitioned_model;
                                          angular_target_iel::Int = 315, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    println("Rank $rank: applying angular refinement at iel=$angular_target_iel " *
            "(only meaningful on the rank owning this element)...")

    mesh = sem.mesh
    if angular_target_iel <= mesh.nelem
        adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
        verify_angular_phi_seam_conformity(mesh, angular_target_iel)
    end

    # ── Same-rank adjacency + rank-cut identification (ORDINARY, conforming
    # version — 442's face itself is spatially conforming, only its
    # ANGULAR mesh differs from its neighbor). ────────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved rank-cut faces")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts present — aborting")
        return nothing
    end

    # ── Checksum-based verification (as designed) ────────────────────────────
    function local_connijk_sum(em)
        s = 0.0
        for e = 1:em.extra_nelem
            n = em.extra_nop[e] + 1
            for a = 1:n, b = 1:n
                s += em.extra_connijk[e,a,b]
            end
        end
        return s
    end

    function local_coords_sum(em)
    	s = 0.0
    	for p = 1:em.extra_npoin
            s += em.extra_coords[1,p] + em.extra_coords[2,p]
    	end
        return s
    end

    my_summaries = Dict{Tuple{Int,Int}, Tuple{Int,Int,Float64,Float64}}()
    for (iel,fid) in keys(neighbor_rank)
        em = mesh.extra_mesh[iel]
        s_coords = local_coords_sum(em)
        s_connijk = local_connijk_sum(em)
        my_summaries[(iel,fid)] = (em.extra_nelem, em.extra_npoin, s_coords, s_connijk)
    end

    send_keys = collect(keys(neighbor_rank))
    send_targets = [Int(neighbor_rank[k]) for k in send_keys]
    my_gels = [TInt(mesh.el2gel[k[1]]) for k in send_keys]

    nelem_vals   = Float64[my_summaries[k][1] for k in send_keys]
    npoin_vals   = Float64[my_summaries[k][2] for k in send_keys]
    coords_vals  = Float64[my_summaries[k][3] for k in send_keys]
    connijk_vals = Float64[my_summaries[k][4] for k in send_keys]

    gel_cache = SendReceiveCache(comm, Float64.(my_gels), send_targets)
    send_and_receive!(gel_cache, Float64.(my_gels))
    nelem_cache = SendReceiveCache(comm, nelem_vals, send_targets)
    send_and_receive!(nelem_cache, nelem_vals)
    npoin_cache = SendReceiveCache(comm, npoin_vals, send_targets)
    send_and_receive!(npoin_cache, npoin_vals)
    coords_cache = SendReceiveCache(comm, coords_vals, send_targets)
    send_and_receive!(coords_cache, coords_vals)
    connijk_cache = SendReceiveCache(comm, connijk_vals, send_targets)
    send_and_receive!(connijk_cache, connijk_vals)

    true_remote_summary = Dict{TInt, Tuple{Int,Int,Float64,Float64}}()
    for m = 1:length(gel_cache.combined_recv_data)
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        true_remote_summary[rgel] = (Int(round(Int, nelem_cache.combined_recv_data[m])),
                                       Int(round(Int, npoin_cache.combined_recv_data[m])),
                                       coords_cache.combined_recv_data[m],
                                       connijk_cache.combined_recv_data[m])
    end

    # ── Run the ACTUAL full exchange under test ──────────────────────────────
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)

    n_checked = 0
    n_mismatch = 0
    for (rgel, recon) in remote_angular_mesh
        n_checked += 1
        haskey(true_remote_summary, rgel) || (println("Rank $rank: ⚠ no true summary for gel=$rgel"); continue)
        true_nelem, true_npoin, true_coords_sum, true_connijk_sum = true_remote_summary[rgel]

        recon_coords_sum = sum(recon.extra_coords)
        recon_connijk_sum = 0.0
        for e = 1:recon.extra_nelem
            n = recon.extra_nop[e] + 1
            for a=1:n, b=1:n
                recon_connijk_sum += recon.extra_connijk[e][a,b]
            end
        end

        ok = (recon.extra_nelem == true_nelem) && (recon.extra_npoin == true_npoin) &&
             isapprox(recon_coords_sum, true_coords_sum; atol=1e-6) &&
             isapprox(recon_connijk_sum, true_connijk_sum; atol=1e-6)
        if !ok
            n_mismatch += 1
            println("Rank $rank: ✗ MISMATCH for remote gel=$rgel: " *
                    "recon=($(recon.extra_nelem),$(recon.extra_npoin),$recon_coords_sum,$recon_connijk_sum) " *
                    "true=($true_nelem,$true_npoin,$true_coords_sum,$true_connijk_sum)")
        end
    end
    println("Rank $rank: checked $n_checked reconstructed remote meshes, $n_mismatch mismatches  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))

    return remote_angular_mesh
end

function test_multirank_schedule_spatial_ncf_full(inputs, nparts, distribute, sem, partitioned_model;
                                                    spatial_target_iel::Int = 415,
                                                    Ωx::Float64 = 0.3, Ωy::Float64 = 0.5, Ωz::Float64 = 0.8,
                                                    max_rounds::Int = 200)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    println("Rank $rank: applying spatial refinement at iel=$spatial_target_iel " *
            "(only on rank 0, per test design)...")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Same-rank spatial NCF tables ──────────────────────────────────────────
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
    end
    println("Rank $rank: spatial NCF: $(mesh.num_ncf) same-rank entries")

    # ── Rank-cut identification, EXCLUDING same-rank NCF faces from the
    # suspected-cut list (the earlier fix) ───────────────────────────────────
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved =
        resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved rank-cut faces")
    if n_unresolved > 0
        println("Rank $rank: ⚠ unresolved rank-cuts — aborting")
        return nothing
    end

    # ── The actual scheduling test ─────────────────────────────────────────
    ngl = mesh.ngl
    nelem = mesh.nelem
    is_inflow = Vector{Bool}(undef, 6); is_active = Vector{Bool}(undef, 6)
    processed = falses(nelem)
    round_solved = fill(-1, nelem)

    flat_iel   = Int[]
    flat_fid   = Int[]
    flat_targets = Int[]
    flat_my_gel = TInt[]
    flat_expected_remote_gel = TInt[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        for remote_gel in neighbor_gel_group[(iel,fid)]
            push!(flat_iel, iel); push!(flat_fid, fid)
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_expected_remote_gel, remote_gel)
        end
    end

    n_flat = length(flat_iel)
    round_cache = SendReceiveCache(comm, zeros(Int, n_flat), flat_targets)
    mygel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    send_and_receive!(mygel_cache, Float64.(flat_my_gel))

    gel_to_flat_positions = Dict{TInt,Vector{Int}}()
    for i = 1:n_flat
        push!(get!(gel_to_flat_positions, flat_expected_remote_gel[i], Int[]), i)
    end

    n_recv = length(mygel_cache.combined_recv_data)
    recv_perm = Vector{Vector{Int}}(undef, n_recv)
    for m = 1:n_recv
        recv_gel = TInt(round(Int, mygel_cache.combined_recv_data[m]))
        recv_perm[m] = get(gel_to_flat_positions, recv_gel, Int[])
    end

    remote_round_flat = fill(-1, n_flat)
    progress_history = Int[]

    round_iter = 0
    while true
        round_iter += 1
        if round_iter > max_rounds
            rank == 0 && println("STALL: $(progress_history[max(1,end-9):end])")
            error("test_multirank_schedule_spatial_ncf_full: exceeded max_rounds=$max_rounds")
        end

        for iel = 1:nelem
            processed[iel] && continue
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)
            ready = true
            for fid = 1:6
                (is_inflow[fid] && is_active[fid]) || continue
                jel = neighbor_elem[iel][fid]

                if jel != zero(TInt)
                    processed[jel] || (ready = false; break)
                    continue
                end

                if haskey(ncf_as_child, (iel,fid))
                    ncf_idx = ncf_as_child[(iel,fid)]
                    p_iel = Int(mesh.pip[ncf_idx])
                    processed[p_iel] || (ready = false)
                elseif haskey(ncf_as_parent, (iel,fid))
                    for ncf_idx in ncf_as_parent[(iel,fid)]
                        c_iel = Int(mesh.cip[ncf_idx])
                        processed[c_iel] || (ready = false)
                    end
                elseif haskey(neighbor_rank, (iel,fid))
                    for i = 1:n_flat
                        (flat_iel[i] == iel && flat_fid[i] == fid) || continue
                        remote_round_flat[i] == -1 && (ready = false)
                    end
                end
                ready || break
            end
            if ready
                processed[iel] = true
                round_solved[iel] = round_iter
            end
        end

        send_vals = [round_solved[flat_iel[i]] for i = 1:n_flat]
        send_and_receive!(round_cache, send_vals)
        for m = 1:length(round_cache.combined_recv_data)
            for pos in recv_perm[m]
                remote_round_flat[pos] = round_cache.combined_recv_data[m]
            end
        end

        n_unprocessed_total = MPI.Allreduce(count(!,processed), +, comm)
        push!(progress_history, n_unprocessed_total)
        println("Rank $rank, round $round_iter: unprocessed=$n_unprocessed_total")
        n_unprocessed_total == 0 && break
    end

    rank == 0 && println("Converged in $round_iter rounds")
    return round_iter, processed
end

"""
Key CHANGE: angular_deps_cr is now keyed by (iel,fid,d), not (iel,d) —
a single direction can have GENUINELY DIFFERENT, independent angular-NCF
relationships on different faces simultaneously (one inflow face's
remote element need not be the same as another inflow face's).
"""
#=function build_angular_dependency_graph_crossrank(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                                    neighbor_gel::Dict{Tuple{TInt,TInt},TInt},
                                                    remote_angular_mesh::Dict{TInt,NamedTuple};
                                                    tol::Float64 = 1e-6)
    angular_deps_cr = Dict{Tuple{Int,Int,Int}, Tuple{TInt,Vector{Int}}}()   # (iel,fid,d) now

    for ((iel,fid), target_rank) in neighbor_rank
        haskey(neighbor_gel, (iel,fid)) || continue
        remote_gel = neighbor_gel[(iel,fid)]
        haskey(remote_angular_mesh, remote_gel) || continue
        rem = remote_angular_mesh[remote_gel]
        em = mesh.extra_mesh[iel]

        if em.extra_nelem == rem.extra_nelem
            same = true
            for p = 1:min(em.extra_npoin, rem.extra_npoin)
                if !isapprox(em.extra_coords[1,p], rem.extra_coords[1,p]; atol=tol) ||
                   !isapprox(em.extra_coords[2,p], rem.extra_coords[2,p]; atol=tol)
                    same = false; break
                end
            end
            (same && em.extra_npoin == rem.extra_npoin) && continue
        end

        rem_bounds = [get_angular_element_bounds_ragged(rem.extra_coords, rem.extra_connijk, rem.extra_nop, e)
                      for e = 1:rem.extra_nelem]
        for ip_ext_i = 1:em.extra_npoin
            θi = em.extra_coords[1, ip_ext_i]
            ϕi = em.extra_coords[2, ip_ext_i]

            matched = 0
            for p = 1:rem.extra_npoin
                if isapprox(θi, rem.extra_coords[1,p]; atol=tol) && isapprox(ϕi, rem.extra_coords[2,p]; atol=tol)
                    matched = p; break
                end
            end
            if matched != 0
                angular_deps_cr[(iel, fid, ip_ext_i)] = (remote_gel, [matched])
                continue
            end

            for e_r = 1:rem.extra_nelem
                θmin,θmax,ϕmin,ϕmax = rem_bounds[e_r]
                near_seam = isapprox(ϕmax, 2π; atol=1e-6)
                ϕi_fix = (near_seam && ϕi < 1e-8) ? 2π : ϕi
                inside = (θi >= θmin-tol && θi <= θmax+tol && ϕi_fix >= ϕmin-tol && ϕi_fix <= ϕmax+tol)
                inside || continue
                n = rem.extra_nop[e_r]+1
                rem_nodes = [rem.extra_connijk[e_r][a,b] for a=1:n, b=1:n] |> vec
                angular_deps_cr[(iel, fid, ip_ext_i)] = (remote_gel, rem_nodes)
                break
            end
        end
    end

    return angular_deps_cr
end=#

function build_angular_dependency_graph_crossrank(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                                    neighbor_gel::Dict{Tuple{TInt,TInt},TInt},
                                                    remote_angular_mesh::Dict{TInt,NamedTuple};
                                                    tol::Float64 = 1e-6,
                                                    trace_iel::Int = 0, trace_fid::Int = 0, trace_d::Int = 0)
    rank = MPI.Comm_rank(get_mpi_comm())
    angular_deps_cr = Dict{Tuple{Int,Int,Int}, Tuple{TInt,Vector{Int}}}()

    for ((iel,fid), target_rank) in neighbor_rank
        haskey(neighbor_gel, (iel,fid)) || continue
        remote_gel = neighbor_gel[(iel,fid)]
        haskey(remote_angular_mesh, remote_gel) || continue
        rem = remote_angular_mesh[remote_gel]
        em = mesh.extra_mesh[iel]

        is_trace_face = (iel==trace_iel && fid==trace_fid)
        if is_trace_face
            println("Rank $rank: [TRACE] entering (iel=$iel,fid=$fid) -> remote_gel=$remote_gel " *
                    "em.extra_nelem=$(em.extra_nelem) em.extra_npoin=$(em.extra_npoin) " *
                    "rem.extra_nelem=$(rem.extra_nelem) rem.extra_npoin=$(rem.extra_npoin)")
        end

        if em.extra_nelem == rem.extra_nelem
            same = true
            for p = 1:min(em.extra_npoin, rem.extra_npoin)
                if !isapprox(em.extra_coords[1,p], rem.extra_coords[1,p]; atol=tol) ||
                   !isapprox(em.extra_coords[2,p], rem.extra_coords[2,p]; atol=tol)
                    same = false; break
                end
            end
            if same && em.extra_npoin == rem.extra_npoin
                is_trace_face && println("Rank $rank: [TRACE] SKIPPED as identical mesh (nelem AND npoin match, all coords match)")
                continue
            end
        end

        rem_bounds = [get_angular_element_bounds_ragged(rem.extra_coords, rem.extra_connijk, rem.extra_nop, e)
                      for e = 1:rem.extra_nelem]
        for ip_ext_i = 1:em.extra_npoin
            θi = em.extra_coords[1, ip_ext_i]
            ϕi = em.extra_coords[2, ip_ext_i]

            is_trace_dir = is_trace_face && (ip_ext_i == trace_d)
            is_trace_dir && println("Rank $rank: [TRACE] processing ip_ext_i=$ip_ext_i: θi=$θi, ϕi=$ϕi")

            matched = 0
            for p = 1:rem.extra_npoin
                if isapprox(θi, rem.extra_coords[1,p]; atol=tol) && isapprox(ϕi, rem.extra_coords[2,p]; atol=tol)
                    matched = p; break
                end
            end
            if matched != 0
                is_trace_dir && println("Rank $rank: [TRACE] EXACT MATCH found at remote p=$matched " *
                                          "(rem.extra_coords=($(rem.extra_coords[1,matched]),$(rem.extra_coords[2,matched])))")
                angular_deps_cr[(iel, fid, ip_ext_i)] = (remote_gel, [matched])
                continue
            end
            is_trace_dir && println("Rank $rank: [TRACE] no exact match — checking containment against $(rem.extra_nelem) remote elements")

            found_containing = false
            for e_r = 1:rem.extra_nelem
                θmin,θmax,ϕmin,ϕmax = rem_bounds[e_r]
                near_seam = isapprox(ϕmax, 2π; atol=1e-6)
                ϕi_fix = (near_seam && ϕi < 1e-8) ? 2π : ϕi
                inside = (θi >= θmin-tol && θi <= θmax+tol && ϕi_fix >= ϕmin-tol && ϕi_fix <= ϕmax+tol)
                if is_trace_dir
                    println("Rank $rank: [TRACE]   e_r=$e_r bounds=(θ:$θmin-$θmax, ϕ:$ϕmin-$ϕmax) " *
                             "near_seam=$near_seam ϕi_fix=$ϕi_fix inside=$inside")
                end
                inside || continue
                found_containing = true
                n = rem.extra_nop[e_r]+1
                rem_nodes = [rem.extra_connijk[e_r][a,b] for a=1:n, b=1:n] |> vec
                angular_deps_cr[(iel, fid, ip_ext_i)] = (remote_gel, rem_nodes)
                is_trace_dir && println("Rank $rank: [TRACE] MATCHED containing e_r=$e_r, rem_nodes=$rem_nodes")
                break
            end
            is_trace_dir && !found_containing && println("Rank $rank: [TRACE] ⚠⚠⚠ NO containing element found for " *
                                                          "(iel=$iel,fid=$fid,d=$ip_ext_i) — NO entry created, " *
                                                          "will fall through to ORDINARY fallback incorrectly!")
        end
    end

    return angular_deps_cr
end


function test_multirank_schedule_angular_ncf_full(inputs, nparts, distribute, sem, partitioned_model;
                                                    angular_target_iel::Int = 442, angular_target_e_ext::Int = 4,
                                                    Ωx::Float64 = 0.3, Ωy::Float64 = 0.5, Ωz::Float64 = 0.8,
                                                    max_rounds::Int = 200)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    metrics = sem.metrics
    if angular_target_iel <= mesh.nelem
        adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
        verify_angular_phi_seam_conformity(mesh, angular_target_iel)
    end
    println("Rank $rank: mesh has $(mesh.nelem) elements")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    n_unresolved > 0 && (println("⚠ aborting"); return nothing)

    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    println("Rank $rank: $(length(angular_deps_cr)) cross-rank angular NCF (iel,fid, d) entries")

    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    round_solved = [fill(-1, n_dirs[iel]) for iel = 1:mesh.nelem]
    total_pairs = sum(n_dirs)

    flat_targets = Int[]; flat_my_gel = TInt[]; flat_my_pos = Int[]
    flat_src_iel = Int[]; flat_src_d = Int[]
    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_my_pos, d)
            push!(flat_src_iel, iel)
            push!(flat_src_d, d)
        end
    end
    n_flat = length(flat_targets)

    round_cache = SendReceiveCache(comm, zeros(Int, n_flat), flat_targets)
    gel_cache   = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    pos_cache   = SendReceiveCache(comm, Float64.(flat_my_pos), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(pos_cache, Float64.(flat_my_pos))

    wishlist = Dict{Tuple{TInt,Int}, Vector{Tuple{Int,Int}}}()
    for ((iel,fid,d), (remote_gel, remote_positions)) in angular_deps_cr
        for rp in remote_positions
            push!(get!(wishlist, (remote_gel,rp), Tuple{Int,Int}[]), (iel,d))
        end
    end
    n_recv = length(gel_cache.combined_recv_data)
    recv_match = Vector{Vector{Tuple{Int,Int}}}(undef, n_recv)
    for m = 1:n_recv
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        rpos = Int(round(Int, pos_cache.combined_recv_data[m]))
        recv_match[m] = get(wishlist, (rgel,rpos), Tuple{Int,Int}[])
    end

    remote_round_by_pair = Dict{Tuple{Int,Int}, Int}()
    for ((iel,fid,d), _) in angular_deps_cr
        remote_round_by_pair[(iel,d)] = -1
    end

    n_solved_total = 0
    iter_round = 0
    progress_history = Int[]
    just_resolved_stall = false

    while n_solved_total < total_pairs
        iter_round += 1

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωxd, Ωyd, Ωzd, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel = neighbor_elem[iel][fid]
                    if jel != zero(TInt)
                        key = (iel,d)
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                for dn in dep_nodes
                                    solved[jel][dn] || (ready = false)
                                end
                            end
                        else
                            (d <= n_dirs[jel] && !solved[jel][d]) && (ready = false)
                        end
                        continue
                    end
                    if haskey(neighbor_rank, (iel,fid))
                        if haskey(angular_deps_cr, (iel,fid,d))
                            remote_round_by_pair[(iel,d)] == -1 && (ready = false)
                        end
                        # else: ordinary conforming cross-rank exact-match
                        # direction — no explicit gate here (flagged
                        # elsewhere as a known incompleteness; not the
                        # focus of this fix).
                    end
                    ready || break
                end
                if ready
                    solved[iel][d] = true
                    round_solved[iel][d] = iter_round
                    n_newly_solved += 1
                end
            end
        end
        n_solved_total = sum(count(identity, solved[iel]) for iel = 1:mesh.nelem)

        send_vals = [round_solved[flat_src_iel[i]][flat_src_d[i]] for i = 1:n_flat]
        send_and_receive!(round_cache, send_vals)
        for m = 1:length(round_cache.combined_recv_data)
            for (iel,d) in recv_match[m]
                remote_round_by_pair[(iel,d)] = round_cache.combined_recv_data[m]
            end
        end

        n_unprocessed_total = MPI.Allreduce(total_pairs - n_solved_total, +, comm)
        push!(progress_history, n_unprocessed_total)
        println("Rank $rank, round $iter_round: solved_this_round=$n_newly_solved, global_unprocessed=$n_unprocessed_total")
        n_unprocessed_total == 0 && break

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        if n_newly_solved_global == 0
            if just_resolved_stall
                debug_angular_ncf_crossrank_stall(mesh, metrics, neighbor_elem, neighbor_rank,
                                                   angular_deps, angular_deps_cr,
                                                   solved, n_dirs, remote_round_by_pair, rank)
                error("test_multirank_schedule_angular_ncf_full: made NO progress even immediately " *
                      "after distributed Gauss-Seidel resolution — genuine deadlock.")
            end
            resolve_stall_distributed_gauss_seidel!(mesh, metrics, neighbor_elem, neighbor_rank,
                                                     angular_deps, angular_deps_cr,
                                                     solved, round_solved, n_dirs, iter_round, comm)
            just_resolved_stall = true
            n_solved_total = sum(count(identity, solved[iel]) for iel = 1:mesh.nelem)
        else
            just_resolved_stall = false
        end

        if iter_round > max_rounds
            error("test_multirank_schedule_angular_ncf_full: exceeded max_rounds=$max_rounds " *
                  "(outer safety cap, not the stall-detection mechanism itself)")
        end
    end

    rank == 0 && println("Converged in $iter_round rounds")
    return iter_round, solved
end

"""
Adapter: computes bounds directly from the ragged (Vector{Matrix{Int}})
reconstructed representation, WITHOUT needing get_angular_element_bounds'
flat-3D indexing convention at all — same logic, different indexing.
"""
function get_angular_element_bounds_ragged(coords::Matrix{Float64}, connijk_ragged::Vector{Matrix{Int}},
                                            nop::Vector{Int}, e_ext::Int)
    n = nop[e_ext] + 1
    grid = connijk_ragged[e_ext]
    θs = Float64[]; ϕs_raw = Float64[]
    for i = 1:n, j = 1:n
        ip = grid[i,j]
        push!(θs, coords[1,ip])
        push!(ϕs_raw, coords[2,ip])
    end

    # ── SEAM FIX: if this element's own nodes contain BOTH values near 0
    # AND values well above 0 (e.g. > π), this element straddles the
    # φ=0/2π seam — unwrap the near-zero nodes to 2π before taking bounds,
    # the SAME correction already applied to query points, now applied to
    # the element's OWN bounds computation. ─────────────────────────────────
    has_near_zero = any(ϕ -> ϕ < 1e-6, ϕs_raw)
    has_large = any(ϕ -> ϕ > π, ϕs_raw)
    ϕs = if has_near_zero && has_large
        [ϕ < 1e-6 ? 2π : ϕ for ϕ in ϕs_raw]
    else
        ϕs_raw
    end

    return (minimum(θs), maximum(θs), minimum(ϕs), maximum(ϕs))
end


function debug_angular_ncf_crossrank_stall(mesh, metrics, neighbor_elem, neighbor_rank,
                                            angular_deps, angular_deps_cr,
                                            solved, n_dirs, remote_round_by_pair, rank)
    println("Rank $rank: STUCK pairs and their FULL per-face breakdown:")
    ngl = mesh.ngl
    n_shown = 0
    for iel = 1:mesh.nelem
        for d = 1:n_dirs[iel]
            solved[iel][d] && continue
            n_shown >= 10 && break

            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            println("  STUCK: (iel=$iel, d=$d), Ω=($Ωx,$Ωy,$Ωz)")
            for fid = 1:6
                (is_inflow[fid] && is_active[fid]) || continue
                jel = neighbor_elem[iel][fid]
                if jel != zero(TInt)
                    key = (iel,d)
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            unsolved = [dn for dn in dep_nodes if !solved[dep_jel][dn]]
                            println("    fid=$fid: SAME-RANK angular_deps -> jel=$jel, unsolved=$unsolved")
                        end
                    else
                        ok = d <= n_dirs[jel] && solved[jel][d]
                        println("    fid=$fid: SAME-RANK FALLBACK -> jel=$jel, d=$d exists&solved? $ok " *
                                "(n_dirs[jel]=$(n_dirs[jel]))")
                    end
                    continue
                end
                # jel == 0
                if haskey(neighbor_rank, (iel,fid))
                    if haskey(angular_deps_cr, (iel,fid,d))
                        r = remote_round_by_pair[(iel,d)]
                        println("    fid=$fid: CROSS-RANK angular_deps_cr -> target_rank=$(neighbor_rank[(iel,fid)]), " *
                                "remote_round=$r " * (r == -1 ? "(BLOCKED)" : "(should be ready)"))
                    else
                        println("    fid=$fid: CROSS-RANK, target_rank=$(neighbor_rank[(iel,fid)]), " *
                                "NO angular_deps_cr entry — falls to UNHANDLED no-op branch " *
                                "(flagged as incomplete in the readiness loop)")
                    end
                else
                    println("    fid=$fid: TRUE BOUNDARY (jel==0, no neighbor_rank entry) — should be ready")
                end
            end
            n_shown += 1
        end
        n_shown >= 10 && break
    end
end

"""
    resolve_stall_distributed_gauss_seidel!(...)

At a GLOBAL stall (confirmed via Allreduce — nothing progresses on ANY
rank), treats the ENTIRE current residual of unsolved (iel,d) pairs
(across BOTH ranks) as one candidate pool, and resolves it via
distributed Gauss-Seidel: same-rank dependencies read directly from local
storage (which may itself be mid-iteration), cross-rank dependencies
exchanged via MPI each sweep, convergence checked via a GLOBAL
Allreduce(max_delta, MAX) — not per-rank, since one rank converging while
its remote data is still changing would be a false convergence.

Scheduling-only version: "solving" a candidate just assigns it the
CURRENT outer iter_round as its round number (matching the placeholder
convention used elsewhere in this test), and "value" is that round
number — real physics would instead solve the local DG system, mirroring
resolve_angular_cycle!'s approach.
"""
function resolve_stall_distributed_gauss_seidel!(mesh, metrics, neighbor_elem, neighbor_rank,
                                                   angular_deps, angular_deps_cr,
                                                   solved, round_solved, n_dirs, iter_round, comm;
                                                   max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)

    # ── Candidate pool: EVERYTHING still unsolved on THIS rank. ──────────────
    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed Gauss-Seidel over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total across ranks)")

    # ── Build the cross-rank exchange list for candidates with an
    # angular_deps_cr entry — same self-describing (gel,pos) tagging as
    # the scheduling exchange, but now carrying an ITERATE VALUE
    # (round-number placeholder here) instead of a one-shot round number. ──
    flat_targets = Int[]; flat_my_gel = TInt[]; flat_my_pos = Int[]
    flat_src_iel = Int[]; flat_src_d = Int[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_my_pos, d)
            push!(flat_src_iel, iel)
            push!(flat_src_d, d)
        end
    end
    n_flat = length(flat_targets)
    val_cache = SendReceiveCache(comm, zeros(Int, n_flat), flat_targets)
    gel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    pos_cache = SendReceiveCache(comm, Float64.(flat_my_pos), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(pos_cache, Float64.(flat_my_pos))

    wishlist = Dict{Tuple{TInt,Int}, Vector{Tuple{Int,Int}}}()
    for (iel,d) in candidates
        haskey(angular_deps_cr, (iel,fid,d)) || continue
        remote_gel, remote_positions = angular_deps_cr[(iel,fid,d)]
        for rp in remote_positions
            push!(get!(wishlist, (remote_gel,rp), Tuple{Int,Int}[]), (iel,d))
        end
    end
    n_recv = length(gel_cache.combined_recv_data)
    recv_match = Vector{Vector{Tuple{Int,Int}}}(undef, n_recv)
    for m = 1:n_recv
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        rpos = Int(round(Int, pos_cache.combined_recv_data[m]))
        recv_match[m] = get(wishlist, (rgel,rpos), Tuple{Int,Int}[])
    end

    # ── Initialize: local candidates get round_solved = 0 (placeholder
    # starting guess); remote-dependency values start unresolved (-1). ──────
    for (iel,d) in candidates
        round_solved[iel][d] = 0
    end
    remote_val_by_pair = Dict{Tuple{Int,Int,Int}, Int}(k => -1 for k in keys(
        Dict((iel,fid,d) => true for (iel,fid,d) in candidates if haskey(angular_deps_cr,(iel,fid,d)))))

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0

        for (iel,d) in candidates
            old_val = round_solved[iel][d]

            # NEW value = 1 + max(dependency values seen so far) — same
            # "propagate a round-like counter" placeholder convention as
            # the rest of this scheduling-only test; real physics would
            # solve the local DG system here instead.
            dep_vals = Int[]
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, mesh.ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)
            for fid = 1:6
                (is_inflow[fid] && is_active[fid]) || continue
                jel = neighbor_elem[iel][fid]
                if jel != zero(TInt)
                    key = (iel,d)
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            for dn in dep_nodes
                                push!(dep_vals, round_solved[dep_jel][dn])
                            end
                        end
                    elseif d <= n_dirs[jel]
                        push!(dep_vals, round_solved[jel][d])
                    end
                elseif haskey(angular_deps_cr, (iel,fid,d))
                    push!(dep_vals, get(remote_val_by_pair, (iel,d), -1))
                end
            end
            # With a genuinely convergent placeholder — bounded averaging, not
	    # unbounded "1 + max" propagation:
	    new_val = isempty(dep_vals) ? 0 : round(Int, sum(max.(dep_vals, 0)) / length(dep_vals))
	    round_solved[iel][d] = new_val
            max_delta_local = max(max_delta_local, abs(new_val - old_val))
        end

        send_vals = [round_solved[flat_src_iel[i]][flat_src_d[i]] for i = 1:n_flat]
        send_and_receive!(val_cache, send_vals)
        for m = 1:length(val_cache.combined_recv_data)
            for (iel,d) in recv_match[m]
                remote_val_by_pair[(iel,d)] = val_cache.combined_recv_data[m]
            end
        end

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global == 0 && break
        gs_iter == max_gs_iter && @warn "distributed Gauss-Seidel did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

"""
    exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind, neighbor_gel_group, comm)

For every rank-cut face, exchanges the PHYSICAL COORDINATES (x,y,z) of
the RELEVANT remote face's nodes — the parent's face-node grid (:i_am_fine)
or the 4 children's face-node grids (:i_am_coarse) — building a
coord_cache::Dict{Int,NTuple{3,Float64}} keyed by GLOBAL node id, reusable
by build_spatial_interpolation_matrices_with_coord_cache exactly as
designed for this purpose in the CG-era code.

Geometry-only, direction-independent — built ONCE, reused across every
direction and every sweep call (mirrors ncf_L1L2's one-time construction
in the single-rank case).
"""
function exchange_ncf_face_coordinates(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                        neighbor_kind::Dict{Tuple{TInt,TInt},Symbol},
                                        neighbor_gel_group::Dict{Tuple{TInt,TInt},Vector{TInt}},
                                        comm)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl

    flat_targets = Int[]
    flat_my_gel  = TInt[]
    flat_pos     = Int[]
    flat_x = Float64[]; flat_y = Float64[]; flat_z = Float64[]; flat_gip = Float64[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        pos = 0
        for k in krange, j in jrange, i in irange
            pos += 1
            ip = mesh.connijk[iel,i,j,k]
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_pos, pos)
            push!(flat_x, mesh.x[ip]); push!(flat_y, mesh.y[ip]); push!(flat_z, mesh.z[ip])
            push!(flat_gip, Float64(mesh.ip2gip[ip]))
        end
    end

    n_flat = length(flat_targets)
    gel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    pos_cache = SendReceiveCache(comm, Float64.(flat_pos), flat_targets)
    x_cache   = SendReceiveCache(comm, flat_x, flat_targets)
    y_cache   = SendReceiveCache(comm, flat_y, flat_targets)
    z_cache   = SendReceiveCache(comm, flat_z, flat_targets)
    gip_cache = SendReceiveCache(comm, flat_gip, flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(pos_cache, Float64.(flat_pos))
    send_and_receive!(x_cache, flat_x)
    send_and_receive!(y_cache, flat_y)
    send_and_receive!(z_cache, flat_z)
    send_and_receive!(gip_cache, flat_gip)

    coord_cache = Dict{Int, NTuple{3,Float64}}()
    gelpos_to_gip = Dict{Tuple{TInt,Int}, Int}()   # NEW — reusable by the value exchange

    needed_gels = Set{TInt}()
    for ((iel,fid), kind) in neighbor_kind
        for g in neighbor_gel_group[(iel,fid)]
            push!(needed_gels, g)
        end
    end

    n_recv = length(gel_cache.combined_recv_data)
    for m = 1:n_recv
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        rgel in needed_gels || continue
        rpos = Int(round(Int, pos_cache.combined_recv_data[m]))
        gip = Int(round(Int, gip_cache.combined_recv_data[m]))
        coord_cache[gip] = (x_cache.combined_recv_data[m], y_cache.combined_recv_data[m], z_cache.combined_recv_data[m])
        gelpos_to_gip[(rgel,rpos)] = gip   # NEW
    end

    return coord_cache, gelpos_to_gip   # CHANGED — now returns both
end


"""
    setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model; spatial_target_iel)

Shared setup for every spatial-NCF-across-MPI test: refines
`spatial_target_iel` (rank 0 only), builds same-rank adjacency and
spatial-NCF tables, identifies and resolves rank-cut faces via the
NCF-aware resolver, EXCLUDING same-rank NCF faces from the suspected-cut
list (the fix established earlier in this effort).

Returns everything downstream tests need: mesh, metrics, neighbor_elem,
neighbor_fid, ncf_as_child, ncf_as_parent, neighbor_rank, neighbor_kind,
neighbor_gel_group.
"""
function setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                     spatial_target_iel::Int = 303)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
    end
    println("Rank $rank: spatial NCF: $(mesh.num_ncf) same-rank entries")

    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved =
        resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        error("Rank $rank: unresolved rank-cut faces — aborting")
    end

    return (mesh=mesh, metrics=metrics, sem=sem, partitioned_model=partitioned_model,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid,
            ncf_as_child=ncf_as_child, ncf_as_parent=ncf_as_parent,
            neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
            neighbor_gel_group=neighbor_gel_group)
end

function test_exchange_ncf_face_coordinates_full(inputs, nparts, distribute, sem, partitioned_model;
                                                   spatial_target_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind, neighbor_gel_group, comm)
    println("Rank $rank: coord_cache has $(length(coord_cache)) entries")

    all_needed_gels = Set{TInt}()
    for ((iel,fid), kind) in neighbor_kind
        for g in neighbor_gel_group[(iel,fid)]
            push!(all_needed_gels, g)
        end
    end

    flat_targets = Int[]; flat_gel = TInt[]; flat_gip = Float64[]
    flat_x = Float64[]; flat_y = Float64[]; flat_z = Float64[]
    ngl = mesh.ngl
    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for k in krange, j in jrange, i in irange
            ip = mesh.connijk[iel,i,j,k]
            push!(flat_targets, Int(target_rank)); push!(flat_gel, my_gel)
            push!(flat_gip, Float64(mesh.ip2gip[ip]))
            push!(flat_x, mesh.x[ip]); push!(flat_y, mesh.y[ip]); push!(flat_z, mesh.z[ip])
        end
    end
    gel_cache2 = SendReceiveCache(comm, Float64.(flat_gel), flat_targets)
    gip_cache2 = SendReceiveCache(comm, flat_gip, flat_targets)
    x_cache2 = SendReceiveCache(comm, flat_x, flat_targets)
    y_cache2 = SendReceiveCache(comm, flat_y, flat_targets)
    z_cache2 = SendReceiveCache(comm, flat_z, flat_targets)
    send_and_receive!(gel_cache2, Float64.(flat_gel))
    send_and_receive!(gip_cache2, flat_gip)
    send_and_receive!(x_cache2, flat_x)
    send_and_receive!(y_cache2, flat_y)
    send_and_receive!(z_cache2, flat_z)

    n_checked = 0
    n_mismatch = 0
    for m = 1:length(gel_cache2.combined_recv_data)
        rgel = TInt(round(Int, gel_cache2.combined_recv_data[m]))
        rgel in all_needed_gels || continue
        gip = Int(round(Int, gip_cache2.combined_recv_data[m]))
        true_xyz = (x_cache2.combined_recv_data[m], y_cache2.combined_recv_data[m], z_cache2.combined_recv_data[m])
        n_checked += 1
        if !haskey(coord_cache, gip)
            n_mismatch += 1
            println("Rank $rank: ✗ gip=$gip MISSING from coord_cache")
            continue
        end
        recon_xyz = coord_cache[gip]
        if !isapprox(collect(recon_xyz), collect(true_xyz); atol=1e-8)
            n_mismatch += 1
            println("Rank $rank: ✗ gip=$gip MISMATCH: recon=$recon_xyz true=$true_xyz")
        end
    end
    println("Rank $rank: checked $n_checked coordinate entries, $n_mismatch mismatches  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))
    return coord_cache
end

"""
    setup_angular_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model; angular_target_iel, angular_target_e_ext)

Shared setup for every angular-NCF-across-MPI test: refines
`angular_target_iel` angularly (with its φ-seam partner, conforming),
identifies and resolves ORDINARY (spatially conforming) rank-cut faces,
exchanges each face's full remote angular mesh, and builds the cross-rank
angular dependency graph — mirroring setup_spatial_ncf_mpi_test's role for
the spatial case, but for the angular mechanism instead.

Per the one-buffer-element policy already established: assumes every
rank-cut face touched here is spatially conforming (angular refinement
only, no spatial NCF at this same interface).
"""
function setup_angular_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                     angular_target_iel::Int = 442, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    metrics = sem.metrics
    if angular_target_iel <= mesh.nelem
        adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
        verify_angular_phi_seam_conformity(mesh, angular_target_iel)
    end
    println("Rank $rank: mesh has $(mesh.nelem) elements")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        error("Rank $rank: unresolved rank-cut faces — aborting")
    end

    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    println("Rank $rank: $(length(angular_deps_cr)) cross-rank angular NCF (iel,d) entries")

    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    return (mesh=mesh, metrics=metrics, sem=sem, partitioned_model=partitioned_model,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid,
            neighbor_rank=neighbor_rank, neighbor_gel=neighbor_gel,
            remote_angular_mesh=remote_angular_mesh,
            angular_deps=angular_deps, angular_deps_cr=angular_deps_cr)
end

"""
    build_ncf_value_exchange(mesh, neighbor_rank, comm)

Reusable, geometry-only exchange structure for per-direction DG VALUE
exchange — same (gel,position) tagging as exchange_ncf_face_coordinates,
built ONCE and reused across every direction/round.
"""
function build_ncf_value_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]
    flat_my_gel  = TInt[]
    flat_pos     = Int[]
    flat_local_ip = Int[]   # MY OWN local ip at each flat position — for reading x_dg when SENDING

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        pos = 0
        for k in krange, j in jrange, i in irange
            pos += 1
            ip = mesh.connijk[iel,i,j,k]
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_pos, pos)
            push!(flat_local_ip, ip)
        end
    end

    n_flat = length(flat_targets)
    gel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    pos_cache = SendReceiveCache(comm, Float64.(flat_pos), flat_targets)
    val_cache = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(pos_cache, Float64.(flat_pos))

    return (val_cache=val_cache, gel_cache=gel_cache, pos_cache=pos_cache,
            flat_local_ip=flat_local_ip, n_flat=n_flat)
end

"""
    exchange_ncf_values!(exch, gelpos_to_gip, send_values) -> Dict{Int,Float64}

Sends `send_values` (length n_flat, one per MY OWN face node, SAME order
as exch.flat_local_ip) and returns a Dict{Int,Float64} keyed by REMOTE
global ip (same key space as coord_cache) — using the ALREADY-VALIDATED
gelpos_to_gip mapping from exchange_ncf_face_coordinates to correctly
identify which remote gip each incoming (gel,pos)-tagged value belongs to.
"""
function exchange_ncf_values!(exch, gelpos_to_gip::Dict{Tuple{TInt,Int},Int}, send_values::Vector{Float64})
    send_and_receive!(exch.val_cache, send_values)
    value_by_gip = Dict{Int,Float64}()
    n_recv = length(exch.val_cache.combined_recv_data)
    for m = 1:n_recv
        rgel = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rpos = Int(round(Int, exch.pos_cache.combined_recv_data[m]))
        gip = get(gelpos_to_gip, (rgel,rpos), 0)
        gip == 0 && continue   # not a gip we actually need — discard
        value_by_gip[gip] = exch.val_cache.combined_recv_data[m]
    end
    return value_by_gip
end

function test_exchange_ncf_values_full(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group


    test_ipc_ordering_vs_face_ip_range(mesh)

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    exch = build_ncf_value_exchange(mesh, neighbor_rank, comm)

    # Fabricated payload: each of MY OWN face nodes sends its OWN gip as
    # the "value" — a correct round-trip means value_by_gip[gip] == gip.
    send_values = Float64[Float64(mesh.ip2gip[ip]) for ip in exch.flat_local_ip]
    value_by_gip = exchange_ncf_values!(exch, gelpos_to_gip, send_values)

    n_checked = 0
    n_mismatch = 0
    for (gip, val) in value_by_gip
        n_checked += 1
        if !isapprox(val, Float64(gip); atol=1e-8)
            n_mismatch += 1
            println("Rank $rank: ✗ gip=$gip got value=$val (expected $gip)")
        end
    end
    println("Rank $rank: checked $n_checked value entries, $n_mismatch mismatches  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))
    println("Rank $rank: value_by_gip has $(length(value_by_gip)) entries, " *
            "coord_cache has $(length(coord_cache)) entries " *
            (length(value_by_gip) == length(coord_cache) ? "(match ✓)" : "(MISMATCH ✗)"))

    return coord_cache, gelpos_to_gip, exch, value_by_gip
end

function test_ipc_ordering_vs_face_ip_range(mesh)
    ngl = mesh.ngl
    ncf_idx = 1
    child_iel = Int(mesh.cip[ncf_idx])
    child_fid = Int(mesh.lfid[ncf_idx])

    ipc_list_order = Int.(mesh.IPc_list[:, ncf_idx])

    irange, jrange, krange = face_ip_range(child_fid, ngl)
    face_ip_range_order = Int[]
    for k in krange, j in jrange, i in irange
        push!(face_ip_range_order, mesh.connijk[child_iel, i, j, k])
    end

    println("NCF #$ncf_idx: child_iel=$child_iel, child_fid=$child_fid")
    println("mesh.IPc_list order:       $ipc_list_order")
    println("face_ip_range order:       $face_ip_range_order")
    println("Identical order? $(ipc_list_order == face_ip_range_order)")
    println("Same SET of nodes (order-independent)? $(Set(ipc_list_order) == Set(face_ip_range_order))")

    if ipc_list_order != face_ip_range_order && Set(ipc_list_order) == Set(face_ip_range_order)
        # Show the permutation mapping, if one exists cleanly
        perm = [findfirst(==(n), face_ip_range_order) for n in ipc_list_order]
        println("Permutation (IPc_list position -> face_ip_range position): $perm")
    end

    return ipc_list_order, face_ip_range_order
end

"""
    build_dg_ncf_parent_to_child_remote(mesh, coord_cache, value_by_gip,
                                         child_iel, remote_parent_gel, fid, d,
                                         local_gel_for_parent_gips, ngl)

Cross-rank analogue of build_dg_ncf_parent_to_child_local: the PARENT
lives on a REMOTE rank. `coord_cache`/`value_by_gip` (keyed by GLOBAL ip)
supply the remote parent's face-node coordinates and solved values;
`face_ip_range`'s ordering is used directly (CONFIRMED identical to
mesh.IPc_list/IPp_list's native ordering).

Since the REMOTE parent's own local (i,j,k)/IPp-style node list isn't
directly available (only its face-node GLOBAL ids, via the exchange), we
build L1/L2 using build_spatial_interpolation_matrices_with_coord_cache,
which was DESIGNED for exactly this — parent coordinates via a cache
keyed by global ip, child coordinates from local mesh data.
"""
function build_dg_ncf_parent_to_child_remote(mesh, coord_cache::Dict{Int,NTuple{3,Float64}},
                                              value_by_gip::Dict{Int,Float64},
                                              child_iel::Int, fid::Int,
                                              parent_global_gips::Vector{Int}, ngl::Int)
    IPc_nodes = Int[]
    ngl2 = ngl*ngl
    irange, jrange, krange = face_ip_range(fid, ngl)
    for k in krange, j in jrange, i in irange
        push!(IPc_nodes, mesh.connijk[child_iel,i,j,k])
    end

    L1, L2 = build_spatial_interpolation_matrices_with_coord_cache(
        fid, IPc_nodes, parent_global_gips, mesh, coord_cache, ngl)

    # Gather the remote parent's face-node VALUES, in the SAME order as
    # parent_global_gips (which itself must match face_ip_range's order,
    # per the confirmed identity above).
    parent_grid = zeros(ngl, ngl)
    for idx = 1:length(parent_global_gips)
        gip = parent_global_gips[idx]
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        parent_grid[p1,p2] = get(value_by_gip, gip, NaN)
    end

    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        ci = div(node_idx-1, ngl) + 1
        cj = mod(node_idx-1, ngl) + 1
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += w_1[p1]*w_2[p2]*parent_grid[p1,p2]
        end
        child_vals[node_idx] = val
    end
    return child_vals
end

function test_dg_ncf_parent_to_child_remote(inputs, nparts, distribute, sem, partitioned_model;
                                             spatial_target_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group
    ngl = mesh.ngl

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    exch = build_ncf_value_exchange(mesh, neighbor_rank, comm)

    # ── Well-scaled fabricated field — same discipline as every prior
    # interpolation test. Every rank computes it at its OWN face nodes
    # (ground truth), then we exchange and check the RECEIVED values. ───────
    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    send_values = Float64[field(mesh.x[ip], mesh.y[ip], mesh.z[ip]) for ip in exch.flat_local_ip]
    value_by_gip = exchange_ncf_values!(exch, gelpos_to_gip, send_values)

    # ── Find a face where I AM THE CHILD (:i_am_fine) — this rank needs
    # to have at least one such face for the test to run. ───────────────────
    test_iel = 0; test_fid = 0; remote_gel = zero(TInt)
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_fine || continue
        test_iel = iel; test_fid = fid
        remote_gel = neighbor_gel_group[(iel,fid)][1]
        break
    end
    if test_iel == 0
        println("Rank $rank: no :i_am_fine face found — nothing to test on this rank")
        return nothing
    end
    println("Rank $rank: testing (iel=$test_iel, fid=$test_fid), remote parent gel=$remote_gel")

    ngl2 = ngl*ngl
    parent_global_gips = [gelpos_to_gip[(remote_gel, pos)] for pos = 1:ngl2]

    child_vals = build_dg_ncf_parent_to_child_remote(mesh, coord_cache, value_by_gip,
                                                       test_iel, test_fid, parent_global_gips, ngl)

    # ── Compare against EXACT field at the child's own face nodes ───────────
    irange, jrange, krange = face_ip_range(test_fid, ngl)
    IPc_nodes = Int[]
    for k in krange, j in jrange, i in irange
        push!(IPc_nodes, mesh.connijk[test_iel,i,j,k])
    end

    max_err = 0.0
    for idx = 1:length(IPc_nodes)
        node = IPc_nodes[idx]
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        max_err = max(max_err, abs(child_vals[idx] - exact))
    end
    println("Rank $rank: parent->child REMOTE interpolation max abs error = $max_err")
    return max_err
end

"""
    build_dg_ncf_child_to_parent_remote(mesh, coord_cache, value_by_gip,
                                         parent_iel, fid, remote_children_gels, ngl)

Cross-rank analogue of build_dg_ncf_child_to_parent_v5_local: the (up to)
4 CHILDREN live on a REMOTE rank. For each of the parent's own face
nodes, finds which ONE remote child's face-node coordinate range (via
coord_cache) covers that physical point, then evaluates via direct 2D
Lagrange interpolation from that child's remote face-node VALUES
(value_by_gip) — same "one node, one covering child, direct evaluation"
principle as the local version, no projection.
"""
function build_dg_ncf_child_to_parent_remote(mesh, coord_cache::Dict{Int,NTuple{3,Float64}},
                                              value_by_gip::Dict{Int,Float64},
                                              parent_iel::Int, fid::Int,
                                              remote_children_gels::Vector{TInt},
                                              gelpos_to_gip::Dict{Tuple{TInt,Int}, Int}, ngl::Int)
    ngl2 = ngl*ngl
    irange, jrange, krange = face_ip_range(fid, ngl)
    IPp_nodes = Int[]
    for k in krange, j in jrange, i in irange
        push!(IPp_nodes, mesh.connijk[parent_iel,i,j,k])
    end

    parent_result = fill(NaN, ngl, ngl)

    for remote_gel in remote_children_gels
        child_global_gips = [gelpos_to_gip[(remote_gel,pos)] for pos = 1:ngl2]

        # Empirically determine which physical axes vary (SAME technique
        # as the local version's axis_of_max_spread — no assumption about
        # x/y/z from fid at all).
        child_coords = [coord_cache[gip] for gip in child_global_gips]
        xs = [c[1] for c in child_coords]; ys = [c[2] for c in child_coords]; zs = [c[3] for c in child_coords]

        child_line_1_idx = [(i-1)*ngl+1 for i = 1:ngl]   # SAME face_ip_range-confirmed convention
        child_line_2_idx = [j for j = 1:ngl]

        xs1=[xs[p] for p in child_line_2_idx]; ys1=[ys[p] for p in child_line_2_idx]; zs1=[zs[p] for p in child_line_2_idx]
        axis_c1 = axis_of_max_spread((xs1,ys1,zs1))
        xs2=[xs[p] for p in child_line_1_idx]; ys2=[ys[p] for p in child_line_1_idx]; zs2=[zs[p] for p in child_line_1_idx]
        axis_c2 = axis_of_max_spread((xs2,ys2,zs2))
        coord_of(p, axis) = axis==1 ? xs[p] : (axis==2 ? ys[p] : zs[p])

        c1_coords = [coord_of(p,axis_c1) for p in child_line_2_idx]
        c2_coords = [coord_of(p,axis_c2) for p in child_line_1_idx]
        c1_min,c1_max = extrema(c1_coords); c2_min,c2_max = extrema(c2_coords)
        tol1 = 1e-6*max(c1_max-c1_min,1.0); tol2 = 1e-6*max(c2_max-c2_min,1.0)

        child_grid = zeros(ngl,ngl)
        for idx = 1:ngl2
            c1 = ((idx-1)%ngl)+1; c2 = ((idx-1)÷ngl)+1
            gip = child_global_gips[idx]
            child_grid[c1,c2] = get(value_by_gip, gip, NaN)
        end

        for idx = 1:length(IPp_nodes)
            node = IPp_nodes[idx]
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            pc1 = coord_of_local(mesh, node, axis_c1)
            pc2 = coord_of_local(mesh, node, axis_c2)
            (pc1 < c1_min-tol1 || pc1 > c1_max+tol1 || pc2 < c2_min-tol2 || pc2 > c2_max+tol2) && continue

            L_c1 = build_1d_lagrange_matrix(c1_coords, [pc1])
            L_c2 = build_1d_lagrange_matrix(c2_coords, [pc2])
            val = 0.0
            for a=1:ngl, b=1:ngl
                val += L_c1[1,a]*L_c2[1,b]*child_grid[a,b]
            end
            parent_result[p1,p2] = val
        end
    end
    return parent_result
end

@inline coord_of_local(mesh, node, axis) = axis==1 ? mesh.x[node] : (axis==2 ? mesh.y[node] : mesh.z[node])

function test_dg_ncf_child_to_parent_remote(inputs, nparts, distribute, sem, partitioned_model;
                                             spatial_target_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group
    ngl = mesh.ngl

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    exch = build_ncf_value_exchange(mesh, neighbor_rank, comm)

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    send_values = Float64[field(mesh.x[ip], mesh.y[ip], mesh.z[ip]) for ip in exch.flat_local_ip]
    value_by_gip = exchange_ncf_values!(exch, gelpos_to_gip, send_values)

    # ── Find a face where I AM THE COARSE side ───────────────────────────────
    test_iel = 0; test_fid = 0; remote_children_gels = TInt[]
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_coarse || continue
        test_iel = iel; test_fid = fid
        remote_children_gels = neighbor_gel_group[(iel,fid)]
        break
    end
    if test_iel == 0
        println("Rank $rank: no :i_am_coarse face found — nothing to test on this rank")
        return nothing
    end
    println("Rank $rank: testing (iel=$test_iel, fid=$test_fid), remote children gels=$remote_children_gels")

    parent_computed = build_dg_ncf_child_to_parent_remote(mesh, coord_cache, value_by_gip,
                                                           test_iel, test_fid, remote_children_gels,
                                                           gelpos_to_gip, ngl)

    irange, jrange, krange = face_ip_range(test_fid, ngl)
    IPp_nodes = Int[]
    for k in krange, j in jrange, i in irange
        push!(IPp_nodes, mesh.connijk[test_iel,i,j,k])
    end

    n_nan = count(isnan, parent_computed)
    n_nan > 0 && println("Rank $rank: ⚠ $n_nan parent face nodes never claimed by any remote child")

    max_err = 0.0
    for idx = 1:length(IPp_nodes)
        node = IPp_nodes[idx]
        p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        max_err = max(max_err, abs(parent_computed[p1,p2] - exact))
    end
    println("Rank $rank: child->parent REMOTE interpolation max abs error = $max_err")
    return max_err
end

function build_ncf_value_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, connijk_dg, comm)
    ngl = mesh.ngl
    flat_targets = Int[]; flat_my_gel = TInt[]; flat_pos = Int[]
    flat_owning_iel = Int[]; flat_dg_ip = Int[]   # NEW

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        pos = 0
        for k in krange, j in jrange, i in irange
            pos += 1
            push!(flat_targets, Int(target_rank))
            push!(flat_my_gel, my_gel)
            push!(flat_pos, pos)
            push!(flat_owning_iel, iel)
            push!(flat_dg_ip, connijk_dg[iel,i,j,k])
        end
    end

    n_flat = length(flat_targets)
    gel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    pos_cache = SendReceiveCache(comm, Float64.(flat_pos), flat_targets)
    val_cache = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(pos_cache, Float64.(flat_pos))

    return (val_cache=val_cache, gel_cache=gel_cache, pos_cache=pos_cache,
            flat_owning_iel=flat_owning_iel, flat_dg_ip=flat_dg_ip, n_flat=n_flat)
end

function run_dg_sweep_multirank_spatial_ncf!(x_dg::Vector{Float64}, mesh, metrics, ψ, dψ, ω,
                                              Ωx::Float64, Ωy::Float64, Ωz::Float64,
                                              neighbor_elem, neighbor_rank, neighbor_gel, exchange, value_cache,
                                              ncf_as_child, ncf_L1L2, ncf_as_parent,
                                              neighbor_kind, neighbor_gel_group,
                                              ncf_coord_cache, ncf_gelpos_to_gip, ncf_exch,
                                              connijk_dg, face_neighbor_dg,
                                              κ, σ, intΦ_eff, bc_fn, comm::MPI.Comm;
                                              extra_source = nothing, max_rounds::Int = 200)
    rank  = MPI.Comm_rank(comm)
    ngl   = mesh.ngl
    nelem = mesh.nelem
    ngl2  = ngl*ngl

    is_inflow_all = Vector{Vector{Bool}}(undef, nelem)
    is_active_all = Vector{Vector{Bool}}(undef, nelem)
    for iel = 1:nelem
        ii = Vector{Bool}(undef,6); ia = Vector{Bool}(undef,6)
        classify_element_faces!(ii, ia, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        is_inflow_all[iel] = ii; is_active_all[iel] = ia
    end

    processed = falses(nelem)

    remote_valid = Dict{Tuple{Int,Int}, Vector{Bool}}()
    remote_value = Dict{Tuple{Int,Int}, Vector{Float64}}()
    remote_gel   = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ((iel,fid),_) in neighbor_rank
        n = length(exchange.face_local_order[(TInt(iel),TInt(fid))])
        remote_valid[(iel,fid)] = falses(n)
        remote_value[(iel,fid)] = zeros(Float64, n)
        remote_gel[(iel,fid)]   = fill(-1, n)
    end

    ncf_remote_valid = Dict{Int,Bool}()
    ncf_remote_value = Dict{Int,Float64}()
    for ((iel,fid), kind) in neighbor_kind
        for g in neighbor_gel_group[(iel,fid)]
            for pos = 1:ngl2
                gip = ncf_gelpos_to_gip[(g,pos)]
                ncf_remote_valid[gip] = false
                ncf_remote_value[gip] = NaN
            end
        end
    end

    A_vol = zeros(ngl^3, ngl^3)
    RHS_local = zeros(ngl^3)

    iter_round = 0
    while true
        iter_round += 1
        iter_round > max_rounds && error("run_dg_sweep_multirank_spatial_ncf!: exceeded max_rounds=$max_rounds")

        for iel = 1:nelem
            processed[iel] && continue
            ii = is_inflow_all[iel]; ia = is_active_all[iel]

            # ── Readiness — CORRECTED structure: same-rank NCF and
            # cross-rank (conforming or NCF) BOTH live under jel==0,
            # matching the serial run_dg_sweep!'s structure exactly. ────────
            ready = true
            for fid = 1:6
                (ia[fid] && ii[fid]) || continue
                jel = neighbor_elem[iel][fid]

                if jel == zero(TInt)
                    if haskey(ncf_as_child, (iel,fid))
                        ncf_idx = ncf_as_child[(iel,fid)]
                        parent_iel = Int(mesh.pip[ncf_idx])
                        processed[parent_iel] || (ready = false; break)
                    elseif haskey(ncf_as_parent, (iel,fid))
                        all_children_done = true
                        for ncf_idx in ncf_as_parent[(iel,fid)]
                            c_iel = Int(mesh.cip[ncf_idx])
                            processed[c_iel] || (all_children_done = false; break)
                        end
                        all_children_done || (ready = false; break)
                    elseif haskey(neighbor_kind, (iel,fid))
                        kind = neighbor_kind[(iel,fid)]
                        if kind == :conforming
                            all(remote_valid[(iel,fid)]) || (ready = false; break)
                        elseif kind == :i_am_fine
                            remote_gel_1 = neighbor_gel_group[(iel,fid)][1]
                            ok = true
                            for pos = 1:ngl2
                                gip = ncf_gelpos_to_gip[(remote_gel_1,pos)]
                                ncf_remote_valid[gip] || (ok = false; break)
                            end
                            ok || (ready = false; break)
                        elseif kind == :i_am_coarse
                            ok = true
                            for g in neighbor_gel_group[(iel,fid)]
                                for pos = 1:ngl2
                                    gip = ncf_gelpos_to_gip[(g,pos)]
                                    ncf_remote_valid[gip] || (ok = false; break)
                                end
                                ok || break
                            end
                            ok || (ready = false; break)
                        end
                    end
                    # else: true physical boundary — always ready
                else
                    processed[jel] || (ready = false; break)
                end
            end
            ready || continue

            fill!(A_vol, 0.0); fill!(RHS_local, 0.0)
            build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            extra_source !== nothing &&
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source, mesh.connijk)

            for fid = 1:6
                ia[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)

                I_upwind = (i,j,k) -> 0.0
                if ii[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            L1,L2 = ncf_L1L2[ncf_idx]
                            IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                            child_vals = build_dg_ncf_parent_to_child(mesh, connijk_dg, x_dg,
                                                                        iel, parent_iel, fid,
                                                                        IPc_nodes, IPp_nodes, L1, L2, ngl)
                            I_upwind = (i,j,k) -> child_vals[face_flat_index(fid, ngl, i, j, k)]
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                       for idx in ncf_group]
                            parent_grid = build_dg_ncf_child_to_parent_v5(mesh, connijk_dg, x_dg, iel, fid, entries, ngl)
                            I_upwind = (i,j,k) -> begin
                                node_idx = face_flat_index(fid, ngl, i, j, k)
                                p1 = ((node_idx-1) % ngl) + 1
                                p2 = ((node_idx-1) ÷ ngl) + 1
                                parent_grid[p1,p2]
                            end
                        elseif haskey(neighbor_kind, (iel,fid))
                            kind = neighbor_kind[(iel,fid)]
                            if kind == :conforming
                                rvals = remote_value[(iel,fid)]
                                canonical_order = exchange.face_local_order[(TInt(iel),TInt(fid))]
                                ijk_to_p = Dict{Tuple{Int,Int,Int},Int}(canonical_order[pp] => pp for pp = 1:length(canonical_order))
                                I_upwind = (i,j,k) -> rvals[ijk_to_p[(i,j,k)]]
                            elseif kind == :i_am_fine
                                remote_gel_1 = neighbor_gel_group[(iel,fid)][1]
                                parent_global_gips = [ncf_gelpos_to_gip[(remote_gel_1,pos)] for pos=1:ngl2]
                                child_vals = build_dg_ncf_parent_to_child_remote(mesh, ncf_coord_cache,
                                    ncf_remote_value, iel, fid, parent_global_gips, ngl)
                                I_upwind = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                            elseif kind == :i_am_coarse
                                remote_gels = neighbor_gel_group[(iel,fid)]
                                parent_grid = build_dg_ncf_child_to_parent_remote(mesh, ncf_coord_cache,
                                    ncf_remote_value, iel, fid, remote_gels, ncf_gelpos_to_gip, ngl)
                                I_upwind = (i,j,k) -> begin
                                    node_idx = face_flat_index(fid,ngl,i,j,k)
                                    p1 = ((node_idx-1)%ngl)+1; p2 = ((node_idx-1)÷ngl)+1
                                    parent_grid[p1,p2]
                                end
                            end
                        else
                            I_upwind = (i,j,k) -> bc_fn(mesh.x[mesh.connijk[iel,i,j,k]],
                                                         mesh.y[mesh.connijk[iel,i,j,k]],
                                                         mesh.z[mesh.connijk[iel,i,j,k]])
                        end
                    else
                        I_upwind = (i,j,k) -> x_dg[face_neighbor_dg[face_flat_index(fid, ngl, i, j, k), fid, iel]]
                    end
                end

                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx, Ωy, Ωz,
                                      ii[fid], ia[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_local = A_vol \ RHS_local
            for k=1:ngl, j=1:ngl, i=1:ngl
                x_dg[connijk_dg[iel,i,j,k]] = x_local[lidx(i,j,k,ngl)]
            end
            processed[iel] = true
        end

        send_vals     = Vector{Float64}(undef, length(exchange.send_keys))
        send_gel_vals = Vector{Int}(undef, length(exchange.send_keys))
        for (idx,(iel,fid,pos)) in enumerate(exchange.send_keys)
            i,j,k = exchange.face_local_order[(TInt(iel),TInt(fid))][pos]
            send_vals[idx] = processed[iel] ? x_dg[connijk_dg[iel,i,j,k]] : NaN
            send_gel_vals[idx] = Int(mesh.el2gel[iel])
        end
        send_and_receive!(value_cache, send_vals)
        send_and_receive!(exchange.cache_gel, send_gel_vals)
        for m = 1:length(value_cache.combined_recv_data)
            val = value_cache.combined_recv_data[m]
            isnan(val) && continue
            recv_gel = exchange.cache_gel.combined_recv_data[m]
            for pos in exchange.recv_perm[m]
                iel, fid, p = exchange.send_keys[pos]
                remote_gel[(iel,fid)][p]   = recv_gel
                remote_value[(iel,fid)][p] = val
                remote_valid[(iel,fid)][p] = true
            end
        end

        ncf_send_vals = [processed[ncf_exch.flat_owning_iel[k]] ? x_dg[ncf_exch.flat_dg_ip[k]] : NaN
                          for k = 1:ncf_exch.n_flat]
        send_and_receive!(ncf_exch.val_cache, ncf_send_vals)
        for m = 1:length(ncf_exch.val_cache.combined_recv_data)
            val = ncf_exch.val_cache.combined_recv_data[m]
            isnan(val) && continue
            rgel = TInt(round(Int, ncf_exch.gel_cache.combined_recv_data[m]))
            rpos = Int(round(Int, ncf_exch.pos_cache.combined_recv_data[m]))
            gip = get(ncf_gelpos_to_gip, (rgel,rpos), 0)
            gip == 0 && continue
            ncf_remote_value[gip] = val
            ncf_remote_valid[gip] = true
        end

        n_unprocessed_total = MPI.Allreduce(count(!,processed), +, comm)
        n_unprocessed_total == 0 && break
    end

    return x_dg
end

function test_dg_sweep_multirank_spatial_ncf_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                          spatial_target_iel::Int = 215,
                                                          Ωx::Float64 = 0.3, Ωy::Float64 = 0.5, Ωz::Float64 = 0.8,
                                                          κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh    = setup.mesh
    metrics = setup.metrics
    sem2    = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl

    neighbor_elem = setup.neighbor_elem
    neighbor_fid  = setup.neighbor_fid
    ncf_as_child  = setup.ncf_as_child
    ncf_as_parent = setup.ncf_as_parent
    neighbor_rank_ncf = setup.neighbor_rank      # from resolve_rank_cut_neighbors_ncf
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group

    println("Rank $rank: building DG storage + L1/L2 for same-rank NCF...")
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
    end

    # ── Build a SEPARATE neighbor_rank for the ORDINARY conforming
    # rank-cut faces only, needed by the ALREADY-VALIDATED conforming
    # exchange machinery (build_rank_cut_exchange expects a plain
    # (iel,fid)->rank dict, one entry per conforming face). ─────────────────
    

    # ── Build BOTH inputs build_rank_cut_exchange actually requires ──────────
    conforming_rank = Dict{Tuple{TInt,TInt}, TInt}()
    conforming_gel  = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
    	conforming_rank[(TInt(iel),TInt(fid))] = neighbor_rank_ncf[(iel,fid)]
        conforming_gel[(TInt(iel),TInt(fid))]  = neighbor_gel_group[(iel,fid)][1]
    end
    println("Rank $rank: $(length(conforming_rank)) ordinary conforming rank-cut faces")

    exchange = build_rank_cut_exchange(mesh, conforming_rank, conforming_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        error("Rank $rank: conforming exchange setup failed")
    end
    value_cache = build_rank_cut_value_cache(exchange, comm)
    
    println("Rank $rank: building NCF coordinate/value exchange...")
    ncf_coord_cache, ncf_gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank_ncf, neighbor_kind,
                                                                         neighbor_gel_group, comm)
    ncf_exch = build_ncf_value_exchange(mesh, neighbor_rank_ncf, connijk_dg, comm)

    # ── Exact exponential solution, same as every prior exponential test ────
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

    x_dg = fill(NaN, npoin_dg)
    println("Rank $rank: running combined multirank spatial-NCF sweep...")
    run_dg_sweep_multirank_spatial_ncf!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz,
                                     neighbor_elem, conforming_rank, conforming_gel, exchange, value_cache,
                                     ncf_as_child, ncf_L1L2, ncf_as_parent,
                                     neighbor_kind, neighbor_gel_group,
                                     ncf_coord_cache, ncf_gelpos_to_gip, ncf_exch,
                                     connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn, comm)

    n_nan = count(isnan, x_dg)
    n_nan_total = MPI.Allreduce(n_nan, +, comm)
    if n_nan_total > 0
        println("Rank $rank: ⚠ $n_nan unwritten DG entries locally ($n_nan_total total)")
    end

    max_err_local = 0.0
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        gip = mesh.connijk[iel,i,j,k]
        exact = bc_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        v = x_dg[connijk_dg[iel,i,j,k]]
        max_err_local = max(max_err_local, abs(v - exact))
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)

    if rank == 0
        println("Multirank spatial-NCF exponential test — max abs error (all ranks): $max_err_global")
    end
    return max_err_global
end

function test_trace_element_415(inputs, nparts, distribute, sem, partitioned_model;
                                 spatial_target_iel::Int = 415,
                                 Ωx::Float64 = 1.0, Ωy::Float64 = 0.0, Ωz::Float64 = 0.0,
                                 κ_val::Float64 = 0.0001, I0::Float64 = 4.2,
                                 trace_iel::Int = 415)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    # ══════════════════════════════════════════════════════════════════════
    # PATH A: via setup_spatial_ncf_mpi_test (used by the multirank sweep)
    # ══════════════════════════════════════════════════════════════════════
    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh    = setup.mesh
    metrics = setup.metrics
    sem2    = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl

    neighbor_elem_A = setup.neighbor_elem
    neighbor_fid_A  = setup.neighbor_fid
    ncf_as_child_A  = setup.ncf_as_child
    ncf_as_parent_A = setup.ncf_as_parent
    neighbor_rank_ncf = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem_A, neighbor_fid_A)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem_A, neighbor_fid_A, face_correspondence)

    ncf_L1L2_A = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2_A[ncf_idx] = (L1, L2)
    end

    # ══════════════════════════════════════════════════════════════════════
    # PATH B: EXACT replica of the known-good standalone test's construction
    # (test_dg_sweep_spatial_ncf_all_directions), built INDEPENDENTLY from
    # the SAME mesh — same loop, same order, same everything, just a
    # SEPARATE set of Dicts, to check whether the RESULT differs at all.
    # ══════════════════════════════════════════════════════════════════════
    neighbor_elem_B, neighbor_fid_B = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child_B = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2_B     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent_B = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children_B = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child_B[(child_iel, child_fid)] = ncf_idx
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2_B[ncf_idx] = (L1, L2)
        push!(get!(ncf_as_parent_B, (parent_iel,parent_fid), Int[]), ncf_idx)
        push!(get!(ncf_parent_to_children_B, (parent_iel,parent_fid), Int[]), child_iel)
    end
    ncf_child_to_parent_elemid_B = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child_B)

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON 1: do neighbor_elem_A and neighbor_elem_B agree?
    # ══════════════════════════════════════════════════════════════════════
    println("Rank $rank: === COMPARISON 1: neighbor_elem A vs B ===")
    n_diff_ne = 0
    for iel = 1:mesh.nelem
        if neighbor_elem_A[iel] != neighbor_elem_B[iel]
            n_diff_ne += 1
            n_diff_ne <= 5 && println("  iel=$iel: A=$(neighbor_elem_A[iel]) B=$(neighbor_elem_B[iel])")
        end
    end
    println("  Total differing elements: $n_diff_ne / $(mesh.nelem)")

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON 2: do ncf_as_child_A and ncf_as_child_B agree (same keys,
    # same ncf_idx values)?
    # ══════════════════════════════════════════════════════════════════════
    println("Rank $rank: === COMPARISON 2: ncf_as_child A vs B ===")
    keys_A = Set(keys(ncf_as_child_A)); keys_B = Set(keys(ncf_as_child_B))
    println("  Keys only in A: $(setdiff(keys_A, keys_B))")
    println("  Keys only in B: $(setdiff(keys_B, keys_A))")
    for k in intersect(keys_A, keys_B)
        if ncf_as_child_A[k] != ncf_as_child_B[k]
            println("  key=$k: A_ncf_idx=$(ncf_as_child_A[k]) B_ncf_idx=$(ncf_as_child_B[k])")
        end
    end

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON 3: do ncf_as_parent_A and ncf_as_parent_B agree?
    # ══════════════════════════════════════════════════════════════════════
    println("Rank $rank: === COMPARISON 3: ncf_as_parent A vs B ===")
    keysP_A = Set(keys(ncf_as_parent_A)); keysP_B = Set(keys(ncf_as_parent_B))
    println("  Keys only in A: $(setdiff(keysP_A, keysP_B))")
    println("  Keys only in B: $(setdiff(keysP_B, keysP_A))")
    for k in intersect(keysP_A, keysP_B)
        sA = Set(ncf_as_parent_A[k]); sB = Set(ncf_as_parent_B[k])
        if sA != sB
            println("  key=$k: A=$(ncf_as_parent_A[k]) B=$(ncf_as_parent_B[k])")
        end
    end

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON 4: do ncf_L1L2_A and ncf_L1L2_B agree, entry by entry?
    # ══════════════════════════════════════════════════════════════════════
    println("Rank $rank: === COMPARISON 4: ncf_L1L2 A vs B ===")
    n_diff_L = 0
    for ncf_idx = 1:mesh.num_ncf
        L1A, L2A = ncf_L1L2_A[ncf_idx]
        L1B, L2B = ncf_L1L2_B[ncf_idx]
        if !isapprox(L1A, L1B; atol=1e-10) || !isapprox(L2A, L2B; atol=1e-10)
            n_diff_L += 1
            n_diff_L <= 3 && println("  ncf_idx=$ncf_idx: L1 diff=$(maximum(abs.(L1A.-L1B))), L2 diff=$(maximum(abs.(L2A.-L2B)))")
        end
    end
    println("  Total differing ncf_L1L2 entries: $n_diff_L / $(mesh.num_ncf)")

    # ══════════════════════════════════════════════════════════════════════
    # COMPARISON 5: run BOTH serial sweeps (using set A's tables, and set
    # B's tables) with build_conforming_sweep_order, compare results.
    # ══════════════════════════════════════════════════════════════════════
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

    ncf_child_to_parent_elemid_A = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child_A)
    ncf_parent_to_children_A = Dict{Tuple{Int,Int}, Vector{Int}}()
    for (key, group) in ncf_as_parent_A
        ncf_parent_to_children_A[key] = [Int(mesh.cip[idx]) for idx in group]
    end

    order_A, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem_A,
                                               metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                               metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                               metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                               ncf_child_to_parent = ncf_child_to_parent_elemid_A,
                                               ncf_parent_to_children = ncf_parent_to_children_A)
    x_dg_A = fill(NaN, npoin_dg)
    run_dg_sweep!(x_dg_A, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order_A, neighbor_elem_A,
                  connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                  ncf_as_child=ncf_as_child_A, ncf_as_parent=ncf_as_parent_A, ncf_L1L2=ncf_L1L2_A)

    order_B, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem_B,
                                               metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                               metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                               metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                               ncf_child_to_parent = ncf_child_to_parent_elemid_B,
                                               ncf_parent_to_children = ncf_parent_to_children_B)
    x_dg_B = fill(NaN, npoin_dg)
    run_dg_sweep!(x_dg_B, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order_B, neighbor_elem_B,
                  connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                  ncf_as_child=ncf_as_child_B, ncf_as_parent=ncf_as_parent_B, ncf_L1L2=ncf_L1L2_B)

    println("Rank $rank: === COMPARISON 5: serial sweep using A-tables vs B-tables ===")
    for (i,j,k) in [(1,1,1), (ngl,1,1), (1,ngl,1), (1,1,ngl), (ngl,ngl,ngl)]
        gip = mesh.connijk[trace_iel,i,j,k]
        exact = bc_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        vA = x_dg_A[connijk_dg[trace_iel,i,j,k]]
        vB = x_dg_B[connijk_dg[trace_iel,i,j,k]]
        println("  (i,j,k)=($i,$j,$k): A=$vA B=$vB exact=$exact  A_err=$(abs(vA-exact)) B_err=$(abs(vB-exact))")
    end

    max_err_A = 0.0; max_err_B = 0.0
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        exact = bc_fn(mesh.x[g], mesh.y[g], mesh.z[g])
        max_err_A = max(max_err_A, abs(x_dg_A[connijk_dg[iel,i,j,k]] - exact))
        max_err_B = max(max_err_B, abs(x_dg_B[connijk_dg[iel,i,j,k]] - exact))
    end
    println("  GLOBAL max error using A-tables: $max_err_A")
    println("  GLOBAL max error using B-tables: $max_err_B")

    return (neighbor_elem_A=neighbor_elem_A, neighbor_elem_B=neighbor_elem_B,
            ncf_as_child_A=ncf_as_child_A, ncf_as_child_B=ncf_as_child_B,
            ncf_as_parent_A=ncf_as_parent_A, ncf_as_parent_B=ncf_as_parent_B,
            ncf_L1L2_A=ncf_L1L2_A, ncf_L1L2_B=ncf_L1L2_B,
            x_dg_A=x_dg_A, x_dg_B=x_dg_B, max_err_A=max_err_A, max_err_B=max_err_B)
end

function test_dg_sweep_spatial_ncf_all_directions_415(inputs, nparts, distribute, sem, partitioned_model;
                                                       spatial_target_iel::Int = 415,
                                                       κ_val=0.0001, I0=4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics; ngl = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω
    extra_npoin = mesh.extra_mesh.extra_npoin

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    mesh.num_ncf == 0 && (println("⚠ no NCF"); return nothing)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end
    ncf_child_to_parent_elemid = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    worst_err = 0.0; worst_ip_ext = 0
    all_errs = Float64[]
    n_cycle_failures = 0

    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]; ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

        local order
        try
            order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                                      metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                                      ncf_child_to_parent = ncf_child_to_parent_elemid,
                                                      ncf_parent_to_children = ncf_parent_to_children)
        catch e
            n_cycle_failures += 1
            println("ip_ext=$ip_ext (Ω=($Ωx,$Ωy,$Ωz)): ✗ CYCLE/ORDERING FAILURE — $(sprint(showerror, e))")
            continue
        end

        x_dg = fill(NaN, npoin_dg)
        run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                      connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                      ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2, ncf_as_parent=ncf_as_parent)

        n_nan = count(isnan, x_dg)
        if n_nan > 0
            println("ip_ext=$ip_ext: ⚠ $n_nan unwritten entries — skipping error calc")
            continue
        end

        max_err = 0.0
        for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
            g = mesh.connijk[iel,i,j,k]
            I_exact = bc_fn(mesh.x[g],mesh.y[g],mesh.z[g])
            max_err = max(max_err, abs(x_dg[connijk_dg[iel,i,j,k]] - I_exact))
        end
        push!(all_errs, max_err)
        if max_err > worst_err
            worst_err = max_err; worst_ip_ext = ip_ext
        end
    end

    println("─"^60)
    println("Cycle/ordering failures: $n_cycle_failures / $extra_npoin")
    if !isempty(all_errs)
        println("Worst direction: ip_ext=$worst_ip_ext, max abs error = $worst_err")
        println("Mean error: $(sum(all_errs)/length(all_errs))")
        println("Median error: $(sort(all_errs)[length(all_errs)÷2 + 1])")
    end
    return all_errs
end

function test_dg_sweep_multirank_spatial_ncf_all_directions(inputs, nparts, distribute, sem, partitioned_model;
                                                             spatial_target_iel::Int = 415,
                                                             κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh    = setup.mesh
    metrics = setup.metrics
    sem2    = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl
    extra_npoin = mesh.extra_mesh.extra_npoin

    neighbor_elem = setup.neighbor_elem
    neighbor_fid  = setup.neighbor_fid
    ncf_as_child  = setup.ncf_as_child
    ncf_as_parent = setup.ncf_as_parent
    neighbor_rank_ncf = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group

    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    # ── Geometry-only setup, built ONCE, reused across every direction ──────
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
    end

    conforming_rank = Dict{Tuple{TInt,TInt}, TInt}()
    conforming_gel  = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        conforming_rank[(TInt(iel),TInt(fid))] = neighbor_rank_ncf[(iel,fid)]
        conforming_gel[(TInt(iel),TInt(fid))]  = neighbor_gel_group[(iel,fid)][1]
    end
    exchange = build_rank_cut_exchange(mesh, conforming_rank, conforming_gel, comm)
    if exchange.n_unmatched > 0 || exchange.n_mismatch > 0
        error("Rank $rank: conforming exchange setup failed")
    end
    value_cache = build_rank_cut_value_cache(exchange, comm)

    ncf_coord_cache, ncf_gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank_ncf, neighbor_kind,
                                                                         neighbor_gel_group, comm)
    ncf_exch = build_ncf_value_exchange(mesh, neighbor_rank_ncf, connijk_dg, comm)

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    worst_err = 0.0; worst_ip_ext = 0
    all_errs = Float64[]
    n_nan_failures = 0

    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]; ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

        x_dg = fill(NaN, npoin_dg)
        run_dg_sweep_multirank_spatial_ncf!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz,
                                             neighbor_elem, conforming_rank, conforming_gel, exchange, value_cache,
                                             ncf_as_child, ncf_L1L2, ncf_as_parent,
                                             neighbor_kind, neighbor_gel_group,
                                             ncf_coord_cache, ncf_gelpos_to_gip, ncf_exch,
                                             connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn, comm)

        n_nan_local = count(isnan, x_dg)
        n_nan_total = MPI.Allreduce(n_nan_local, +, comm)
        if n_nan_total > 0
            n_nan_failures += 1
            println("Rank $rank, ip_ext=$ip_ext: ⚠ $n_nan_local unwritten locally ($n_nan_total total) — skipping")
            continue
        end

        max_err_local = 0.0
        for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
            g = mesh.connijk[iel,i,j,k]
            I_exact = bc_fn(mesh.x[g], mesh.y[g], mesh.z[g])
            max_err_local = max(max_err_local, abs(x_dg[connijk_dg[iel,i,j,k]] - I_exact))
        end
        max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)

        push!(all_errs, max_err_global)
        if max_err_global > worst_err
            worst_err = max_err_global; worst_ip_ext = ip_ext
        end
        rank == 0 && println("ip_ext=$ip_ext: max_err=$max_err_global")
    end

    if rank == 0
        println("─"^60)
        println("NaN failures: $n_nan_failures / $extra_npoin")
        if !isempty(all_errs)
            println("Worst direction: ip_ext=$worst_ip_ext, max abs error = $worst_err")
            println("Mean error: $(sum(all_errs)/length(all_errs))")
            println("Median error: $(sort(all_errs)[length(all_errs)÷2 + 1])")
        end
    end
    return all_errs
end

#=function build_ncf_angular_value_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                           n_dirs::Vector{Int}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]
    flat_my_gel  = TInt[]
    flat_d       = Int[]
    flat_spos    = Int[]
    flat_owning_iel = Int[]
    flat_i = Int[]; flat_j = Int[]; flat_k = Int[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for d = 1:n_dirs[iel]
            spos = 0
            for k in krange, j in jrange, i in irange
                spos += 1
                push!(flat_targets, Int(target_rank))
                push!(flat_my_gel, my_gel)
                push!(flat_d, d)
                push!(flat_spos, spos)
                push!(flat_owning_iel, iel)
                push!(flat_i, i); push!(flat_j, j); push!(flat_k, k)
            end
        end
    end

    n_flat = length(flat_targets)
    gel_cache  = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    d_cache    = SendReceiveCache(comm, Float64.(flat_d), flat_targets)
    spos_cache = SendReceiveCache(comm, Float64.(flat_spos), flat_targets)
    val_cache  = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(d_cache, Float64.(flat_d))
    send_and_receive!(spos_cache, Float64.(flat_spos))

    return (val_cache=val_cache, gel_cache=gel_cache, d_cache=d_cache, spos_cache=spos_cache,
            flat_owning_iel=flat_owning_iel, flat_d=flat_d,   # NOW stored for sending
            flat_i=flat_i, flat_j=flat_j, flat_k=flat_k, n_flat=n_flat)
end=#

function build_ncf_angular_value_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                           n_dirs::Vector{Int}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]
    flat_my_gel  = TInt[]
    flat_d       = Int[]
    flat_fid     = Int[]    # NEW — the face this position is relative to
    flat_spos    = Int[]
    flat_owning_iel = Int[]
    flat_i = Int[]; flat_j = Int[]; flat_k = Int[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for d = 1:n_dirs[iel]
            spos = 0
            for k in krange, j in jrange, i in irange
                spos += 1
                push!(flat_targets, Int(target_rank))
                push!(flat_my_gel, my_gel)
                push!(flat_d, d)
                push!(flat_fid, fid)
                push!(flat_spos, spos)
                push!(flat_owning_iel, iel)
                push!(flat_i, i); push!(flat_j, j); push!(flat_k, k)
            end
        end
    end

    n_flat = length(flat_targets)
    gel_cache  = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    d_cache    = SendReceiveCache(comm, Float64.(flat_d), flat_targets)
    fid_cache  = SendReceiveCache(comm, Float64.(flat_fid), flat_targets)   # NEW
    spos_cache = SendReceiveCache(comm, Float64.(flat_spos), flat_targets)
    val_cache  = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(d_cache, Float64.(flat_d))
    send_and_receive!(fid_cache, Float64.(flat_fid))   # NEW
    send_and_receive!(spos_cache, Float64.(flat_spos))

    return (val_cache=val_cache, gel_cache=gel_cache, d_cache=d_cache, fid_cache=fid_cache,
            spos_cache=spos_cache, flat_owning_iel=flat_owning_iel, flat_d=flat_d,
            flat_i=flat_i, flat_j=flat_j, flat_k=flat_k, n_flat=n_flat)
end

function exchange_ncf_angular_values!(exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]
        d   = exch.flat_d[m]
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        send_vals[m] = haskey(x_dg, (iel,d)) ? x_dg[(iel,d)][lidx(i,j,k,ngl)] : NaN
    end
    send_and_receive!(exch.val_cache, send_vals)

    received = Dict{Tuple{TInt,Int,Int,Int}, Float64}()   # KEY: (gel,d,fid,spos) — 4-tuple now
    for m = 1:length(exch.val_cache.combined_recv_data)
        val = exch.val_cache.combined_recv_data[m]
        isnan(val) && continue
        rgel  = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd    = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rfid  = Int(round(Int, exch.fid_cache.combined_recv_data[m]))   # NEW
        rspos = Int(round(Int, exch.spos_cache.combined_recv_data[m]))
        received[(rgel,rd,rfid,rspos)] = val
    end
    return received
end

#=function exchange_ncf_angular_values!(exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]
        d   = exch.flat_d[m]
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        send_vals[m] = haskey(x_dg, (iel,d)) ? x_dg[(iel,d)][lidx(i,j,k,ngl)] : NaN
    end
    send_and_receive!(exch.val_cache, send_vals)

    received = Dict{Tuple{TInt,Int,Int}, Float64}()
    for m = 1:length(exch.val_cache.combined_recv_data)
        val = exch.val_cache.combined_recv_data[m]
        isnan(val) && continue
        rgel  = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd    = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rspos = Int(round(Int, exch.spos_cache.combined_recv_data[m]))
        received[(rgel,rd,rspos)] = val
    end
    return received
end=#

function test_ncf_angular_value_exchange(inputs, nparts, distribute, sem, partitioned_model;
                                          angular_target_iel::Int = 442, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_angular_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        angular_target_iel=angular_target_iel,
                                        angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    ngl = mesh.ngl

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    exch = build_ncf_angular_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    println("Rank $rank: built exchange with $(exch.n_flat) flat entries")

    # ── Fabricated x_dg: every (iel,d) at every local node encodes its own
    # (my_gel, d, spatial_pos) into the value — but spatial_pos depends on
    # WHICH face we're on, and a node can be on MULTIPLE faces (edges/
    # corners), so instead we encode based on (gel, d) ONLY at the node's
    # OWN (i,j,k), verified per flat entry using the SAME (i,j,k) the
    # exchange itself will read — i.e., value = gel*1e6 + d*1e3 + LOCAL
    # flat-index encoded via i,j,k directly (not spos, to avoid ambiguity
    # at shared edge/corner nodes touching multiple faces). ──────────────────
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        my_gel = Int(mesh.el2gel[iel])
        for d = 1:n_dirs[iel]
            vals = zeros(ngl^3)
            for k=1:ngl, j=1:ngl, i=1:ngl
                # Encode (my_gel, d, i, j, k) uniquely and DEcodably.
                vals[lidx(i,j,k,ngl)] = my_gel*1e8 + d*1e5 + i*1e3 + j*1e1 + k
            end
            x_dg[(iel,d)] = vals
        end
    end

    received = exchange_ncf_angular_values!(exch, x_dg, ngl)
    println("Rank $rank: received $(length(received)) entries")

    # ── Verify: for each received (rgel,rd,rspos)->val, decode val and
    # check it matches what we'd expect IF we knew the remote element's
    # OWN (i,j,k) at that (fid,spos) — but we don't have that directly, so
    # instead verify SELF-CONSISTENCY: does the SAME (rgel,rd) always
    # decode to a gel matching rgel and a d matching rd? ─────────────────────
    n_checked = 0
    n_mismatch = 0
    for ((rgel,rd,rspos), val) in received
        n_checked += 1
        decoded_gel = Int(floor(val / 1e8))
        decoded_d   = Int(floor((val - decoded_gel*1e8) / 1e5))
        if decoded_gel != rgel || decoded_d != rd
            n_mismatch += 1
            println("Rank $rank: ✗ MISMATCH (rgel=$rgel,rd=$rd,rspos=$rspos): " *
                    "decoded_gel=$decoded_gel decoded_d=$decoded_d")
        end
    end
    println("Rank $rank: checked $n_checked entries, $n_mismatch mismatches  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))

    return exch, received
end

"""
    build_angular_ncf_upwind_unified_remote(iel, fid, remote_gel, d, src_e_ext,
                                             mesh, em_remote_placeholder, received, ngl)

Cross-rank analogue of build_angular_ncf_upwind_unified: the source
element's angular mesh lives on a REMOTE rank (already reconstructed via
remote_angular_mesh, from exchange_angular_mesh_across_rank_cut), and its
SOLVED VALUES arrive via `received::Dict{(gel,d,spos)}` from
exchange_ncf_angular_values!, keyed by spatial position (spos), matching
face_ip_range's canonical order (same confirmed-safe convention as the
spatial case).
"""
function build_angular_ncf_upwind_unified_remote(mesh, iel::Int, remote_gel::TInt, d::Int,
                                                  src_e_ext::Int, em_remote, received::Dict{Tuple{TInt,Int,Int},Float64})
    em_iel = mesh.extra_mesh[iel]
    θd = em_iel.extra_coords[1,d]; ϕd = em_iel.extra_coords[2,d]

    nsrc = em_remote.extra_nop[src_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds_ragged(em_remote.extra_coords, em_remote.extra_connijk, em_remote.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd

    θ_src = [em_remote.extra_coords[1, em_remote.extra_connijk[src_e_ext][a,1]] for a=1:nsrc]
    ϕ_src = Float64[]
    for b=1:nsrc
        raw = em_remote.extra_coords[2, em_remote.extra_connijk[src_e_ext][1,b]]
        push!(ϕ_src, (near_seam && raw<1e-8) ? 2π : raw)
    end
    Lθ = build_1d_lagrange_matrix(θ_src, [θd]); Lϕ = build_1d_lagrange_matrix(ϕ_src, [ϕd_fix])
    src_ips = [em_remote.extra_connijk[src_e_ext][a,b] for a=1:nsrc, b=1:nsrc]

    return (i,j,k) -> begin
        gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
        val = 0.0
        for a=1:nsrc, b=1:nsrc
            v = get(received, (remote_gel, src_ips[a,b], gid), NaN)
            val += Lθ[1,a]*Lϕ[1,b]*v
        end
        val
    end
end

function find_containing_angular_element_ragged(em_ragged, node_set::Set{Int})
    found_e = 0
    for e = 1:em_ragged.extra_nelem
        n = em_ragged.extra_nop[e] + 1
        elem_nodes = Set{Int}()
        grid = em_ragged.extra_connijk[e]
        for a = 1:n, b = 1:n
            push!(elem_nodes, grid[a,b])
        end
        elem_nodes == node_set && (found_e = e)
    end
    return found_e
end

function test_angular_ncf_upwind_unified_remote(inputs, nparts, distribute, sem, partitioned_model;
                                                 angular_target_iel::Int = 442, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_angular_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        angular_target_iel=angular_target_iel,
                                        angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    angular_deps_cr = setup.angular_deps_cr
    remote_angular_mesh = setup.remote_angular_mesh
    ngl = mesh.ngl

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    exch = build_ncf_angular_value_exchange(mesh, neighbor_rank, n_dirs, comm)

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    Pspatial(t) = sum(0.3^k*t^k for k=0:ngl-1)
    Pang(t) = 0.3 + 0.4*t

    # ── Find a GROUP-case cross-rank entry FIRST, so we know src_e_ext's
    # true node count (nsrc) BEFORE fabricating data — Pang's degree must
    # match nsrc-1, exactly as the validated single-rank test did. ──────────
    test_key = nothing
    for ((iel,fid,d), (remote_gel, remote_nodes)) in angular_deps_cr
        length(remote_nodes) > 1 || continue
        test_key = (iel, d, remote_gel, remote_nodes)
        break
    end
    if test_key === nothing
        println("Rank $rank: no group-case entry found")
        return nothing
    end
    test_iel, test_d, remote_gel, remote_nodes = test_key

    test_fid = 0
    for fid = 1:6
        haskey(neighbor_rank, (test_iel,fid)) || continue
        test_fid = fid; break
    end
    println("Rank $rank: testing (iel=$test_iel, d=$test_d, fid=$test_fid), remote_gel=$remote_gel")

    em_remote = remote_angular_mesh[remote_gel]
    src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
    if src_e_ext == 0
        println("Rank $rank: ⚠ could not locate src_e_ext")
        return nothing
    end

    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds_ragged(em_remote.extra_coords, em_remote.extra_connijk,
                                                              em_remote.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ
    nsrc = em_remote.extra_nop[src_e_ext] + 1
    field(x,y,z,θ,ϕ) = Pspatial((x-x0)/Lx) + Pspatial((y-y0)/Ly) + Pspatial((z-z0)/Lz) +
                    Pang(θ/π) + Pang(ϕ/(2π))

    # ── Fabricate x_dg ONLY for the remote element's directions actually
    # needed (remote_nodes) — using THIS SAME field, evaluated at the
    # REMOTE element's OWN (θ,ϕ) for each of its directions. But x_dg here
    # is OUR OWN storage (we're sending OUR values to be received on the
    # OTHER rank) — meaning we need to fabricate x_dg for OUR OWN elements
    # that will be READ by the remote side as ITS dependency. Since we are
    # testing OURSELVES receiving remote_gel's data, we fabricate as if WE
    # WERE remote_gel's owning rank — i.e. populate x_dg for test_iel's
    # OWN local elements using field(), and TRUST that whichever rank
    # OWNS remote_gel does the SAME (this test only meaningfully checks
    # what happens ON THIS RANK, receiving from the OTHER rank's ACTUAL
    # fabricated send — so BOTH ranks must run this SAME test, fabricating
    # with the SAME field formula, for the comparison to be valid). ─────────
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θd = em.extra_coords[1,d]; ϕd = em.extra_coords[2,d]
            vals = zeros(ngl^3)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                vals[lidx(i,j,k,ngl)] = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θd, fix_ϕ(ϕd))
            end
            x_dg[(iel,d)] = vals
        end
    end

    received = exchange_ncf_angular_values!(exch, x_dg, ngl)

    closure = build_angular_ncf_upwind_unified_remote(test_iel, test_fid, remote_gel, test_d,
                                                        src_e_ext, mesh, em_remote, received, ngl)

    θd_test = mesh.extra_mesh[test_iel].extra_coords[1,test_d]
    ϕd_test = mesh.extra_mesh[test_iel].extra_coords[2,test_d]

    irange, jrange, krange = face_ip_range(test_fid, ngl)
    max_err = 0.0
    for (i,j,k) in Iterators.product(irange,jrange,krange)
        computed = closure(i,j,k)
        gip = mesh.connijk[test_iel,i,j,k]
        exact = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θd_test, fix_ϕ(ϕd_test))
        max_err = max(max_err, abs(computed-exact))
    end
    println("Rank $rank: remote angular upwind closure max abs error = $max_err")
    return max_err
end

function resolve_stall_distributed_gauss_seidel_physics!(mesh, metrics, ψ, dψ, ω,
                                                          neighbor_elem, face_correspondence, neighbor_gel,
                                                          angular_deps, angular_deps_cr, remote_angular_mesh,
                                                          ncf_exch, x_dg, solved, n_dirs,
                                                          κ, σ, intΦ_eff, κ_val, I0, comm;
                                                          max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-physics over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total)")

    for (iel,d) in candidates
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    received = exchange_ncf_angular_values!(ncf_exch, x_dg, ngl)

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0.0

        for (iel,d) in candidates
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
                I_upwind = (i,j,k) -> 0.0

                if is_inflow[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        if haskey(angular_deps_cr, (iel,fid,d))
                            remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                            if length(remote_nodes) == 1
                                rn = remote_nodes[1]
                                fully_avail = all(haskey(received, (remote_gel, rn, spos)) for spos = 1:ngl2)
                                if fully_avail
                                    rn_c = rn
                                    I_upwind = (i,j,k) -> begin
                                        spos = face_flat_index(fid, ngl, i, j, k)
                                        received[(remote_gel, rn_c, spos)]
                                    end
                                end
                            else
                                fully_avail = all(haskey(received, (remote_gel, rn, spos))
                                                   for rn in remote_nodes for spos = 1:ngl2)
                                if fully_avail
                                    em_remote = remote_angular_mesh[remote_gel]
                                    src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                    if src_e_ext != 0
                                        I_upwind = build_angular_ncf_upwind_unified_remote(
                                            iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received, ngl)
                                    end
                                end
                            end
                        elseif haskey(neighbor_gel, (iel,fid))
                            # ── NEW: ordinary conforming cross-rank face ────────
                            remote_gel_ord = neighbor_gel[(iel,fid)]
                            fully_avail = all(haskey(received, (remote_gel_ord, d, fid, spos)) for spos = 1:ngl2)
                            if fully_avail
                                gel_c = remote_gel_ord; d_c = d
                                I_upwind = (i,j,k) -> begin
                                    spos = face_flat_index(fid, ngl, i, j, k)
                                    received[(gel_c, d_c, spos)]
                                end
                            end
                        else
                            I_upwind = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                            end
                        end
                    else
                        key = (iel,d)
                        handled = false
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                handled = true
                                if length(dep_nodes) == 1
                                    dj = dep_nodes[1]; fid_c=fid
                                    if haskey(x_dg,(jel,dj))
                                        I_upwind = (i,j,k) -> begin
                                            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                            x_dg[(jel,dj)][lidx(corr[1],corr[2],corr[3],ngl)]
                                        end
                                    end
                                else
                                    if all(haskey(x_dg,(jel,dn)) for dn in dep_nodes)
                                        em_jel = mesh.extra_mesh[jel]
                                        src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                        I_upwind = build_angular_ncf_upwind_unified(
                                            iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    end
                                end
                            end
                        end
                        if !handled && haskey(x_dg,(jel,d))
                            d_c=d; fid_c=fid
                            I_upwind = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end

        received = exchange_ncf_angular_values!(ncf_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-physics iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-physics did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

function run_dg_sweep_multirank_angular_ncf!(mesh, metrics, ψ, dψ, ω,
                                              neighbor_elem, neighbor_fid, face_correspondence,
                                              neighbor_rank, neighbor_gel, angular_deps, angular_deps_cr,
                                              remote_angular_mesh, ncf_exch,
                                              κ, σ, intΦ_eff, κ_val::Float64, I0::Float64, comm::MPI.Comm;
                                              max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    received = Dict{Tuple{TInt,Int,Int}, Float64}()

    while n_solved_total < total_pairs
        iter_round += 1

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]

                    if jel_raw == zero(TInt)
                        if haskey(neighbor_rank, (iel,fid))
                            key = (iel,fid,d)
                            if haskey(angular_deps_cr, key)
                                remote_gel, remote_nodes = angular_deps_cr[key]

                                if length(remote_nodes) == 1
                                    rn = remote_nodes[1]
                                    fully_avail = all(haskey(received, (remote_gel, rn, spos)) for spos = 1:ngl2)
                                    if !fully_avail
                                        ready = false; break
                                    else
                                        rn_c = rn
                                        upwind_fns[fid] = (i,j,k) -> begin
                                            spos = face_flat_index(fid, ngl, i, j, k)
                                            received[(remote_gel, rn_c, spos)]
                                        end
                                    end
                                else
                                    fully_avail = all(haskey(received, (remote_gel, rn, spos))
                                                       for rn in remote_nodes for spos = 1:ngl2)
                                    if !fully_avail
                                        ready = false; break
                                    else
                                        em_remote = remote_angular_mesh[remote_gel]
                                        src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                        if src_e_ext == 0
                                            ready = false; break
                                        else
                                            upwind_fns[fid] = build_angular_ncf_upwind_unified_remote(
                                                iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received, ngl)
                                        end
                                    end
                                end
                            elseif haskey(neighbor_gel, (iel,fid))
                                # ── Ordinary conforming cross-rank face —
                                # direction d means the SAME thing on both
                                # sides (angularly conforming), plain lookup. ──
                                remote_gel_ord = neighbor_gel[(iel,fid)]
                                fully_avail = all(haskey(received, (remote_gel_ord, d, fid, spos)) for spos = 1:ngl2)
                                if !fully_avail
                                    ready = false; break
                                else
                                    gel_c = remote_gel_ord; d_c = d
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        spos = face_flat_index(fid, ngl, i, j, k)
                                        received[(gel_c, d_c, spos)]
                                    end
                                end
                            else
                                ready = false; break   # rank-cut face with no correspondence info at all — genuine gap, block rather than silently zero
                            end
                        else
                            bc = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1])+Ωy*(y-r_bc[2])+Ωz*(z-r_bc[3])))
                            upwind_fns[fid] = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                            end
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    key = (iel,d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]; ready=false
                                else
                                    dj_c=dj; fid_c=fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved; ready=false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c=d; fid_c=fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
        received = exchange_ncf_angular_values!(ncf_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        if n_newly_solved_global == 0
            if just_resolved_stall
                error("run_dg_sweep_multirank_angular_ncf!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_physics!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence, neighbor_gel, angular_deps, angular_deps_cr, remote_angular_mesh,
                ncf_exch, x_dg, solved, n_dirs, κ, σ, intΦ_eff, κ_val, I0, comm;
                max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
            received = exchange_ncf_angular_values!(ncf_exch, x_dg, ngl)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs
end

function test_dg_sweep_multirank_angular_ncf_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                          angular_target_iel::Int = 442, angular_target_e_ext::Int = 4,
                                                          κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_angular_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        angular_target_iel=angular_target_iel,
                                        angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh; metrics = setup.metrics
    sem2 = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl

    neighbor_elem = setup.neighbor_elem
    neighbor_fid  = setup.neighbor_fid
    neighbor_rank = setup.neighbor_rank
    neighbor_gel  = setup.neighbor_gel
    angular_deps  = setup.angular_deps
    angular_deps_cr = setup.angular_deps_cr
    remote_angular_mesh = setup.remote_angular_mesh

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    ncf_exch = build_ncf_angular_value_exchange(mesh, neighbor_rank, n_dirs, comm)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    x_dg, n_dirs_ret = run_dg_sweep_multirank_angular_ncf!(mesh, metrics, ψ, dψ, ω,
                                                            neighbor_elem, neighbor_fid, face_correspondence,
                                                            neighbor_rank, neighbor_gel, angular_deps, angular_deps_cr,
                                                            remote_angular_mesh, ncf_exch,
                                                            κ, σz, intΦ_eff, κ_val, I0, comm)

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                max_err_local = max(max_err_local, abs(x_dg[(iel,d)][lidx(i,j,k,ngl)] - exact))
            end
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Multirank angular-NCF exponential test — max abs error: $max_err_global")
    return max_err_global
end


function debug_angular_deps_cr_correspondence(mesh, angular_deps_cr, remote_angular_mesh; n_show::Int=10, tol::Float64=1e-6)
    rank = MPI.Comm_rank(get_mpi_comm())
    ngl = mesh.ngl
    n_checked = 0
    n_bad = 0

    for ((iel,fid,d), (remote_gel, remote_nodes)) in angular_deps_cr
        n_checked >= n_show && break
        n_checked += 1

        em_iel = mesh.extra_mesh[iel]
        θd = em_iel.extra_coords[1,d]; ϕd = em_iel.extra_coords[2,d]

        if !haskey(remote_angular_mesh, remote_gel)
            println("Rank $rank: ✗ (iel=$iel,d=$d) claims remote_gel=$remote_gel, NOT in remote_angular_mesh")
            n_bad += 1
            continue
        end
        em_remote = remote_angular_mesh[remote_gel]

        if length(remote_nodes) == 1
            # SINGLETON — the remote node's own (θ,ϕ) should EXACTLY match (θd,ϕd)
            rn = remote_nodes[1]
            if rn < 1 || rn > em_remote.extra_npoin
                println("Rank $rank: ✗ (iel=$iel,d=$d) SINGLETON rn=$rn OUT OF BOUNDS " *
                        "(remote npoin=$(em_remote.extra_npoin))")
                n_bad += 1
                continue
            end
            θr = em_remote.extra_coords[1,rn]; ϕr = em_remote.extra_coords[2,rn]
            θ_ok = isapprox(θd, θr; atol=tol)
            ϕ_ok = isapprox(ϕd, ϕr; atol=tol) || isapprox(abs(ϕd-ϕr), 2π; atol=tol)  # seam-tolerant
            status = (θ_ok && ϕ_ok) ? "✓" : "✗ MISMATCH"
            (θ_ok && ϕ_ok) || (n_bad += 1)
            println("Rank $rank: $status (iel=$iel,d=$d) SINGLETON: mine=(θ=$θd,ϕ=$ϕd) " *
                    "remote(rn=$rn)=(θ=$θr,ϕ=$ϕr)")
        else
            # GROUP — my (θd,ϕd) should fall INSIDE the containing remote
            # element's bounds; ALSO verify remote_nodes forms a genuine
            # single element's full node set.
            src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
            if src_e_ext == 0
                println("Rank $rank: ✗ (iel=$iel,d=$d) GROUP remote_nodes=$remote_nodes " *
                        "does NOT match any single remote element's node set")
                n_bad += 1
                continue
            end
            θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds_ragged(em_remote.extra_coords,
                                                                       em_remote.extra_connijk, em_remote.extra_nop, src_e_ext)
            near_seam = isapprox(ϕmax, 2π; atol=1e-6)
            ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd
            inside = (θd >= θmin-tol && θd <= θmax+tol && ϕd_fix >= ϕmin-tol && ϕd_fix <= ϕmax+tol)
            status = inside ? "✓" : "✗ OUTSIDE BOUNDS"
            inside || (n_bad += 1)
            println("Rank $rank: $status (iel=$iel,d=$d) GROUP: mine=(θ=$θd,ϕ=$ϕd) " *
                    "remote src_e_ext=$src_e_ext bounds=θ∈[$θmin,$θmax] ϕ∈[$ϕmin,$ϕmax]")
        end
    end

    println("Rank $rank: checked $n_checked entries, $n_bad bad  " * (n_bad==0 ? "✓ ALL PASS" : "✗ FOUND ISSUES"))
    return n_checked, n_bad
end

function debug_compare_pre_stall_values(mesh, x_dg_multirank, n_dirs, κ_val, I0, comm)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    max_err_pre_stall = 0.0
    n_checked = 0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            haskey(x_dg_multirank, (iel,d)) || continue
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                v = x_dg_multirank[(iel,d)][lidx(i,j,k,ngl)]
                err = abs(v-exact)
                max_err_pre_stall = max(max_err_pre_stall, err)
                n_checked += 1
            end
        end
    end
    println("Rank $rank: checked $n_checked node-direction values, max error = $max_err_pre_stall")
end

function debug_trace_cross_rank_exchange_soundness(mesh, x_dg, ncf_exch, received, solved, n_dirs,
                                                     angular_deps_cr, κ_val, I0, comm; n_show::Int=5)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    # ── STEP 1: for a few of MY OWN solved (iel,d) pairs that get SENT
    # across the rank-cut boundary, check what value I'm ABOUT TO SEND is
    # actually correct (matches exact solution, if it's a true-boundary-
    # fed direction — otherwise just report the value for inspection). ──────
    n_sent_checked = 0
    for m = 1:ncf_exch.n_flat
        n_sent_checked >= n_show && break
        iel = ncf_exch.flat_owning_iel[m]; d = ncf_exch.flat_d[m]
        solved[iel][d] || continue
        i,j,k = ncf_exch.flat_i[m], ncf_exch.flat_j[m], ncf_exch.flat_k[m]
        gip = mesh.connijk[iel,i,j,k]
        em = mesh.extra_mesh[iel]
        θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
        sent_val = x_dg[(iel,d)][lidx(i,j,k,ngl)]
        my_gel = Int(mesh.el2gel[iel])
        println("Rank $rank SEND-CHECK: (iel=$iel,d=$d,gel=$my_gel) at (i,j,k)=($i,$j,$k) " *
                "sending=$sent_val exact=$exact err=$(abs(sent_val-exact))")
        n_sent_checked += 1
    end

    # ── STEP 2: for a few RECEIVED entries, report the raw value and try to
    # cross-check it against the exact solution IF we can independently
    # compute the physical (x,y,z) for that remote (gel,d,spos) — which we
    # generally CANNOT without the remote element's mesh, so instead we
    # report the raw received value paired with its tag for manual
    # inspection against the SEND-CHECK output from the OTHER rank. ─────────
    n_recv_checked = 0
    for (key, val) in received
        n_recv_checked >= n_show && break
        rgel, rd, rspos = key
        println("Rank $rank RECV-CHECK: remote(gel=$rgel,d=$rd,spos=$rspos) received_value=$val")
        n_recv_checked += 1
    end

    # ── STEP 3: for a cross-rank-dependent (iel,d) that's ABOUT TO become
    # solvable THIS round, show exactly which received[...] entries its
    # I_upwind will pull from, and what value results — BEFORE the actual
    # solve happens (this must be called from INSIDE the sweep, right
    # before build_dg_flux_terms! for a genuinely-ready cross-rank face). ───
    n_consume_checked = 0
    for ((iel,fid,d), (remote_gel, remote_nodes)) in angular_deps_cr
        n_consume_checked >= n_show && break
        solved[iel][d] && continue  # only check NOT-YET-solved (about to be)
        if length(remote_nodes) == 1
            rn = remote_nodes[1]
            avail = all(haskey(received,(remote_gel,rn,fid,spos)) for spos=1:ngl*ngl)
            avail || continue
            sample_val = received[(remote_gel,rn,1)]
            println("Rank $rank CONSUME-CHECK: (iel=$iel,d=$d) about to use SINGLETON " *
                    "remote(gel=$remote_gel,d=$rn) sample_value(spos=1)=$sample_val")
        else
            avail = all(haskey(received,(remote_gel,rn,fid,spos)) for rn in remote_nodes for spos=1:ngl*ngl)
            avail || continue
            sample_vals = [received[(remote_gel,rn,1)] for rn in remote_nodes[1:min(3,end)]]
            println("Rank $rank CONSUME-CHECK: (iel=$iel,d=$d) about to use GROUP " *
                    "remote_gel=$remote_gel sample_values(spos=1, first 3 nodes)=$sample_vals")
        end
        n_consume_checked += 1
    end
end

"""
    build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received,
                                                child_iel, fid, d, parent_global_gips, ngl)

Per-direction extension of build_dg_ncf_parent_to_child_remote: identical
geometry/interpolation-weight construction (coord_cache, unchanged), but
reads the remote parent's face-node VALUES for THIS SPECIFIC direction `d`
from `received::Dict{(gel,d,spos)}` (the universal per-direction exchange)
instead of `value_by_gip::Dict{Int,Float64}` (which had no direction
dimension at all — single-direction sweeps only).
"""
function build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache::Dict{Int,NTuple{3,Float64}},
                                                     received::Dict{Tuple{TInt,Int,Int},Float64},
                                                     child_iel::Int, fid::Int, d::Int,
                                                     remote_gel::TInt, parent_global_gips::Vector{Int}, ngl::Int)
    IPc_nodes = Int[]
    irange, jrange, krange = face_ip_range(fid, ngl)
    for k in krange, j in jrange, i in irange
        push!(IPc_nodes, mesh.connijk[child_iel,i,j,k])
    end

    L1, L2 = build_spatial_interpolation_matrices_with_coord_cache(
        fid, IPc_nodes, parent_global_gips, mesh, coord_cache, ngl)

    ngl2 = ngl*ngl
    parent_grid = zeros(ngl, ngl)
    for idx = 1:ngl2
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        # KEY CHANGE: value now comes from received[(remote_gel, d, spos)],
        # NOT value_by_gip[gip] — same spos convention (idx here IS the
        # canonical face_ip_range position, matching the exchange's tagging).
        parent_grid[p1,p2] = get(received, (remote_gel, d, idx), NaN)
    end

    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        ci = div(node_idx-1, ngl) + 1
        cj = mod(node_idx-1, ngl) + 1
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += w_1[p1]*w_2[p2]*parent_grid[p1,p2]
        end
        child_vals[node_idx] = val
    end
    return child_vals
end

"""
    build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received,
                                                parent_iel, fid, d,
                                                remote_children_gels, gelpos_to_gip, ngl)

Per-direction extension of build_dg_ncf_child_to_parent_remote — same
per-node covering-child containment logic (unchanged, using coord_cache),
but values come from received[(gel,d,spos)] for the given direction `d`.
"""
function build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache::Dict{Int,NTuple{3,Float64}},
                                                     received::Dict{Tuple{TInt,Int,Int},Float64},
                                                     parent_iel::Int, fid::Int, d::Int,
                                                     remote_children_gels::Vector{TInt},
                                                     gelpos_to_gip::Dict{Tuple{TInt,Int}, Int}, ngl::Int)
    ngl2 = ngl*ngl
    irange, jrange, krange = face_ip_range(fid, ngl)
    IPp_nodes = Int[]
    for k in krange, j in jrange, i in irange
        push!(IPp_nodes, mesh.connijk[parent_iel,i,j,k])
    end

    parent_result = fill(NaN, ngl, ngl)

    for remote_gel in remote_children_gels
        child_global_gips = [gelpos_to_gip[(remote_gel,pos)] for pos = 1:ngl2]

        child_coords = [coord_cache[gip] for gip in child_global_gips]
        xs = [c[1] for c in child_coords]; ys = [c[2] for c in child_coords]; zs = [c[3] for c in child_coords]

        child_line_1_idx = [(i-1)*ngl+1 for i = 1:ngl]
        child_line_2_idx = [j for j = 1:ngl]

        xs1=[xs[p] for p in child_line_2_idx]; ys1=[ys[p] for p in child_line_2_idx]; zs1=[zs[p] for p in child_line_2_idx]
        axis_c1 = axis_of_max_spread((xs1,ys1,zs1))
        xs2=[xs[p] for p in child_line_1_idx]; ys2=[ys[p] for p in child_line_1_idx]; zs2=[zs[p] for p in child_line_1_idx]
        axis_c2 = axis_of_max_spread((xs2,ys2,zs2))
        coord_of(p, axis) = axis==1 ? xs[p] : (axis==2 ? ys[p] : zs[p])

        c1_coords = [coord_of(p,axis_c1) for p in child_line_2_idx]
        c2_coords = [coord_of(p,axis_c2) for p in child_line_1_idx]
        c1_min,c1_max = extrema(c1_coords); c2_min,c2_max = extrema(c2_coords)
        tol1 = 1e-6*max(c1_max-c1_min,1.0); tol2 = 1e-6*max(c2_max-c2_min,1.0)

        # KEY CHANGE: gather values for THIS direction d, from received.
        child_grid = zeros(ngl,ngl)
        for idx = 1:ngl2
            c1 = ((idx-1)%ngl)+1; c2 = ((idx-1)÷ngl)+1
            child_grid[c1,c2] = get(received, (remote_gel, d, idx), NaN)
        end

        for idx = 1:length(IPp_nodes)
            node = IPp_nodes[idx]
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            pc1 = coord_of_local(mesh, node, axis_c1)
            pc2 = coord_of_local(mesh, node, axis_c2)
            (pc1 < c1_min-tol1 || pc1 > c1_max+tol1 || pc2 < c2_min-tol2 || pc2 > c2_max+tol2) && continue

            L_c1 = build_1d_lagrange_matrix(c1_coords, [pc1])
            L_c2 = build_1d_lagrange_matrix(c2_coords, [pc2])
            val = 0.0
            for a=1:ngl, b=1:ngl
                val += L_c1[1,a]*L_c2[1,b]*child_grid[a,b]
            end
            parent_result[p1,p2] = val
        end
    end
    return parent_result
end

"""
    build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)

Per-direction value exchange for the SPATIAL-NCF remote closures — same
(gel, d, spos) self-describing tagging as the angular version, but built
and used SEPARATELY (per the decision to keep spatial and angular
cross-rank mechanisms independent for this first working version, with
unification deferred to a later optimization pass).
"""
function build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                                   n_dirs::Vector{Int}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]
    flat_my_gel  = TInt[]
    flat_d       = Int[]
    flat_spos    = Int[]
    flat_owning_iel = Int[]
    flat_i = Int[]; flat_j = Int[]; flat_k = Int[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for d = 1:n_dirs[iel]
            spos = 0
            for k in krange, j in jrange, i in irange
                spos += 1
                push!(flat_targets, Int(target_rank))
                push!(flat_my_gel, my_gel)
                push!(flat_d, d)
                push!(flat_spos, spos)
                push!(flat_owning_iel, iel)
                push!(flat_i, i); push!(flat_j, j); push!(flat_k, k)
            end
        end
    end

    n_flat = length(flat_targets)
    gel_cache  = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    d_cache    = SendReceiveCache(comm, Float64.(flat_d), flat_targets)
    spos_cache = SendReceiveCache(comm, Float64.(flat_spos), flat_targets)
    val_cache  = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(d_cache, Float64.(flat_d))
    send_and_receive!(spos_cache, Float64.(flat_spos))

    return (val_cache=val_cache, gel_cache=gel_cache, d_cache=d_cache, spos_cache=spos_cache,
            flat_owning_iel=flat_owning_iel, flat_d=flat_d,
            flat_i=flat_i, flat_j=flat_j, flat_k=flat_k, n_flat=n_flat)
end

"""
    exchange_ncf_spatial_perdir_values!(exch, x_dg, ngl) -> Dict{Tuple{TInt,Int,Int},Float64}

Same mechanics as exchange_ncf_angular_values!, kept as a SEPARATE
function (not shared) per the decision to avoid unifying the two
mechanisms until a later optimization pass.
"""
function exchange_ncf_spatial_perdir_values!(exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]
        d   = exch.flat_d[m]
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        send_vals[m] = haskey(x_dg, (iel,d)) ? x_dg[(iel,d)][lidx(i,j,k,ngl)] : NaN
    end
    send_and_receive!(exch.val_cache, send_vals)

    received = Dict{Tuple{TInt,Int,Int}, Float64}()
    for m = 1:length(exch.val_cache.combined_recv_data)
        val = exch.val_cache.combined_recv_data[m]
        isnan(val) && continue
        rgel  = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd    = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rspos = Int(round(Int, exch.spos_cache.combined_recv_data[m]))
        received[(rgel,rd,rspos)] = val
    end
    return received
end

function test_ncf_spatial_perdir_closures_remote(inputs, nparts, distribute, sem, partitioned_model;
                                                  spatial_target_iel::Int = 303)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_spatial_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel=spatial_target_iel)
    mesh = setup.mesh
    neighbor_rank = setup.neighbor_rank
    neighbor_kind = setup.neighbor_kind
    neighbor_gel_group = setup.neighbor_gel_group
    ngl = mesh.ngl

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)

    # ── Fake per-direction structure: pretend every element has, say, 5
    # directions, just to exercise the per-direction plumbing (this test
    # is purely about the SPATIAL interpolation correctness across
    # directions, not real angular meshes at all). ─────────────────────────
    n_dirs = fill(5, mesh.nelem)
    exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)

    # ── GLOBAL, rank-independent field — same discipline as the angular
    # closure fix: no locally-discovered normalization. ─────────────────────
    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z,d) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz) + 0.1*d   # small, EXACT per-direction offset

    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        vals = zeros(ngl^3)
        for k=1:ngl,j=1:ngl,i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            vals[lidx(i,j,k,ngl)] = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], d)
        end
        x_dg[(iel,d)] = vals
    end

    received = exchange_ncf_spatial_perdir_values!(exch, x_dg, ngl)

    # ── Test PARENT->CHILD (:i_am_fine side) ─────────────────────────────────
    test_iel = 0; test_fid = 0; remote_gel_fine = zero(TInt)
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_fine || continue
        test_iel = iel; test_fid = fid
        remote_gel_fine = neighbor_gel_group[(iel,fid)][1]
        break
    end
    if test_iel != 0
        test_d = 3
        ngl2 = ngl*ngl
        parent_global_gips = [gelpos_to_gip[(remote_gel_fine,pos)] for pos = 1:ngl2]
        child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received,
            test_iel, test_fid, test_d, remote_gel_fine, parent_global_gips, ngl)

        irange,jrange,krange = face_ip_range(test_fid, ngl)
        IPc_nodes = Int[]
        for k in krange, j in jrange, i in irange
            push!(IPc_nodes, mesh.connijk[test_iel,i,j,k])
        end
        max_err = 0.0
        for idx = 1:length(IPc_nodes)
            node = IPc_nodes[idx]
            exact = field(mesh.x[node], mesh.y[node], mesh.z[node], test_d)
            max_err = max(max_err, abs(child_vals[idx]-exact))
        end
        println("Rank $rank: parent->child PER-DIR remote max abs error = $max_err")
    else
        println("Rank $rank: no :i_am_fine face — skipping parent->child test")
    end

    # ── Test CHILD->PARENT (:i_am_coarse side) ───────────────────────────────
    test_iel2 = 0; test_fid2 = 0; remote_gels_coarse = TInt[]
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_coarse || continue
        test_iel2 = iel; test_fid2 = fid
        remote_gels_coarse = neighbor_gel_group[(iel,fid)]
        break
    end
    if test_iel2 != 0
        test_d2 = 3
        parent_computed = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received,
            test_iel2, test_fid2, test_d2, remote_gels_coarse, gelpos_to_gip, ngl)

        irange,jrange,krange = face_ip_range(test_fid2, ngl)
        IPp_nodes = Int[]
        for k in krange, j in jrange, i in irange
            push!(IPp_nodes, mesh.connijk[test_iel2,i,j,k])
        end
        max_err2 = 0.0
        for idx = 1:length(IPp_nodes)
            node = IPp_nodes[idx]
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            exact = field(mesh.x[node], mesh.y[node], mesh.z[node], test_d2)
            max_err2 = max(max_err2, abs(parent_computed[p1,p2]-exact))
        end
        println("Rank $rank: child->parent PER-DIR remote max abs error = $max_err2")
    else
        println("Rank $rank: no :i_am_coarse face — skipping child->parent test")
    end
end

"""
    setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                 spatial_target_iel, angular_target_iel, angular_target_e_ext)

Merged setup for the combined spatial+angular NCF-across-MPI sweep: both
refinements applied (far apart, per the one-buffer-element policy), all
four correspondence/exchange structures built.
"""
function setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                      spatial_target_iel::Int = 415,
                                      angular_target_iel::Int = 442, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    # ── Spatial refinement (rank 0 only, via sem_setup) ──────────────────────
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    # ── Angular refinement (only meaningful on the rank owning the target) ──
    if angular_target_iel <= mesh.nelem
        adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
        verify_angular_phi_seam_conformity(mesh, angular_target_iel)
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Same-rank spatial-NCF tables ─────────────────────────────────────────
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ngl = mesh.ngl
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
    end
    println("Rank $rank: spatial NCF: $(mesh.num_ncf) same-rank entries")

    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    # ── Rank-cut identification, excluding same-rank NCF faces ──────────────
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end
    suspected_cuts = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in tagged_face_keys && continue
            push!(suspected_cuts, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_cuts)) suspected rank-cut faces")

    # ── Cross-rank correspondence: NCF-aware (handles BOTH spatially
    # conforming AND non-conforming rank-cut faces in one pass). ────────────
    neighbor_rank, neighbor_kind, neighbor_gel_group, n_unresolved =
        resolve_rank_cut_neighbors_ncf(mesh, suspected_cuts, get_mpi_comm())
    println("Rank $rank: $(length(neighbor_rank)) resolved, $n_unresolved unresolved")
    if n_unresolved > 0
        error("Rank $rank: unresolved rank-cut faces — aborting")
    end

    # Plain (iel,fid)->remote_gel, for ordinary conforming faces only (needed
    # by BOTH the ordinary-cross-rank branch and the coordinate exchange).
    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    # ── Cross-rank spatial-NCF: coordinate exchange (geometry only) ─────────
    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)

    # ── Cross-rank angular-NCF: full angular mesh exchange + dependency graph ─
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    println("Rank $rank: $(length(angular_deps_cr)) cross-rank angular NCF entries")

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    return (mesh=mesh, metrics=metrics, sem=sem, partitioned_model=partitioned_model,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid, face_correspondence=face_correspondence,
            ncf_as_child=ncf_as_child, ncf_as_parent=ncf_as_parent, ncf_L1L2=ncf_L1L2,
            neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind, neighbor_gel_group=neighbor_gel_group,
            neighbor_gel=neighbor_gel, coord_cache=coord_cache, gelpos_to_gip=gelpos_to_gip,
            remote_angular_mesh=remote_angular_mesh, angular_deps=angular_deps, angular_deps_cr=angular_deps_cr)
end

function test_setup_combined_ncf_mpi(inputs, nparts, distribute, sem, partitioned_model;
                                      spatial_target_iel::Int = 303,
                                      angular_target_iel::Int = 330, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())

    setup = setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                         spatial_target_iel=spatial_target_iel,
                                         angular_target_iel=angular_target_iel,
                                         angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh

    println("Rank $rank: === SETUP VALIDATION ===")

    # ── 1. Basic counts ──────────────────────────────────────────────────────
    println("  mesh.nelem=$(mesh.nelem), num_ncf=$(mesh.num_ncf)")
    println("  neighbor_rank entries: $(length(setup.neighbor_rank))")
    println("  neighbor_kind breakdown:")
    for kind in (:conforming, :i_am_coarse, :i_am_fine)
        n = count(==(kind), values(setup.neighbor_kind))
        println("    $kind: $n")
    end
    println("  neighbor_gel (ordinary conforming) entries: $(length(setup.neighbor_gel))")
    println("  ncf_as_child (same-rank spatial) entries: $(length(setup.ncf_as_child))")
    println("  ncf_as_parent (same-rank spatial) entries: $(length(setup.ncf_as_parent))")
    println("  angular_deps_cr (cross-rank angular) entries: $(length(setup.angular_deps_cr))")
    println("  angular_deps (same-rank angular) entries: $(length(setup.angular_deps))")
    println("  coord_cache entries: $(length(setup.coord_cache))")
    println("  remote_angular_mesh entries: $(length(setup.remote_angular_mesh))")

    # ── 2. Consistency: every neighbor_kind key should ALSO appear in
    # neighbor_rank (they're built from the same source, should never diverge). ──
    n_kind_not_in_rank = count(k -> !haskey(setup.neighbor_rank, k), keys(setup.neighbor_kind))
    println("  neighbor_kind keys missing from neighbor_rank: $n_kind_not_in_rank  " *
            (n_kind_not_in_rank==0 ? "✓" : "✗"))

    # ── 3. Every :conforming face should have an entry in neighbor_gel ───────
    n_conforming = count(==(:conforming), values(setup.neighbor_kind))
    n_conforming_missing_gel = 0
    for ((iel,fid), kind) in setup.neighbor_kind
        kind == :conforming || continue
        haskey(setup.neighbor_gel, (TInt(iel),TInt(fid))) || (n_conforming_missing_gel += 1)
    end
    println("  :conforming faces missing from neighbor_gel: $n_conforming_missing_gel / $n_conforming  " *
            (n_conforming_missing_gel==0 ? "✓" : "✗"))

    # ── 4. Sanity: same-rank NCF and cross-rank sets should be DISJOINT
    # (a face can't be BOTH same-rank NCF and cross-rank — same face-key
    # space (iel,fid), should never collide). ────────────────────────────────
    same_rank_ncf_keys = union(Set(keys(setup.ncf_as_child)),
                                Set(k for (k,v) in setup.ncf_as_parent))
    cross_rank_keys = Set(keys(setup.neighbor_rank))
    n_overlap = length(intersect(same_rank_ncf_keys, cross_rank_keys))
    println("  Overlap between same-rank-NCF and cross-rank face keys: $n_overlap  " *
            (n_overlap==0 ? "✓" : "✗ SHOULD BE ZERO"))

    # ── 5. Every cross-rank angular NCF (iel,d) key's iel should ALSO be a
    # key (for SOME fid) in neighbor_rank (it's a cross-rank face). ─────────
    n_dep_cr_orphan = 0
    for (iel,fid,d) in keys(setup.angular_deps_cr)
        any(haskey(setup.neighbor_rank,(iel,fid)) for fid=1:6) || (n_dep_cr_orphan += 1)
    end
    println("  angular_deps_cr entries whose iel has NO neighbor_rank face: $n_dep_cr_orphan  " *
            (n_dep_cr_orphan==0 ? "✓" : "✗"))

    return setup
end

function resolve_stall_distributed_gauss_seidel_combined!(mesh, metrics, ψ, dψ, ω,
                                                            neighbor_elem, face_correspondence,
                                                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                            coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                            angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                            x_dg, solved, n_dirs,
                                                            κ, σ, intΦ_eff, κ_val, I0, comm;
                                                            max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-combined over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total)")

    for (iel,d) in candidates
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
    received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0.0

        for (iel,d) in candidates
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
                I_upwind = (i,j,k) -> 0.0

                if is_inflow[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            if haskey(x_dg,(parent_iel,d))
                                L1,L2 = ncf_L1L2[ncf_idx]
                                IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                                child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d,
                                                                                 IPc_nodes, IPp_nodes, L1, L2, ngl)
                                I_upwind = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                            end
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            children = [Int(mesh.cip[idx]) for idx in ncf_group]
                            if all(c -> haskey(x_dg,(c,d)), children)
                                entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                           for idx in ncf_group]
                                parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
                                I_upwind = (i,j,k) -> begin
                                    ni = face_flat_index(fid,ngl,i,j,k)
                                    p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                    parent_grid[p1,p2]
                                end
                            end
                        elseif haskey(neighbor_rank, (iel,fid))
                            kind = neighbor_kind[(iel,fid)]

                            if kind == :i_am_fine
                                remote_gel_f = neighbor_gel_group[(iel,fid)][1]
                                if all(haskey(received_sp,(remote_gel_f,d,spos)) for spos=1:ngl2)
                                    parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
                                    child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gel_f, parent_global_gips, ngl)
                                    I_upwind = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                                end
                            elseif kind == :i_am_coarse
                                remote_gels_c = neighbor_gel_group[(iel,fid)]
                                if all(haskey(received_sp,(g,d,spos)) for g in remote_gels_c for spos=1:ngl2)
                                    parent_grid = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gels_c, gelpos_to_gip, ngl)
                                    I_upwind = (i,j,k) -> begin
                                        ni=face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                        parent_grid[p1,p2]
                                    end
                                end
                            elseif kind == :conforming
                                # ── SAME FIX AS THE MAIN SWEEP: check
                                # angular_deps_cr FIRST, per direction, before
                                # trusting the spatially-conforming ordinary path. ──
                                if haskey(angular_deps_cr, (iel,fid,d))
                                    remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                                    if length(remote_nodes) == 1
                                        rn = remote_nodes[1]
                                        if all(haskey(received_ang,(remote_gel,rn,fid,spos)) for spos=1:ngl2)
                                            rn_c=rn
                                            fid_cc = fid
					    I_upwind = (i,j,k) -> received_ang[(remote_gel,rn_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
                                        end
                                    else
                                        if all(haskey(received_ang,(remote_gel,rn,fid,spos)) for rn in remote_nodes for spos=1:ngl2)
                                            em_remote = remote_angular_mesh[remote_gel]
                                            src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                            if src_e_ext != 0
                                                I_upwind = build_angular_ncf_upwind_unified_remote(
                                                    iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received_ang, ngl)
                                            end
                                        end
                                    end
                                else
                                    remote_gel_ord = neighbor_gel[(iel,fid)]
                                    if all(haskey(received_ang,(remote_gel_ord,d,fid,spos)) for spos=1:ngl2)
                                        gel_c=remote_gel_ord; d_c=d
                                        fid_cc = fid
					I_upwind = (i,j,k) -> received_ang[(gel_c,d_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
                                    end
                                end
                            end
                        else
                            I_upwind = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                            end
                        end
                    else
                        key = (iel,d)
                        handled = false
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                handled = true
                                if length(dep_nodes) == 1
                                    dj = dep_nodes[1]; fid_c=fid
                                    if haskey(x_dg,(jel,dj))
                                        I_upwind = (i,j,k) -> begin
                                            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                            x_dg[(jel,dj)][lidx(corr[1],corr[2],corr[3],ngl)]
                                        end
                                    end
                                else
                                    if all(haskey(x_dg,(jel,dn)) for dn in dep_nodes)
                                        em_jel = mesh.extra_mesh[jel]
                                        src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                        I_upwind = build_angular_ncf_upwind_unified(
                                            iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    end
                                end
                            end
                        end
                        if !handled && haskey(x_dg,(jel,d))
                            d_c=d; fid_c=fid
                            I_upwind = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end

        received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-combined iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-combined did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

function run_dg_sweep_multirank_combined_ncf_instrumented!(mesh, metrics, ψ, dψ, ω,
                                                            neighbor_elem, neighbor_fid, face_correspondence,
                                                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                            coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                            angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                            κ, σ, intΦ_eff, κ_val::Float64, I0::Float64, comm::MPI.Comm;
                                                            max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9,
                                                            err_threshold::Float64 = 1.0)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    exact_at(iel,Ωx,Ωy,Ωz,i,j,k) = begin
        gip = mesh.connijk[iel,i,j,k]
        I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
    end

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    received_ang = Dict{Tuple{TInt,Int,Int,Int}, Float64}()
    received_sp  = Dict{Tuple{TInt,Int,Int}, Float64}()

    while n_solved_total < total_pairs
        iter_round += 1
        n_newly_solved = 0

        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()
                upwind_source = Dict{Int,String}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]

                    if jel_raw == zero(TInt)
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            if d > n_dirs[parent_iel] || !solved[parent_iel][d]
                                ready = false; break
                            else
                                L1,L2 = ncf_L1L2[ncf_idx]
                                IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                                child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d,
                                                                                 IPc_nodes, IPp_nodes, L1, L2, ngl)
                                upwind_fns[fid] = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                                upwind_source[fid] = "SAME-RANK-SPATIAL-NCF-CHILD"
                            end
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            children = [Int(mesh.cip[idx]) for idx in ncf_group]
                            if any(c -> d > n_dirs[c] || !solved[c][d], children)
                                ready = false; break
                            else
                                entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                           for idx in ncf_group]
                                parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
                                upwind_fns[fid] = (i,j,k) -> begin
                                    ni = face_flat_index(fid,ngl,i,j,k)
                                    p1 = ((ni-1)%ngl)+1; p2 = ((ni-1)÷ngl)+1
                                    parent_grid[p1,p2]
                                end
                                upwind_source[fid] = "SAME-RANK-SPATIAL-NCF-PARENT"
                            end
                        elseif haskey(neighbor_rank, (iel,fid))
                            kind = neighbor_kind[(iel,fid)]

                            if kind == :i_am_fine
                                remote_gel_f = neighbor_gel_group[(iel,fid)][1]
                                fully_avail = all(haskey(received_sp, (remote_gel_f, d, spos)) for spos=1:ngl2)
                                if !fully_avail
                                    ready = false; break
                                else
                                    parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
                                    child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gel_f, parent_global_gips, ngl)
                                    upwind_fns[fid] = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                                    upwind_source[fid] = "CROSS-RANK-SPATIAL-NCF-FINE"
                                end
                            elseif kind == :i_am_coarse
                                remote_gels_c = neighbor_gel_group[(iel,fid)]
                                fully_avail = all(haskey(received_sp,(g,d,spos)) for g in remote_gels_c for spos=1:ngl2)
                                if !fully_avail
                                    ready = false; break
                                else
                                    parent_grid = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gels_c, gelpos_to_gip, ngl)
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        ni = face_flat_index(fid,ngl,i,j,k)
                                        p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                        parent_grid[p1,p2]
                                    end
                                    upwind_source[fid] = "CROSS-RANK-SPATIAL-NCF-COARSE"
                                end
                            elseif kind == :conforming
                                # SPATIALLY conforming — MUST ALSO check ANGULAR
                                # correspondence PER-DIRECTION before trusting
                                # the ordinary same-index lookup (a remote
                                # element can be PARTIALLY angularly refined).
                                if haskey(angular_deps_cr, (iel,fid,d))
                                    remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                                    if length(remote_nodes) == 1
                                        rn = remote_nodes[1]
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos)) for spos=1:ngl2)
					if !fully_avail
                                            ready = false; break
                                        else
                                            rn_c=rn
                                            upwind_fns[fid] = (i,j,k) -> begin
                                                spos = face_flat_index(fid,ngl,i,j,k)
                                                received_ang[(remote_gel,rn_c,fid,spos)]
                                            end
                                            upwind_source[fid] = "CROSS-RANK-ANGULAR-NCF-SINGLETON"
                                        end
                                    else
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos)) for rn in remote_nodes for spos=1:ngl2)
					if !fully_avail
                                            ready = false; break
                                        else
                                            em_remote = remote_angular_mesh[remote_gel]
                                            src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                            if src_e_ext == 0
                                                ready = false; break
                                            else
                                                upwind_fns[fid] = build_angular_ncf_upwind_unified_remote(
                                                    iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received_ang, ngl)
                                                upwind_source[fid] = "CROSS-RANK-ANGULAR-NCF-GROUP"
                                            end
                                        end
                                    end
                                else
                                    # Genuinely conforming in BOTH space and
                                    # angle for THIS direction — safe ordinary lookup.
                                    remote_gel_ord = neighbor_gel[(iel,fid)]
                                    fully_avail = all(haskey(received_ang, (remote_gel_ord, d, fid, spos)) for spos=1:ngl2)
				    if !fully_avail
                                        ready = false; break
                                    else
                                        gel_c=remote_gel_ord; d_c=d
                                        upwind_fns[fid] = (i,j,k) -> begin
                                            spos = face_flat_index(fid,ngl,i,j,k)
                                            received_ang[(gel_c,d_c,fid,spos)]
                                        end
                                        upwind_source[fid] = "CROSS-RANK-ORDINARY"
                                    end
                                end
                            else
                                error("run_dg_sweep_multirank_combined_ncf!: face (iel=$iel,fid=$fid) has " *
                                      "neighbor_kind=$kind — not one of {:conforming,:i_am_fine,:i_am_coarse}.")
                            end
                        end
                        haskey(upwind_fns,fid) || begin
                            bc = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1])+Ωy*(y-r_bc[2])+Ωz*(z-r_bc[3])))
                            upwind_fns[fid] = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                            end
                            upwind_source[fid] = "TRUE-BOUNDARY"
                        end
                        continue
                    end
                    jel = Int(jel_raw)
                    key = (iel,d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]; ready=false
                                else
                                    dj_c=dj; fid_c=fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                    upwind_source[fid] = "SAME-RANK-ANGULAR-SINGLETON"
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved; ready=false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    upwind_source[fid] = "SAME-RANK-ANGULAR-GROUP"
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c=d; fid_c=fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                            upwind_source[fid] = "SAME-RANK-ORDINARY"
                        end
                    end
                    ready || break
                end
                ready || continue

                for fid = 1:6
                    haskey(upwind_fns,fid) || continue
                    irange,jrange,krange = face_ip_range(fid, ngl)
                    for (i,j,k) in Iterators.product(irange,jrange,krange)
                        val = upwind_fns[fid](i,j,k)
                        ex = exact_at(iel,Ωx,Ωy,Ωz,i,j,k)
                        if abs(val-ex) > err_threshold
                            println("Rank $rank ⚠⚠⚠ LARGE I_upwind ERROR: (iel=$iel,d=$d,fid=$fid,i=$i,j=$j,k=$k) " *
                                    "source=$(upwind_source[fid]) val=$val exact=$ex")
                        end
                    end
                end

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    if is_active[fid] && is_inflow[fid] && !haskey(upwind_fns,fid)
                        println("Rank $rank ⚠⚠⚠ INCOMPLETE upwind_fns: (iel=$iel,d=$d,fid=$fid) — SILENT ZERO FALLBACK")
                    end
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1

                for k=1:ngl,j=1:ngl,i=1:ngl
                    ex = exact_at(iel,Ωx,Ωy,Ωz,i,j,k)
                    val = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                    if abs(val-ex) > err_threshold
                        println("Rank $rank ⚠⚠⚠ LARGE SOLVE ERROR: (iel=$iel,d=$d,i=$i,j=$j,k=$k) val=$val exact=$ex")
                    end
                end
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        if n_newly_solved_global == 0
            if just_resolved_stall
                error("run_dg_sweep_multirank_combined_ncf!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_combined!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                coord_cache, gelpos_to_gip, spatial_perdir_exch,
                angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                x_dg, solved, n_dirs, κ, σ, intΦ_eff, κ_val, I0, comm;
                max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
            received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
            received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs
end

function test_dg_sweep_multirank_combined_ncf_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                            spatial_target_iel::Int = 415,
                                                            angular_target_iel::Int = 442, angular_target_e_ext::Int = 4,
                                                            κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                         spatial_target_iel=spatial_target_iel,
                                         angular_target_iel=angular_target_iel,
                                         angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh; metrics = setup.metrics
    sem2 = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    angular_exch = build_ncf_angular_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    x_dg, n_dirs_ret = run_dg_sweep_multirank_combined_ncf_instrumented!(
        mesh, metrics, ψ, dψ, ω,
        setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence,
        setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group,
        setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
        setup.coord_cache, setup.gelpos_to_gip, spatial_perdir_exch,
        setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh, angular_exch,
        κ, σz, intΦ_eff, κ_val, I0, comm)

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                max_err_local = max(max_err_local, abs(x_dg[(iel,d)][lidx(i,j,k,ngl)] - exact))
            end
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Multirank COMBINED NCF exponential test — max abs error: $max_err_global")
    return max_err_global
end

function test_check_exchange_crosstalk(mesh, neighbor_rank, neighbor_gel, angular_exch, spatial_perdir_exch,
                                        x_dg, ngl, comm; iel_check::Int=35, d_check::Int=31)
    rank = MPI.Comm_rank(comm)

    received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
    received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

    # Find the specific (fid, remote_gel) for iel_check that's cross-rank
    test_fid = 0; remote_gel = zero(TInt)
    for fid = 1:6
        haskey(neighbor_rank, (iel_check,fid)) || continue
        test_fid = fid
        remote_gel = neighbor_gel[(iel_check,fid)]
        break
    end
    if test_fid == 0
        println("Rank $rank: iel=$iel_check has no cross-rank face — cannot test")
        return nothing
    end
    println("Rank $rank: checking (iel=$iel_check, d=$d_check, fid=$test_fid), remote_gel=$remote_gel")

    ngl2 = ngl*ngl
    n_mismatch = 0
    n_checked = 0
    for spos = 1:ngl2
        key = (remote_gel, d_check, spos)
        has_ang = haskey(received_ang, key)
        has_sp  = haskey(received_sp, key)
        val_ang = has_ang ? received_ang[key] : nothing
        val_sp  = has_sp  ? received_sp[key]  : nothing
        n_checked += 1
        if has_ang && has_sp
            if !isapprox(val_ang, val_sp; atol=1e-10)
                n_mismatch += 1
                println("Rank $rank: ✗ MISMATCH at spos=$spos: received_ang=$val_ang received_sp=$val_sp")
            end
        end
        println("Rank $rank: spos=$spos: has_ang=$has_ang (val=$val_ang), has_sp=$has_sp (val=$val_sp)")
    end
    println("Rank $rank: checked $n_checked spos values, $n_mismatch mismatches between the two exchanges " *
            "(where BOTH had an entry for the same key)")

    return received_ang, received_sp
end

function debug_check_direction_match(mesh, remote_angular_mesh, iel_check::Int, d_check::Int, remote_gel::TInt)
    rank = MPI.Comm_rank(get_mpi_comm())
    em_local = mesh.extra_mesh[iel_check]
    θ_local = em_local.extra_coords[1, d_check]
    ϕ_local = em_local.extra_coords[2, d_check]
    println("Rank $rank: LOCAL (iel=$iel_check, d=$d_check): (θ,ϕ)=($θ_local,$ϕ_local)")

    if haskey(remote_angular_mesh, remote_gel)
        em_remote = remote_angular_mesh[remote_gel]
        if d_check <= em_remote.extra_npoin
            θ_remote = em_remote.extra_coords[1, d_check]
            ϕ_remote = em_remote.extra_coords[2, d_check]
            println("Rank $rank: REMOTE (gel=$remote_gel, d=$d_check): (θ,ϕ)=($θ_remote,$ϕ_remote)")
            match = isapprox(θ_local,θ_remote;atol=1e-6) && isapprox(ϕ_local,ϕ_remote;atol=1e-6)
            println("Rank $rank: MATCH? $match")
        else
            println("Rank $rank: REMOTE gel=$remote_gel only has $(em_remote.extra_npoin) directions — d=$d_check out of range")
        end
    else
        println("Rank $rank: remote_gel=$remote_gel not found in remote_angular_mesh " *
                "(this element's angular mesh was never exchanged for THIS purpose — see note)")
    end
end

function run_dg_sweep_multirank_combined_ncf!(mesh, metrics, ψ, dψ, ω,
                                               neighbor_elem, neighbor_fid, face_correspondence,
                                               neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                               ncf_as_child, ncf_L1L2, ncf_as_parent,
                                               coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                               angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                               κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                               extra_source_all = nothing,
                                               max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    received_ang = Dict{Tuple{TInt,Int,Int,Int}, Float64}()
    received_sp  = Dict{Tuple{TInt,Int,Int}, Float64}()

    while n_solved_total < total_pairs
        iter_round += 1
        n_newly_solved = 0

        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]

                    if jel_raw == zero(TInt)
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            if d > n_dirs[parent_iel] || !solved[parent_iel][d]
                                ready = false; break
                            else
                                L1,L2 = ncf_L1L2[ncf_idx]
                                IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                                child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d,
                                                                                 IPc_nodes, IPp_nodes, L1, L2, ngl)
                                upwind_fns[fid] = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                            end
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            children = [Int(mesh.cip[idx]) for idx in ncf_group]
                            if any(c -> d > n_dirs[c] || !solved[c][d], children)
                                ready = false; break
                            else
                                entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                           for idx in ncf_group]
                                parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
                                upwind_fns[fid] = (i,j,k) -> begin
                                    ni = face_flat_index(fid,ngl,i,j,k)
                                    p1 = ((ni-1)%ngl)+1; p2 = ((ni-1)÷ngl)+1
                                    parent_grid[p1,p2]
                                end
                            end
                        elseif haskey(neighbor_rank, (iel,fid))
                            kind = neighbor_kind[(iel,fid)]
                            if kind == :i_am_fine
                                remote_gel_f = neighbor_gel_group[(iel,fid)][1]
                                fully_avail = all(haskey(received_sp, (remote_gel_f, d, spos)) for spos=1:ngl2)
                                if !fully_avail
                                    ready = false; break
                                else
                                    parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
                                    child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gel_f, parent_global_gips, ngl)
                                    upwind_fns[fid] = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                                end
                            elseif kind == :i_am_coarse
                                remote_gels_c = neighbor_gel_group[(iel,fid)]
                                fully_avail = all(haskey(received_sp,(g,d,spos)) for g in remote_gels_c for spos=1:ngl2)
                                if !fully_avail
                                    ready = false; break
                                else
                                    parent_grid = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gels_c, gelpos_to_gip, ngl)
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        ni = face_flat_index(fid,ngl,i,j,k)
                                        p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                        parent_grid[p1,p2]
                                    end
                                end
                            elseif kind == :conforming
                                if haskey(angular_deps_cr, (iel,fid,d))
                                    remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                                    if length(remote_nodes) == 1
                                        rn = remote_nodes[1]
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos)) for spos=1:ngl2)
					if !fully_avail
                                            ready = false; break
                                        else
                                            rn_c=rn
                                            fid_cc = fid
					    upwind_fns[fid] = (i,j,k) -> received_ang[(remote_gel,rn_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
					end
                                    else
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos)) for rn in remote_nodes for spos=1:ngl2)
					if !fully_avail
                                            ready = false; break
                                        else
                                            em_remote = remote_angular_mesh[remote_gel]
                                            src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                            if src_e_ext == 0
                                                ready = false; break
                                            else
                                                upwind_fns[fid] = build_angular_ncf_upwind_unified_remote(
                                                    iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received_ang, ngl)
                                            end
                                        end
                                    end
                                else
                                    remote_gel_ord = neighbor_gel[(iel,fid)]
                                    fully_avail = all(haskey(received_ang, (remote_gel_ord, d, fid, spos)) for spos=1:ngl2)
				    if !fully_avail
                                        ready = false; break
                                    else
                                        gel_c=remote_gel_ord; d_c=d
                                        fid_cc = fid
					upwind_fns[fid] = (i,j,k) -> received_ang[(gel_c,d_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
				    end
                                end
                            else
                                error("run_dg_sweep_multirank_combined_ncf!: bad neighbor_kind=$kind")
                            end
                        end
                        haskey(upwind_fns,fid) || begin
                            bc = bc_fn_all[(iel,d)]
                            upwind_fns[fid] = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                            end
                        end
                        continue
                    end
                    jel = Int(jel_raw)
                    key = (iel,d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]; ready=false
                                else
                                    dj_c=dj; fid_c=fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved; ready=false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c=d; fid_c=fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        if n_newly_solved_global == 0
            if just_resolved_stall
                error("run_dg_sweep_multirank_combined_ncf!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_combined!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                coord_cache, gelpos_to_gip, spatial_perdir_exch,
                angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                x_dg, solved, n_dirs, κ, σ, intΦ_eff, bc_fn_all, comm;
                extra_source_all=extra_source_all, max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
            received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
            received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs
end 

function resolve_stall_distributed_gauss_seidel_combined!(mesh, metrics, ψ, dψ, ω,
                                                            neighbor_elem, face_correspondence,
                                                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                            coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                            angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                            x_dg, solved, n_dirs,
                                                            κ, σ, intΦ_eff, bc_fn_all, comm;
                                                            extra_source_all = nothing,
                                                            max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-combined over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total)")

    for (iel,d) in candidates
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
    received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0.0

        for (iel,d) in candidates
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
            end

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
                I_upwind = (i,j,k) -> 0.0

                if is_inflow[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            if haskey(x_dg,(parent_iel,d))
                                L1,L2 = ncf_L1L2[ncf_idx]
                                IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                                child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d,
                                                                                 IPc_nodes, IPp_nodes, L1, L2, ngl)
                                I_upwind = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                            end
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            children = [Int(mesh.cip[idx]) for idx in ncf_group]
                            if all(c -> haskey(x_dg,(c,d)), children)
                                entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                           for idx in ncf_group]
                                parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
                                I_upwind = (i,j,k) -> begin
                                    ni = face_flat_index(fid,ngl,i,j,k)
                                    p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                    parent_grid[p1,p2]
                                end
                            end
                        elseif haskey(neighbor_rank, (iel,fid))
                            kind = neighbor_kind[(iel,fid)]
                            if kind == :i_am_fine
                                remote_gel_f = neighbor_gel_group[(iel,fid)][1]
                                if all(haskey(received_sp,(remote_gel_f,d,spos)) for spos=1:ngl2)
                                    parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
                                    child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gel_f, parent_global_gips, ngl)
                                    I_upwind = (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]
                                end
                            elseif kind == :i_am_coarse
                                remote_gels_c = neighbor_gel_group[(iel,fid)]
                                if all(haskey(received_sp,(g,d,spos)) for g in remote_gels_c for spos=1:ngl2)
                                    parent_grid = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, received_sp,
                                        iel, fid, d, remote_gels_c, gelpos_to_gip, ngl)
                                    I_upwind = (i,j,k) -> begin
                                        ni=face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                                        parent_grid[p1,p2]
                                    end
                                end
                            elseif kind == :conforming
                                if haskey(angular_deps_cr, (iel,fid,d))
                                    remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                                    if length(remote_nodes) == 1
                                        rn = remote_nodes[1]
                                        if all(haskey(received_ang,(remote_gel,rn,fid,spos)) for spos=1:ngl2)
                                            rn_c=rn
					    fid_cc = fid
                                            I_upwind = (i,j,k) -> received_ang[(remote_gel,rn_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
                                        end
                                    else
                                        if all(haskey(received_ang,(remote_gel,rn,fid,spos)) for rn in remote_nodes for spos=1:ngl2)
                                            em_remote = remote_angular_mesh[remote_gel]
                                            src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                                            if src_e_ext != 0
                                                I_upwind = build_angular_ncf_upwind_unified_remote(
                                                    iel, fid, remote_gel, d, src_e_ext, mesh, em_remote, received_ang, ngl)
                                            end
                                        end
                                    end
                                else
                                    remote_gel_ord = neighbor_gel[(iel,fid)]
                                    if all(haskey(received_ang,(remote_gel_ord,d,fid,spos)) for spos=1:ngl2)
                                        gel_c=remote_gel_ord; d_c=d
                                        fid_cc = fid
					I_upwind = (i,j,k) -> received_ang[(gel_c,d_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
                                    end
                                end
                            end
                        else
                            bc = bc_fn_all[(iel,d)]
                            I_upwind = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                            end
                        end
                    else
                        key = (iel,d)
                        handled = false
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                handled = true
                                if length(dep_nodes) == 1
                                    dj = dep_nodes[1]; fid_c=fid
                                    if haskey(x_dg,(jel,dj))
                                        I_upwind = (i,j,k) -> begin
                                            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                            x_dg[(jel,dj)][lidx(corr[1],corr[2],corr[3],ngl)]
                                        end
                                    end
                                else
                                    if all(haskey(x_dg,(jel,dn)) for dn in dep_nodes)
                                        em_jel = mesh.extra_mesh[jel]
                                        src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                        I_upwind = build_angular_ncf_upwind_unified(
                                            iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    end
                                end
                            end
                        end
                        if !handled && haskey(x_dg,(jel,d))
                            d_c=d; fid_c=fid
                            I_upwind = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end

        received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-combined iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-combined did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

function test_manufactured_isotropic_combined_ncf_multirank(inputs, nparts, distribute, sem, partitioned_model;
                                                             spatial_target_iel::Int = 415,
                                                             angular_target_iel::Int = 442, angular_target_e_ext::Int = 4,
                                                             I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                             κ_val=0.001, σ_val=0.0005,
                                                             max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                         spatial_target_iel=spatial_target_iel,
                                         angular_target_iel=angular_target_iel,
                                         angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh; metrics = setup.metrics
    sem2 = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    angular_exch = build_ncf_angular_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
    println("Rank $rank: grouped $(mesh.nelem) elements into $(length(groups)) angular-mesh signatures")

    w_ext_all = Dict{Int, Vector{Float64}}(); W_all = Dict{Int, Float64}()
    coeff_all = Dict{Int, Vector{Float64}}(); C_all = Dict{Int, Matrix{Float64}}()
    Ωm_all = Dict{Int, Vector{NTuple{3,Float64}}}()
    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        W = sum(w_ext)
        w_ext_all[rep] = w_ext; W_all[rep] = W
        C_all[rep] = [σ_val*w_ext[m]/W for kk=1:em.extra_npoin, m=1:em.extra_npoin]
        Ωm_all[rep] = [(sin(em.extra_coords[1,m])*cos(em.extra_coords[2,m]),
                        sin(em.extra_coords[1,m])*sin(em.extra_coords[2,m]),
                        cos(em.extra_coords[1,m])) for m = 1:em.extra_npoin]
        coeff_all[rep] = σ_val .* w_ext ./ W
    end

    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]
        C_iel = C_all[rep]; Ω_m = Ωm_all[rep]
        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)
            q_vec = zeros(npoin_cg)
            for k=1:ngl,j=1:ngl,i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = sum(C_iel[dk,m]*I_manuf(x_,y_,z_, Ω_m[m]...) for m = 1:em.extra_npoin)
                q_vec[g] = Ω_dot_grad + κ_val*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
        end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        I_cg_all[(iel,d)] = zeros(ngl^3)
    end

    for it = 1:max_iter
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]
            rep = representative_of[iel]
            coeff = coeff_all[rep]
            for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_k_local = @view q_k[[mesh.connijk[iel,i,j,k] for k=1:ngl,j=1:ngl,i=1:ngl][:]]
                for m = 1:em.extra_npoin
                    coeff[m] == 0.0 && continue
                    q_k_local .+= coeff[m] .* I_cg_all[(iel,m)]
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end

        x_dg, _ = run_dg_sweep_multirank_combined_ncf!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.coord_cache, setup.gelpos_to_gip, spatial_perdir_exch,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh, angular_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm; extra_source_all=extra_source_all)

        max_delta_local = 0.0
        for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
            max_delta_local = max(max_delta_local, maximum(abs.(x_dg[(iel,d)] .- I_cg_all[(iel,d)])))
            I_cg_all[(iel,d)] = copy(x_dg[(iel,d)])
        end
        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        rank == 0 && println("  iter $it: max_delta=$max_delta_global")
        max_delta_global < tol && (rank==0 && println("Converged after $it iterations"); break)
    end

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
                computed = I_cg_all[(iel,d)][lidx(i,j,k,ngl)]
                max_err_local = max(max_err_local, abs(computed - I_true))
            end
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Multirank COMBINED NCF manufactured solution — max abs error: $max_err_global")
    return max_err_global
end

function test_dg_sweep_multirank_combined_ncf_exponential_instrumented(inputs, nparts, distribute, sem, partitioned_model;
                                                                        spatial_target_iel::Int = 415,
                                                                        angular_target_iel::Int = 442, angular_target_e_ext::Int = 4,
                                                                        κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_combined_ncf_mpi_test(inputs, nparts, distribute, sem, partitioned_model;
                                         spatial_target_iel=spatial_target_iel,
                                         angular_target_iel=angular_target_iel,
                                         angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh; metrics = setup.metrics
    sem2 = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    angular_exch = build_ncf_angular_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    # ── Per-(iel,d) BC — each direction's OWN (θ,ϕ)-derived (Ωx,Ωy,Ωz),
    # matching exactly how the sweep computes Ω internally for that (iel,d). ──
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            bc_fn_all[(iel,d)] = (x,y,z) -> I0*exp(-κ_val*(Ωxd*(x-r_bc[1])+Ωyd*(y-r_bc[2])+Ωzd*(z-r_bc[3])))
        end
    end

    # ── exact_fn takes (x,y,z,Ωxd,Ωyd,Ωzd) — the sweep passes EACH
    # (iel,d)'s OWN freshly-computed Ω here, so this must be a PURE
    # function of those six numbers, no other hidden state. ─────────────────
    exact_fn = (x,y,z,Ωxd,Ωyd,Ωzd) -> I0*exp(-κ_val*(Ωxd*(x-r_bc[1])+Ωyd*(y-r_bc[2])+Ωzd*(z-r_bc[3])))

    x_dg, _, solved = run_dg_sweep_multirank_combined_ncf_periodic_instrumented!(
        mesh, metrics, ψ, dψ, ω,
        setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence,
        Dict{Tuple{Int,Int},Dict{Tuple{Int,Int,Int},Tuple{Int,Int,Int}}}(),
        setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group,
        Dict{Tuple{Int,Int},Bool}(),
        Dict{Tuple{Int,Int},Symbol}(), Dict{Tuple{Int,Int},Vector{Int}}(), Dict{Tuple{Int,Int},Bool}(),
        setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
        setup.coord_cache, setup.gelpos_to_gip, spatial_perdir_exch,
        setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh, angular_exch,
        κ, σz, intΦ_eff, bc_fn_all, comm;
        exact_fn=exact_fn, err_threshold=1.0, block_periodic=false)

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            solved[iel][d] || continue
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωxd*(mesh.x[gip]-r_bc[1])+Ωyd*(mesh.y[gip]-r_bc[2])+Ωzd*(mesh.z[gip]-r_bc[3])))
                max_err_local = max(max_err_local, abs(x_dg[(iel,d)][lidx(i,j,k,ngl)] - exact))
            end
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Instrumented exponential (combined NCF, no periodicity) — max abs error: $max_err_global")
    return max_err_global
end

function debug_cross_rank_ordinary_lookup(mesh, neighbor_gel, received_ang, x_dg, iel::Int, fid::Int, d::Int, ngl::Int,
                                            i::Int, j::Int, k::Int)
    rank = MPI.Comm_rank(get_mpi_comm())
    remote_gel_ord = neighbor_gel[(iel,fid)]
    spos = face_flat_index(fid, ngl, i, j, k)
    val = get(received_ang, (remote_gel_ord, d, fid, spos), NaN)

    println("Rank $rank: DEBUG (iel=$iel,fid=$fid,d=$d,i=$i,j=$j,k=$k): " *
            "neighbor_gel resolves to remote_gel=$remote_gel_ord, spos=$spos, received_ang value=$val")

    # Also print EVERY entry in received_ang matching this remote_gel and
    # direction d, across ALL spos, to see the full picture for this
    # (gel,d) pair.
    matching = [(k2,v) for (k2,v) in received_ang if k2[1]==remote_gel_ord && k2[2]==d]
    println("Rank $rank: all received_ang entries for (remote_gel=$remote_gel_ord,d=$d): $(sort(matching, by=x->x[1][3]))")
end

"""
    exchange_ncf_angular_values_traced!(exch, x_dg, ngl; trace_gel, trace_d, trace_spos)

Same as exchange_ncf_angular_values!, but prints full detail for ANY
flat entry matching the trace target — on the SEND side (what local
node/value is being sent under this tag) and, separately, prints what
gets decoded on the RECEIVE side for the same tag.
"""
function exchange_ncf_angular_values_traced!(mesh, exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int;
                                              trace_gel::Int, trace_d::Int, trace_spos::Int)
    rank = MPI.Comm_rank(get_mpi_comm())
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]
        d   = exch.flat_d[m]
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        val = haskey(x_dg, (iel,d)) ? x_dg[(iel,d)][lidx(i,j,k,ngl)] : NaN
        send_vals[m] = val

        my_gel = Int(mesh.el2gel[iel])
        if my_gel == trace_gel && d == trace_d
            println("Rank $rank: [SEND-TRACE] flat_idx=$m (iel=$iel,gel=$my_gel,d=$d) " *
                    "local (i,j,k)=($i,$j,$k) -> sending val=$val")
        end
    end
    send_and_receive!(exch.val_cache, send_vals)

    received = Dict{Tuple{TInt,Int,Int,Int}, Float64}()
    for m = 1:length(exch.val_cache.combined_recv_data)
        val = exch.val_cache.combined_recv_data[m]
        isnan(val) && continue
        rgel  = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd    = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rspos = Int(round(Int, exch.spos_cache.combined_recv_data[m]))
        rfid  = Int(round(Int, exch.fid_cache.combined_recv_data[m])) 
	received[(rgel,rd,rfid,rspos)] = val
        if rgel == trace_gel && rd == trace_d && rspos == trace_spos
            println("Rank $rank: [RECV-TRACE] decoded (gel=$rgel,d=$rd,spos=$rspos) = $val")
        end
    end
    return received
end

function build_ncf_gid_value_exchange(mesh, face_list::Dict{Tuple{TInt,TInt},TInt}, n_dirs::Vector{Int}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]; flat_my_gel = TInt[]; flat_d = Int[]; flat_gid = Int[]
    flat_owning_iel = Int[]; flat_i = Int[]; flat_j = Int[]; flat_k = Int[]

    for ((iel,fid), target_rank) in face_list
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for d = 1:n_dirs[iel]
            for k in krange, j in jrange, i in irange
                gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
                push!(flat_targets, Int(target_rank)); push!(flat_my_gel, my_gel)
                push!(flat_d, d); push!(flat_gid, gid)
                push!(flat_owning_iel, iel); push!(flat_i, i); push!(flat_j, j); push!(flat_k, k)
            end
        end
    end

    n_flat = length(flat_targets)
    gel_cache = SendReceiveCache(comm, Float64.(flat_my_gel), flat_targets)
    d_cache   = SendReceiveCache(comm, Float64.(flat_d), flat_targets)
    gid_cache = SendReceiveCache(comm, Float64.(flat_gid), flat_targets)
    val_cache = SendReceiveCache(comm, zeros(Float64, n_flat), flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_my_gel))
    send_and_receive!(d_cache, Float64.(flat_d))
    send_and_receive!(gid_cache, Float64.(flat_gid))

    return (val_cache=val_cache, gel_cache=gel_cache, d_cache=d_cache, gid_cache=gid_cache,
            flat_owning_iel=flat_owning_iel, flat_d=flat_d,
            flat_i=flat_i, flat_j=flat_j, flat_k=flat_k, n_flat=n_flat)
end

function exchange_ncf_gid_values!(exch, x_dg::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]; d = exch.flat_d[m]
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        send_vals[m] = haskey(x_dg, (iel,d)) ? x_dg[(iel,d)][lidx(i,j,k,ngl)] : NaN
    end
    send_and_receive!(exch.val_cache, send_vals)

    received = Dict{Tuple{TInt,Int,Int}, Float64}()   # (gel, d, gid) -> val
    for m = 1:length(exch.val_cache.combined_recv_data)
        val = exch.val_cache.combined_recv_data[m]
        isnan(val) && continue
        rgel = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd   = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rgid = Int(round(Int, exch.gid_cache.combined_recv_data[m]))
        received[(rgel,rd,rgid)] = val   # duplicate sends of the SAME physical point carry
                                            # the SAME value → overwrite is always safe
    end
    return received
end

"""
    build_face_upwind(...)

Returns (available::Bool, closure::Union{Function,Nothing}) for ONE face
of ONE (iel,d). `available=false` means not ready yet — callers decide
what that means (block for the ordinary sweep, fall back to a default
for GS). `is_solved(jel,dj)` abstracts "has this (jel,dj) got a usable
value" so the SAME helper serves both the real sweep (checking a genuine
solved[] array) and GS (checking haskey(x_dg,...) on its candidate pool).
"""
#=function build_face_upwind(mesh, iel::Int, fid::Int, d::Int, ngl::Int,
                            neighbor_elem, face_correspondence,
                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                            neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                            periodic_face_correspondence,
                            angular_deps, angular_deps_cr, remote_angular_mesh,
                            coord_cache, gelpos_to_gip, spatial_perdir_received,
                            gid_received, periodic_gid_received,
                            is_solved::Function, x_dg, I_cg_all_prev, bc_fn_all)
    ngl2 = ngl*ngl

    if get(same_rank_is_periodic, (iel,fid), false)
        partner_lel = same_rank_partner[(iel,fid)][1]
        fid_c = fid
        if I_cg_all_prev !== nothing && haskey(I_cg_all_prev, (partner_lel,d))
            partner_vals = I_cg_all_prev[(partner_lel,d)]
            return true, (i,j,k) -> begin
                corr = periodic_face_correspondence[(iel,fid_c)][(i,j,k)]
                partner_vals[lidx(corr[1],corr[2],corr[3],ngl)]
            end
        else
            return true, (i,j,k) -> 0.0
        end
    end
    if get(neighbor_is_periodic, (iel,fid), false)
        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
        return true, (i,j,k) -> begin
            gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
            get(periodic_gid_received, (remote_gel_p, d, gid), 0.0)
        end
    end

    jel_raw = neighbor_elem[iel][fid]

    if jel_raw != zero(TInt)
        jel = Int(jel_raw)
        key = (iel,d)
        if haskey(angular_deps, key)
            for (dep_jel, dep_nodes) in angular_deps[key]
                dep_jel == jel || continue
                if length(dep_nodes) == 1
                    dj = dep_nodes[1]
                    is_solved(jel,dj) || return false, nothing
                    dj_c=dj; fid_c=fid
                    return true, (i,j,k) -> begin
                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                    end
                else
                    all(dn -> is_solved(jel,dn), dep_nodes) || return false, nothing
                    em_jel = mesh.extra_mesh[jel]
                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                    return true, build_angular_ncf_upwind_unified(iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                end
            end
        end
        is_solved(jel,d) || return false, nothing
        d_c=d; fid_c=fid
        return true, (i,j,k) -> begin
            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
            x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
        end
    end

    # jel_raw == 0
    if haskey(ncf_as_child, (iel,fid))
        ncf_idx = ncf_as_child[(iel,fid)]
        parent_iel = Int(mesh.pip[ncf_idx])
        is_solved(parent_iel,d) || return false, nothing
        L1,L2 = ncf_L1L2[ncf_idx]
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d, IPc_nodes, IPp_nodes, L1, L2, ngl)
        return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

    elseif haskey(ncf_as_parent, (iel,fid))
        ncf_group = ncf_as_parent[(iel,fid)]
        children = [Int(mesh.cip[idx]) for idx in ncf_group]
        all(c -> is_solved(c,d), children) || return false, nothing
        entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in ncf_group]
        parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
        return true, (i,j,k) -> begin
            ni = face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
            parent_grid[p1,p2]
        end

    elseif haskey(neighbor_rank, (iel,fid))
        kind = neighbor_kind[(iel,fid)]

        if kind == :i_am_fine
            #=remote_gel_f = neighbor_gel_group[(iel,fid)][1]
            fully = all(haskey(spatial_perdir_received, (remote_gel_f, d, spos)) for spos = 1:ngl2)
            fully || return false, nothing
            parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
            child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, spatial_perdir_received,
                iel, fid, d, remote_gel_f, parent_global_gips, ngl)
            return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]=#
	    remote_gel_f = neighbor_gel_group[(iel,fid)][1]
    fully = all(haskey(spatial_perdir_received, (remote_gel_f, d, spos)) for spos = 1:ngl2)
    fully || return false, nothing
    parent_global_gips = [gelpos_to_gip[(remote_gel_f,pos)] for pos=1:ngl2]
    if iel == 201 && fid == 4 && d == 9   # TEMPORARY debug trigger — matches the failing case
        child_vals = build_dg_ncf_parent_to_child_remote_perdir_debug(mesh, coord_cache, spatial_perdir_received,
            iel, fid, d, remote_gel_f, parent_global_gips, ngl, get_mpi_comm(); target_i=1, target_j=1, target_k=1)
    else
        child_vals = build_dg_ncf_parent_to_child_remote_perdir(mesh, coord_cache, spatial_perdir_received,
            iel, fid, d, remote_gel_f, parent_global_gips, ngl)
    end
    return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

        elseif kind == :i_am_coarse
            remote_gels_c = neighbor_gel_group[(iel,fid)]
            fully = all(haskey(spatial_perdir_received,(g,d,spos)) for g in remote_gels_c for spos=1:ngl2)
            fully || return false, nothing
            parent_grid = build_dg_ncf_child_to_parent_remote_perdir(mesh, coord_cache, spatial_perdir_received,
                iel, fid, d, remote_gels_c, gelpos_to_gip, ngl)
            return true, (i,j,k) -> begin
                ni=face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                parent_grid[p1,p2]
            end

        elseif kind == :conforming
            if haskey(angular_deps_cr, (iel,fid,d))
                remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                if length(remote_nodes) == 1
                    rn = remote_nodes[1]
                    return true, (i,j,k) -> begin
                        gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
                        get(gid_received, (remote_gel, rn, gid), NaN)
                    end
                else
                    em_remote = remote_angular_mesh[remote_gel]
                    src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                    src_e_ext == 0 && return false, nothing
                    return true, build_angular_ncf_upwind_unified_remote(mesh, iel, remote_gel, d, src_e_ext, em_remote, gid_received)
                end
            else
                remote_gel_ord = neighbor_gel[(iel,fid)]
                return true, (i,j,k) -> begin
                    gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
                    get(gid_received, (remote_gel_ord, d, gid), NaN)
                end
            end
        end
    else
        bc = bc_fn_all[(iel,d)]
        return true, (i,j,k) -> begin
            gip = mesh.connijk[iel,i,j,k]
            bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        end
    end
end=#

"""
    build_face_upwind(...)

Returns (available::Bool, closure::Union{Function,Nothing}) for ONE face
of ONE (iel,d). Spatial-NCF cross-rank branches now use the fully
coordinate-driven, gid-based reconstruction (no transmitted position
index at all) via remote_gid_coords + gid_received.
"""
function build_face_upwind(mesh, iel::Int, fid::Int, d::Int, ngl::Int,
                            neighbor_elem, face_correspondence,
                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                            neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                            periodic_face_correspondence,
                            angular_deps, angular_deps_cr, remote_angular_mesh,
                            remote_gid_coords, gid_received, periodic_gid_received,
                            is_solved::Function, x_dg, I_cg_all_prev, bc_fn_all)
    ngl2 = ngl*ngl

    if get(same_rank_is_periodic, (iel,fid), false)
        partner_lel = same_rank_partner[(iel,fid)][1]
        fid_c = fid
        if I_cg_all_prev !== nothing && haskey(I_cg_all_prev, (partner_lel,d))
            partner_vals = I_cg_all_prev[(partner_lel,d)]
            return true, (i,j,k) -> begin
                corr = periodic_face_correspondence[(iel,fid_c)][(i,j,k)]
                partner_vals[lidx(corr[1],corr[2],corr[3],ngl)]
            end
        else
            return true, (i,j,k) -> 0.0
        end
    end
    if get(neighbor_is_periodic, (iel,fid), false)
        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
        return true, (i,j,k) -> begin
            gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
            get(periodic_gid_received, (remote_gel_p, d, gid), 0.0)
        end
    end

    jel_raw = neighbor_elem[iel][fid]

    if jel_raw != zero(TInt)
        jel = Int(jel_raw)
        key = (iel,d)
        if haskey(angular_deps, key)
            for (dep_jel, dep_nodes) in angular_deps[key]
                dep_jel == jel || continue
                if length(dep_nodes) == 1
                    dj = dep_nodes[1]
                    is_solved(jel,dj) || return false, nothing
                    dj_c=dj; fid_c=fid
                    return true, (i,j,k) -> begin
                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                    end
                else
                    all(dn -> is_solved(jel,dn), dep_nodes) || return false, nothing
                    em_jel = mesh.extra_mesh[jel]
                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                    return true, build_angular_ncf_upwind_unified(iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                end
            end
        end
        is_solved(jel,d) || return false, nothing
        d_c=d; fid_c=fid
        return true, (i,j,k) -> begin
            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
            x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
        end
    end

    if haskey(ncf_as_child, (iel,fid))
        ncf_idx = ncf_as_child[(iel,fid)]
        parent_iel = Int(mesh.pip[ncf_idx])
        is_solved(parent_iel,d) || return false, nothing
        L1,L2 = ncf_L1L2[ncf_idx]
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d, IPc_nodes, IPp_nodes, L1, L2, ngl)
        return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

    elseif haskey(ncf_as_parent, (iel,fid))
        ncf_group = ncf_as_parent[(iel,fid)]
        children = [Int(mesh.cip[idx]) for idx in ncf_group]
        all(c -> is_solved(c,d), children) || return false, nothing
        entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in ncf_group]
        parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
        return true, (i,j,k) -> begin
            ni = face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
            parent_grid[p1,p2]
        end

    elseif haskey(neighbor_rank, (iel,fid))
        kind = neighbor_kind[(iel,fid)]

        if kind == :i_am_fine
            remote_gel_f = neighbor_gel_group[(iel,fid)][1]
            haskey(remote_gid_coords, remote_gel_f) || return false, nothing
            child_vals = build_dg_ncf_parent_to_child_remote_by_coords(mesh, iel, fid, d, remote_gel_f,
                remote_gid_coords, gid_received, ngl)
            any(isnan, child_vals) && return false, nothing
            return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

        elseif kind == :i_am_coarse
            remote_gels_c = neighbor_gel_group[(iel,fid)]
            all(g -> haskey(remote_gid_coords, g), remote_gels_c) || return false, nothing
            parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, remote_gels_c,
                remote_gid_coords, gid_received, ngl)
            any(isnan, parent_grid) && return false, nothing
            return true, (i,j,k) -> begin
                ni=face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                parent_grid[p1,p2]
            end

        elseif kind == :conforming
            if haskey(angular_deps_cr, (iel,fid,d))
                remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
                if length(remote_nodes) == 1
                    rn = remote_nodes[1]
                    return true, (i,j,k) -> begin
                        gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
                        get(gid_received, (remote_gel, rn, gid), NaN)
                    end
                else
                    em_remote = remote_angular_mesh[remote_gel]
                    src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
                    src_e_ext == 0 && return false, nothing
                    return true, build_angular_ncf_upwind_unified_remote(mesh, iel, remote_gel, d, src_e_ext, em_remote, gid_received)
                end
            else
                remote_gel_ord = neighbor_gel[(iel,fid)]
                return true, (i,j,k) -> begin
                    gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
                    get(gid_received, (remote_gel_ord, d, gid), NaN)
                end
            end
        end
    else
        bc = bc_fn_all[(iel,d)]
        return true, (i,j,k) -> begin
            gip = mesh.connijk[iel,i,j,k]
            bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        end
    end
end

#=function run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          periodic_face_correspondence,
                                          neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                          same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                          ncf_as_child, ncf_L1L2, ncf_as_parent,
                                          coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                          angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch,
                                          periodic_gid_exch,
                                          κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                          extra_source_all = nothing,
                                          I_cg_all_prev = nothing,
                                          max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    periodic_gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    spatial_perdir_received = Dict{Tuple{TInt,Int,Int}, Float64}()

    is_solved_sweep = (jel,dj) -> solved[jel][dj]

    while true
        iter_round += 1
        n_newly_solved = 0

        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue

                    avail, clo = build_face_upwind(mesh, iel, fid, d, ngl,
                        neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                        neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                        periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                        coord_cache, gelpos_to_gip, spatial_perdir_received,
                        gid_received, periodic_gid_received,
                        is_solved_sweep, x_dg, I_cg_all_prev, bc_fn_all)

                    if !avail
                        ready = false; break
                    end

                    # ── Trial-evaluate at ONE face node to catch NaN
                    # (the lazy-availability signal for gid-based cross-rank
                    # branches) BEFORE trusting this closure. ────────────────
                    irange,jrange,krange = face_ip_range(fid, ngl)
                    trial = clo(first(irange), first(jrange), first(krange))
                    if isnan(trial)
                        ready = false; break
                    end
                    upwind_fns[fid] = clo
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
        spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        periodic_gid_received = I_cg_all_prev === nothing ? periodic_gid_received :
                                 exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        if n_newly_solved_global == 0
            if just_resolved_stall
                error("run_dg_sweep_multirank_unified!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence, periodic_face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                same_rank_kind, same_rank_partner, same_rank_is_periodic,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                coord_cache, gelpos_to_gip, spatial_perdir_exch,
                angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                x_dg, solved, n_dirs, κ, σ, intΦ_eff, bc_fn_all, comm;
                extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev,
                max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs
end=#

function run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          periodic_face_correspondence,
                                          neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                          same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                          ncf_as_child, ncf_L1L2, ncf_as_parent,
                                          remote_gid_coords, spatial_perdir_exch,
                                          angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                                          κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                          extra_source_all = nothing,
                                          I_cg_all_prev = nothing,
                                          max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    periodic_gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    spatial_perdir_received = Dict{Tuple{TInt,Int,Int}, Float64}()

    is_solved_sweep = (jel,dj) -> solved[jel][dj]

    while true
        iter_round += 1
        n_newly_solved = 0

        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue

                    avail, clo = build_face_upwind(mesh, iel, fid, d, ngl,
                        neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                        neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                        periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                        remote_gid_coords, gid_received, periodic_gid_received,
                        is_solved_sweep, x_dg, I_cg_all_prev, bc_fn_all)

                    if !avail
                        ready = false; break
                    end
                    irange,jrange,krange = face_ip_range(fid, ngl)
                    trial = clo(first(irange), first(jrange), first(krange))
                    if isnan(trial)
                        ready = false; break
                    end
                    upwind_fns[fid] = clo
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
        spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        periodic_gid_received = I_cg_all_prev === nothing ? periodic_gid_received :
                                 exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        if n_newly_solved_global == 0
            if just_resolved_stall
                error("run_dg_sweep_multirank_unified!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence, periodic_face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                same_rank_kind, same_rank_partner, same_rank_is_periodic,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                remote_gid_coords, spatial_perdir_exch,
                angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                x_dg, solved, n_dirs, κ, σ, intΦ_eff, bc_fn_all, comm;
                extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev,
                max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs
end

function resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                                                           neighbor_elem, face_correspondence, periodic_face_correspondence,
                                                           neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                                           same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                           remote_gid_coords, spatial_perdir_exch,
                                                           angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                                                           x_dg, solved, n_dirs,
                                                           κ, σ, intΦ_eff, bc_fn_all, comm;
                                                           extra_source_all = nothing,
                                                           I_cg_all_prev = nothing,
                                                           max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-unified over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total)")

    for (iel,d) in candidates
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    is_solved_gs = (jel,dj) -> haskey(x_dg,(jel,dj))
    gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
    spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
    periodic_gid_received = I_cg_all_prev === nothing ? Dict{Tuple{TInt,Int,Int},Float64}() :
                              exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0.0

        for (iel,d) in candidates
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
            end

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)

                I_upwind = (i,j,k) -> 0.0
                if is_inflow[fid]
                    avail, clo = build_face_upwind(mesh, iel, fid, d, ngl,
                        neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                        neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                        periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                        remote_gid_coords, gid_received, periodic_gid_received,
                        is_solved_gs, x_dg, I_cg_all_prev, bc_fn_all)

                    if avail
                        I_upwind = (i,j,k) -> begin
                            v = clo(i,j,k)
                            isnan(v) ? 0.0 : v
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end

        gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
        spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        periodic_gid_received = I_cg_all_prev === nothing ? periodic_gid_received :
                                  exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-unified iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-unified did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

#=function resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                                                           neighbor_elem, face_correspondence, periodic_face_correspondence,
                                                           neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                                           same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                           coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                           angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                                                           x_dg, solved, n_dirs,
                                                           κ, σ, intΦ_eff, bc_fn_all, comm;
                                                           extra_source_all = nothing,
                                                           I_cg_all_prev = nothing,
                                                           max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-unified over $(length(candidates)) local candidates " *
            "($(n_candidates_total) total)")

    for (iel,d) in candidates
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    is_solved_gs = (jel,dj) -> haskey(x_dg,(jel,dj))
    gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
    spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
    periodic_gid_received = I_cg_all_prev === nothing ? Dict{Tuple{TInt,Int,Int},Float64}() :
                              exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

    for gs_iter = 1:max_gs_iter
        max_delta_local = 0.0

        for (iel,d) in candidates
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
            end

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)

                I_upwind = (i,j,k) -> 0.0
                if is_inflow[fid]
                    avail, clo = build_face_upwind(mesh, iel, fid, d, ngl,
                        neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                        neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                        periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                        coord_cache, gelpos_to_gip, spatial_perdir_received,
                        gid_received, periodic_gid_received,
                        is_solved_gs, x_dg, I_cg_all_prev, bc_fn_all)

                    if avail
                        # Wrap: fall back to 0.0 per-point if a specific
                        # sample happens to be NaN (lazy-availability
                        # branches) — GS NEVER blocks, always proceeds. ────────
                        I_upwind = (i,j,k) -> begin
                            v = clo(i,j,k)
                            isnan(v) ? 0.0 : v
                        end
                    end
                    # avail==false (e.g. a genuinely same-rank NOT-YET-solved
                    # same-rank dependency outside the candidate pool) ->
                    # falls through to the default 0.0 closure already set.
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end

        gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
        spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        periodic_gid_received = I_cg_all_prev === nothing ? periodic_gid_received :
                                  exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-unified iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-unified did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end=#



function test_exchange_ncf_gid_values(inputs, nparts, distribute, sem, partitioned_model;
                                       spatial_target_iel::Int = 315,
                                       angular_target_iel::Int = 242, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_combined_ncf_mpi_test_periodic(inputs, nparts, distribute, sem, partitioned_model;
                                                  spatial_target_iel=spatial_target_iel,
                                                  angular_target_iel=angular_target_iel,
                                                  angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    n_periodic_cross = count(identity, values(setup.neighbor_is_periodic))
    println("Rank $rank: cross-rank periodic faces: $n_periodic_cross " *
            "(0 is expected/fine on a fully rigid mesh)")

    # ── ORDINARY + angular-NCF exchange (neighbor_rank as-is; periodic
    # faces, if any, are ALSO in here, but that's fine — they're just
    # extra entries this specific test doesn't care about). ─────────────────
    gid_exch = build_ncf_gid_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    println("Rank $rank: ordinary/angular gid exchange has $(gid_exch.n_flat) flat entries")

    # ── PERIODIC-ONLY exchange, filtered to neighbor_is_periodic — EMPTY
    # tables/dicts if the mesh has no periodic faces at all. ─────────────────
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in setup.neighbor_rank
                                                        if get(setup.neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)
    println("Rank $rank: periodic-only gid exchange has $(periodic_gid_exch.n_flat) flat entries")

    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        my_gel = Int(mesh.el2gel[iel])
        for d = 1:n_dirs[iel]
            x_dg[(iel,d)] = fill(Float64(my_gel*1e6 + d), ngl^3)
        end
    end

    function check(exch, label)
        received = exchange_ncf_gid_values!(exch, x_dg, ngl)
        n_checked = 0; n_mismatch = 0
        for ((rgel,rd,rgid), val) in received
            n_checked += 1
            decoded_gel = Int(floor(val / 1e6))
            decoded_d = Int(round(val - decoded_gel*1e6))
            if decoded_gel != rgel || decoded_d != rd
                n_mismatch += 1
                n_mismatch <= 5 && println("Rank $rank: ✗ [$label] MISMATCH (rgel=$rgel,rd=$rd,rgid=$rgid): decoded=($decoded_gel,$decoded_d)")
            end
        end
        println("Rank $rank: [$label] checked $n_checked entries, $n_mismatch mismatches  " *
                (n_mismatch==0 ? "✓ PASS" : "✗ FAIL"))
        return received
    end

    received_ordinary = check(gid_exch, "ORDINARY+ANGULAR")
    received_periodic = check(periodic_gid_exch, "PERIODIC-ONLY")

    return received_ordinary, received_periodic
end

function test_dg_sweep_multirank_unified_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                       spatial_target_iel::Int = 315,
                                                       angular_target_iel::Int = 242, angular_target_e_ext::Int = 4,
                                                       κ_val::Float64 = 0.0001, I0::Float64 = 4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    setup = setup_combined_ncf_mpi_test_periodic(inputs, nparts, distribute, sem, partitioned_model;
                                                  spatial_target_iel=spatial_target_iel,
                                                  angular_target_iel=angular_target_iel,
                                                  angular_target_e_ext=angular_target_e_ext)
    mesh = setup.mesh; metrics = setup.metrics
    sem2 = setup.sem
    ψ, dψ, ω = sem2.basis.ψ, sem2.basis.dψ, sem2.ω
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    gid_exch = build_ncf_gid_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)

    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in setup.neighbor_rank
                                                        if get(setup.neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            bc_fn_all[(iel,d)] = (x,y,z) -> I0*exp(-κ_val*(Ωxd*(x-r_bc[1])+Ωyd*(y-r_bc[2])+Ωzd*(z-r_bc[3])))
        end
    end

    x_dg, _ = run_dg_sweep_multirank_unified_instrumented!(mesh, metrics, ψ, dψ, ω,
        setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
        setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
        setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
        setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
        setup.coord_cache, setup.gelpos_to_gip, spatial_perdir_exch,
        setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh, gid_exch, periodic_gid_exch,
        κ, σz, intΦ_eff, bc_fn_all, comm;
        I_cg_all_prev=nothing, trace_iel=1)   # no periodicity being exercised in THIS pure-transport test

    max_err_local = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωxd = sin(θ)*cos(ϕ); Ωyd = sin(θ)*sin(ϕ); Ωzd = cos(θ)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωxd*(mesh.x[gip]-r_bc[1])+Ωyd*(mesh.y[gip]-r_bc[2])+Ωzd*(mesh.z[gip]-r_bc[3])))
                max_err_local = max(max_err_local, abs(x_dg[(iel,d)][lidx(i,j,k,ngl)] - exact))
            end
        end
    end
    max_err_global = MPI.Allreduce(max_err_local, MPI.MAX, comm)
    rank == 0 && println("Unified sweep exponential (rigid mesh) — max abs error: $max_err_global")
    return max_err_global
end

function run_dg_sweep_multirank_unified_instrumented!(mesh, metrics, ψ, dψ, ω,
                                                        neighbor_elem, neighbor_fid, face_correspondence,
                                                        periodic_face_correspondence,
                                                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                                        same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                        ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                        coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                        angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                                                        κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                                        extra_source_all = nothing,
                                                        I_cg_all_prev = nothing,
                                                        exact_fn = nothing, err_threshold::Float64 = 1.0,
                                                        trace_iel::Int = 0,
                                                        max_rounds::Int = 500, max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    just_resolved_stall = false
    gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    periodic_gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    spatial_perdir_received = Dict{Tuple{TInt,Int,Int}, Float64}()

    is_solved_sweep = (jel,dj) -> solved[jel][dj]

    while true
        iter_round += 1
        n_newly_solved = 0

        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                upwind_fns = Dict{Int,Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue

                    avail, clo = build_face_upwind(mesh, iel, fid, d, ngl,
                        neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                        neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                        periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                        coord_cache, gelpos_to_gip, spatial_perdir_received,
                        gid_received, periodic_gid_received,
                        is_solved_sweep, x_dg, I_cg_all_prev, bc_fn_all)

                    if !avail
                        ready = false; break
                    end
                    irange,jrange,krange = face_ip_range(fid, ngl)
                    trial = clo(first(irange), first(jrange), first(krange))
                    if isnan(trial)
                        ready = false; break
                    end
                    upwind_fns[fid] = clo
                end
                ready || continue

                if exact_fn !== nothing
                    for fid = 1:6
                        haskey(upwind_fns,fid) || continue
                        irange,jrange,krange = face_ip_range(fid, ngl)
                        for (i,j,k) in Iterators.product(irange,jrange,krange)
                            val = upwind_fns[fid](i,j,k)
                            gip = mesh.connijk[iel,i,j,k]
                            ex = exact_fn(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                            if abs(val-ex) > err_threshold
                                println("Rank $rank ⚠⚠⚠ LARGE I_upwind ERROR: (iel=$iel,d=$d,fid=$fid,i=$i,j=$j,k=$k) val=$val exact=$ex")
                            end
                        end
                    end
                end

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns,fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1

                if exact_fn !== nothing
                    for k=1:ngl,j=1:ngl,i=1:ngl
                        gip = mesh.connijk[iel,i,j,k]
                        ex = exact_fn(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                        val = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                        if abs(val-ex) > err_threshold
                            println("Rank $rank ⚠⚠⚠ LARGE SOLVE ERROR: (iel=$iel,d=$d,i=$i,j=$j,k=$k) val=$val exact=$ex")
                        end
                    end
                end
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
        spatial_perdir_received = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)
        periodic_gid_received = I_cg_all_prev === nothing ? periodic_gid_received :
                                 exchange_ncf_gid_values!(periodic_gid_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

	println("Rank $rank: round $iter_round, solved_this_round=$n_newly_solved, " *
        "total_solved=$n_solved_total/$total_pairs")

        if n_newly_solved_global == 0
	    debug_dump_stuck_elements(mesh, solved, n_dirs, comm)
	    debug_check_spatial_ncf_crossrank_identity(mesh, neighbor_kind, neighbor_gel_group, spatial_perdir_received, ngl, 201, 4, 9, comm)
	    if trace_iel != 0 && trace_iel <= mesh.nelem


                for d = 1:n_dirs[trace_iel]
                    solved[trace_iel][d] && continue
                    em = mesh.extra_mesh[trace_iel]
                    θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                    classify_element_faces!(is_inflow, is_active, trace_iel, Ωx, Ωy, Ωz, ngl,
                                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                             metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    println("Rank $rank: [TRACE going into GS] (iel=$trace_iel,d=$d) UNSOLVED")
                    for fid = 1:6
                        (is_inflow[fid] && is_active[fid]) || continue
                        avail, clo = build_face_upwind(mesh, trace_iel, fid, d, ngl,
                            neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                            neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                            periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                            coord_cache, gelpos_to_gip, spatial_perdir_received,
                            gid_received, periodic_gid_received,
                            is_solved_sweep, x_dg, I_cg_all_prev, bc_fn_all)
                        jel_raw = neighbor_elem[trace_iel][fid]

                        sample_str = "N/A"
                        if avail
                            irange,jrange,krange = face_ip_range(fid, ngl)
                            s = clo(first(irange), first(jrange), first(krange))
                            sample_str = isnan(s) ? "NaN (cross-rank data not yet arrived)" : string(s)
                        end
                        println("    fid=$fid: jel_raw=$jel_raw avail=$avail sample=$sample_str")

                        if jel_raw != zero(TInt)
                            jel = Int(jel_raw)
                            key = (trace_iel, d)
                            found_dep = false
                            if haskey(angular_deps, key)
                                for (dep_jel, dep_nodes) in angular_deps[key]
                                    dep_jel == jel || continue
                                    found_dep = true
                                    println("      -> angular_deps entry on jel=$jel, needs directions=$dep_nodes")
                                    for dn in dep_nodes
                                        println("         jel=$jel, dn=$dn: solved=$(solved[jel][dn])")
                                    end
                                end
                            end
                            if !found_dep
                                println("      -> ORDINARY same-index fallback: jel=$jel needs its OWN d=$d: solved=$(solved[jel][d])")
                            end
                        end
                    end
                end
            end
            if just_resolved_stall
                error("run_dg_sweep_multirank_unified_instrumented!: no progress after GS resolution")
            end
            resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence, periodic_face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                same_rank_kind, same_rank_partner, same_rank_is_periodic,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                coord_cache, gelpos_to_gip, spatial_perdir_exch,
                angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                x_dg, solved, n_dirs, κ, σ, intΦ_eff, bc_fn_all, comm;
                extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev,
                max_gs_iter=max_gs_iter, gs_tol=gs_tol)
            just_resolved_stall = true
            n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs, solved
end

function debug_check_crossrank_angular_lookup(mesh, angular_deps_cr, gid_received, neighbor_gel,
                                                iel::Int, fid::Int, d::Int, ngl::Int, comm)
    rank = MPI.Comm_rank(comm)

    println("Rank $rank: === checking (iel=$iel,fid=$fid,d=$d) ===")

    has_cr = haskey(angular_deps_cr, (iel,fid,d))
    println("Rank $rank: haskey(angular_deps_cr,($iel,$fid,$d)) = $has_cr")

    remote_gel_to_check = nothing
    remote_nodes_to_check = nothing

    if has_cr
        remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
        println("Rank $rank: angular_deps_cr entry -> remote_gel=$remote_gel, remote_nodes=$remote_nodes " *
                "(length=$(length(remote_nodes)))")
        remote_gel_to_check = remote_gel
        remote_nodes_to_check = remote_nodes
    else
        haskey(neighbor_gel, (iel,fid)) || (println("Rank $rank: NO angular_deps_cr AND NO neighbor_gel entry — this face has no cross-rank mapping at all!"); return)
        remote_gel_ord = neighbor_gel[(iel,fid)]
        println("Rank $rank: falls to ORDINARY branch -> neighbor_gel=$remote_gel_ord (expects remote's SAME direction d=$d)")
        remote_gel_to_check = remote_gel_ord
        remote_nodes_to_check = [d]
    end

    # ── Check what gid_received ACTUALLY has for this remote_gel, across
    # ALL directions and gids — a full inventory, not just the one we expect. ──
    matching_gel_entries = [(k,v) for (k,v) in gid_received if k[1] == remote_gel_to_check]
    println("Rank $rank: gid_received has $(length(matching_gel_entries)) TOTAL entries for remote_gel=$remote_gel_to_check " *
            "(across ALL directions/gids)")

    for rn in remote_nodes_to_check
        matching_dir = [(k,v) for (k,v) in gid_received if k[1]==remote_gel_to_check && k[2]==rn]
        println("Rank $rank:   remote direction $rn: $(length(matching_dir)) entries in gid_received")
        if length(matching_dir) > 0 && length(matching_dir) <= 5
            for (k,v) in matching_dir
                println("      $k -> $v")
            end
        end
    end

    # ── ALSO check: does iel's OWN face (fid) even have nodes whose gid
    # we'd be looking up? Show what gids THIS element's face actually has. ────
    irange,jrange,krange = face_ip_range(fid, ngl)
    my_gids = [Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]]) for k in krange for j in jrange for i in irange]
    println("Rank $rank: MY face's own gids (what we'll search gid_received[($remote_gel_to_check,*,gid)] for): " *
            "$(my_gids[1:min(5,length(my_gids))])... ($(length(my_gids)) total)")

    # ── Direct overlap check: of MY gids, how many actually appear as the
    # THIRD element of ANY key for this remote_gel in gid_received? ──────────
    received_gids_for_remote = Set(k[3] for (k,v) in gid_received if k[1]==remote_gel_to_check)
    overlap = intersect(Set(my_gids), received_gids_for_remote)
    println("Rank $rank: overlap between MY face's gids and gids ACTUALLY RECEIVED for remote_gel=$remote_gel_to_check: " *
            "$(length(overlap))/$(length(my_gids))")
    length(overlap) < length(my_gids) && println("Rank $rank: ⚠ MISSING gids: $(setdiff(Set(my_gids), received_gids_for_remote))")
end

function debug_dump_stuck_elements(mesh, solved, n_dirs, comm; max_show::Int=10)
    rank = MPI.Comm_rank(comm)
    n_shown = 0
    stuck_gels = Set{Int}()
    for iel = 1:mesh.nelem
        for d = 1:n_dirs[iel]
            solved[iel][d] && continue
            push!(stuck_gels, Int(mesh.el2gel[iel]))
        end
    end
    println("Rank $rank: $(length(stuck_gels)) DISTINCT stuck elements (out of $(count(x->!all(x),solved)) with any stuck direction)")
    for g in collect(stuck_gels)[1:min(max_show,length(stuck_gels))]
        println("Rank $rank:   stuck gel=$g")
    end
    return stuck_gels
end

function debug_check_ordinary_crossrank_identity(mesh, neighbor_gel, angular_deps_cr, gid_received,
                                                    iel::Int, fid::Int, d::Int, ngl::Int, comm)
    rank = MPI.Comm_rank(comm)
    println("Rank $rank: === checking (iel=$iel,fid=$fid,d=$d) ===")

    has_cr = haskey(angular_deps_cr, (iel,fid,d))
    println("Rank $rank: angular_deps_cr entry exists: $has_cr")
    if has_cr
        remote_gel, remote_nodes = angular_deps_cr[(iel,fid,d)]
        println("Rank $rank: angular_deps_cr -> remote_gel=$remote_gel, remote_nodes=$remote_nodes")
    end

    has_ord = haskey(neighbor_gel, (iel,fid))
    println("Rank $rank: neighbor_gel entry exists: $has_ord")
    if has_ord
        remote_gel_ord = neighbor_gel[(iel,fid)]
        println("Rank $rank: neighbor_gel -> remote_gel=$remote_gel_ord")
        irange,jrange,krange = face_ip_range(fid, ngl)
        gid1 = Int(mesh.ip2gip[mesh.connijk[iel,first(irange),first(jrange),first(krange)]])
        val = get(gid_received, (remote_gel_ord, d, gid1), NaN)
        println("Rank $rank: gid_received[($remote_gel_ord,$d,$gid1)] = $val")
    end

    # ── Also print MY OWN element's position, to sanity-check whether
    # gel=201 sits at a multi-rank corner at all. ────────────────────────────
    println("Rank $rank: my (iel=$iel) corner coords: " *
            "($(mesh.x[mesh.connijk[iel,1,1,1]]), $(mesh.y[mesh.connijk[iel,1,1,1]]), $(mesh.z[mesh.connijk[iel,1,1,1]])) " *
            "to ($(mesh.x[mesh.connijk[iel,ngl,ngl,ngl]]), $(mesh.y[mesh.connijk[iel,ngl,ngl,ngl]]), $(mesh.z[mesh.connijk[iel,ngl,ngl,ngl]]))")
end

function debug_check_spatial_ncf_crossrank_identity(mesh, neighbor_kind, neighbor_gel_group,
                                                      spatial_perdir_received, ngl::Int,
                                                      iel::Int, fid::Int, d::Int, comm)
    rank = MPI.Comm_rank(comm)
    println("Rank $rank: === checking SPATIAL-NCF identity for (iel=$iel,fid=$fid,d=$d) ===")

    haskey(neighbor_kind, (iel,fid)) || (println("Rank $rank: NO neighbor_kind entry at all for this face"); return)
    kind = neighbor_kind[(iel,fid)]
    println("Rank $rank: neighbor_kind[($iel,$fid)] = $kind")

    kind in (:i_am_fine, :i_am_coarse) || (println("Rank $rank: kind is $kind, not spatial-NCF — nothing to trace here"); return)

    remote_gels = neighbor_gel_group[(iel,fid)]
    println("Rank $rank: neighbor_gel_group -> $remote_gels")

    ngl2 = ngl*ngl
    for rg in remote_gels
        n_entries = count(k -> k[1]==rg && k[2]==d, keys(spatial_perdir_received))
        println("Rank $rank:   remote_gel=$rg, direction d=$d: $n_entries/$ngl2 spos entries present in spatial_perdir_received")
        vals = [(k[3],v) for (k,v) in spatial_perdir_received if k[1]==rg && k[2]==d]
        sorted_vals = sort(vals, by=x->x[1])
        if length(sorted_vals) <= 30
            println("Rank $rank:     entries: $sorted_vals")
        else
            println("Rank $rank:     first 10: $(sorted_vals[1:10])")
            println("Rank $rank:     last 10: $(sorted_vals[end-9:end])")
        end

        # ── Check for DUPLICATE-LOOKING position collisions: are there
        # spos values 1..ngl2 that repeat suspiciously, or gaps? ──────────────
        spos_list = [v[1] for v in sorted_vals]
        expected = collect(1:ngl2)
        missing_spos = setdiff(expected, spos_list)
        extra_spos = setdiff(spos_list, expected)
        !isempty(missing_spos) && println("Rank $rank:     ⚠ MISSING spos: $missing_spos")
        !isempty(extra_spos) && println("Rank $rank:     ⚠ UNEXPECTED spos (out of 1:$ngl2 range): $extra_spos")
    end
end

function build_dg_ncf_parent_to_child_remote_perdir_debug(mesh, coord_cache::Dict{Int,NTuple{3,Float64}},
                                                            received::Dict{Tuple{TInt,Int,Int},Float64},
                                                            child_iel::Int, fid::Int, d::Int,
                                                            remote_gel::TInt, parent_global_gips::Vector{Int}, ngl::Int,
                                                            comm; target_i::Int=1, target_j::Int=1, target_k::Int=1)
    rank = MPI.Comm_rank(comm)
    IPc_nodes = Int[]
    IPc_ijk = Tuple{Int,Int,Int}[]
    irange, jrange, krange = face_ip_range(fid, ngl)
    for k in krange, j in jrange, i in irange
        push!(IPc_nodes, mesh.connijk[child_iel,i,j,k])
        push!(IPc_ijk, (i,j,k))
    end

    target_node_idx = findfirst(==( (target_i,target_j,target_k) ), IPc_ijk)
    println("Rank $rank: [NCF-DEBUG] target (i=$target_i,j=$target_j,k=$target_k) -> node_idx=$target_node_idx " *
            "(out of $(length(IPc_nodes)) face nodes), IPc_nodes[target]=$(target_node_idx===nothing ? "N/A" : IPc_nodes[target_node_idx])")

    L1, L2 = build_spatial_interpolation_matrices_with_coord_cache(
        fid, IPc_nodes, parent_global_gips, mesh, coord_cache, ngl)

    println("Rank $rank: [NCF-DEBUG] parent_global_gips = $parent_global_gips")
    println("Rank $rank: [NCF-DEBUG] L1 size=$(size(L1)), L2 size=$(size(L2))")

    ngl2 = ngl*ngl
    parent_grid = zeros(ngl, ngl)
    n_nan = 0
    for idx = 1:ngl2
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        v = get(received, (remote_gel, d, idx), NaN)
        isnan(v) && (n_nan += 1)
        parent_grid[p1,p2] = v
    end
    println("Rank $rank: [NCF-DEBUG] parent_grid (from received): $parent_grid")
    println("Rank $rank: [NCF-DEBUG] $n_nan / $ngl2 entries were NaN (missing)")

    if target_node_idx !== nothing
        ci = div(target_node_idx-1, ngl) + 1
        cj = mod(target_node_idx-1, ngl) + 1
        println("Rank $rank: [NCF-DEBUG] target node_idx=$target_node_idx -> ci=$ci, cj=$cj")
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        println("Rank $rank: [NCF-DEBUG] w_1 (L1[cj=$cj,:]) = $w_1   (sum=$(sum(w_1)))")
        println("Rank $rank: [NCF-DEBUG] w_2 (L2[ci=$ci,:]) = $w_2   (sum=$(sum(w_2)))")
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += w_1[p1]*w_2[p2]*parent_grid[p1,p2]
        end
        println("Rank $rank: [NCF-DEBUG] computed val = $val")
    end

    ngl2 = ngl*ngl
    parent_grid2 = zeros(ngl, ngl)
    for idx = 1:ngl2
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        parent_grid2[p1,p2] = get(received, (remote_gel, d, idx), NaN)
    end
    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        ci = div(node_idx-1, ngl) + 1
        cj = mod(node_idx-1, ngl) + 1
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        v = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            v += w_1[p1]*w_2[p2]*parent_grid2[p1,p2]
        end
        child_vals[node_idx] = v
    end
    return child_vals
end

function exchange_ncf_face_coordinates_v2(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, comm)
    ngl = mesh.ngl
    flat_targets = Int[]; flat_gel = TInt[]; flat_gid = Float64[]
    flat_x = Float64[]; flat_y = Float64[]; flat_z = Float64[]

    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for k in krange, j in jrange, i in irange
            ip = mesh.connijk[iel,i,j,k]
            push!(flat_targets, Int(target_rank)); push!(flat_gel, my_gel)
            push!(flat_gid, Float64(mesh.ip2gip[ip]))
            push!(flat_x, mesh.x[ip]); push!(flat_y, mesh.y[ip]); push!(flat_z, mesh.z[ip])
        end
    end

    gel_cache = SendReceiveCache(comm, Float64.(flat_gel), flat_targets)
    gid_cache = SendReceiveCache(comm, flat_gid, flat_targets)
    x_cache   = SendReceiveCache(comm, flat_x, flat_targets)
    y_cache   = SendReceiveCache(comm, flat_y, flat_targets)
    z_cache   = SendReceiveCache(comm, flat_z, flat_targets)
    send_and_receive!(gel_cache, Float64.(flat_gel))
    send_and_receive!(gid_cache, flat_gid)
    send_and_receive!(x_cache, flat_x)
    send_and_receive!(y_cache, flat_y)
    send_and_receive!(z_cache, flat_z)

    # ── Grouped by owning gel — duplicates (same gid sent via multiple
    # faces) simply overwrite with the SAME coordinate, harmless. ───────────
    remote_gid_coords = Dict{TInt, Dict{Int,NTuple{3,Float64}}}()
    for m = 1:length(gel_cache.combined_recv_data)
        rgel = TInt(round(Int, gel_cache.combined_recv_data[m]))
        gid  = Int(round(Int, gid_cache.combined_recv_data[m]))
        coord = (x_cache.combined_recv_data[m], y_cache.combined_recv_data[m], z_cache.combined_recv_data[m])
        d = get!(remote_gid_coords, rgel, Dict{Int,NTuple{3,Float64}}())
        d[gid] = coord
    end
    return remote_gid_coords
end

function build_dg_ncf_parent_to_child_remote_by_coords(mesh, child_iel::Int, fid::Int, d::Int,
                                                        remote_gel::TInt,
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::Dict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
                                                        tol::Float64 = 1e-4)
    irange,jrange,krange = face_ip_range(fid, ngl)
    IPc_nodes = Int[]; IPc_coords = NTuple{3,Float64}[]
    for k in krange, j in jrange, i in irange
        ip = mesh.connijk[child_iel,i,j,k]
        push!(IPc_nodes, ip)
        push!(IPc_coords, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
    end

    xs = [c[1] for c in IPc_coords]; ys = [c[2] for c in IPc_coords]; zs = [c[3] for c in IPc_coords]
    spread = (maximum(xs)-minimum(xs), maximum(ys)-minimum(ys), maximum(zs)-minimum(zs))
    excl_axis = argmin(spread)
    my_normal = excl_axis==1 ? sum(xs)/length(xs) : (excl_axis==2 ? sum(ys)/length(ys) : sum(zs)/length(zs))
    axes = setdiff(1:3, [excl_axis]); a1, a2 = axes[1], axes[2]
    getc(c,ax) = ax==1 ? c[1] : (ax==2 ? c[2] : c[3])

    haskey(remote_gid_coords, remote_gel) || error("no coordinate data received for remote_gel=$remote_gel")
    candidates = [(gid,c) for (gid,c) in remote_gid_coords[remote_gel] if isapprox(getc(c,excl_axis), my_normal; atol=tol)]
    length(candidates) == ngl*ngl || error("build_dg_ncf_parent_to_child_remote_by_coords: expected $(ngl*ngl) " *
        "matching remote nodes for remote_gel=$remote_gel (child_iel=$child_iel,fid=$fid,d=$d) — found $(length(candidates))")

    parent_a1_vals = sort(unique(round.(getc.(map(x->x[2],candidates), a1); digits=8)))
    parent_a2_vals = sort(unique(round.(getc.(map(x->x[2],candidates), a2); digits=8)))
    (length(parent_a1_vals) == ngl && length(parent_a2_vals) == ngl) ||
        error("degenerate parent grid: $(length(parent_a1_vals))x$(length(parent_a2_vals)) distinct axis values, expected $ngl x $ngl")

    parent_grid_gid = Dict{Tuple{Float64,Float64}, Int}()
    for (gid,c) in candidates
        key = (round(getc(c,a1);digits=8), round(getc(c,a2);digits=8))
        parent_grid_gid[key] = gid
    end

    parent_val_grid = zeros(ngl,ngl)
    for (i1,v1) in enumerate(parent_a1_vals), (i2,v2) in enumerate(parent_a2_vals)
        gid = parent_grid_gid[(v1,v2)]
        parent_val_grid[i1,i2] = get(gid_received, (remote_gel, d, gid), NaN)
    end

    child_a1 = getc.(IPc_coords, a1); child_a2 = getc.(IPc_coords, a2)
    L1 = build_1d_lagrange_matrix(parent_a1_vals, child_a1)   # (length(IPc), ngl)
    L2 = build_1d_lagrange_matrix(parent_a2_vals, child_a2)   # (length(IPc), ngl)

    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += L1[node_idx,p1] * L2[node_idx,p2] * parent_val_grid[p1,p2]
        end
        child_vals[node_idx] = val
    end
    return child_vals
end

function test_build_dg_ncf_parent_to_child_remote_by_coords(inputs, nparts, distribute, sem, partitioned_model;
                                                              spatial_target_iel::Int = 315)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl = mesh.ngl
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ncf_as_child_keys = Set{Tuple{Int,Int}}(); ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end
    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem, fid = 1:6
        neighbor_elem[iel][fid] == zero(TInt) || continue
        (iel,fid) in ncf_as_child_keys && continue
        (iel,fid) in ncf_as_parent_keys && continue
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1=mesh.conn_face_el[1,fe,iel]; n2=mesh.conn_face_el[2,fe,iel]; n3=mesh.conn_face_el[3,fe,iel]; n4=mesh.conn_face_el[4,fe,iel]
        key = Tuple(sort(TInt[n1,n2,n3,n4]))
        key in rigid_face_keys && continue
        push!(suspected_faces, (iel,fid))
    end
    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved")

    n_dirs = fill(1, mesh.nelem)   # single "direction" — pure spatial test
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    remote_gid_coords = exchange_ncf_face_coordinates_v2(mesh, neighbor_rank, comm)

    # ── Known polynomial field: f(x,y,z) = 1 + 2x + 3y + 4z (degree 1,
    # exactly reproducible by ANY order-ngl Lagrange interpolation). ────────
    f(x,y,z) = 1.0 + 2.0*x/1000.0 + 3.0*y/1000.0 + 4.0*z/1000.0   # scaled for reasonable magnitude

    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        vals = zeros(ngl^3)
        for k=1:ngl,j=1:ngl,i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            vals[lidx(i,j,k,ngl)] = f(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        end
        x_dg[(iel,1)] = vals
    end
    gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)

    n_checked = 0; max_err = 0.0
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_fine || continue
        remote_gel_f = neighbor_gel_group[(iel,fid)][1]
        child_vals = build_dg_ncf_parent_to_child_remote_by_coords(mesh, iel, fid, 1, remote_gel_f,
                                                                     remote_gid_coords, gid_received, ngl)
        irange,jrange,krange = face_ip_range(fid, ngl)
        idx = 0
        for k in krange, j in jrange, i in irange
            idx += 1
            gip = mesh.connijk[iel,i,j,k]
            exact = f(mesh.x[gip], mesh.y[gip], mesh.z[gip])
            err = abs(child_vals[idx] - exact)
            max_err = max(max_err, err)
            n_checked += 1
        end
    end
    println("Rank $rank: checked $n_checked fine-child face nodes, max_err=$max_err  " *
            (max_err < 1e-9 ? "✓ PASS" : "✗ FAIL"))
end

function build_dg_ncf_child_to_parent_remote_by_coords(mesh, coarse_iel::Int, fid::Int, d::Int,
                                                        remote_gels::Vector{TInt},
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::Dict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
                                                        tol::Float64 = 1e-4)
    irange,jrange,krange = face_ip_range(fid, ngl)
    coarse_coords = NTuple{3,Float64}[]
    for k in krange, j in jrange, i in irange
        ip = mesh.connijk[coarse_iel,i,j,k]
        push!(coarse_coords, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
    end

    xs = [c[1] for c in coarse_coords]; ys = [c[2] for c in coarse_coords]; zs = [c[3] for c in coarse_coords]
    spread = (maximum(xs)-minimum(xs), maximum(ys)-minimum(ys), maximum(zs)-minimum(zs))
    excl_axis = argmin(spread)
    axes = setdiff(1:3, [excl_axis]); a1, a2 = axes[1], axes[2]
    getc(c,ax) = ax==1 ? c[1] : (ax==2 ? c[2] : c[3])

    # ── For EACH remote_gel (fine child), recover its OWN ngl×ngl grid of
    # (a1,a2) values from its coordinate set, and its bounding box on the
    # shared plane — used to determine which child a coarse point falls
    # inside, and to build that child's OWN 1D Lagrange node lines. ─────────
    child_grids = Dict{TInt, NamedTuple}()
    for rg in remote_gels
        haskey(remote_gid_coords, rg) || error("no coordinate data for remote_gel=$rg")
        coords_this = collect(remote_gid_coords[rg])   # [(gid,coord), ...]
        # normal-axis filter: keep only nodes on the shared plane, matching
        # the coarse face's own normal-axis value (same principle as the
        # fine-to-coarse direction).
        my_normal_val = excl_axis==1 ? sum(xs)/length(xs) : (excl_axis==2 ? sum(ys)/length(ys) : sum(zs)/length(zs))
        on_plane = [(gid,c) for (gid,c) in coords_this if isapprox(getc(c,excl_axis), my_normal_val; atol=tol)]
        length(on_plane) == ngl*ngl || error("build_dg_ncf_child_to_parent_remote_by_coords: remote_gel=$rg " *
            "has $(length(on_plane)) on-plane nodes, expected $(ngl*ngl)")

        a1_vals = sort(unique(round.(getc.(map(x->x[2],on_plane), a1); digits=8)))
        a2_vals = sort(unique(round.(getc.(map(x->x[2],on_plane), a2); digits=8)))
        (length(a1_vals)==ngl && length(a2_vals)==ngl) || error("degenerate child grid for remote_gel=$rg")

        grid_gid = Dict{Tuple{Float64,Float64},Int}()
        for (gid,c) in on_plane
            grid_gid[(round(getc(c,a1);digits=8), round(getc(c,a2);digits=8))] = gid
        end

        bbox = (minimum(a1_vals)-tol, maximum(a1_vals)+tol, minimum(a2_vals)-tol, maximum(a2_vals)+tol)
        child_grids[rg] = (a1_vals=a1_vals, a2_vals=a2_vals, grid_gid=grid_gid, bbox=bbox)
    end

    parent_grid = zeros(ngl,ngl)
    for idx = 1:(ngl*ngl)
        p1 = ((idx-1) % ngl) + 1; p2 = ((idx-1) ÷ ngl) + 1
        cc = coarse_coords[idx]
        cv1 = getc(cc,a1); cv2 = getc(cc,a2)

        owning_rg = nothing
        for (rg, cg) in child_grids
            (cv1 >= cg.bbox[1] && cv1 <= cg.bbox[2] && cv2 >= cg.bbox[3] && cv2 <= cg.bbox[4]) || continue
            owning_rg = rg; break
        end
        owning_rg === nothing && error("build_dg_ncf_child_to_parent_remote_by_coords: coarse node " *
            "(a1=$cv1,a2=$cv2) not contained in ANY child's bounding box")

        cg = child_grids[owning_rg]
        val_grid = zeros(ngl,ngl)
        for (i1,v1) in enumerate(cg.a1_vals), (i2,v2) in enumerate(cg.a2_vals)
            gid = cg.grid_gid[(v1,v2)]
            val_grid[i1,i2] = get(gid_received, (owning_rg, d, gid), NaN)
        end

        # EXACT match short-circuit (shared corner/edge nodes — common case)
        exact_i1 = findfirst(v -> isapprox(v, cv1; atol=tol), cg.a1_vals)
        exact_i2 = findfirst(v -> isapprox(v, cv2; atol=tol), cg.a2_vals)
        if exact_i1 !== nothing && exact_i2 !== nothing
            parent_grid[p1,p2] = val_grid[exact_i1, exact_i2]
        else
            L1 = build_1d_lagrange_matrix(cg.a1_vals, [cv1])
            L2 = build_1d_lagrange_matrix(cg.a2_vals, [cv2])
            v = 0.0
            for q1=1:ngl, q2=1:ngl
                v += L1[1,q1]*L2[1,q2]*val_grid[q1,q2]
            end
            parent_grid[p1,p2] = v
        end
    end
    return parent_grid
end

function test_build_dg_ncf_child_to_parent_remote_by_coords(inputs, nparts, distribute, sem, partitioned_model;
                                                              spatial_target_iel::Int = 315)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl = mesh.ngl
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf), ngl=$ngl")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ncf_as_child_keys = Set{Tuple{Int,Int}}(); ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end
    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem, fid = 1:6
        neighbor_elem[iel][fid] == zero(TInt) || continue
        (iel,fid) in ncf_as_child_keys && continue
        (iel,fid) in ncf_as_parent_keys && continue
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1=mesh.conn_face_el[1,fe,iel]; n2=mesh.conn_face_el[2,fe,iel]; n3=mesh.conn_face_el[3,fe,iel]; n4=mesh.conn_face_el[4,fe,iel]
        key = Tuple(sort(TInt[n1,n2,n3,n4]))
        key in rigid_face_keys && continue
        push!(suspected_faces, (iel,fid))
    end
    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved")

    n_dirs = fill(1, mesh.nelem)
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    remote_gid_coords = exchange_ncf_face_coordinates_v2(mesh, neighbor_rank, comm)

    # ── Degree (ngl-1) polynomial in EACH coordinate — exactly reproducible
    # by an order-(ngl-1) Lagrange basis, sensitive to node selection. ───────
    p = ngl - 1
    f(x,y,z) = sum((x/1000.0)^k for k=0:p) + sum((y/1000.0)^k for k=0:p) + sum((z/1000.0)^k for k=0:p)

    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        vals = zeros(ngl^3)
        for k=1:ngl,j=1:ngl,i=1:ngl
            gip = mesh.connijk[iel,i,j,k]
            vals[lidx(i,j,k,ngl)] = f(mesh.x[gip], mesh.y[gip], mesh.z[gip])
        end
        x_dg[(iel,1)] = vals
    end
    gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)

    n_checked = 0; max_err = 0.0
    for ((iel,fid), kind) in neighbor_kind
        kind == :i_am_coarse || continue
        remote_gels_c = neighbor_gel_group[(iel,fid)]
        parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, 1, remote_gels_c,
                                                                      remote_gid_coords, gid_received, ngl)
        irange,jrange,krange = face_ip_range(fid, ngl)
        idx = 0
        for k in krange, j in jrange, i in irange
            idx += 1
            p1 = ((idx-1) % ngl) + 1; p2 = ((idx-1) ÷ ngl) + 1
            gip = mesh.connijk[iel,i,j,k]
            exact = f(mesh.x[gip], mesh.y[gip], mesh.z[gip])
            err = abs(parent_grid[p1,p2] - exact)
            max_err = max(max_err, err)
            n_checked += 1
        end
    end
    println("Rank $rank: checked $n_checked coarse face nodes, max_err=$max_err  " *
            (max_err < 1e-9 ? "✓ PASS" : "✗ FAIL"))
end

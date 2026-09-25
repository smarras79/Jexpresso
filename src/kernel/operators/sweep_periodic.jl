using Random

"""
    build_tagged_face_keys_periodic_aware(mesh)

Same purpose as the original tagged_face_keys construction, but splits
boundary-tagged faces into TWO sets: genuinely rigid/true-boundary faces
(rigid_face_keys) and periodic-tagged faces (periodic_face_keys) — so
periodic faces are NEVER swallowed as ordinary true boundaries.
"""
function build_tagged_face_keys_periodic_aware(mesh)
    rigid_face_keys = Set{NTuple{4,TInt}}()
    periodic_face_keys = Set{NTuple{4,TInt}}()

    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, mesh.ngl, mesh.ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, mesh.ngl]
        key = Tuple(sort(TInt[ip1, ip2, ip3, ip4]))

        tag = mesh.bdy_face_type[iface_bdy]   # exact field name/type TBD — see caveat
        tag_str = string(tag)
        if occursin("periodic", lowercase(tag_str))
            push!(periodic_face_keys, key)
        else
            push!(rigid_face_keys, key)
        end
    end

    return rigid_face_keys, periodic_face_keys
end

"""
    resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, comm)

Generalizes resolve_rank_cut_neighbors_ncf to ALSO match SAME-RANK
candidates (needed for periodic pairs, which can legitimately be on the
same rank) — same gid/coordinate matching and containment logic,
unchanged, just without the `e.rank != rank` exclusion.

Returns, in ADDITION to the usual cross-rank tables, same-rank tables:
  same_rank_periodic_kind :: Dict{(iel,fid), Symbol}  (:conforming/:i_am_coarse/:i_am_fine)
  same_rank_periodic_partner :: Dict{(iel,fid), Vector{Int}}  (LOCAL element id(s), same rank)
Cross-rank tables (neighbor_rank, neighbor_kind, neighbor_gel_group) are
returned exactly as before, for entries that turned out to be cross-rank.
"""
function resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces::Vector{Tuple{TInt,TInt}}, comm;
                                                     tol::Float64 = 1e-6)
    rank = MPI.Comm_rank(comm)

    local_keys   = NTuple{4,TInt}[]
    local_coords = Vector{NTuple{4,NTuple{3,Float64}}}()
    local_iel = TInt[]; local_fid = TInt[]; local_gel = TInt[]; local_lel = TInt[]  # NEW: local elem id too

    for (iel, fid) in suspected_faces
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
        key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        coords = ((mesh.x[n1],mesh.y[n1],mesh.z[n1]), (mesh.x[n2],mesh.y[n2],mesh.z[n2]),
                  (mesh.x[n3],mesh.y[n3],mesh.z[n3]), (mesh.x[n4],mesh.y[n4],mesh.z[n4]))
        push!(local_keys, key); push!(local_coords, coords)
        push!(local_iel, iel); push!(local_fid, fid); push!(local_gel, mesh.el2gel[iel]); push!(local_lel, iel)
    end

    n_local = TInt(length(local_keys))
    all_n = MPI.Allgather(n_local, comm)
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

    neighbor_rank = Dict{Tuple{TInt,TInt}, TInt}()
    neighbor_kind = Dict{Tuple{TInt,TInt}, Symbol}()
    neighbor_gel_group = Dict{Tuple{TInt,TInt}, Vector{TInt}}()
    same_rank_kind = Dict{Tuple{Int,Int}, Symbol}()
    same_rank_partner = Dict{Tuple{Int,Int}, Vector{Int}}()
    n_unresolved = 0


    # Map global gel -> local element id, IF it's on my own rank
    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)

    for i = eachindex(local_keys)
        my_key = local_keys[i]; my_iel = local_iel[i]; my_fid = local_fid[i]
        my_bbox = bbox_excl_normal(local_coords[i])

        # Exact match — NOW INCLUDING same-rank candidates (the key change).
        exact_matches = [e for e in entries if e.key == my_key && e.gel != local_gel[i]]  # exclude SELF only
        if length(exact_matches) == 1
            m = exact_matches[1]
            if m.rank == rank
                lel = get(gel_to_local, m.gel, 0)
                lel == 0 && (n_unresolved += 1; continue)
                same_rank_kind[(my_iel,my_fid)] = :conforming
                same_rank_partner[(my_iel,my_fid)] = [lel]
            else
                neighbor_rank[(my_iel,my_fid)] = m.rank
                neighbor_kind[(my_iel,my_fid)] = :conforming
                neighbor_gel_group[(my_iel,my_fid)] = [m.gel]
            end
            continue
        end

        candidates_inside = NamedTuple[]; candidates_containing = NamedTuple[]
        for e in entries
            e.gel == local_gel[i] && continue
            ebbox = bbox_excl_normal(e.corners)
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
            all_same_rank = all(c.rank == rank for c in candidates_inside)
            if all_same_rank
                lels = [get(gel_to_local, c.gel, 0) for c in candidates_inside]
                any(==(0), lels) ? (n_unresolved += 1) : begin
                    same_rank_kind[(my_iel,my_fid)] = :i_am_coarse
                    same_rank_partner[(my_iel,my_fid)] = lels
                end
            elseif all(c.rank == candidates_inside[1].rank for c in candidates_inside)
                neighbor_rank[(my_iel,my_fid)] = candidates_inside[1].rank
                neighbor_kind[(my_iel,my_fid)] = :i_am_coarse
                neighbor_gel_group[(my_iel,my_fid)] = [c.gel for c in candidates_inside]
            else
                n_unresolved += 1   # MIXED same-rank/cross-rank children — not yet supported
            end
        elseif length(candidates_containing) == 1
            c = candidates_containing[1]
            if c.rank == rank
                lel = get(gel_to_local, c.gel, 0)
                lel == 0 ? (n_unresolved += 1) : begin
                    same_rank_kind[(my_iel,my_fid)] = :i_am_fine
                    same_rank_partner[(my_iel,my_fid)] = [lel]
                end
            else
                neighbor_rank[(my_iel,my_fid)] = c.rank
                neighbor_kind[(my_iel,my_fid)] = :i_am_fine
                neighbor_gel_group[(my_iel,my_fid)] = [c.gel]
            end
        else
            n_unresolved += 1
        end
    end

    return neighbor_rank, neighbor_kind, neighbor_gel_group, same_rank_kind, same_rank_partner, n_unresolved
end

function test_periodic_detection_and_resolution(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    println("Rank $rank: mesh has $(mesh.nelem) elements")

    # ── Gap A check ──────────────────────────────────────────────────────────
    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    println("Rank $rank: rigid_face_keys=$(length(rigid_face_keys)), periodic_face_keys=$(length(periodic_face_keys))")

    n_periodic_in_rigid = length(intersect(rigid_face_keys, periodic_face_keys))
    println("Rank $rank: periodic keys also in rigid set: $n_periodic_in_rigid  " *
            (n_periodic_in_rigid==0 ? "✓" : "✗ SHOULD BE ZERO"))

    # ── Build suspected faces: same-rank NCF-excluded, RIGID-excluded,
    # but PERIODIC-INCLUDED (the actual fix). ────────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue   # exclude TRUE rigid boundaries only
            push!(suspected_faces, (iel, fid))
        end
    end
    println("Rank $rank: $(length(suspected_faces)) suspected periodic/rank-cut faces (rigid excluded)")

    # ── Gap B: resolve, now INCLUDING same-rank matches ──────────────────────
    neighbor_rank, neighbor_kind, neighbor_gel_group, same_rank_kind, same_rank_partner, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, comm)

    println("Rank $rank: cross-rank matches: $(length(neighbor_rank)), " *
            "same-rank matches: $(length(same_rank_kind)), unresolved: $n_unresolved  " *
            (n_unresolved==0 ? "✓" : "✗"))

    # ── Mutual consistency: for each SAME-RANK :conforming match, does the
    # partner ALSO report back a :conforming match pointing at ME? ──────────
    n_checked = 0; n_bad = 0
    for ((iel,fid), kind) in same_rank_kind
        kind == :conforming || continue
        partner_lel = same_rank_partner[(iel,fid)][1]
        n_checked += 1
        found_back = false
        for fid2 = 1:6
            haskey(same_rank_kind, (partner_lel,fid2)) || continue
            same_rank_kind[(partner_lel,fid2)] == :conforming || continue
            same_rank_partner[(partner_lel,fid2)] == [iel] || continue
            found_back = true; break
        end
        found_back || (n_bad += 1; println("Rank $rank: ✗ (iel=$iel,fid=$fid)->lel=$partner_lel " *
                                            "has NO reciprocal match"))
    end
    println("Rank $rank: mutual-consistency check: $n_checked same-rank conforming pairs, $n_bad inconsistent  " *
            (n_bad==0 ? "✓" : "✗"))

    return neighbor_rank, neighbor_kind, neighbor_gel_group, same_rank_kind, same_rank_partner
end


function test_periodic_ncf_detection(inputs, nparts, distribute, sem, partitioned_model;
                                      spatial_target_iel::Int = 1)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    # ── Refine an element ON THE PERIODIC BOUNDARY itself, so its periodic
    # partner face becomes a genuine coarse/fine NCF relationship. ──────────
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    debug_periodic_ncf_containment(mesh, suspected_faces, get_mpi_comm())

    println("Rank $rank: $(length(suspected_faces)) suspected periodic/rank-cut faces")

    neighbor_rank, neighbor_kind, neighbor_gel_group, same_rank_kind, same_rank_partner, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, comm)

    println("Rank $rank: cross-rank=$(length(neighbor_rank)), same-rank=$(length(same_rank_kind)), " *
            "unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    n_ncf_same = count(v -> v != :conforming, values(same_rank_kind))
    n_ncf_cross = count(v -> v != :conforming, values(neighbor_kind))
    println("Rank $rank: NCF (non-conforming) matches — same-rank: $n_ncf_same, cross-rank: $n_ncf_cross")
    for ((iel,fid),kind) in same_rank_kind
        kind == :conforming && continue
        println("  SAME-RANK NCF: (iel=$iel,fid=$fid) kind=$kind partner=$(same_rank_partner[(iel,fid)])")
    end
    for ((iel,fid),kind) in neighbor_kind
        kind == :conforming && continue
        println("  CROSS-RANK NCF: (iel=$iel,fid=$fid) kind=$kind remote_gels=$(neighbor_gel_group[(iel,fid)])")
    end

    return neighbor_rank, neighbor_kind, neighbor_gel_group, same_rank_kind, same_rank_partner
end

function debug_periodic_ncf_containment(mesh, suspected_faces, comm; target_iel::Int = 1)
    rank = MPI.Comm_rank(comm)

    local_keys   = NTuple{4,TInt}[]
    local_coords = Vector{NTuple{4,NTuple{3,Float64}}}()
    local_iel = TInt[]; local_fid = TInt[]; local_gel = TInt[]

    for (iel, fid) in suspected_faces
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
        key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        coords = ((mesh.x[n1],mesh.y[n1],mesh.z[n1]), (mesh.x[n2],mesh.y[n2],mesh.z[n2]),
                  (mesh.x[n3],mesh.y[n3],mesh.z[n3]), (mesh.x[n4],mesh.y[n4],mesh.z[n4]))
        push!(local_keys, key); push!(local_coords, coords)
        push!(local_iel, iel); push!(local_fid, fid); push!(local_gel, mesh.el2gel[iel])
    end

    n_local = TInt(length(local_keys))
    all_n = MPI.Allgather(n_local, comm)
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

    for i = eachindex(local_keys)
        local_iel[i] == target_iel || continue
        my_bbox = bbox_excl_normal(local_coords[i])
        println("Rank $rank: MY (iel=$target_iel,fid=$(local_fid[i])) bbox_excl_normal=$my_bbox " *
                "corners=$(local_coords[i])")

        for e in entries
            e.gel == local_gel[i] && continue
            ebbox = bbox_excl_normal(e.corners)
            inside = (ebbox[1] >= my_bbox[1]-1e-6 && ebbox[2] <= my_bbox[2]+1e-6 &&
                      ebbox[3] >= my_bbox[3]-1e-6 && ebbox[4] <= my_bbox[4]+1e-6 &&
                      ebbox[5] >= my_bbox[5]-1e-6 && ebbox[6] <= my_bbox[6]+1e-6)
            containing = (my_bbox[1] >= ebbox[1]-1e-6 && my_bbox[2] <= ebbox[2]+1e-6 &&
                          my_bbox[3] >= ebbox[3]-1e-6 && my_bbox[4] <= ebbox[4]+1e-6 &&
                          my_bbox[5] >= ebbox[5]-1e-6 && my_bbox[6] <= ebbox[6]+1e-6)
            (inside || containing) || continue
            println("  candidate (iel=$(e.iel),fid=$(e.fid),rank=$(e.rank)): ebbox=$ebbox " *
                    "inside=$inside containing=$containing corners=$(e.corners)")
        end
    end
end

"""
    build_periodic_equivalence_classes(mesh, same_rank_kind, same_rank_partner,
                                        neighbor_rank, neighbor_kind, neighbor_gel_group, comm)

Computes the TRANSITIVE CLOSURE of periodic pairing (same-rank AND
cross-rank) into equivalence classes: if A~B and B~C then A,B,C all
belong to the same class. Each class must be entirely refined or
entirely left alone as a unit.

Returns: Dict{TInt, Set{TInt}} mapping EVERY periodically-connected
element's global id (gel) to the FULL SET of gels in its class
(including itself). Elements with no periodic pairing at all are
absent from this dict (their "class" is just themselves, trivially).
"""
function build_periodic_equivalence_classes(mesh, same_rank_kind::Dict{Tuple{Int,Int},Symbol},
                                             same_rank_partner::Dict{Tuple{Int,Int},Vector{Int}},
                                             neighbor_rank::Dict{Tuple{TInt,TInt},TInt},
                                             neighbor_kind::Dict{Tuple{TInt,TInt},Symbol},
                                             neighbor_gel_group::Dict{Tuple{TInt,TInt},Vector{TInt}},
                                             comm)
    rank = MPI.Comm_rank(comm)

    # ── Build a LOCAL adjacency list (gel -> [periodic partner gels]),
    # from BOTH same-rank and cross-rank periodic pairs (kind doesn't
    # matter for equivalence-class purposes — :conforming since we're
    # enforcing angular conformity across periodic boundaries, per the
    # simplification just agreed on). ────────────────────────────────────────
    local_adj = Dict{TInt, Vector{TInt}}()
    function add_edge!(g1::TInt, g2::TInt)
        push!(get!(local_adj, g1, TInt[]), g2)
        push!(get!(local_adj, g2, TInt[]), g1)
    end

    for ((iel,fid), kind) in same_rank_kind
        my_gel = TInt(mesh.el2gel[iel])
        for lel in same_rank_partner[(iel,fid)]
            add_edge!(my_gel, TInt(mesh.el2gel[lel]))
        end
    end
    for ((iel,fid), kind) in neighbor_kind
        my_gel = TInt(mesh.el2gel[iel])
        for rgel in neighbor_gel_group[(iel,fid)]
            add_edge!(my_gel, rgel)
        end
    end

    # ── Exchange ALL local edges globally (small — one entry per periodic
    # face pair, not per node) so every rank can compute the SAME global
    # transitive closure independently — avoiding any need for iterative
    # cross-rank consensus-building. ─────────────────────────────────────────
    local_edges = Tuple{TInt,TInt}[]
    for (g1, partners) in local_adj
        for g2 in partners
            g1 < g2 && push!(local_edges, (g1,g2))   # dedupe direction
        end
    end
    n_local = length(local_edges)
    all_n = MPI.Allgather(TInt(n_local), comm)
    flat_local = Float64[]
    for (g1,g2) in local_edges
        append!(flat_local, Float64[g1,g2])
    end
    counts = all_n .* TInt(2)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    global_edges = Tuple{TInt,TInt}[]
    for m = 1:(length(all_flat)÷2)
        push!(global_edges, (TInt(all_flat[2m-1]), TInt(all_flat[2m])))
    end

    # ── Union-Find over the GLOBAL edge list — every rank computes this
    # identically, no coordination needed beyond the one Allgatherv above. ──
    parent = Dict{TInt,TInt}()
    find(x) = begin
        haskey(parent, x) || (parent[x] = x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    union!(a,b) = begin
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end

    for (g1,g2) in global_edges
        union!(g1,g2)
    end

    # ── Build final classes: root -> Set of all members ──────────────────────
    all_gels = Set{TInt}()
    for (g1,g2) in global_edges
        push!(all_gels, g1); push!(all_gels, g2)
    end
    classes_by_root = Dict{TInt, Set{TInt}}()
    for g in all_gels
        r = find(g)
        push!(get!(classes_by_root, r, Set{TInt}()), g)
    end

    gel_to_class = Dict{TInt, Set{TInt}}()
    for (root, members) in classes_by_root
        for g in members
            gel_to_class[g] = members
        end
    end

    return gel_to_class
end

function test_periodic_equivalence_classes(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                        neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)




    # ── Classify EVERY local element by its class SIZE (1 = not periodic
    # at all — i.e., absent from gel_to_class; 2 = single periodic;
    # 4 = double periodic corner). ────────────────────────────────────────────
    size_counts = Dict{Int,Int}()
    for iel = 1:mesh.nelem
        gel = TInt(mesh.el2gel[iel])
        sz = haskey(gel_to_class, gel) ? length(gel_to_class[gel]) : 1
        size_counts[sz] = get(size_counts, sz, 0) + 1
    end
    println("Rank $rank: class-size distribution: $size_counts")

    # ── Mutual consistency: pick one 4-member class (if any exist locally)
    # and confirm EVERY member (even remote ones) agrees on the SAME set. ────
    checked_a_quad = false
    for (gel, members) in gel_to_class
        length(members) == 4 || continue
        checked_a_quad = true
        println("Rank $rank: sample QUAD class containing gel=$gel: members=$(sort(collect(members)))")
        break
    end
    checked_a_quad || println("Rank $rank: no quad (double-periodic) class found locally")

    return gel_to_class, size_counts
end

function debug_periodic_class_merge(mesh, same_rank_kind, same_rank_partner,
                                     neighbor_rank, neighbor_kind, neighbor_gel_group, comm;
                                     suspect_gels::Vector{Int} = [401, 410, 491, 500, 501, 510, 591, 600])
    rank = MPI.Comm_rank(comm)

    local_adj = Dict{TInt, Vector{TInt}}()
    edge_source = Dict{Tuple{TInt,TInt}, String}()   # NEW — record WHERE each edge came from

    function add_edge!(g1::TInt, g2::TInt, source::String)
        push!(get!(local_adj, g1, TInt[]), g2)
        push!(get!(local_adj, g2, TInt[]), g1)
        key = g1 < g2 ? (g1,g2) : (g2,g1)
        edge_source[key] = source
    end

    for ((iel,fid), kind) in same_rank_kind
        my_gel = TInt(mesh.el2gel[iel])
        for lel in same_rank_partner[(iel,fid)]
            add_edge!(my_gel, TInt(mesh.el2gel[lel]), "SAME-RANK (iel=$iel,fid=$fid)->lel=$lel")
        end
    end
    for ((iel,fid), kind) in neighbor_kind
        my_gel = TInt(mesh.el2gel[iel])
        for rgel in neighbor_gel_group[(iel,fid)]
            add_edge!(my_gel, rgel, "CROSS-RANK (iel=$iel,fid=$fid)->rgel=$rgel")
        end
    end

    local_edges = Tuple{TInt,TInt}[]
    for (g1, partners) in local_adj
        for g2 in partners
            g1 < g2 && push!(local_edges, (g1,g2))
        end
    end
    n_local = length(local_edges)
    all_n = MPI.Allgather(TInt(n_local), comm)
    flat_local = Float64[]
    for (g1,g2) in local_edges
        append!(flat_local, Float64[g1,g2])
    end
    counts = all_n .* TInt(2)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)
    global_edges = Tuple{TInt,TInt}[]
    for m = 1:(length(all_flat)÷2)
        push!(global_edges, (TInt(all_flat[2m-1]), TInt(all_flat[2m])))
    end

    # ── Print EVERY edge touching any suspect gel, with its LOCAL source
    # label if this rank knows it (the edge may have been contributed by
    # the OTHER rank, in which case we only have the (g1,g2) pair itself). ──
    suspect_set = Set(TInt.(suspect_gels))
    println("Rank $rank: edges touching suspect gels $suspect_gels:")
    for (g1,g2) in global_edges
        (g1 in suspect_set || g2 in suspect_set) || continue
        key = g1 < g2 ? (g1,g2) : (g2,g1)
        src = get(edge_source, key, "(from OTHER rank — no local source label)")
        println("  ($g1,$g2): $src")
    end
end

function build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                 neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)
    rank = MPI.Comm_rank(comm)
    local_adj = Dict{TInt, Vector{TInt}}()
    function add_edge!(g1::TInt, g2::TInt)
        push!(get!(local_adj, g1, TInt[]), g2)
        push!(get!(local_adj, g2, TInt[]), g1)
    end

    for ((iel,fid), kind) in same_rank_kind
        get(same_rank_is_periodic, (iel,fid), false) || continue   # NEW filter
        my_gel = TInt(mesh.el2gel[iel])
        for lel in same_rank_partner[(iel,fid)]
            add_edge!(my_gel, TInt(mesh.el2gel[lel]))
        end
    end
    for ((iel,fid), kind) in neighbor_kind
        get(neighbor_is_periodic, (iel,fid), false) || continue   # NEW filter
        my_gel = TInt(mesh.el2gel[iel])
        for rgel in neighbor_gel_group[(iel,fid)]
            add_edge!(my_gel, rgel)
        end
    end

    local_edges = Tuple{TInt,TInt}[]
    for (g1, partners) in local_adj
        for g2 in partners
            g1 < g2 && push!(local_edges, (g1,g2))
        end
    end
    n_local = length(local_edges)
    all_n = MPI.Allgather(TInt(n_local), comm)
    flat_local = Float64[]
    for (g1,g2) in local_edges
        append!(flat_local, Float64[g1,g2])
    end
    counts = all_n .* TInt(2)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)
    global_edges = Tuple{TInt,TInt}[]
    for m = 1:(length(all_flat)÷2)
        push!(global_edges, (TInt(all_flat[2m-1]), TInt(all_flat[2m])))
    end

    parent = Dict{TInt,TInt}()
    find(x) = begin
        haskey(parent, x) || (parent[x] = x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    union!(a,b) = begin
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    for (g1,g2) in global_edges
        union!(g1,g2)
    end

    all_gels = Set{TInt}()
    for (g1,g2) in global_edges
        push!(all_gels, g1); push!(all_gels, g2)
    end
    classes_by_root = Dict{TInt, Set{TInt}}()
    for g in all_gels
        r = find(g)
        push!(get!(classes_by_root, r, Set{TInt}()), g)
    end
    gel_to_class = Dict{TInt, Set{TInt}}()
    for (root, members) in classes_by_root
        for g in members
            gel_to_class[g] = members
        end
    end
    return gel_to_class
end

"""
    build_angular_refinement_list(mesh, candidates, spatial_ncf_same_rank_keys, spatial_ncf_cross_rank_keys,
                                   gel_to_class, comm; buffer_hops::Int = 1)

Input: candidates :: Vector{Tuple{Int,Int}} — (spatial_iel, angular_e_ext)
pairs THIS RANK wants to angularly refine.

1. One-buffer-element filter: removes candidates whose spatial_iel is
   itself spatially-NCF-refined, OR is within `buffer_hops` same-rank
   graph-adjacency of a spatially-NCF-refined element (cross-rank
   adjacency intentionally NOT traversed here — the one-buffer-element
   policy has always been about LOCAL geometric proximity, and a
   same-rank hop-based check is what's already used elsewhere in this
   effort for "far enough apart").
2. Periodic propagation: expands each SURVIVING candidate to its FULL
   periodic equivalence class (via gel_to_class), assigning the SAME
   angular_e_ext to every member — including cross-rank members, which
   get communicated back to their owning rank.
3. Post-propagation re-check: for each EXPANDED class, if ANY member
   fails the one-buffer-element check, the ENTIRE class is dropped.

Returns: Dict{Int, Vector{Tuple{Int,Int}}} mapping RANK -> the (LOCAL
spatial_iel, angular_e_ext) pairs that rank should call
adapt_single_angular_element_conforming! on. (Includes this rank's own
entry too, for uniform handling.)
"""
#=function build_angular_refinement_list(mesh, candidates::Vector{Tuple{Int,Int}},
                                        neighbor_elem, ncf_as_child_keys::Set{Tuple{Int,Int}},
                                        ncf_as_parent_keys::Set{Tuple{Int,Int}},
                                        gel_to_class::Dict{TInt,Set{TInt}}, comm;
                                        neighbor_rank = nothing, neighbor_kind = nothing, neighbor_gel_group = nothing,
                                        buffer_hops::Int = 1)
    rank = MPI.Comm_rank(comm)

    spatial_ncf_elements = Set{Int}()
    for (iel,fid) in ncf_as_child_keys; push!(spatial_ncf_elements, iel); end
    for (iel,fid) in ncf_as_parent_keys; push!(spatial_ncf_elements, iel); end

    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)
    for iel in collect(spatial_ncf_elements)
        g = TInt(mesh.el2gel[iel])
        haskey(gel_to_class, g) || continue
        for partner_g in gel_to_class[g]
            haskey(gel_to_local, partner_g) || continue
            push!(spatial_ncf_elements, gel_to_local[partner_g])
        end
    end

    # ── Cross-rank spatial-NCF-adjacent elements ALSO count as "touching"
    # spatial refinement — the buffer check must be partition-independent,
    # so an element with a genuine cross-rank spatial-NCF face (:i_am_fine
    # or :i_am_coarse) must ALSO be added to spatial_ncf_elements, the
    # SAME way a same-rank spatial-NCF neighbor already is. ──────────────────
    if neighbor_rank !== nothing
        for ((iel,fid), target_rank) in neighbor_rank
            kind = get(neighbor_kind, (iel,fid), :conforming)
            (kind == :i_am_fine || kind == :i_am_coarse) || continue
            push!(spatial_ncf_elements, iel)
        end
    end

    function bfs_within_hops(start_iel, max_hops)
        dist = Dict{Int,Int}(start_iel => 0)
        frontier = [start_iel]
        depth = 0
        while !isempty(frontier) && depth < max_hops
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

    function violates_buffer(iel::Int)
        iel in spatial_ncf_elements && return true
        nearby = bfs_within_hops(iel, buffer_hops)
        any(e -> e in spatial_ncf_elements, keys(nearby))
    end

    step1_survivors = Tuple{Int,Int}[]
    for (iel, e_ext) in candidates
        if violates_buffer(iel)
            println("Rank $rank: candidate (iel=$iel,e_ext=$e_ext) REMOVED — violates one-buffer-element (step 1)")
        else
            push!(step1_survivors, (iel, e_ext))
        end
    end

    expanded_gel_pairs = Dict{TInt, Vector{Int}}()
    for (iel, e_ext) in step1_survivors
        my_gel = TInt(mesh.el2gel[iel])
        class = get(gel_to_class, my_gel, Set([my_gel]))
        for g in class
            list = get!(expanded_gel_pairs, g, Int[])
            e_ext in list || push!(list, e_ext)
        end
    end

    flat_pairs = Tuple{TInt,Int}[]
    for (g, exts) in expanded_gel_pairs, e in exts
        push!(flat_pairs, (g,e))
    end
    n_local = length(flat_pairs)
    all_n = MPI.Allgather(TInt(n_local), comm)
    flat_local = Float64[]
    for (g,e) in flat_pairs
        append!(flat_local, Float64[g,e])
    end
    counts = all_n .* TInt(2)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    global_gel_pairs = Dict{TInt, Vector{Int}}()
    for m = 1:(length(all_flat)÷2)
        g = TInt(all_flat[2m-1]); e = Int(all_flat[2m])
        list = get!(global_gel_pairs, g, Int[])
        e in list || push!(list, e)
    end

    local_violation_gels = TInt[]
    for g in keys(global_gel_pairs)
        haskey(gel_to_local, g) || continue
        lel = gel_to_local[g]
        violates_buffer(lel) && push!(local_violation_gels, g)
    end

    n_local_v = length(local_violation_gels)
    all_nv = MPI.Allgather(TInt(n_local_v), comm)
    flat_v = Float64.(local_violation_gels)
    all_vf = MPI.Allgatherv(flat_v, Cint.(all_nv), comm)
    violating_gels = Set{TInt}(TInt.(all_vf))

    final_gel_pairs = Dict{TInt, Vector{Int}}()
    for (g, exts) in global_gel_pairs
        class = get(gel_to_class, g, Set([g]))
        any(m -> m in violating_gels, class) && continue
        final_gel_pairs[g] = exts
    end

    my_list = Tuple{Int,Int}[]
    for (g, exts) in final_gel_pairs
        haskey(gel_to_local, g) || continue
        lel = gel_to_local[g]
        for e in exts
            push!(my_list, (lel, e))
        end
    end

    n_dropped_by_step3 = length(global_gel_pairs) - length(final_gel_pairs)
    println("Rank $rank: final angular refinement list has $(length(my_list)) local entries " *
            "($(n_dropped_by_step3) gels dropped globally by step 3)")

    return my_list
end=#

function build_angular_refinement_list(mesh, candidates::Vector{Tuple{Int,Int}},
                                        neighbor_elem, ncf_as_child_keys::Set{Tuple{Int,Int}},
                                        ncf_as_parent_keys::Set{Tuple{Int,Int}},
                                        gel_to_class::Dict{TInt,Set{TInt}}, comm;
                                        neighbor_rank = nothing, neighbor_kind = nothing, neighbor_gel_group = nothing,
                                        buffer_hops::Int = 1)
    rank = MPI.Comm_rank(comm)

    spatial_ncf_elements = Set{Int}()
    for (iel,fid) in ncf_as_child_keys; push!(spatial_ncf_elements, iel); end
    for (iel,fid) in ncf_as_parent_keys; push!(spatial_ncf_elements, iel); end

    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)
    #=for iel in collect(spatial_ncf_elements)
        g = TInt(mesh.el2gel[iel])
        haskey(gel_to_class, g) || continue
        for partner_g in gel_to_class[g]
            haskey(gel_to_local, partner_g) || continue
            push!(spatial_ncf_elements, gel_to_local[partner_g])
        end
    end=#

    # ── GLOBAL, rank-independent adjacency graph (same-rank + cross-rank
    # edges combined) and GLOBAL refined-gel set — built ONCE, used for
    # EVERY buffer check below, so the BFS can cross rank boundaries any
    # number of times within buffer_hops, not just one hop. ─────────────────
    global_adj, global_refined_gels = build_global_element_adjacency_and_refined_set(
    	mesh, neighbor_elem, spatial_ncf_elements,
    	neighbor_rank === nothing ? Dict{Tuple{TInt,TInt},TInt}() : neighbor_rank,
    	neighbor_gel_group === nothing ? Dict{Tuple{TInt,TInt},Vector{TInt}}() : neighbor_gel_group,
    	gel_to_class, comm)

    violates(iel::Int) = violates_buffer_global(mesh, iel, global_adj, global_refined_gels, buffer_hops)

    step1_survivors = Tuple{Int,Int}[]
    for (iel, e_ext) in candidates
        if violates(iel)
            println("Rank $rank: candidate (iel=$iel,e_ext=$e_ext) REMOVED — violates one-buffer-element (step 1)")
        else
            push!(step1_survivors, (iel, e_ext))
        end
    end

    expanded_gel_pairs = Dict{TInt, Vector{Int}}()
    for (iel, e_ext) in step1_survivors
        my_gel = TInt(mesh.el2gel[iel])
        class = get(gel_to_class, my_gel, Set([my_gel]))
        for g in class
            list = get!(expanded_gel_pairs, g, Int[])
            e_ext in list || push!(list, e_ext)
        end
    end

    flat_pairs = Tuple{TInt,Int}[]
    for (g, exts) in expanded_gel_pairs, e in exts
        push!(flat_pairs, (g,e))
    end
    n_local = length(flat_pairs)
    all_n = MPI.Allgather(TInt(n_local), comm)
    flat_local = Float64[]
    for (g,e) in flat_pairs
        append!(flat_local, Float64[g,e])
    end
    counts = all_n .* TInt(2)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    global_gel_pairs = Dict{TInt, Vector{Int}}()
    for m = 1:(length(all_flat)÷2)
        g = TInt(all_flat[2m-1]); e = Int(all_flat[2m])
        list = get!(global_gel_pairs, g, Int[])
        e in list || push!(list, e)
    end

    local_violation_gels = TInt[]
    for g in keys(global_gel_pairs)
        haskey(gel_to_local, g) || continue
        lel = gel_to_local[g]
        violates(lel) && push!(local_violation_gels, g)
    end

    n_local_v = length(local_violation_gels)
    all_nv = MPI.Allgather(TInt(n_local_v), comm)
    flat_v = Float64.(local_violation_gels)
    all_vf = MPI.Allgatherv(flat_v, Cint.(all_nv), comm)
    violating_gels = Set{TInt}(TInt.(all_vf))

    final_gel_pairs = Dict{TInt, Vector{Int}}()
    for (g, exts) in global_gel_pairs
        class = get(gel_to_class, g, Set([g]))
        any(m -> m in violating_gels, class) && continue
        final_gel_pairs[g] = exts
    end

    my_list = Tuple{Int,Int}[]
    for (g, exts) in final_gel_pairs
        haskey(gel_to_local, g) || continue
        lel = gel_to_local[g]
        for e in exts
            push!(my_list, (lel, e))
        end
    end

    n_dropped_by_step3 = length(global_gel_pairs) - length(final_gel_pairs)
    println("Rank $rank: final angular refinement list has $(length(my_list)) local entries " *
            "($(n_dropped_by_step3) gels dropped globally by step 3)")

    return my_list
end

"""
Finds a candidate (iel, e_ext) such that iel ITSELF passes the step-1
buffer check (not near any spatial-NCF element), but its periodic
partner (from gel_to_class) DOES violate the buffer check — i.e., a
genuine step-3 drop case.
"""
function find_step3_test_candidate(mesh, neighbor_elem, spatial_ncf_elements::Set{Int},
                                    gel_to_class::Dict{TInt,Set{TInt}}, buffer_hops::Int=1)
    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)

    function bfs_within_hops(start_iel, max_hops)
        dist = Dict{Int,Int}(start_iel => 0)
        frontier = [start_iel]
        depth = 0
        while !isempty(frontier) && depth < max_hops
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
    violates(iel) = iel in spatial_ncf_elements || any(e -> e in spatial_ncf_elements, keys(bfs_within_hops(iel, buffer_hops)))

    for (g, class) in gel_to_class
        haskey(gel_to_local, g) || continue
        my_lel = gel_to_local[g]
        violates(my_lel) && continue   # candidate itself must PASS step 1
        for partner_g in class
            partner_g == g && continue
            haskey(gel_to_local, partner_g) || continue   # (only checking same-rank partners for simplicity)
            partner_lel = gel_to_local[partner_g]
            if violates(partner_lel)
                return (my_lel, g, partner_lel, partner_g)   # found: candidate, its partner-that-violates
            end
        end
    end
    return nothing
end

function test_build_angular_refinement_list(inputs, nparts, distribute, sem, partitioned_model;
                                              spatial_target_iel::Int = 1)
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
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)



    spatial_ncf_elements = union(Set(k[1] for k in ncf_as_child_keys), Set(k[1] for k in ncf_as_parent_keys))

    # ── Discover a genuine step-3 test case automatically ────────────────────
    step3_case = find_step3_test_candidate(mesh, neighbor_elem, spatial_ncf_elements, gel_to_class)
    if step3_case !== nothing
        cand_lel, cand_gel, viol_lel, viol_gel = step3_case
        println("Rank $rank: discovered step-3 case: candidate iel=$cand_lel (gel=$cand_gel), " *
                "periodic partner iel=$viol_lel (gel=$viol_gel) VIOLATES buffer")
    else
        println("Rank $rank: no step-3 case found locally")
    end

    # ── Build candidate list: element 1 (should be REMOVED, step 1),
    # element 500 (should PROPAGATE cleanly, step 2), plus the discovered
    # step-3 candidate if found on this rank. ────────────────────────────────
    candidates = Tuple{Int,Int}[]
    if 1 <= mesh.nelem
        push!(candidates, (1, 4))
    end
    if 500 <= mesh.nelem
        push!(candidates, (500, 4))
    end
    if step3_case !== nothing
        push!(candidates, (step3_case[1], 5))   # different e_ext, to distinguish in output
    end
    println("Rank $rank: testing candidates: $candidates")

    
    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                    ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                    neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                    neighbor_gel_group=neighbor_gel_group)
    
    println("Rank $rank: RESULT: $result")

    return result, gel_to_class, step3_case
end

function test_build_angular_refinement_list_step3(inputs, nparts, distribute, sem, partitioned_model;
                                                    spatial_target_iel_rank0::Int = 1,
                                                    spatial_target_iel_rank1::Int = 91)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0
        if rank == 0
            adapt_flags[spatial_target_iel_rank0] = 1
        elseif rank == 1
            adapt_flags[spatial_target_iel_rank1] = 1
        end
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    spatial_ncf_elements = union(Set(k[1] for k in ncf_as_child_keys), Set(k[1] for k in ncf_as_parent_keys))
    println("Rank $rank: spatial_ncf_elements = $spatial_ncf_elements")

    candidates = Tuple{Int,Int}[]
    if 1 <= mesh.nelem
        push!(candidates, (1, 4))
    end
    if 500 <= mesh.nelem
        push!(candidates, (500, 4))
    end
    println("Rank $rank: testing candidates: $candidates")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                    ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                    neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                    neighbor_gel_group=neighbor_gel_group)
    
    println("Rank $rank: RESULT: $result")

    return result
end

function build_periodic_value_exchange(mesh, neighbor_rank, neighbor_is_periodic, n_dirs, comm)
    ngl = mesh.ngl
    flat_targets = Int[]; flat_my_gel = TInt[]; flat_d = Int[]; flat_spos = Int[]
    flat_owning_iel = Int[]; flat_i = Int[]; flat_j = Int[]; flat_k = Int[]

    for ((iel,fid), target_rank) in neighbor_rank
        get(neighbor_is_periodic, (iel,fid), false) || continue
        my_gel = TInt(mesh.el2gel[iel])
        irange, jrange, krange = face_ip_range(fid, ngl)
        for d = 1:n_dirs[iel]
            spos = 0
            for k in krange, j in jrange, i in irange
                spos += 1
                push!(flat_targets, Int(target_rank)); push!(flat_my_gel, my_gel)
                push!(flat_d, d); push!(flat_spos, spos)
                push!(flat_owning_iel, iel); push!(flat_i, i); push!(flat_j, j); push!(flat_k, k)
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
            flat_owning_iel=flat_owning_iel, flat_d=flat_d, flat_i=flat_i, flat_j=flat_j, flat_k=flat_k, n_flat=n_flat)
end

function exchange_periodic_values!(exch, I_cg_all::Dict{Tuple{Int,Int},Vector{Float64}}, ngl::Int)
    send_vals = Vector{Float64}(undef, exch.n_flat)
    for m = 1:exch.n_flat
        iel = exch.flat_owning_iel[m]
        d   = exch.flat_d[m]  # NOTE: flat_d here reused for direction — matches angular exch's naming
        i,j,k = exch.flat_i[m], exch.flat_j[m], exch.flat_k[m]
        send_vals[m] = haskey(I_cg_all, (iel,d)) ? I_cg_all[(iel,d)][lidx(i,j,k,ngl)] : 0.0
    end
    send_and_receive!(exch.val_cache, send_vals)
    received = Dict{Tuple{TInt,Int,Int}, Float64}()
    for m = 1:length(exch.val_cache.combined_recv_data)
        rgel  = TInt(round(Int, exch.gel_cache.combined_recv_data[m]))
        rd    = Int(round(Int, exch.d_cache.combined_recv_data[m]))
        rspos = Int(round(Int, exch.spos_cache.combined_recv_data[m]))
        received[(rgel,rd,rspos)] = exch.val_cache.combined_recv_data[m]
    end
    return received
end

function test_periodic_value_exchange(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    nsuspect =0
    nrigid = 0
    nperiodic = 0
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
	    key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            nsuspect += 1 
	    if (key in rigid_face_keys)  nrigid += 1 end
	    if (key in periodic_face_keys) nperiodic +=1 end
        end
    end

    @info rank, nsuspect, nrigid, nperiodic, nsuspect-nrigid-nperiodic



    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    println("Rank $rank: neighbor_rank total entries: $(length(neighbor_rank))")
println("Rank $rank: neighbor_is_periodic true count: $(count(identity, values(neighbor_is_periodic)))")
println("Rank $rank: neighbor_is_periodic total entries: $(length(neighbor_is_periodic))")
println("Rank $rank: periodic_face_keys total entries: $(length(periodic_face_keys))")
println("Rank $rank: suspected_faces total entries: $(length(suspected_faces))")
println("Rank $rank: same_rank_is_periodic total entries: $(length(same_rank_is_periodic))")
@info rank, maximum(mesh.x), maximum(mesh.y), maximum(mesh.z), minimum(mesh.x), minimum(mesh.y), minimum(mesh.z)

    # ── Fake per-direction structure (pure exchange-mechanism test) ──────────
    n_dirs = fill(5, mesh.nelem)
    exch = build_periodic_value_exchange(mesh, neighbor_rank, neighbor_is_periodic, n_dirs, comm)
    println("Rank $rank: periodic exchange has $(exch.n_flat) flat entries")

    # ── Fabricated I_cg_all: encode (my_gel, d) into the value itself ────────
    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        my_gel = Int(mesh.el2gel[iel])
        for d = 1:n_dirs[iel]
            vals = fill(Float64(my_gel*1e6 + d), ngl^3)   # constant per (iel,d), encodes identity
            I_cg_all[(iel,d)] = vals
        end
    end

    received = exchange_periodic_values!(exch, I_cg_all, ngl)
    println("Rank $rank: received $(length(received)) periodic entries")

    n_checked = 0; n_mismatch = 0
    for ((rgel,rd,rspos), val) in received
        n_checked += 1
        decoded_gel = Int(floor(val / 1e6))
        decoded_d = Int(round(val - decoded_gel*1e6))
        if decoded_gel != rgel || decoded_d != rd
            n_mismatch += 1
            println("Rank $rank: ✗ MISMATCH (rgel=$rgel,rd=$rd,rspos=$rspos): decoded=($decoded_gel,$decoded_d)")
        end
    end
    println("Rank $rank: checked $n_checked entries, $n_mismatch mismatches  " *
            (n_mismatch==0 ? "✓ PASS" : "✗ FAIL"))

    return received
end

function test_periodic_value_exchange_both(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    mesh = sem.mesh
    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child_keys = Set{Tuple{Int,Int}}()
    ncf_as_parent_keys = Set{Tuple{Int,Int}}()
    for ncf_idx = 1:mesh.num_ncf
        push!(ncf_as_child_keys, (Int(mesh.cip[ncf_idx]), Int(mesh.lfid[ncf_idx])))
        push!(ncf_as_parent_keys, (Int(mesh.pip[ncf_idx]), OPP_FACE[Int(mesh.lfid[ncf_idx])]))
    end

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    println("Rank $rank: same-rank periodic entries: $(count(identity,values(same_rank_is_periodic)))")
    println("Rank $rank: cross-rank periodic entries: $(count(identity,values(neighbor_is_periodic)))")

    n_dirs = fill(5, mesh.nelem)

    # ── Fabricated I_cg_all: encode (my_gel, d) into the value itself ────────
    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        my_gel = Int(mesh.el2gel[iel])
        for d = 1:n_dirs[iel]
            I_cg_all[(iel,d)] = fill(Float64(my_gel*1e6 + d), ngl^3)
        end
    end

    # ── CROSS-RANK check (via the exchange) ──────────────────────────────────
    exch = build_periodic_value_exchange(mesh, neighbor_rank, neighbor_is_periodic, n_dirs, comm)
    println("Rank $rank: cross-rank periodic exchange has $(exch.n_flat) flat entries")
    received = exchange_periodic_values!(exch, I_cg_all, ngl)
    n_checked_cr = 0; n_mismatch_cr = 0
    for ((rgel,rd,rspos), val) in received
        n_checked_cr += 1
        decoded_gel = Int(floor(val / 1e6)); decoded_d = Int(round(val - decoded_gel*1e6))
        (decoded_gel != rgel || decoded_d != rd) && (n_mismatch_cr += 1)
    end
    println("Rank $rank: CROSS-RANK: checked $n_checked_cr, $n_mismatch_cr mismatches  " *
            (n_mismatch_cr==0 ? "✓ PASS" : "✗ FAIL"))

    # ── SAME-RANK check (direct, no MPI needed — both sides in local memory) ─
    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)
    n_checked_sr = 0; n_mismatch_sr = 0
    for ((iel,fid), kind) in same_rank_kind
        get(same_rank_is_periodic, (iel,fid), false) || continue
        partner_lel = same_rank_partner[(iel,fid)][1]
        for d = 1:n_dirs[iel]
            n_checked_sr += 1
            expected_val = Float64(Int(mesh.el2gel[partner_lel])*1e6 + d)
            actual_val = I_cg_all[(partner_lel,d)][1]   # direct local lookup — this IS what the sweep would read
            actual_val != expected_val && (n_mismatch_sr += 1)
        end
    end
    println("Rank $rank: SAME-RANK: checked $n_checked_sr, $n_mismatch_sr mismatches  " *
            (n_mismatch_sr==0 ? "✓ PASS" : "✗ FAIL"))

    return received, I_cg_all
end

function run_dg_sweep_multirank_combined_ncf_periodic!(mesh, metrics, ψ, dψ, ω,
                                                        neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
                                                        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                        neighbor_is_periodic,
                                                        same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                        ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                        coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                        angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                        κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                                        extra_source_all = nothing,
                                                        I_cg_all_prev = nothing,       # PREVIOUS outer iter's solution
                                                        periodic_received = nothing,    # cross-rank periodic snapshot
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

    while true #n_solved_total < total_pairs
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

                    # ── PERIODIC CHECK — FIRST, before ordinary dispatch.
                    # Always ready: reads a FROZEN prior-outer-iteration
                    # snapshot, never blocks. ──────────────────────────────────
                    if get(same_rank_is_periodic, (iel,fid), false)
    			partner_lel = same_rank_partner[(iel,fid)][1]
    			fid_c = fid
    			if I_cg_all_prev !== nothing && haskey(I_cg_all_prev, (partner_lel,d))
        		    partner_vals = I_cg_all_prev[(partner_lel,d)]
        		    upwind_fns[fid] = (i,j,k) -> begin
            			corr = periodic_face_correspondence[(iel,fid_c)][(i,j,k)]
            			partner_vals[lidx(corr[1],corr[2],corr[3],ngl)]
        		    end
    			else
        		    upwind_fns[fid] = (i,j,k) -> 0.0
    			end
    			continue
		    end
		    if get(neighbor_is_periodic, (iel,fid), false)
                        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
                        fid_c = fid
                        if periodic_received !== nothing
                            upwind_fns[fid] = (i,j,k) -> begin
                                spos = face_flat_index(fid_c, ngl, i, j, k)
                                get(periodic_received, (remote_gel_p, d, spos), 0.0)
                            end
                        else
                            upwind_fns[fid] = (i,j,k) -> 0.0
                        end
                        continue
                    end

                    # ── ORDINARY / NCF dispatch — UNCHANGED from before ────────
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
					    fid_cc=fid
                                            upwind_fns[fid] = (i,j,k) -> received_ang[(remote_gel,rn_c,fid_cc,face_flat_index(fid_cc,ngl,i,j,k))]
                                        end
                                    else
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos))
                                                           for rn in remote_nodes for spos=1:ngl2)
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

	n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)
	if iter_round % 1 == 0   # print EVERY round for now, given the uncertainty
    		println("Rank $rank: round $iter_round, solved $n_solved_total/$total_pairs " *
            	"(newly_solved this round: $n_newly_solved)")
	end

        received_ang = exchange_ncf_angular_values!(angular_exch, x_dg, ngl)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break

        # Fix inside run_dg_sweep_multirank_combined_ncf_periodic!'s stall branch:
	if n_newly_solved_global == 0
    	    just_resolved_stall && error("no progress after GS resolution")
    	    resolve_stall_distributed_gauss_seidel_combined_periodic!(mesh, metrics, ψ, dψ, ω,
        	neighbor_elem, face_correspondence, periodic_face_correspondence,
        	neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
        	same_rank_kind, same_rank_partner, same_rank_is_periodic,
        	ncf_as_child, ncf_L1L2, ncf_as_parent,
        	coord_cache, gelpos_to_gip, spatial_perdir_exch,
        	angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
        	x_dg, solved, n_dirs, κ, σ, intΦ_eff, bc_fn_all, comm;
        	extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev, periodic_received=periodic_received,
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

function resolve_stall_distributed_gauss_seidel_combined_periodic!(mesh, metrics, ψ, dψ, ω,
                                                                     neighbor_elem, face_correspondence,
                                                                     periodic_face_correspondence,
								     neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                                     neighbor_is_periodic,
                                                                     same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                                     ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                                     coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                                     angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                                     x_dg, solved, n_dirs,
                                                                     κ, σ, intΦ_eff, bc_fn_all, comm;
                                                                     extra_source_all = nothing,
                                                                     I_cg_all_prev = nothing,
                                                                     periodic_received = nothing,
                                                                     max_gs_iter::Int = 100, gs_tol::Float64 = 1e-9)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    n_candidates_total = MPI.Allreduce(length(candidates), +, comm)
    println("Rank $rank: distributed GS-combined-periodic over $(length(candidates)) local candidates " *
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
                    # ── PERIODIC CHECK — FIRST, same as the main sweep ────────
                    
		    if get(same_rank_is_periodic, (iel,fid), false)
    			partner_lel = same_rank_partner[(iel,fid)][1]
    			fid_c = fid
    			if I_cg_all_prev !== nothing && haskey(I_cg_all_prev, (partner_lel,d))
        		    partner_vals = I_cg_all_prev[(partner_lel,d)]
        		    I_upwind = (i,j,k) -> begin
            		        corr = periodic_face_correspondence[(iel,fid_c)][(i,j,k)]
            			partner_vals[lidx(corr[1],corr[2],corr[3],ngl)]
        		    end
    		        end
		    elseif get(neighbor_is_periodic, (iel,fid), false)
                        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
                        fid_c = fid
                        if periodic_received !== nothing
                            I_upwind = (i,j,k) -> begin
                                spos = face_flat_index(fid_c, ngl, i, j, k)
                                get(periodic_received, (remote_gel_p, d, spos), 0.0)
                            end
                        end
                    else
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
                                        if all(haskey(received_ang,(remote_gel_ord,d,spos)) for spos=1:ngl2)
                                            gel_c=remote_gel_ord; d_c=d
					    fid_cc=fid
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
        println("Rank $rank: GS-combined-periodic iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-combined-periodic did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

function setup_combined_ncf_mpi_test_periodic(inputs, nparts, distribute, sem, partitioned_model;
                                               spatial_target_iel::Int = 1,
                                               angular_target_iel::Int = 442, angular_target_e_ext::Int = 4)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0 && spatial_target_iel > 0 && inputs[:ladapt] 
        adapt_flags[spatial_target_iel] = 1

    end
   
    if sem.mesh.nelem > 0 && spatial_target_iel > 0 && inputs[:ladapt]
        sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)

    end

 
    mesh = sem.mesh; metrics = sem.metrics

    if angular_target_iel <= mesh.nelem && angular_target_iel > 0 && rank == 0
	@info "target gel", mesh.el2gel[angular_target_iel] 
	adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
        verify_angular_phi_seam_conformity(mesh, angular_target_iel)
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: $(length(neighbor_rank)) resolved cross-rank, $(length(same_rank_kind)) same-rank, unresolved=$n_unresolved")
    n_unresolved > 0 && error("unresolved faces")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh; trace_iel=1, trace_fid=3, trace_d=106)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    debug_dump_all_angular_deps_for_element(mesh, angular_deps, angular_deps_cr, 511, comm)

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)

    return (mesh=mesh, metrics=metrics, sem=sem, partitioned_model=partitioned_model,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid, face_correspondence=face_correspondence,
            ncf_as_child=ncf_as_child, ncf_as_parent=ncf_as_parent, ncf_L1L2=ncf_L1L2,
            neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind, neighbor_gel_group=neighbor_gel_group,
            neighbor_gel=neighbor_gel, neighbor_is_periodic=neighbor_is_periodic, periodic_face_correspondence=periodic_face_correspondence,
            same_rank_kind=same_rank_kind, same_rank_partner=same_rank_partner, same_rank_is_periodic=same_rank_is_periodic,
            coord_cache=coord_cache, gelpos_to_gip=gelpos_to_gip,
            remote_angular_mesh=remote_angular_mesh, angular_deps=angular_deps, angular_deps_cr=angular_deps_cr)
end

function test_manufactured_isotropic_combined_ncf_multirank_periodic(inputs, nparts, distribute, sem, partitioned_model;
                                                                       spatial_target_iel::Int = 0,   # 0 = NO spatial refinement
                                                                       angular_target_iel::Int = 0,   # 0 = NO angular refinement
                                                                       angular_target_e_ext::Int = 4,
                                                                       I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                                       κ_val=0.001, σ_val=0.0005,
                                                                       max_iter=500, tol=1e-9)
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

    angular_exch = build_ncf_angular_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, setup.neighbor_rank, n_dirs, comm)
    periodic_exch = build_periodic_value_exchange(mesh, setup.neighbor_rank, setup.neighbor_is_periodic, n_dirs, comm)

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
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
            bc_fn_all[(iel,dk)] = (x,y,z) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)   # only used at TRUE boundaries now
        end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        I_cg_all[(iel,d)] = zeros(ngl^3)
    end
    periodic_received = nothing


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

        # ── Refresh periodic cross-rank snapshot from the PREVIOUS iteration's
        # converged I_cg_all — same-rank reads I_cg_all DIRECTLY inside the
        # sweep, no exchange needed. ─────────────────────────────────────────
        periodic_received = exchange_periodic_values!(periodic_exch, I_cg_all, ngl)

        x_dg, _ = run_dg_sweep_multirank_combined_ncf_periodic!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group,
            setup.neighbor_is_periodic, setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.coord_cache, setup.gelpos_to_gip, spatial_perdir_exch,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh, angular_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all, periodic_received=periodic_received)

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
    rank == 0 && println("Multirank COMBINED NCF PERIODIC manufactured solution — max abs error: $max_err_global")
    return max_err_global
end

"""
    build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)

Same-rank periodic analogue of face_correspondence's construction: maps
THIS element's local (i,j,k) at a periodic face to its periodic partner's
local (i,j,k) — matched via SHARED GLOBAL ID (mesh.ip2gip), since a
periodic pair's two faces never share a literal connijk node identity
(they're physically distinct locations), only a periodic gid pairing.
"""
function build_periodic_face_correspondence_same_rank(mesh, same_rank_kind::Dict{Tuple{Int,Int},Symbol},
                                                        same_rank_partner::Dict{Tuple{Int,Int},Vector{Int}},
                                                        ngl::Int)
    periodic_face_correspondence = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int,Int},Tuple{Int,Int,Int}}}()

    for ((iel,fid), kind) in same_rank_kind
        kind == :conforming || continue   # periodicity is always :conforming, per our policy
        jel = same_rank_partner[(iel,fid)][1]

        # Find jel's own local fid for this SAME periodic relationship —
        # by checking which of jel's faces maps back to iel.
        jfid = 0
        for fid2 = 1:6
            haskey(same_rank_kind, (jel,fid2)) || continue
            same_rank_kind[(jel,fid2)] == :conforming || continue
            same_rank_partner[(jel,fid2)] == [iel] || continue
            jfid = fid2; break
        end
        jfid == 0 && continue   # no reciprocal found — shouldn't happen given earlier mutual-consistency checks

        irange_i,jrange_i,krange_i = face_ip_range(fid, ngl)
        irange_j,jrange_j,krange_j = face_ip_range(jfid, ngl)

        jel_gid_to_ijk = Dict{Int,Tuple{Int,Int,Int}}()
        for k in krange_j, j in jrange_j, i in irange_j
            g = mesh.ip2gip[mesh.connijk[jel,i,j,k]]   # KEY CHANGE: gid, not raw connijk
            jel_gid_to_ijk[g] = (i,j,k)
        end

        lmap = Dict{Tuple{Int,Int,Int},Tuple{Int,Int,Int}}()
        for k in krange_i, j in jrange_i, i in irange_i
            g = mesh.ip2gip[mesh.connijk[iel,i,j,k]]   # KEY CHANGE: gid, not raw connijk
            haskey(jel_gid_to_ijk, g) && (lmap[(i,j,k)] = jel_gid_to_ijk[g])
        end
        periodic_face_correspondence[(iel,fid)] = lmap
    end

    return periodic_face_correspondence
end

function test_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)

    n_faces_checked = 0
    n_incomplete = 0
    n_non_bijective = 0
    for ((iel,fid), kind) in same_rank_kind
        kind == :conforming || continue
        haskey(periodic_face_correspondence, (iel,fid)) || (n_incomplete += 1; continue)
        lmap = periodic_face_correspondence[(iel,fid)]
        irange,jrange,krange = face_ip_range(fid, ngl)
        expected_n = length(irange)*length(jrange)*length(krange)
        n_faces_checked += 1
        if length(lmap) != expected_n
            n_incomplete += 1
            println("  (iel=$iel,fid=$fid): expected $expected_n mapped nodes, got $(length(lmap))")
        end
        targets = collect(values(lmap))
        length(unique(targets)) != length(targets) && (n_non_bijective += 1)
    end
    println("Checked $n_faces_checked periodic same-rank faces: $n_incomplete incomplete, $n_non_bijective non-bijective  " *
            ((n_incomplete==0 && n_non_bijective==0) ? "✓ PASS" : "✗ FAIL"))
    return periodic_face_correspondence
end

function generate_stress_test_targets(nelem_base::Int; seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4)
    rng = MersenneTwister(seed)
    n_spatial = round(Int, frac * nelem_base)
    n_angular = round(Int, frac * nelem_base)


    target_gels = sort(randperm(rng, nelem_base)[1:n_spatial+n_angular])
	
    spatial_target_gels = target_gels[1:n_spatial]

    # Continue drawing from the SAME rng stream — independent draw over the
    # FULL range (overlaps with spatial_target_gels allowed; rare at this
    # scale, and harmless since they're trivially rejected by the
    # one-buffer-element check anyway).
    angular_target_gels = target_gels[n_spatial:n_spatial+n_angular] #sort(randperm(rng, nelem_base)[1:n_angular])
    angular_target_e_exts = [rand(rng, 1:max_e_ext) for _ in angular_target_gels]

    @info frac, length(target_gels), length(spatial_target_gels), length(angular_target_gels)

    return spatial_target_gels, angular_target_gels, angular_target_e_exts
end

function test_stress_combined_ncf_periodic_multirank(inputs, nparts, distribute, sem, partitioned_model;
                                                       seed::Int=42, frac::Float64=0.010, max_e_ext::Int=4,
                                                       I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                       κ_val=0.001, σ_val=0.0005,
                                                       max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    nelem_base = sem.mesh.gnelem   # NOTE: exact field name TBD — see caveat below
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(nelem_base; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    # ── Apply spatial refinement: translate GLOBAL target gels to LOCAL
    # adapt_flags on THIS rank's pre-refinement mesh (gels == local ids
    # here, since refinement hasn't happened yet — the base mesh's el2gel
    # should be the identity locally-to-globally at this stage... — NEED
    # TO CONFIRM this assumption holds for the BASE, unrefined partition). ──
    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    # ── Build the full correspondence/exchange stack (spatial NCF +
    # periodic), same as setup_combined_ncf_mpi_test_periodic but WITHOUT
    # angular refinement yet — that comes next, via build_angular_refinement_list. ──
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end


    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    debug_multi_face_angular_ncf_elements(mesh, neighbor_kind, comm)

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)

    # ── Periodic equivalence classes, needed for build_angular_refinement_list ──
    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    # ── Translate angular target gels -> LOCAL candidates on THIS rank ────
    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        idx = findfirst(==(gel), angular_target_gels)
        idx === nothing && continue
        push!(local_angular_candidates, (iel, angular_target_e_exts[idx]))
    end
    println("Rank $rank: $(length(local_angular_candidates)) local angular candidates before buffer/periodic checks")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                    ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                    neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                    neighbor_gel_group=neighbor_gel_group)
    
    println("Rank $rank: $(length(final_angular_list)) angular refinements to APPLY after buffer/periodic checks")

    for (iel, e_ext) in final_angular_list
        adapt_single_angular_element_conforming!(mesh, iel, e_ext)
        verify_angular_phi_seam_conformity(mesh, iel)
    end

    # ── Cross-rank angular/spatial exchange setup — SAME as every other
    # combined test, built AFTER all refinement (spatial + angular) is done. ──
    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    angular_exch = build_ncf_angular_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    periodic_exch = build_periodic_value_exchange(mesh, neighbor_rank, neighbor_is_periodic, n_dirs, comm)

    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    # ── Manufactured solution outer loop — IDENTICAL to the smaller test ────
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
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
    periodic_received = nothing

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

        periodic_received = exchange_periodic_values!(periodic_exch, I_cg_all, ngl)

        #=if it == 1
    	    x_dg, _, _ = run_dg_sweep_multirank_combined_ncf_periodic_instrumented!(mesh, metrics, ψ, dψ, ω,
        	neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
        	neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
        	neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
        	ncf_as_child, ncf_L1L2, ncf_as_parent,
        	coord_cache, gelpos_to_gip, spatial_perdir_exch,
        	angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
        	κ, σz, intΦ_eff, bc_fn_all, comm;
        	extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all, periodic_received=periodic_received,
        	exact_fn=I_manuf, err_threshold=1.0)
	else=#


	    x_dg, _ = run_dg_sweep_multirank_combined_ncf_periodic!(mesh, metrics, ψ, dψ, ω,
            	neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
            	neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
            	neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
            	ncf_as_child, ncf_L1L2, ncf_as_parent,
            	coord_cache, gelpos_to_gip, spatial_perdir_exch,
            	angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
            	κ, σz, intΦ_eff, bc_fn_all, comm;
            	extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all, periodic_received=periodic_received)
	#end

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
    rank == 0 && println("STRESS TEST combined NCF+periodic multirank — max abs error: $max_err_global")
    return max_err_global
end

function debug_unresolved_face(mesh, suspected_faces, comm; target_iel::Int, target_fid::Int)
    rank = MPI.Comm_rank(comm)

    local_keys = NTuple{4,TInt}[]; local_iel = TInt[]; local_fid = TInt[]; local_gel = TInt[]
    for (iel, fid) in suspected_faces
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
        key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        push!(local_keys, key); push!(local_iel, iel); push!(local_fid, fid); push!(local_gel, mesh.el2gel[iel])
    end

    n_local = TInt(length(local_keys))
    all_n = MPI.Allgather(n_local, comm)
    flat_local = Float64[]
    for i = eachindex(local_keys)
        k = local_keys[i]
        append!(flat_local, Float64[k[1],k[2],k[3],k[4], rank, local_iel[i], local_fid[i], local_gel[i]])
    end
    counts = all_n .* TInt(8)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)
    n_total = length(all_flat) ÷ 8

    my_idx = findfirst(i -> local_iel[i]==target_iel && local_fid[i]==target_fid, eachindex(local_iel))
    my_idx === nothing && (println("Rank $rank: target not in MY suspected_faces"); return)
    my_key = local_keys[my_idx]
    println("Rank $rank: target (iel=$target_iel,fid=$target_fid) key=$my_key gel=$(local_gel[my_idx])")

    println("Rank $rank: searching ALL $n_total gathered entries for this key:")
    for m = 1:n_total
        b = (m-1)*8
        key = (TInt(all_flat[b+1]),TInt(all_flat[b+2]),TInt(all_flat[b+3]),TInt(all_flat[b+4]))
        key == my_key || continue
        println("  MATCH FOUND: rank=$(Int(all_flat[b+5])) iel=$(Int(all_flat[b+6])) fid=$(Int(all_flat[b+7])) gel=$(Int(all_flat[b+8]))")
    end
end

function normal_axis_of_face(corners)
    xs = [c[1] for c in corners]; ys = [c[2] for c in corners]; zs = [c[3] for c in corners]
    spread = (maximum(xs)-minimum(xs), maximum(ys)-minimum(ys), maximum(zs)-minimum(zs))
    return argmin(spread)
end

"""
Returns (in-plane bbox with excluded axis set to -Inf/Inf for the
CONTAINMENT check, excluded axis index, mean coordinate along the
excluded/normal axis) — the normal value is used SEPARATELY to
disambiguate real neighbors (same normal value) from unrelated faces
that coincidentally share the same in-plane footprint.
"""
function bbox_and_normal(corners)
    excl = normal_axis_of_face(corners)
    xs = [c[1] for c in corners]; ys = [c[2] for c in corners]; zs = [c[3] for c in corners]
    xr = excl==1 ? (-Inf,Inf) : (minimum(xs),maximum(xs))
    yr = excl==2 ? (-Inf,Inf) : (minimum(ys),maximum(ys))
    zr = excl==3 ? (-Inf,Inf) : (minimum(zs),maximum(zs))
    normal_val = excl==1 ? sum(xs)/4 : (excl==2 ? sum(ys)/4 : sum(zs)/4)
    return (xr[1],xr[2], yr[1],yr[2], zr[1],zr[2]), excl, normal_val
end

function resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces::Vector{Tuple{TInt,TInt}},
                                                        periodic_face_keys::Set{NTuple{4,TInt}},
                                                        comm; tol::Float64 = 1e-6)
    rank = MPI.Comm_rank(comm)

    local_keys   = NTuple{4,TInt}[]
    local_coords = Vector{NTuple{4,NTuple{3,Float64}}}()
    local_iel = TInt[]; local_fid = TInt[]; local_gel = TInt[]
    local_is_periodic = Bool[]

    for (iel, fid) in suspected_faces
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
        raw_key = Tuple(sort(TInt[n1, n2, n3, n4]))
        gid_key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        coords = ((mesh.x[n1],mesh.y[n1],mesh.z[n1]), (mesh.x[n2],mesh.y[n2],mesh.z[n2]),
                  (mesh.x[n3],mesh.y[n3],mesh.z[n3]), (mesh.x[n4],mesh.y[n4],mesh.z[n4]))
        push!(local_keys, gid_key); push!(local_coords, coords)
        push!(local_iel, iel); push!(local_fid, fid); push!(local_gel, mesh.el2gel[iel])
        push!(local_is_periodic, raw_key in periodic_face_keys)
    end

    n_local = TInt(length(local_keys))
    all_n = MPI.Allgather(n_local, comm)
    flat_local = Float64[]
    for i = eachindex(local_keys)
        k = local_keys[i]; c = local_coords[i]
        append!(flat_local, Float64[k[1],k[2],k[3],k[4]])
        for corner in c
            append!(flat_local, Float64[corner[1],corner[2],corner[3]])
        end
        append!(flat_local, Float64[rank, local_iel[i], local_fid[i], local_gel[i],
                                     local_is_periodic[i] ? 1.0 : 0.0])
    end
    counts = all_n .* TInt(21)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    n_total = length(all_flat) ÷ 21
    entries = Vector{NamedTuple}(undef, n_total)
    for m = 1:n_total
        b = (m-1)*21
        key = (TInt(all_flat[b+1]),TInt(all_flat[b+2]),TInt(all_flat[b+3]),TInt(all_flat[b+4]))
        corners = ((all_flat[b+5],all_flat[b+6],all_flat[b+7]), (all_flat[b+8],all_flat[b+9],all_flat[b+10]),
                   (all_flat[b+11],all_flat[b+12],all_flat[b+13]), (all_flat[b+14],all_flat[b+15],all_flat[b+16]))
        entries[m] = (key=key, corners=corners, rank=TInt(all_flat[b+17]), iel=TInt(all_flat[b+18]),
                      fid=TInt(all_flat[b+19]), gel=TInt(all_flat[b+20]), is_periodic=(all_flat[b+21] > 0.5))
    end

    neighbor_rank = Dict{Tuple{TInt,TInt}, TInt}()
    neighbor_kind = Dict{Tuple{TInt,TInt}, Symbol}()
    neighbor_gel_group = Dict{Tuple{TInt,TInt}, Vector{TInt}}()
    neighbor_is_periodic = Dict{Tuple{TInt,TInt}, Bool}()
    same_rank_kind = Dict{Tuple{Int,Int}, Symbol}()
    same_rank_partner = Dict{Tuple{Int,Int}, Vector{Int}}()
    same_rank_is_periodic = Dict{Tuple{Int,Int}, Bool}()
    n_unresolved = 0

    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)

    for i = eachindex(local_keys)
        my_key = local_keys[i]; my_iel = local_iel[i]; my_fid = local_fid[i]
        my_is_periodic = local_is_periodic[i]
        my_bbox, my_excl, my_normal = bbox_and_normal(local_coords[i])

        exact_matches = [e for e in entries if e.key == my_key && e.gel != local_gel[i]]
        if length(exact_matches) == 1
            m = exact_matches[1]
            if m.rank == rank
                lel = get(gel_to_local, m.gel, 0)
                lel == 0 && (n_unresolved += 1; continue)
                same_rank_kind[(my_iel,my_fid)] = :conforming
                same_rank_partner[(my_iel,my_fid)] = [lel]
                same_rank_is_periodic[(my_iel,my_fid)] = my_is_periodic
            else
                neighbor_rank[(my_iel,my_fid)] = m.rank
                neighbor_kind[(my_iel,my_fid)] = :conforming
                neighbor_gel_group[(my_iel,my_fid)] = [m.gel]
                neighbor_is_periodic[(my_iel,my_fid)] = my_is_periodic
            end
            continue
        end

        # ── Containment probes — NOW with normal-axis disambiguation:
        # a candidate only counts if (a) in-plane bbox contains/is-contained,
        # AND (b) EITHER it's a genuine periodic partner (normal axis allowed
        # to differ by domain length) OR its normal-axis value actually
        # COINCIDES with mine (confirming same physical plane — rules out
        # unrelated, coincidentally-aligned faces elsewhere in the mesh). ────
        candidates_inside = NamedTuple[]; candidates_containing = NamedTuple[]
        for e in entries
            e.gel == local_gel[i] && continue
            ebbox, e_excl, e_normal = bbox_and_normal(e.corners)
            e_excl == my_excl || continue   # must share the SAME normal axis to even be comparable

            # Disambiguation: for a NON-periodic candidate, require the
            # normal-axis values to coincide (same physical plane).
            plausible_plane = my_is_periodic || (e_normal !== nothing && isapprox(my_normal, e_normal; atol=tol*100))
            plausible_plane || continue

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
            all_same_rank = all(c.rank == rank for c in candidates_inside)
            if all_same_rank
                lels = [get(gel_to_local, c.gel, 0) for c in candidates_inside]
                any(==(0), lels) ? (n_unresolved += 1) : begin
                    same_rank_kind[(my_iel,my_fid)] = :i_am_coarse
                    same_rank_partner[(my_iel,my_fid)] = lels
                    same_rank_is_periodic[(my_iel,my_fid)] = my_is_periodic
                end
            elseif all(c.rank == candidates_inside[1].rank for c in candidates_inside)
                neighbor_rank[(my_iel,my_fid)] = candidates_inside[1].rank
                neighbor_kind[(my_iel,my_fid)] = :i_am_coarse
                neighbor_gel_group[(my_iel,my_fid)] = [c.gel for c in candidates_inside]
                neighbor_is_periodic[(my_iel,my_fid)] = my_is_periodic
            else
                n_unresolved += 1
            end
        elseif length(candidates_containing) == 1
            c = candidates_containing[1]
            if c.rank == rank
                lel = get(gel_to_local, c.gel, 0)
                lel == 0 ? (n_unresolved += 1) : begin
                    same_rank_kind[(my_iel,my_fid)] = :i_am_fine
                    same_rank_partner[(my_iel,my_fid)] = [lel]
                    same_rank_is_periodic[(my_iel,my_fid)] = my_is_periodic
                end
            else
                neighbor_rank[(my_iel,my_fid)] = c.rank
                neighbor_kind[(my_iel,my_fid)] = :i_am_fine
                neighbor_gel_group[(my_iel,my_fid)] = [c.gel]
                neighbor_is_periodic[(my_iel,my_fid)] = my_is_periodic
            end
        else
            n_unresolved += 1
        end
    end

    return neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
           same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved
end

function debug_multi_face_angular_ncf_elements(mesh, neighbor_kind, comm)
    rank = MPI.Comm_rank(comm)

    # Group faces by their OWNING local element, counting how many
    # NON-conforming (genuine angular-NCF-relevant) faces each element has.
    per_element_faces = Dict{Int, Vector{Int}}()
    for ((iel,fid), kind) in neighbor_kind
        push!(get!(per_element_faces, iel, Int[]), fid)
    end

    n_multi = 0
    for (iel, fids) in per_element_faces
        length(fids) > 1 || continue
        n_multi += 1
        n_multi <= 10 && println("Rank $rank: iel=$iel has $(length(fids)) cross-rank faces: fids=$fids " *
                                   "kinds=$([neighbor_kind[(iel,f)] for f in fids])")
    end
    println("Rank $rank: $n_multi local elements have MORE THAN ONE cross-rank face " *
            "(out of $(length(per_element_faces)) total cross-rank-touching elements)")

    return n_multi
end


function run_dg_sweep_multirank_combined_ncf_periodic_instrumented!(mesh, metrics, ψ, dψ, ω,
                                                                      neighbor_elem, neighbor_fid, face_correspondence,
                                                                      periodic_face_correspondence,
                                                                      neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                                                                      neighbor_is_periodic,
                                                                      same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                                      ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                                      coord_cache, gelpos_to_gip, spatial_perdir_exch,
                                                                      angular_deps, angular_deps_cr, remote_angular_mesh, angular_exch,
                                                                      κ, σ, intΦ_eff, bc_fn_all, comm::MPI.Comm;
                                                                      extra_source_all = nothing,
                                                                      I_cg_all_prev = nothing,
                                                                      periodic_received = nothing,
                                                                      exact_fn = nothing,
                                                                      err_threshold::Float64 = 1.0,
                                                                      block_periodic::Bool = true,
                                                                      max_rounds::Int = 500)
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl
    ngl2 = ngl*ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    iter_round = 0
    received_ang = Dict{Tuple{TInt,Int,Int,Int}, Float64}()
    received_sp  = Dict{Tuple{TInt,Int,Int}, Float64}()

    find_owning_rank_and_local_iel(mesh, 678, comm)

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
                upwind_source = Dict{Int,String}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue

                    if get(same_rank_is_periodic, (iel,fid), false)
                        if block_periodic
                            ready = false; break
                        end
                        partner_lel = same_rank_partner[(iel,fid)][1]
                        fid_c = fid
                        if I_cg_all_prev !== nothing && haskey(I_cg_all_prev, (partner_lel,d))
                            partner_vals = I_cg_all_prev[(partner_lel,d)]
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = periodic_face_correspondence[(iel,fid_c)][(i,j,k)]
                                partner_vals[lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        else
                            upwind_fns[fid] = (i,j,k) -> 0.0
                        end
                        upwind_source[fid] = "PERIODIC-SAME-RANK"
                        continue
                    end
                    if get(neighbor_is_periodic, (iel,fid), false)
                        if block_periodic
                            ready = false; break
                        end
                        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
                        fid_c = fid
                        if periodic_received !== nothing
                            upwind_fns[fid] = (i,j,k) -> begin
                                spos = face_flat_index(fid_c, ngl, i, j, k)
                                get(periodic_received, (remote_gel_p, d, spos), 0.0)
                            end
                        else
                            upwind_fns[fid] = (i,j,k) -> 0.0
                        end
                        upwind_source[fid] = "PERIODIC-CROSS-RANK"
                        continue
                    end

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
                                            upwind_source[fid] = "CROSS-RANK-ANGULAR-NCF-SINGLETON"
                                        end
                                    else
                                        fully_avail = all(haskey(received_ang,(remote_gel,rn,fid,spos))
                                                           for rn in remote_nodes for spos=1:ngl2)
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
                                    remote_gel_ord = neighbor_gel[(iel,fid)]
                                    fully_avail = all(haskey(received_ang, (remote_gel_ord, d, fid, spos)) for spos=1:ngl2)
                                    if !fully_avail
                                        ready = false; break
                                    else
                                        gel_c=remote_gel_ord; d_c=d
					fid_cc=fid
                                        upwind_fns[fid] = (i,j,k) -> received_ang[(gel_c,d_c,fid_cc,face_flat_index(fid,ngl,i,j,k))]
                                        upwind_source[fid] = "CROSS-RANK-ORDINARY"
                                    end
                                end
                            end
                        end
                        haskey(upwind_fns,fid) || begin
                            bc = bc_fn_all[(iel,d)]
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

                # ── COMPLETENESS CHECK — same diagnostic that found the
                # ORIGINAL missing-branch bug in this whole effort. ─────────
                missing_fids = Int[]
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    haskey(upwind_fns, fid) || push!(missing_fids, fid)
                end
                if !isempty(missing_fids)
                    println("Rank $rank ⚠⚠⚠ INCOMPLETE upwind_fns for (iel=$iel,d=$d): " *
                             "is_active=$is_active is_inflow=$is_inflow " *
                             "MISSING closures for active+inflow fid(s)=$missing_fids " *
                             "(SILENT ZERO fallback in build_dg_flux_terms!)")
                end

                if exact_fn !== nothing
                    for fid = 1:6
                        haskey(upwind_fns,fid) || continue
                        irange,jrange,krange = face_ip_range(fid, ngl)
                        for (i,j,k) in Iterators.product(irange,jrange,krange)
                            val = upwind_fns[fid](i,j,k)
                            gip = mesh.connijk[iel,i,j,k]
                            ex = exact_fn(mesh.x[gip],mesh.y[gip],mesh.z[gip],Ωx,Ωy,Ωz)
                            if abs(val-ex) > err_threshold
                                if (iel == 328 && d == 29 && fid ==6 && i == 5 && j==4 && k==1)
				    debug_cross_rank_ordinary_lookup(mesh, neighbor_gel, received_ang, x_dg, iel, fid, d, ngl,
                                            i, j, k)
				end
				println("Rank $rank ⚠⚠⚠ LARGE I_upwind ERROR: (iel=$iel,d=$d,fid=$fid,i=$i,j=$j,k=$k) " *
                                        "source=$(upwind_source[fid]) val=$val exact=$ex")
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
                            println("Rank $rank ⚠⚠⚠ LARGE SOLVE ERROR: (iel=$iel,d=$d,i=$i,j=$j,k=$k) " *
                                    "val=$val exact=$ex  [sources=$upwind_source]")
                        end
                    end
                end
            end
        end
        n_solved_total = sum(count(identity,solved[iel]) for iel=1:mesh.nelem)

        received_ang = exchange_ncf_angular_values_traced!(mesh,angular_exch, x_dg, ngl; trace_gel=678, trace_d=29, trace_spos=4)
        received_sp  = exchange_ncf_spatial_perdir_values!(spatial_perdir_exch, x_dg, ngl)

        n_unprocessed_total = MPI.Allreduce(total_pairs-n_solved_total, +, comm)
        n_unprocessed_total == 0 && break
        n_newly_solved_global = MPI.Allreduce(n_newly_solved, +, comm)
        if n_newly_solved_global == 0
            println("Rank $rank: STALLED at $(total_pairs-n_solved_total) unsolved — stopping instrumented run HERE")
            break
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end
    return x_dg, n_dirs, solved
end

function find_owning_rank_and_local_iel(mesh, target_gel::Int, comm)
    rank = MPI.Comm_rank(comm)
    local_iel = 0
    for iel = 1:mesh.nelem
        if Int(mesh.el2gel[iel]) == target_gel
            local_iel = iel
            break
        end
    end
    found = local_iel != 0
    println("Rank $rank: gel=$target_gel " * (found ? "FOUND as local iel=$local_iel" : "not found locally"))
    return found ? (rank, local_iel) : nothing
end

function debug_dump_all_angular_deps_for_element(mesh, angular_deps, angular_deps_cr, target_gel::Int, comm)
    rank = MPI.Comm_rank(comm)
    local_iel = 0
    for iel = 1:mesh.nelem
        Int(mesh.el2gel[iel]) == target_gel && (local_iel = iel; break)
    end
    if local_iel == 0
        println("Rank $rank: gel=$target_gel not owned here")
        return nothing
    end

    n_dirs = mesh.extra_mesh[local_iel].extra_npoin
    println("Rank $rank: === ALL angular deps for gel=$target_gel (local iel=$local_iel), $n_dirs directions ===")
    for d = 1:n_dirs
        key = (local_iel, d)
        entries_found = String[]

        if haskey(angular_deps, key)
            for (dep_jel, dep_nodes) in angular_deps[key]
                dep_gel = Int(mesh.el2gel[dep_jel])
                push!(entries_found, "SAME-RANK dep_gel=$dep_gel nodes=$(sort(dep_nodes))")
            end
        end

        for fid = 1:6
            haskey(angular_deps_cr, (local_iel,fid,d)) || continue
            remote_gel, remote_nodes = angular_deps_cr[(local_iel,fid,d)]
            push!(entries_found, "CROSS-RANK(fid=$fid) dep_gel=$remote_gel nodes=$(sort(remote_nodes))")
        end

        if isempty(entries_found)
            println("  d=$d -> NO ENTRY (implicit same-index fallback)")
        else
            for e in entries_found
                println("  d=$d -> $e")
            end
        end
    end
end

function debug_gs_candidates_summary(mesh, solved, n_dirs, comm; max_show::Int=20)
    rank = MPI.Comm_rank(comm)
    n_shown = 0
    for iel = 1:mesh.nelem
        for d = 1:n_dirs[iel]
            solved[iel][d] && continue
            n_shown >= max_show && return
            gel = Int(mesh.el2gel[iel])
            println("Rank $rank: STUCK (iel=$iel,gel=$gel,d=$d)")
            n_shown += 1
        end
    end
end

function test_stress_combined_ncf_unified_multirank(inputs, nparts, distribute, sem, partitioned_model;
                                                     seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                     I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                     κ_val=0.001, σ_val=0.0005,
                                                     max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    nelem_base = sem.mesh.gnelem
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(nelem_base; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        idx = findfirst(==(gel), angular_target_gels)
        idx === nothing && continue
        push!(local_angular_candidates, (iel, angular_target_e_exts[idx]))
    end
    println("Rank $rank: $(length(local_angular_candidates)) local angular candidates before buffer/periodic checks")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                    ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                    neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                    neighbor_gel_group=neighbor_gel_group)   
 
    println("Rank $rank: $(length(final_angular_list)) angular refinements to APPLY after buffer/periodic checks")

    for (iel, e_ext) in final_angular_list
        adapt_single_angular_element_conforming!(mesh, iel, e_ext)
        verify_angular_phi_seam_conformity(mesh, iel)
    end

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)

    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
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

        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
            same_rank_kind, same_rank_partner, same_rank_is_periodic,
            ncf_as_child, ncf_L1L2, ncf_as_parent,
            coord_cache, gelpos_to_gip, spatial_perdir_exch,
            angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)

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
    rank == 0 && println("STRESS TEST unified sweep — max abs error: $max_err_global")
    return max_err_global
end

#=function test_stress_combined_ncf_unified_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                        seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                        κ_val::Float64=0.0001, I0::Float64=4.2,
                                                        use_instrumentation::Bool=false, err_threshold::Float64=1.0)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    nelem_base = sem.mesh.gnelem
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(nelem_base; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        idx = findfirst(==(gel), angular_target_gels)
        idx === nothing && continue
        push!(local_angular_candidates, (iel, angular_target_e_exts[idx]))
    end
    println("Rank $rank: $(length(local_angular_candidates)) local angular candidates before buffer/periodic checks")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                    ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                    neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                    neighbor_gel_group=neighbor_gel_group)
    
    println("Rank $rank: $(length(final_angular_list)) angular refinements to APPLY after buffer/periodic checks")

    for (iel, e_ext) in final_angular_list
        adapt_single_angular_element_conforming!(mesh, iel, e_ext)
        verify_angular_phi_seam_conformity(mesh, iel)
    end

    coord_cache, gelpos_to_gip = exchange_ncf_face_coordinates(mesh, neighbor_rank, neighbor_kind,
                                                                 neighbor_gel_group, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)

    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

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
    exact_fn = (x,y,z,Ωxd,Ωyd,Ωzd) -> I0*exp(-κ_val*(Ωxd*(x-r_bc[1])+Ωyd*(y-r_bc[2])+Ωzd*(z-r_bc[3])))

    if use_instrumentation
        x_dg, _, solved = run_dg_sweep_multirank_unified_instrumented!(mesh, metrics, ψ, dψ, ω,
            neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
            same_rank_kind, same_rank_partner, same_rank_is_periodic,
            ncf_as_child, ncf_L1L2, ncf_as_parent,
            coord_cache, gelpos_to_gip, spatial_perdir_exch,
            angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm;
            exact_fn=exact_fn, err_threshold=err_threshold, trace_iel=0)
    else
        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
            same_rank_kind, same_rank_partner, same_rank_is_periodic,
            ncf_as_child, ncf_L1L2, ncf_as_parent,
            coord_cache, gelpos_to_gip, spatial_perdir_exch,
            angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm)
    end

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
    rank == 0 && println("STRESS TEST unified sweep exponential — max abs error: $max_err_global")
    return max_err_global
end=#

function test_stress_combined_ncf_unified_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                                        seed::Int=42, frac::Float64=0.010, max_e_ext::Int=4,
                                                        κ_val::Float64=0.0001, I0::Float64=4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    nelem_base = sem.mesh.gnelem
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(nelem_base; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        idx = findfirst(==(gel), angular_target_gels)
        idx === nothing && continue
        push!(local_angular_candidates, (iel, angular_target_e_exts[idx]))
    end
    println("Rank $rank: $(length(local_angular_candidates)) local angular candidates before buffer/periodic checks")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                        ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                        neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                        neighbor_gel_group=neighbor_gel_group)
    println("Rank $rank: $(length(final_angular_list)) angular refinements to APPLY after buffer/periodic checks")

    for (iel, e_ext) in final_angular_list
        adapt_single_angular_element_conforming!(mesh, iel, e_ext)
        verify_angular_phi_seam_conformity(mesh, iel)
    end

    # ── COORDINATE-DRIVEN spatial-NCF cross-rank exchange (gid-based,
    # no transmitted position index) — REPLACES exchange_ncf_face_coordinates. ──
    remote_gid_coords = exchange_ncf_face_coordinates_v2(mesh, neighbor_rank, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)   # unused now, kept for later cleanup
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)

    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

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

    x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
        neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic,
        ncf_as_child, ncf_L1L2, ncf_as_parent,
        remote_gid_coords, spatial_perdir_exch,
        angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
        κ, σz, intΦ_eff, bc_fn_all, comm)

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
    rank == 0 && println("STRESS TEST unified sweep exponential — max abs error: $max_err_global")
    return max_err_global
end

function build_global_element_adjacency_and_refined_set(mesh, neighbor_elem, spatial_ncf_elements_local::Set{Int},
                                                          neighbor_rank, neighbor_gel_group,
                                                          gel_to_class::Dict{TInt,Set{TInt}}, comm)
    rank = MPI.Comm_rank(comm)

    local_edges = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        my_gel = TInt(mesh.el2gel[iel])
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            jgel = TInt(mesh.el2gel[Int(jel)])
            my_gel < jgel && push!(local_edges, (my_gel, jgel))
        end
    end
    for ((iel,fid), target_rank) in neighbor_rank
        my_gel = TInt(mesh.el2gel[iel])
        for rgel in neighbor_gel_group[(iel,fid)]
            g1,g2 = my_gel < rgel ? (my_gel,rgel) : (rgel,my_gel)
            push!(local_edges, (g1,g2))
        end
    end
    local_edges = unique(local_edges)

    n_local = length(local_edges)
    all_n = MPI.Allgather(Int32(n_local), comm)
    flat_local = Float64[]
    for (g1,g2) in local_edges
        push!(flat_local, Float64(g1)); push!(flat_local, Float64(g2))
    end
    counts = Cint.(all_n .* 2)
    all_flat = MPI.Allgatherv(flat_local, counts, comm)
    global_edges = Set{Tuple{TInt,TInt}}()
    for m = 1:(length(all_flat)÷2)
        g1 = TInt(round(Int, all_flat[2m-1])); g2 = TInt(round(Int, all_flat[2m]))
        push!(global_edges, g1 < g2 ? (g1,g2) : (g2,g1))
    end

    global_adj = Dict{TInt, Vector{TInt}}()
    for (g1,g2) in global_edges
        push!(get!(global_adj, g1, TInt[]), g2)
        push!(get!(global_adj, g2, TInt[]), g1)
    end

    local_refined_gels = Float64[]
    for iel in spatial_ncf_elements_local
        push!(local_refined_gels, Float64(mesh.el2gel[iel]))
    end
    n_local_r = length(local_refined_gels)
    all_nr = MPI.Allgather(Int32(n_local_r), comm)
    counts_r = Cint.(all_nr)
    all_rf = MPI.Allgatherv(local_refined_gels, counts_r, comm)
    global_refined_gels_raw = Set{TInt}(TInt(round(Int, v)) for v in all_rf)

    # ── NOW GLOBAL: expand via periodic equivalence classes, closing the
    # set under periodicity — a refined gel's ENTIRE periodic class
    # (wherever its members happen to live) counts as buffer-violating,
    # regardless of which rank owns which member. ───────────────────────────
    global_refined_gels = Set{TInt}()
    for g in global_refined_gels_raw
        push!(global_refined_gels, g)
        haskey(gel_to_class, g) || continue
        for partner_g in gel_to_class[g]
            push!(global_refined_gels, partner_g)
        end
    end

    return global_adj, global_refined_gels
end

function violates_buffer_global(mesh, iel::Int, global_adj, global_refined_gels, buffer_hops::Int)
    my_gel = TInt(mesh.el2gel[iel])
    my_gel in global_refined_gels && return true

    dist = Dict{TInt,Int}(my_gel => 0)
    frontier = [my_gel]
    depth = 0
    while !isempty(frontier) && depth < buffer_hops
        depth += 1
        next_frontier = TInt[]
        for g in frontier
            for ng in get(global_adj, g, TInt[])
                haskey(dist, ng) && continue
                dist[ng] = depth
                ng in global_refined_gels && return true
                push!(next_frontier, ng)
            end
        end
        frontier = next_frontier
    end
    return false
end

function test_stress_combined_ncf_unified_manufactured(inputs, nparts, distribute, sem, partitioned_model;
                                                        seed::Int=42, frac::Float64=0.01, max_e_ext::Int=4,
                                                        I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                        κ_val=0.001, σ_val=0.0005,
                                                        max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    nelem_base = sem.mesh.gnelem
    spatial_target_gels, angular_target_gels, angular_target_e_exts =
        generate_stress_test_targets(nelem_base; seed=seed, frac=frac, max_e_ext=max_e_ext)
    rank == 0 && println("Generated $(length(spatial_target_gels)) spatial targets, " *
                          "$(length(angular_target_gels)) angular targets")

    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh now has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    ngl = mesh.ngl

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
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
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf_v2(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    gel_to_class = build_periodic_equivalence_classes_v2(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        idx = findfirst(==(gel), angular_target_gels)
        idx === nothing && continue
        push!(local_angular_candidates, (iel, angular_target_e_exts[idx]))
    end
    println("Rank $rank: $(length(local_angular_candidates)) local angular candidates before buffer/periodic checks")

    final_angular_list = build_angular_refinement_list(mesh, local_angular_candidates, neighbor_elem,
                                                        ncf_as_child_keys, ncf_as_parent_keys, gel_to_class, comm;
                                                        neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind,
                                                        neighbor_gel_group=neighbor_gel_group)
    println("Rank $rank: $(length(final_angular_list)) angular refinements to APPLY after buffer/periodic checks")

    for (iel, e_ext) in final_angular_list
        adapt_single_angular_element_conforming!(mesh, iel, e_ext)
        verify_angular_phi_seam_conformity(mesh, iel)
    end

    remote_gid_coords = exchange_ncf_face_coordinates_v2(mesh, neighbor_rank, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    spatial_perdir_exch = build_ncf_spatial_perdir_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)

    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    groups, representative_of = group_elements_by_angular_mesh(mesh)
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

        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
            same_rank_kind, same_rank_partner, same_rank_is_periodic,
            ncf_as_child, ncf_L1L2, ncf_as_parent,
            remote_gid_coords, spatial_perdir_exch,
            angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
            κ, σz, intΦ_eff, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all)

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
    rank == 0 && println("STRESS TEST unified sweep manufactured — max abs error: $max_err_global")
    return max_err_global
end

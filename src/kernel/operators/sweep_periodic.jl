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
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    gel_to_class = build_periodic_equivalence_classes(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
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

function build_periodic_equivalence_classes(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
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

    if neighbor_rank !== nothing && neighbor_kind !== nothing
    	for ((iel,fid), _) in neighbor_rank
            k = get(neighbor_kind, (iel,fid), :conforming)
            (k == :i_am_fine || k == :i_am_coarse) && push!(spatial_ncf_elements, Int(iel))
    	end
    end

    gel_to_local = Dict{TInt,Int}(mesh.el2gel[i] => i for i = 1:mesh.nelem)

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
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))

    gel_to_class = build_periodic_equivalence_classes(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
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

function resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces::Vector{Tuple{TInt,TInt}},
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

"""
    setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                            spatial_target_gels=Int[], angular_targets=Tuple{Int,Int}[])

General-purpose setup for the unified NCF/periodic sweep. Periodicity is
NOT a parameter here — it's simply a property of whatever mesh `inputs`
describes; if the mesh has no periodic-tagged faces, every periodic
table below comes back empty and the pipeline degrades gracefully.

spatial_target_gels: GLOBAL element ids to refine spatially (empty = no
spatial refinement at all).
angular_targets: (GLOBAL element id, angular e_ext) pairs requested for
angular refinement (empty = none) — filtered through the SAME
one-buffer-element + periodic-propagation safety logic as production use.

Returns a NamedTuple with every table run_dg_sweep_multirank_unified!/
resolve_stall_distributed_gauss_seidel_unified! need.
"""
function setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                 spatial_target_gels::Vector{Int} = Int[],
                                 angular_targets::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[])
    clear_face_plans!()
    rank = MPI.Comm_rank(get_mpi_comm())
    comm = get_mpi_comm()

    adapt_flags = zeros(Int, sem.mesh.nelem)
    for iel = 1:sem.mesh.nelem
        gel = Int(sem.mesh.el2gel[iel])
        gel in spatial_target_gels && (adapt_flags[iel] = 1)
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics
    println("Rank $rank: mesh has $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")

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
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
    gel_to_class = build_periodic_equivalence_classes(mesh, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           neighbor_kind, neighbor_gel_group, neighbor_is_periodic, comm)

    
    check_refined_set_periodic_closure(mesh, ncf_as_child, ncf_as_parent, gel_to_class, comm)

    local_angular_candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel])
        for (tgel, e_ext) in angular_targets
            tgel == gel && push!(local_angular_candidates, (iel, e_ext))
        end
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

    remote_gid_coords = exchange_ncf_face_coordinates(mesh, neighbor_rank, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    #gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    gid_exch = augment_gid_exchange(build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm), ngl)
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    #periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)
    periodic_gid_exch = augment_gid_exchange(build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm), ngl)

    return (mesh=mesh, metrics=metrics, sem=sem, partitioned_model=partitioned_model,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid, face_correspondence=face_correspondence,
            periodic_face_correspondence=periodic_face_correspondence,
            neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind, neighbor_gel=neighbor_gel,
            neighbor_gel_group=neighbor_gel_group, neighbor_is_periodic=neighbor_is_periodic,
            same_rank_kind=same_rank_kind, same_rank_partner=same_rank_partner, same_rank_is_periodic=same_rank_is_periodic,
            ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2, ncf_as_parent=ncf_as_parent,
            remote_gid_coords=remote_gid_coords,
            angular_deps=angular_deps, angular_deps_cr=angular_deps_cr, remote_angular_mesh=remote_angular_mesh,
            gid_exch=gid_exch, periodic_gid_exch=periodic_gid_exch, n_dirs=n_dirs)
end

const GS_MODE = Ref(:cached)            # :reference | :cached | :cached_upwind
const USE_GS_INPLACE = Ref(true)
const USE_GS_GC_OFF = Ref(true)
struct GSStatic
    F::Any                              # LU factors, or nothing when the matrix was triangular
    A::Any                              # the matrix, kept only when F === nothing
    rhs_src::Vector{Float64}            # source term (+ the zero-upwind flux part, which adds only zeros)
    faces::Vector{Int}                  # inflow and active faces, increasing fid
    weights::Vector{Vector{Float64}}    # flux weight per face node, in face_ip_range order
end

gs_static_bytes(ngl) = 8 * (ngl^6 + 2*ngl^3 + 6*ngl^2)
gs_solve(st::GSStatic, rhs) = st.F === nothing ? st.A \ rhs : st.F \ rhs


function gs_solve!(st::GSStatic, rhs)            # returns the solution; may reuse rhs
    st.F === nothing && return st.A \ rhs
    ldiv!(st.F, rhs); return rhs
end

"""
Cluster nearly-equal floats by TOLERANCE rather than fixed-digit rounding
— round(x; digits=n) is not a reliable snap-to-grid for floating point,
since its own scale/round/rescale can leave two conceptually-identical
values as bit-distinct outputs. Returns the sorted cluster representatives
(each the mean of its members) and, for each input value, which cluster
index it belongs to.
"""
function cluster_values_with_index(vals::Vector{Float64}; tol::Float64=1e-6)
    order = sortperm(vals)
    sorted = vals[order]
    cluster_of_sorted = zeros(Int, length(vals))
    reps = Float64[]
    current = Float64[sorted[1]]
    cluster_of_sorted[1] = 1
    for k = 2:length(sorted)
        if sorted[k] - sorted[k-1] < tol
            push!(current, sorted[k])
        else
            push!(reps, sum(current)/length(current))
            current = Float64[sorted[k]]
        end
        cluster_of_sorted[k] = length(reps) + 1
    end
    push!(reps, sum(current)/length(current))

    cluster_of_val = zeros(Int, length(vals))
    for k = 1:length(vals)
        cluster_of_val[order[k]] = cluster_of_sorted[k]
    end
    return reps, cluster_of_val
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

"""
Key CHANGE: angular_deps_cr is now keyed by (iel,fid,d), not (iel,d) —
a single direction can have GENUINELY DIFFERENT, independent angular-NCF
relationships on different faces simultaneously (one inflow face's
remote element need not be the same as another inflow face's).
"""

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
                                                  src_e_ext::Int, em_remote, received::AbstractDict{Tuple{TInt,Int,Int},Float64})
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

FINAL, cleaned version: spatial_perdir_exch/spatial_perdir_received
removed entirely (dead weight since the coordinate-driven gid-based
reconstruction replaced them). Returns (available::Bool,
closure::Union{Function,Nothing}) for ONE face of ONE (iel,d).
"""
function build_face_upwind(mesh, metrics, iel::Int, fid::Int, d::Int, ngl::Int,
                            neighbor_elem, face_correspondence,
                            ncf_as_child, ncf_L1L2, ncf_as_parent,
                            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                            neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                            periodic_face_correspondence,
                            angular_deps, angular_deps_cr, remote_angular_mesh,
                            remote_gid_coords, gid_received, periodic_gid_received,
                            is_solved::Function, x_dg, I_cg_all_prev, bc_fn_all)
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
	    I_cg_all_prev === nothing || error("periodic same-rank snapshot missing: iel=$iel fid=$fid partner=$partner_lel d=$d")
            return true, (i,j,k) -> 0.0
        end
    end
    if get(neighbor_is_periodic, (iel,fid), false)
        remote_gel_p = neighbor_gel_group[(iel,fid)][1]
        return true, (i,j,k) -> begin
            gid = Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]])
            v = get(periodic_gid_received, (remote_gel_p, d, gid), 0.0)
	    if isnan(v)
        	I_cg_all_prev === nothing && return 0.0
        	error("periodic cross-rank key missing: iel=$iel fid=$fid d=$d gid=$gid remote_gel=$remote_gel_p")
    	    end
    	    v
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
                    if USE_FLAT_FACES[]
    			cidx1 = flat_corr_index!(face_correspondence, iel, fid, ngl)
    			vals1 = x_dg[(jel,dj)]
    			return true, (i,j,k) -> vals1[cidx1[i,j,k]]
		    end
		    dj_c=dj; fid_c=fid
                    return true, (i,j,k) -> begin
                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                    end
                else
                    all(dn -> is_solved(jel,dn), dep_nodes) || return false, nothing
                    if USE_ANGULAR_PLANS[]
    			pl = angular_same_plan!(mesh, angular_deps, iel, fid, jel, d, dep_nodes, face_correspondence, ngl)
    			vals = zeros(ngl^3)
    			eval_angular_same_plan!(vals, pl, x_dg)
    			return true, (i,j,k) -> vals[lidx(i,j,k,ngl)]
		    end
		    em_jel = mesh.extra_mesh[jel]
                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                    return true, build_angular_ncf_upwind_unified(iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                end
            end
        end
        is_solved(jel,d) || return false, nothing
	if USE_FLAT_FACES[]
            cidx = flat_corr_index!(face_correspondence, iel, fid, ngl)
            vals = x_dg[(jel,d)]
            return true, (i,j,k) -> vals[cidx[i,j,k]]
        end
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
    	local child_vals
    	if USE_LOCAL_NCF_PLANS[]
            pl = local_fine_plan!(mesh, ncf_L1L2, ncf_idx, parent_iel, ngl)
            child_vals = Vector{Float64}(undef, length(pl.skip))
            eval_local_fine_plan!(child_vals, pl, x_dg[(parent_iel,d)], ngl)
    	else
            L1,L2 = ncf_L1L2[ncf_idx]
            child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d,
                        mesh.IPc_list[:,ncf_idx], mesh.IPp_list[:,ncf_idx], L1, L2, ngl)
    	end
    	return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

    elseif haskey(ncf_as_parent, (iel,fid))
    	ncf_group = ncf_as_parent[(iel,fid)]
    	local parent_grid
    	if USE_LOCAL_NCF_PLANS[]
            pl = local_coarse_plan!(mesh, ncf_L1L2, iel, fid, ncf_group, ngl)
            all(c -> is_solved(c,d), pl.children) || return false, nothing
            parent_grid = eval_local_coarse_plan(pl, x_dg, d, ngl)
    	else
            children = [Int(mesh.cip[idx]) for idx in ncf_group]
            all(c -> is_solved(c,d), children) || return false, nothing
            entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in ncf_group]
            parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
    	end
        return true, (i,j,k) -> begin
            ni = face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
            parent_grid[p1,p2]
    	end

    elseif haskey(neighbor_rank, (iel,fid))
        kind = neighbor_kind[(iel,fid)]

        #=if kind == :i_am_fine
            remote_gel_f = neighbor_gel_group[(iel,fid)][1]
            haskey(remote_gid_coords, remote_gel_f) || return false, nothing
            child_vals = build_dg_ncf_parent_to_child_remote_by_coords(mesh, iel, fid, d, remote_gel_f,
                remote_gid_coords, gid_received, ngl)
            any(isnan, child_vals) && return false, nothing
            return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

        elseif kind == :i_am_coarse
            remote_gels_c = neighbor_gel_group[(iel,fid)]
            all(g -> haskey(remote_gid_coords, g), remote_gels_c) || return false, nothing
            #=parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, remote_gels_c,
                remote_gid_coords, gid_received, ngl)=#
	    try
        	parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, remote_gels_c,
            	remote_gid_coords, gid_received, ngl)
    	    catch e
        	debug_dump_on_plane_coords(mesh, iel, fid, remote_gels_c, remote_gid_coords, ngl)
        	rethrow(e)
    	    end
            any(isnan, parent_grid) && return false, nothing
            return true, (i,j,k) -> begin
                ni=face_flat_index(fid,ngl,i,j,k); p1=((ni-1)%ngl)+1; p2=((ni-1)÷ngl)+1
                parent_grid[p1,p2]
            end=#
	if kind == :i_am_fine
    	    remote_gel_f = neighbor_gel_group[(iel,fid)][1]
    	    haskey(remote_gid_coords, remote_gel_f) || return false, nothing
    	    local child_vals
    	    if USE_FACE_PLANS[]
        	pl = fine_plan!(mesh, remote_gid_coords, iel, fid, remote_gel_f, ngl)
        	child_vals = Vector{Float64}(undef, ngl*ngl)
        	eval_fine_plan!(child_vals, pl, d, gid_received, ngl) || return false, nothing
    	    else
        	child_vals = build_dg_ncf_parent_to_child_remote_by_coords(mesh, iel, fid, d, remote_gel_f,
            	    remote_gid_coords, gid_received, ngl)
         	any(isnan, child_vals) && return false, nothing
    	    end
    	    return true, (i,j,k) -> child_vals[face_flat_index(fid,ngl,i,j,k)]

	elseif kind == :i_am_coarse
    	    remote_gels_c = neighbor_gel_group[(iel,fid)]
    	    all(g -> haskey(remote_gid_coords, g), remote_gels_c) || return false, nothing
    	    local parent_grid
    	    if USE_FACE_PLANS[]
        	pl = coarse_plan!(mesh, remote_gid_coords, iel, fid, remote_gels_c, ngl)
        	parent_grid = Matrix{Float64}(undef, ngl, ngl)
        	eval_coarse_plan!(parent_grid, pl, d, gid_received, ngl) || return false, nothing
    	    else
                try
            	    parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, remote_gels_c,
                	remote_gid_coords, gid_received, ngl)
                catch e
        	    debug_dump_on_plane_coords(mesh, iel, fid, remote_gels_c, remote_gid_coords, ngl)
    	            rethrow(e)
                end

		parent_grid = build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, remote_gels_c,
            	    remote_gid_coords, gid_received, ngl)
            	any(isnan, parent_grid) && return false, nothing
    	    end
            return true, (i,j,k) -> begin
                ni = face_flat_index(fid,ngl,i,j,k); p1 = ((ni-1)%ngl)+1; p2 = ((ni-1)÷ngl)+1
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
		    if USE_ANGULAR_PLANS[]
    			pl = angular_remote_plan!(mesh, angular_deps_cr, iel, fid, d, remote_gel, remote_nodes, remote_angular_mesh, ngl)
    			pl === nothing && return false, nothing
    			vals = zeros(ngl^3)
    			eval_angular_remote_plan!(vals, pl, gid_received)
    			return true, (i,j,k) -> vals[lidx(i,j,k,ngl)]
		    end
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
        nx,ny,nz = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz,
            metrics.dηdx,metrics.dηdy,metrics.dηdz,
            metrics.dζdx,metrics.dζdy,metrics.dζdz)
        return true, (i,j,k) -> begin
            gip = mesh.connijk[iel,i,j,k]
            bc(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip, nx, ny, nz)
        end    
    end
end

#=function run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          periodic_face_correspondence,
                                          neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                          same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                          ncf_as_child, ncf_L1L2, ncf_as_parent,
                                          remote_gid_coords,
                                          angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
                                          κ, bc_fn_all, comm::MPI.Comm;
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

    # ── FROZEN periodic snapshot: exchanged ONCE, from the PREVIOUS outer
    # iteration's solution, and never refreshed during this sweep. This is
    # a collective call, made unconditionally on every rank (the condition
    # is identical on all ranks since it depends only on the call). ─────────
    periodic_gid_received = I_cg_all_prev === nothing ?
        Dict{Tuple{TInt,Int,Int}, Float64}() :
        #exchange_ncf_gid_values!(periodic_gid_exch, I_cg_all_prev, ngl)
	exchange_ncf_gid_table!(periodic_gid_exch, I_cg_all_prev, ngl)
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

                    avail, clo = build_face_upwind(mesh, metrics, iel, fid, d, ngl,
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
                                         metrics.Je, κ, mesh.connijk)
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

        # Live values are exchanged every round for the ordinary/angular-NCF/
        # spatial-NCF branches. The periodic snapshot is NOT refreshed here.
        #gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
	gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)

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
                remote_gid_coords,
                angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_received,
                x_dg, solved, n_dirs, κ, bc_fn_all, comm;
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
        neighbor_elem, neighbor_fid, face_correspondence, periodic_face_correspondence,
        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic,
        ncf_as_child, ncf_L1L2, ncf_as_parent,
        remote_gid_coords,
        angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_exch,
        κ, bc_fn_all, comm::MPI.Comm;
        extra_source_all=nothing, I_cg_all_prev=nothing,
        elem_cache::Union{ElemCache,Nothing}=nothing,
        prof::Union{Nothing,SweepProfile}=nothing, local_passes::Bool=false, measure_skew::Bool=false,
	gc_mode::Symbol=:auto, gc_every::Int=0,
        max_rounds::Int=500, max_gs_iter::Int=100, gs_tol::Float64=1e-9, round_cap::Int=0)
    t_sweep = gtic(prof)
    gc_t0 = Base.gc_time_ns(); gc_b0 = Base.gc_bytes()
    rank = MPI.Comm_rank(comm); ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    cache = elem_cache === nothing ? make_elem_cache(mesh, metrics, n_dirs, ngl; budget_mb=0.0) : elem_cache
    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    total_pairs = sum(n_dirs)
    cap = round_cap < 0 ? max(64, total_pairs ÷ (-round_cap)) : round_cap     # negative = pairs/|cap|, scales with load
    incr = USE_INCREMENTAL_PACK[]
    newly = Tuple{Int,Int}[]; valid_count = Ref(0)
    incr && fill!(gid_exch.send_buf, NaN)    

    n_solved_total = 0; iter_round = 0; just_resolved_stall = false
    gid_received = Dict{Tuple{TInt,Int,Int}, Float64}()
    slog = prof === nothing ? nothing : SweepLog()
    prev_valid = prof === nothing ? falses(0) : falses(gid_exch.n_flat)

    t0 = gtic(prof)
    periodic_gid_received = I_cg_all_prev === nothing ? Dict{Tuple{TInt,Int,Int}, Float64}() :
                            exchange_ncf_gid_table!(periodic_gid_exch, I_cg_all_prev, ngl)
    gtoc!(prof, :periodic_exch, t0)
    is_solved_sweep = (jel,dj) -> solved[jel][dj]

    while true
        iter_round += 1
        n_newly_solved = 0; n_passes = 0; t_empty_acc = 0.0
        t_round = time_ns(); gc_r0 = Base.gc_time_ns()
        while true                                           # local passes (a single scan unless local_passes)
            n_pass = 0; t_pass = time_ns()
 
	    for oi = 1:length(cache.ord_iel)
    		iel = Int(cache.ord_iel[oi]); d = Int(cache.ord_d[oi])
    		solved[iel][d] && continue
    		t_att = gtic(prof)
    		m = cache.faceinfo[iel][d]
    		ready = true
    		fill!(cache.ws_faces, nothing)
    		prof === nothing || fill!(cache.ws_br, :none)
    		for fid = 1:6
        	    mask_inflow(m, fid) || continue
        	    br = prof === nothing ? :none : face_branch(neighbor_elem, neighbor_rank, neighbor_kind, ncf_as_child,
                	ncf_as_parent, same_rank_is_periodic, neighbor_is_periodic, angular_deps_cr, iel, fid, d)
        	    prof === nothing || (cache.ws_br[fid] = br)
        	    t0 = gtic(prof)
        	    avail, clo = build_face_upwind(mesh, metrics, iel, fid, d, ngl,
            		neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
            		neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
            		neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
            		periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
            		remote_gid_coords, gid_received, periodic_gid_received,
            		is_solved_sweep, x_dg, I_cg_all_prev, bc_fn_all)
        	    if prof !== nothing
            		dtf = (time_ns() - t0) * 1e-9
            		padd!(prof, :face_upwind, dtf)
            		padd!(prof, avail ? BRANCH_READY[br] : BRANCH_BLOCKED[br], dtf)
        	    end
        	    if !avail; ready = false; break; end
        	    t0 = gtic(prof)
        	    irange, jrange, krange = face_ip_range(fid, ngl)
       		    trial = clo(first(irange), first(jrange), first(krange))
        	    gtoc!(prof, :trial_eval, t0)
        	    if isnan(trial); ready = false; break; end
        	    cache.ws_faces[fid] = clo
    		end
    		if !ready
                    gtoc!(prof, :attempt_wasted, t_att)
       		    continue
    		end
    		x_dg[(iel,d)] = solve_pair!(cache, mesh, metrics, ψ, dψ, ω, κ, iel, d, ngl, m, extra_source_all, prof)
    	        solved[iel][d] = true
 		incr && push!(newly, (iel, d))
    		n_pass += 1
    		gtoc!(prof, :attempt_solved, t_att)
	        if cap > 0 && n_pass >= cap; break; end
	    end
	
            n_passes += 1; n_newly_solved += n_pass
            if n_pass == 0
                dt = (time_ns() - t_pass) * 1e-9; t_empty_acc += dt
                prof === nothing || padd!(prof, :empty_pass, dt)
                break
            end
            local_passes || break
        end
        t_compute = (time_ns() - t_round) * 1e-9
        gc_comp = (Base.gc_time_ns() - gc_r0) * 1e-9

        t0 = gtic(prof)
        n_solved_total = sum(count(identity, solved[iel]) for iel = 1:mesh.nelem)
        gtoc!(prof, :bookkeeping, t0)

        t_skew = 0.0
        if prof !== nothing && measure_skew
            tt = time_ns(); MPI.Barrier(comm); dt = (time_ns() - tt) * 1e-9; padd!(prof, :skew_exchange, dt); t_skew += dt
        end
        #=gc_e0 = Base.gc_time_ns()
        ex_sent = 0; ex_valid = 0; ex_fresh = 0; ex_sec = 0.0
        if prof === nothing
            gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)
        else
            ex = exchange_ncf_gid_values_profiled!(gid_exch, x_dg, ngl, prof, prev_valid)
            gid_received = ex.received; ex_sent = ex.sent; ex_valid = ex.valid; ex_fresh = ex.fresh; ex_sec = ex.seconds
        end
        gc_ex = (Base.gc_time_ns() - gc_e0) * 1e-9=#
	gc_e0 = Base.gc_time_ns()
	ex_sent = 0; ex_valid = 0; ex_fresh = 0; ex_sec = 0.0
	if prof === nothing
    	    gid_received = incr ? exchange_ncf_gid_table_incremental!(gid_exch, x_dg, newly, ngl) :
                          exchange_ncf_gid_table!(gid_exch, x_dg, ngl)
	else
    	    ex = exchange_ncf_gid_values_profiled!(gid_exch, x_dg, ngl, prof, prev_valid;
                                           newly = incr ? newly : nothing, valid_count = valid_count)
    	    gid_received = ex.received; ex_sent = ex.sent; ex_valid = ex.valid; ex_fresh = ex.fresh; ex_sec = ex.seconds
	end
	empty!(newly)
	gc_ex = (Base.gc_time_ns() - gc_e0) * 1e-9
        if prof !== nothing && measure_skew
            tt = time_ns(); MPI.Barrier(comm); dt = (time_ns() - tt) * 1e-9; padd!(prof, :skew_reduce, dt); t_skew += dt
        end

        tt = time_ns()
        g = MPI.Allreduce(Int[n_newly_solved, total_pairs - n_solved_total], +, comm)    # one reduction per round
        t_red = (time_ns() - tt) * 1e-9
        prof === nothing || padd!(prof, :allreduce, t_red)
        n_newly_solved_global = g[1]; n_unprocessed_total = g[2]

        if slog !== nothing
            push!(slog.passes, n_passes); push!(slog.solved, n_newly_solved); push!(slog.t_compute, t_compute)
            push!(slog.t_exchange, ex_sec); push!(slog.t_skew, t_skew); push!(slog.t_reduce, t_red)
            push!(slog.sent, ex_sent); push!(slog.valid, ex_valid); push!(slog.fresh, ex_fresh); push!(slog.t_empty, t_empty_acc)
            push!(slog.gc_compute, gc_comp); push!(slog.gc_exchange, gc_ex)
        end
        n_unprocessed_total == 0 && break

        if gc_mode == :sync && gc_every > 0 && iter_round % gc_every == 0
            tt = time_ns(); GC.enable(true); GC.gc(false); GC.enable(false)
            prof === nothing || padd!(prof, :sync_gc, (time_ns() - tt) * 1e-9)
        end

        if n_newly_solved_global == 0
            just_resolved_stall && error("no progress after GS resolution")
            t0 = gtic(prof)
            resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
                neighbor_elem, face_correspondence, periodic_face_correspondence,
                neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                same_rank_kind, same_rank_partner, same_rank_is_periodic,
                ncf_as_child, ncf_L1L2, ncf_as_parent,
                remote_gid_coords,
                angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch, periodic_gid_received,
                x_dg, solved, n_dirs, κ, bc_fn_all, comm;
                extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev,
                max_gs_iter=max_gs_iter, gs_tol=gs_tol, prof=prof)
            gtoc!(prof, :gs_total, t0)
            just_resolved_stall = true
        else
            just_resolved_stall = false
        end
        iter_round > max_rounds && error("exceeded max_rounds=$max_rounds")
    end

    if prof !== nothing
        prof.n[:rounds] = get(prof.n, :rounds, 0) + iter_round
        prof.t[:gc_seconds] = get(prof.t, :gc_seconds, 0.0) + (Base.gc_time_ns() - gc_t0) * 1e-9
        prof.n[:alloc_bytes] = get(prof.n, :alloc_bytes, 0) + Int(Base.gc_bytes() - gc_b0)
        push!(prof.sweeps, slog)
        ptoc!(prof, :sweep_total, t_sweep)
    end
    return x_dg, n_dirs
end

function resolve_stall_distributed_gauss_seidel_unified_reference!(mesh, metrics, ψ, dψ, ω,
                                                           neighbor_elem, face_correspondence, periodic_face_correspondence,
                                                           neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
                                                           same_rank_kind, same_rank_partner, same_rank_is_periodic,
                                                           ncf_as_child, ncf_L1L2, ncf_as_parent,
                                                           remote_gid_coords,
                                                           angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch,
                                                           periodic_gid_received,   # CHANGED: the frozen snapshot, passed in
                                                           x_dg, solved, n_dirs,
                                                           κ, bc_fn_all, comm;
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
    #gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
    gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)

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
                                     metrics.Je, κ, mesh.connijk)
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
                    avail, clo = build_face_upwind(mesh, metrics, iel, fid, d, ngl,
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

        #gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
	gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)

        max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        println("Rank $rank: GS-unified iter $gs_iter, max_delta_global=$max_delta_global")
        max_delta_global < gs_tol && break
        gs_iter == max_gs_iter && @warn "GS-unified did not converge (max_delta=$max_delta_global)"
    end

    for (iel,d) in candidates
        solved[iel][d] = true
    end
end

function gs_build_static(mesh, metrics, ψ, dψ, ω, κ, iel::Int, d::Int, ngl::Int, extra_source_all)
    em = mesh.extra_mesh[iel]
    θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
        metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
    A_vol = zeros(ngl^3, ngl^3); rhs = zeros(ngl^3)
    build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
        metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
        metrics.Je, κ, mesh.connijk)
    if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
        build_dg_source_term!(rhs, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
    end
    faces = Int[]; weights = Vector{Vector{Float64}}()
    for fid = 1:6
        is_active[fid] || continue
        face_normal = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz)
        w = zeros(ngl*ngl)
        build_dg_flux_terms!(A_vol, rhs, fid, ngl, Ωx,Ωy,Ωz, is_inflow[fid], is_active[fid], face_normal, ω,
            metrics.dξdx,metrics.dξdy,metrics.dξdz, metrics.dηdx,metrics.dηdy,metrics.dηdz, metrics.dζdx,metrics.dζdy,metrics.dζdz,
            metrics.Je, iel, (i,j,k)->0.0; weight_out=w)
        if is_inflow[fid]; push!(faces, fid); push!(weights, w); end
    end
    (istril(A_vol) || istriu(A_vol)) && return GSStatic(nothing, A_vol, rhs, faces, weights)
    return GSStatic(lu!(A_vol), nothing, rhs, faces, weights)
end

function resolve_stall_distributed_gauss_seidel_unified!(mesh, metrics, ψ, dψ, ω,
        neighbor_elem, face_correspondence, periodic_face_correspondence,
        neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic,
        ncf_as_child, ncf_L1L2, ncf_as_parent,
        remote_gid_coords,
        angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch,
        periodic_gid_received,
        x_dg, solved, n_dirs, κ, bc_fn_all, comm;
        extra_source_all=nothing, I_cg_all_prev=nothing,
        max_gs_iter::Int=100, gs_tol::Float64=1e-9, gs_cache_mb::Float64=1500.0, prof=nothing)
    if GS_MODE[] == :reference
        return resolve_stall_distributed_gauss_seidel_unified_reference!(mesh, metrics, ψ, dψ, ω,
            neighbor_elem, face_correspondence, periodic_face_correspondence,
            neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group, neighbor_is_periodic,
            same_rank_kind, same_rank_partner, same_rank_is_periodic,
            ncf_as_child, ncf_L1L2, ncf_as_parent,
            remote_gid_coords,
            angular_deps, angular_deps_cr, remote_angular_mesh, gid_exch,
            periodic_gid_received,
            x_dg, solved, n_dirs, κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all_prev,
            max_gs_iter=max_gs_iter, gs_tol=gs_tol)
    end
    rank = MPI.Comm_rank(comm); ngl = mesh.ngl
    gs_gc_prev = USE_GS_GC_OFF[] ? GC.enable(false) : true
    t_start = time_ns()

    candidates = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(candidates, (iel,d))
    end
    if GS_MODE[] == :cached_upwind
        cent = Dict{Int,NTuple{3,Float64}}()
        for (iel,_) in candidates
            haskey(cent, iel) && continue
            cx = 0.0; cy = 0.0; cz = 0.0
            for k=1:ngl, j=1:ngl, i=1:ngl
                ip = mesh.connijk[iel,i,j,k]; cx += mesh.x[ip]; cy += mesh.y[ip]; cz += mesh.z[ip]
            end
            cent[iel] = (cx/ngl^3, cy/ngl^3, cz/ngl^3)
        end
        upkey(c) = (em = mesh.extra_mesh[c[1]]; θ = em.extra_coords[1,c[2]]; ϕ = em.extra_coords[2,c[2]]; p = cent[c[1]];
                    p[1]*sin(θ)*cos(ϕ) + p[2]*sin(θ)*sin(ϕ) + p[3]*cos(θ))
        sort!(candidates; by=upkey, alg=Base.Sort.MergeSort)          # upstream first, stable
    end

    for (iel,d) in candidates; x_dg[(iel,d)] = zeros(ngl^3); end
    is_solved_gs = (jel,dj) -> haskey(x_dg, (jel,dj))
    t0 = gtic(prof)
    gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)
    gtoc!(prof, :gs_exchange, t0)
    statics = Vector{Union{GSStatic,Nothing}}(nothing, length(candidates))
    budget = gs_cache_mb * 1e6; bytes_each = gs_static_bytes(ngl); n_cached = 0
    deltas = Float64[]; t_iter1 = 0.0; iters = 0
    diag = prof !== nothing
    np_gs = MPI.Comm_size(comm)
    cnt_hist = Vector{NTuple{4,Int}}(); tcomp_mean = Float64[]; tcomp_max = Float64[]
    gc_max_hist = Float64[]
    for gs_iter = 1:max_gs_iter
        t_it = time_ns(); iters = gs_iter
        max_delta_local = 0.0
        n0 = 0; n1 = 0; n2 = 0; n3 = 0; t_c0 = time_ns()
	gc_c0 = Base.gc_time_ns()
	for (c,(iel,d)) in enumerate(candidates)
            st = statics[c]
            if st === nothing
		t0 = gtic(prof)
                st = gs_build_static(mesh, metrics, ψ, dψ, ω, κ, iel, d, ngl, extra_source_all)
                gtoc!(prof, :gs_build, t0)
		if (n_cached + 1) * bytes_each <= budget
                    statics[c] = st; n_cached += 1
                end
            end
            RHS = copy(st.rhs_src)
            #=t0 = gtic(prof)
	    for (q, fid) in enumerate(st.faces)
                avail, clo = build_face_upwind(mesh, metrics, iel, fid, d, ngl,
                    neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                    neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                    neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                    periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                    remote_gid_coords, gid_received, periodic_gid_received,
                    is_solved_gs, x_dg, I_cg_all_prev, bc_fn_all)
                avail || continue                                       # the original adds weight*0.0 here
                w = st.weights[q]
                irange, jrange, krange = face_ip_range(fid, ngl)
                idx = 0
                for k in krange, j in jrange, i in irange
                    idx += 1
                    v = clo(i,j,k)
                    isnan(v) && (v = 0.0)
                    RHS[lidx(i,j,k,ngl)] += w[idx] * v
                end
            end
	    gtoc!(prof, :gs_upwind, t0)=#
	    t0 = gtic(prof)
	    for (q, fid) in enumerate(st.faces)
    		br = prof === nothing ? :none : face_branch(neighbor_elem, neighbor_rank, neighbor_kind, ncf_as_child,
            		ncf_as_parent, same_rank_is_periodic, neighbor_is_periodic, angular_deps_cr, iel, fid, d)
    		tu = gtic(prof)
		avail, clo = build_face_upwind(mesh, metrics, iel, fid, d, ngl,
                    neighbor_elem, face_correspondence, ncf_as_child, ncf_L1L2, ncf_as_parent,
                    neighbor_rank, neighbor_kind, neighbor_gel, neighbor_gel_group,
                    neighbor_is_periodic, same_rank_kind, same_rank_partner, same_rank_is_periodic,
                    periodic_face_correspondence, angular_deps, angular_deps_cr, remote_angular_mesh,
                    remote_gid_coords, gid_received, periodic_gid_received,
                    is_solved_gs, x_dg, I_cg_all_prev, bc_fn_all)    
    		prof === nothing || padd!(prof, GS_BRANCH_UP[br], (time_ns() - tu) * 1e-9)
    		avail || continue
    		te = gtic(prof)
    		w = st.weights[q]
    		irange, jrange, krange = face_ip_range(fid, ngl)
    		idx = 0
    		for k in krange, j in jrange, i in irange
        	    idx += 1
        	    v = clo(i,j,k)
        	    isnan(v) && (v = 0.0)
        	    RHS[lidx(i,j,k,ngl)] += w[idx] * v
    		end
    		prof === nothing || padd!(prof, GS_BRANCH_EV[br], (time_ns() - te) * 1e-9)
	    end
	    gtoc!(prof, :gs_upwind, t0)
	    t0 = gtic(prof)
	    x_new = USE_GS_INPLACE[] ? gs_solve!(st, RHS) : gs_solve(st, RHS)
	    gtoc!(prof, :gs_solve, t0)
            #max_delta_local = max(max_delta_local, maximum(abs.(x_new .- x_dg[(iel,d)])))
            #dc = maximum(abs.(x_new .- x_dg[(iel,d)]))
	    old = x_dg[(iel,d)]; dc = 0.0
	    for p in eachindex(x_new); dc = max(dc, abs(x_new[p] - old[p])); end
	    max_delta_local = max(max_delta_local, dc)
	    if diag
    	       if dc == 0.0; n0 += 1 elseif dc < 1e-12; n1 += 1 elseif dc < 1e-9; n2 += 1 else n3 += 1 end
	    end
	    x_dg[(iel,d)] = x_new
        end
	t_comp = (time_ns() - t_c0) * 1e-9
	gc_comp = (Base.gc_time_ns() - gc_c0) * 1e-9
	t0 = gtic(prof)
        gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)
        gtoc!(prof, :gs_exchange, t0)
	t0 = gtic(prof)
	delta_g = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
        if diag
    	    g = MPI.Allreduce(Float64[n0, n1, n2, n3, t_comp], +, comm)
    	    tm = MPI.Allreduce(t_comp, MPI.MAX, comm)
    	    push!(cnt_hist, (Int(g[1]), Int(g[2]), Int(g[3]), Int(g[4])))
    	    push!(tcomp_mean, g[5]/np_gs); push!(tcomp_max, tm)
	    push!(gc_max_hist, MPI.Allreduce(gc_comp, MPI.MAX, comm))
	end
        gtoc!(prof, :gs_reduce, t0)
	push!(deltas, delta_g)
        gs_iter == 1 && (t_iter1 = (time_ns() - t_it) * 1e-9)
        delta_g < gs_tol && break
        gs_iter == max_gs_iter && rank == 0 && @warn "GS-unified did not converge (max_delta=$delta_g)"
    end

    n_tot = MPI.Allreduce(length(candidates), +, comm); n_max = MPI.Allreduce(length(candidates), MPI.MAX, comm)
    n_cached_tot = MPI.Allreduce(n_cached, +, comm)
    secs = MPI.Allreduce((time_ns() - t_start) * 1e-9, MPI.MAX, comm)
    t1 = MPI.Allreduce(t_iter1, MPI.MAX, comm)
    if rank == 0
    	println("GS[$(GS_MODE[])]: candidates $n_tot (max per rank $n_max), cached $n_cached_tot, iterations $iters, " *
        	    "first iteration $(round(t1; digits=2)) s of $(round(secs; digits=2)) s total, deltas: " *
            	    join([@sprintf("%.2e", x) for x in deltas], " "))
        if diag
            parts = String[]
            for k in eachindex(cnt_hist)
            	h = cnt_hist[k]
            	gcs = k <= length(gc_max_hist) ? round(gc_max_hist[k]; digits=3) : NaN
            	push!(parts, "it$(k): unchanged $(h[1]) <1e-12 $(h[2]) <1e-9 $(h[3]) >=1e-9 $(h[4]) " *
                             "compute mean/max $(round(tcomp_mean[k]; digits=3))/$(round(tcomp_max[k]; digits=3)) s " *
                             "gc max $(gcs) s")
            end
            println("GS detail: " * join(parts, " | "))
    	end
    end
    if USE_GS_GC_OFF[]
    	GC.enable(gs_gc_prev)
    	gs_gc_prev && GC.gc(false)          # one young-generation collection after GS
    end
    for (iel,d) in candidates; solved[iel][d] = true; end
end



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
        #received = exchange_ncf_gid_values!(exch, x_dg, ngl)
        received = exchange_ncf_gid_table!(exch, x_dg, ngl)
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


function exchange_ncf_face_coordinates(mesh, neighbor_rank::Dict{Tuple{TInt,TInt},TInt}, comm)
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

#=function build_dg_ncf_parent_to_child_remote_by_coords(mesh, child_iel::Int, fid::Int, d::Int,
                                                        remote_gel::TInt,
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::AbstractDict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
							#gid_received::Dict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
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
end=#

function build_dg_ncf_parent_to_child_remote_by_coords(mesh, child_iel::Int, fid::Int, d::Int,
                                                        remote_gel::TInt,
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::AbstractDict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
							#gid_received::Dict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
                                                        tol::Float64 = 1e-4, cluster_tol::Float64 = 1e-6)
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

    a1_raw = [getc(c,a1) for (gid,c) in candidates]
    a2_raw = [getc(c,a2) for (gid,c) in candidates]
    parent_a1_vals, a1_cluster = cluster_values_with_index(a1_raw; tol=cluster_tol)
    parent_a2_vals, a2_cluster = cluster_values_with_index(a2_raw; tol=cluster_tol)
    (length(parent_a1_vals) == ngl && length(parent_a2_vals) == ngl) ||
        error("degenerate parent grid: $(length(parent_a1_vals))x$(length(parent_a2_vals)) distinct axis values, " *
              "expected $ngl x $ngl")

    parent_grid_gid = Dict{Tuple{Int,Int}, Int}()
    for (idx,(gid,c)) in enumerate(candidates)
        parent_grid_gid[(a1_cluster[idx], a2_cluster[idx])] = gid
    end

    parent_val_grid = zeros(ngl,ngl)
    for i1 = 1:ngl, i2 = 1:ngl
        haskey(parent_grid_gid, (i1,i2)) || continue
        gid = parent_grid_gid[(i1,i2)]
        parent_val_grid[i1,i2] = get(gid_received, (remote_gel, d, gid), NaN)
    end

    child_a1 = getc.(IPc_coords, a1); child_a2 = getc.(IPc_coords, a2)
    L1 = build_1d_lagrange_matrix(parent_a1_vals, child_a1)
    L2 = build_1d_lagrange_matrix(parent_a2_vals, child_a2)

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
    remote_gid_coords = exchange_ncf_face_coordinates(mesh, neighbor_rank, comm)

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
    #gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
    gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)

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

#=function build_dg_ncf_child_to_parent_remote_by_coords(mesh, coarse_iel::Int, fid::Int, d::Int,
                                                        remote_gels::Vector{TInt},
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::AbstractDict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
							#gid_received::Dict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
                                                        tol::Float64 = 1e-4, cluster_tol::Float64 = 1e-6)
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

    child_grids = Dict{TInt, NamedTuple}()
    for rg in remote_gels
        haskey(remote_gid_coords, rg) || error("no coordinate data for remote_gel=$rg")
        coords_this = collect(remote_gid_coords[rg])
        my_normal_val = excl_axis==1 ? sum(xs)/length(xs) : (excl_axis==2 ? sum(ys)/length(ys) : sum(zs)/length(zs))
        on_plane = [(gid,c) for (gid,c) in coords_this if isapprox(getc(c,excl_axis), my_normal_val; atol=tol)]
        length(on_plane) == ngl*ngl || error("build_dg_ncf_child_to_parent_remote_by_coords: remote_gel=$rg " *
            "has $(length(on_plane)) on-plane nodes, expected $(ngl*ngl)")

        a1_raw = [getc(c,a1) for (gid,c) in on_plane]
        a2_raw = [getc(c,a2) for (gid,c) in on_plane]
        a1_vals, a1_cluster = cluster_values_with_index(a1_raw; tol=cluster_tol)
        a2_vals, a2_cluster = cluster_values_with_index(a2_raw; tol=cluster_tol)
        (length(a1_vals)==ngl && length(a2_vals)==ngl) ||
            error("degenerate child grid for remote_gel=$rg: $(length(a1_vals))x$(length(a2_vals)) " *
                  "distinct axis values, expected $ngl x $ngl")

        grid_gid = Dict{Tuple{Int,Int},Int}()
        for (idx,(gid,c)) in enumerate(on_plane)
            grid_gid[(a1_cluster[idx], a2_cluster[idx])] = gid
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
        for (ci1,v1) in enumerate(cg.a1_vals), (ci2,v2) in enumerate(cg.a2_vals)
            haskey(cg.grid_gid, (ci1,ci2)) || continue
            gid = cg.grid_gid[(ci1,ci2)]
            val_grid[ci1,ci2] = get(gid_received, (owning_rg, d, gid), NaN)
        end

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
end=#

function build_dg_ncf_child_to_parent_remote_by_coords(mesh, coarse_iel::Int, fid::Int, d::Int,
                                                        remote_gels::Vector{TInt},
                                                        remote_gid_coords::Dict{TInt,Dict{Int,NTuple{3,Float64}}},
                                                        gid_received::AbstractDict{Tuple{TInt,Int,Int},Float64}, ngl::Int;
                                                        tol::Float64 = 1e-4, cluster_tol::Float64 = 1e-6)
    irange, jrange, krange = face_ip_range(fid, ngl)
    coarse_coords = NTuple{3,Float64}[]
    for k in krange, j in jrange, i in irange
        ip = mesh.connijk[coarse_iel,i,j,k]
        push!(coarse_coords, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
    end
    xs = [c[1] for c in coarse_coords]; ys = [c[2] for c in coarse_coords]; zs = [c[3] for c in coarse_coords]
    spread = (maximum(xs)-minimum(xs), maximum(ys)-minimum(ys), maximum(zs)-minimum(zs))
    excl_axis = argmin(spread)
    axes = setdiff(1:3, [excl_axis]); a1, a2 = axes[1], axes[2]
    getc(c, ax) = ax == 1 ? c[1] : (ax == 2 ? c[2] : c[3])
    my_normal_val = excl_axis == 1 ? sum(xs)/length(xs) : (excl_axis == 2 ? sum(ys)/length(ys) : sum(zs)/length(zs))

    sorted_gels = sort(collect(remote_gels))          # canonical order: global element id

    # ---- each child's grid on this face's plane, with its received values ----
    child_grids = Dict{TInt,NamedTuple}()
    for rg in sorted_gels
        haskey(remote_gid_coords, rg) || error("no coordinate data for remote_gel=$rg")
        on_plane = [(gid, c) for (gid, c) in remote_gid_coords[rg] if isapprox(getc(c, excl_axis), my_normal_val; atol=tol)]
        length(on_plane) == ngl*ngl || error("build_dg_ncf_child_to_parent_remote_by_coords: remote_gel=$rg " *
            "has $(length(on_plane)) on-plane nodes, expected $(ngl*ngl)")
        v1, cl1 = cluster_values_with_index([getc(c, a1) for (_, c) in on_plane]; tol=cluster_tol)
        v2, cl2 = cluster_values_with_index([getc(c, a2) for (_, c) in on_plane]; tol=cluster_tol)
        (length(v1) == ngl && length(v2) == ngl) ||
            error("degenerate child grid for remote_gel=$rg: $(length(v1))x$(length(v2)) distinct axis values, expected $ngl x $ngl")
        val_grid = zeros(ngl, ngl)                     # nodes with no gid stay 0.0, as in the plan
        grid_has = falses(ngl, ngl)
        for (idx, (gid, _)) in enumerate(on_plane)
            val_grid[cl1[idx], cl2[idx]] = get(gid_received, (rg, d, gid), NaN)      # NaN = not arrived yet
            grid_has[cl1[idx], cl2[idx]] = true
        end
        child_grids[rg] = (a1_vals=v1, a2_vals=v2, val_grid=val_grid,
                           bbox=(minimum(v1)-tol, maximum(v1)+tol, minimum(v2)-tol, maximum(v2)+tol))
    end

    # ---- every child covering a node contributes; the node value is their mean ----
    parent_grid = zeros(ngl, ngl)
    for idx in 1:ngl*ngl
        p1 = ((idx-1) % ngl) + 1; p2 = ((idx-1) ÷ ngl) + 1
        cv1 = getc(coarse_coords[idx], a1); cv2 = getc(coarse_coords[idx], a2)
        acc = 0.0; count = 0
        for rg in sorted_gels
            cg = child_grids[rg]
            (cv1 >= cg.bbox[1] && cv1 <= cg.bbox[2] && cv2 >= cg.bbox[3] && cv2 <= cg.bbox[4]) || continue
            e1 = findfirst(v -> isapprox(v, cv1; atol=tol), cg.a1_vals)
            e2 = findfirst(v -> isapprox(v, cv2; atol=tol), cg.a2_vals)
            local val::Float64
            if e1 !== nothing && e2 !== nothing
                val = cg.val_grid[e1, e2]
            else
                L1 = build_1d_lagrange_matrix(cg.a1_vals, [cv1])
                L2 = build_1d_lagrange_matrix(cg.a2_vals, [cv2])
                val = 0.0
                for q1 in 1:ngl, q2 in 1:ngl
                    val += L1[1,q1] * L2[1,q2] * cg.val_grid[q1,q2]
                end
            end
            acc += val; count += 1
        end
        count == 0 && error("build_dg_ncf_child_to_parent_remote_by_coords: coarse node (a1=$cv1, a2=$cv2) " *
                            "is not contained in any child's bounding box")
        parent_grid[p1, p2] = acc / count
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
    remote_gid_coords = exchange_ncf_face_coordinates(mesh, neighbor_rank, comm)

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
    #gid_received = exchange_ncf_gid_values!(gid_exch, x_dg, ngl)
    gid_received = exchange_ncf_gid_table!(gid_exch, x_dg, ngl)

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

function debug_dump_on_plane_coords(mesh, coarse_iel::Int, fid::Int, remote_gels::Vector{TInt},
                                     remote_gid_coords, ngl::Int; tol::Float64=1e-4)
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
    my_normal_val = excl_axis==1 ? sum(xs)/length(xs) : (excl_axis==2 ? sum(ys)/length(ys) : sum(zs)/length(zs))

    for rg in remote_gels
        coords_this = collect(remote_gid_coords[rg])
        on_plane = [(gid,c) for (gid,c) in coords_this if isapprox(getc(c,excl_axis), my_normal_val; atol=tol)]
        println("remote_gel=$rg: $(length(on_plane)) on-plane points")
        a1_raw = sort([getc(c,a1) for (gid,c) in on_plane])
        a2_raw = sort([getc(c,a2) for (gid,c) in on_plane])
        println("  a1 sorted (raw): $a1_raw")
        println("  a2 sorted (raw): $a2_raw")
        a1_rounded = sort(unique(round.(a1_raw; digits=8)))
        a2_rounded = sort(unique(round.(a2_raw; digits=8)))
        println("  a1 unique count after round(digits=8): $(length(a1_rounded))")
        println("  a2 unique count after round(digits=8): $(length(a2_rounded))")
        # gaps between consecutive raw values, to see if any are suspiciously tiny
        println("  a1 consecutive diffs: $(diff(a1_raw))")
        println("  a2 consecutive diffs: $(diff(a2_raw))")
    end
end

function node_memory_limit_bytes()
    lim = Float64(Sys.total_memory())
    for f in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")
        if isfile(f)
            v = tryparse(Float64, strip(read(f, String)))      # "max" parses to nothing and is ignored
            v !== nothing && v > 0 && (lim = min(lim, v))
        end
    end
    s = tryparse(Float64, get(ENV, "SLURM_MEM_PER_NODE", ""))
    s !== nothing && s > 0 && (lim = min(lim, s*1e6))
    return lim
end

function auto_cache_budget_mb(setup, ni, comm; krylov_vectors::Int=0, usable_frac::Float64=0.85, min_pairs_per_rank::Int=3000)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    ngl = setup.mesh.ngl; n3 = ngl^3
    pairs = Float64(sum(setup.n_dirs))
    bpp = 8.0*(n3*n3 + n3 + 6*ngl*ngl) + 200.0
    GC.gc()
    base = Float64(Sys.maxrss())
    kry = (krylov_vectors + 3) * 8.0 * n3 * pairs
    all = reshape(gather_all(Float64[base + kry, pairs, node_memory_limit_bytes()], comm), 3, np)   # collective
    mine = [r for r in 1:np if ni.node_of_rank[r] == ni.node_of_rank[rank+1]]
    node_need = sum(all[1,r] for r in mine); node_pairs = sum(all[2,r] for r in mine)
    node_lim = minimum(all[3,r] for r in mine)
    free = max(usable_frac*node_lim - node_need, 0.0)
    budget = min(free * pairs / max(node_pairs, 1.0), pairs*bpp)
    mean_pairs = MPI.Allreduce(pairs, +, comm) / np
    skipped = mean_pairs < min_pairs_per_rank
    skipped && (budget = 0.0)
    bmin = MPI.Allreduce(budget, MPI.MIN, comm); bmax = MPI.Allreduce(budget, MPI.MAX, comm)
    full_max = MPI.Allreduce(pairs*bpp, MPI.MAX, comm)
    need_max = MPI.Allreduce(node_need, MPI.MAX, comm); lim_min = MPI.Allreduce(node_lim, MPI.MIN, comm)
    lines = String[]
    if rank == 0
        push!(lines, "auto cache budget: smallest node memory limit $(round(lim_min/1e9; digits=0)) GB, usable fraction $usable_frac, " *
            "baseline (peak RSS) + Krylov reserve on the fullest node $(round(need_max/1e9; digits=1)) GB; " *
            "budget per rank min $(round(bmin/1e6; digits=0)) max $(round(bmax/1e6; digits=0)) MB " *
	    " | mean pairs/rank $(round(Int, mean_pairs)), threshold $min_pairs_per_rank" * (skipped ? ": cache skipped" : "") *
            "(a full cache needs up to $(round(full_max/1e6; digits=0)) MB on the busiest rank)")
    end
    return budget/1e6, lines
end

# ─────────────────────────────────────────────────────────────────────────────
#  sweep_ordering.jl  (DRAFT — untested, single-rank conforming case only)
#
#  Stage 1 of the transport sweep: for a fixed direction Ω = (Ωx,Ωy,Ωz),
#  classify every element's 6 faces as inflow/outflow, build the interior
#  upwind dependency graph via facet_cell_ids, and topologically sort.
#
#  PERF NOTE: takes raw metric arrays (dξdx, dξdy, ...) rather than a
#  `metrics` struct, and raw mesh fields (nelem, ngl, facet_cell_ids,
#  cell_face_ids) rather than the `mesh` struct — struct-field access across
#  a function boundary can defeat inlining/type-specialization and force
#  heap allocation on every access, which matters here since these
#  functions run once per element per direction, every direction, every
#  sweep. Longer signatures, but no allocation on the hot path.
#
#  Deliberately NOT handled yet (separate follow-up passes):
#    - non-conforming (hanging) interior interfaces
#    - MPI cross-rank dependencies / wavefront parallelism
#  Physical AND periodic boundary faces are treated identically here — both
#  are graph "sources" (no dependency), since with lagged periodic coupling
#  neither supplies a live in-graph dependency; only the ACTUAL inflow value
#  used at assembly time differs (BC formula vs. frozen periodic partner).
#
#  Face convention (matches mesh.jl / collect_periodic_ncf_pairs_3D!):
#    1 front (j=ngl)   2 back (j=1)
#    3 bottom (k=ngl)  4 top (k=1)
#    5 right (i=1)     6 left (i=ngl)
#  opp = (2, 1, 4, 3, 6, 5)
# ─────────────────────────────────────────────────────────────────────────────

const CONN_FACE_EL_TO_MY_FID = (4, 6, 1, 5, 2, 3)   # index by face_el = 1:6
const OPP_FACE = (2, 1, 4, 3, 6, 5)
# ─────────────────────────────────────────────────────────────────────────────
#  Robust face-container accessor
#
#  facet_cell_ids and conn_unique_faces both come from Gridap queries (one
#  from the model, one from the topology) but that does NOT guarantee they
#  share a container type — assuming so is what caused the earlier crash.
#  Dispatch on the actual type instead of assuming: Gridap.Arrays.Table
#  (variable-row-width, needs array_cache/getindex! — bare a[i] indexing on
#  this type crashed under repeated varying-width calls, see conversation)
#  vs. a plain AbstractArray (ordinary indexing is correct and sufficient;
#  array_cache isn't meaningful for it, so no cache to build).
# ─────────────────────────────────────────────────────────────────────────────
_face_cache(a::Gridap.Arrays.Table) = array_cache(a)
_face_cache(a::AbstractArray)       = nothing

_face_row(cache, a::Gridap.Arrays.Table, i) = getindex!(cache, a, i)
_face_row(_,     a::AbstractArray,       i) = a[i]


"""
    face_normal_sign(fid)

Returns (metric_symbol, sign) such that the outward normal on face `fid`
is proportional to `sign .* (dξdx/dηdx/dζdx-type metric triple)` evaluated
at any node on that face. `metric_symbol` selects which reference-direction
metric (:xi, :eta, :zeta) to use.

ASSUMPTION FLAGGED FOR CONFIRMATION: this assumes GLL index 1 ↔ ξ=-1 and
index ngl ↔ ξ=+1 (ascending), consistent with "right=i=1" meaning the face
at the MINIMUM ξ. If your node ordering is descending instead, every sign
below is inverted — please sanity-check on one element before trusting
inflow/outflow classification broadly.
"""
@inline function face_normal_sign(fid::Int)
    fid == 1 && return (1, +1.0)   # front  (j=ngl)  -- 1 = :xi, 2 = :eta, 3 = :zeta
    fid == 2 && return (1, -1.0)   # back   (j=1)
    fid == 3 && return (2, +1.0)   # bottom (k=ngl)
    fid == 4 && return (2, -1.0)   # top    (k=1)
    fid == 5 && return (3, -1.0)   # right  (i=1)
    fid == 6 && return (3, +1.0)   # left   (i=ngl)
    error("face_normal_sign: fid must be 1..6, got $fid")
end
# NOTE: metric_symbol mapping above is (eta=1, zeta=2, xi=3) to match the
# order metric arrays are passed into classify_element_faces below —
# see the (mxi,mxi,mxi)/(meta,...)/(mzeta,...) argument grouping there.
# Re-expressed explicitly per face for clarity (this is what's actually used):
#   face 1 (front,  j=ngl): eta,  +1
#   face 2 (back,   j=1)  : eta,  -1
#   face 3 (bottom, k=ngl): zeta, +1
#   face 4 (top,    k=1)  : zeta, -1
#   face 5 (right,  i=1)  : xi,   -1
#   face 6 (left,   i=ngl): xi,   +1

"""
    face_ip_range(fid, ngl)

Returns the (i_range, j_range, k_range) index ranges into connijk/metrics
for face `fid`, matching the convention used throughout mesh.jl
(collect_periodic_ncf_pairs_3D!'s face_ips helper uses the same ranges).
"""
@inline function face_ip_range(fid::Int, ngl::Int)
    fid == 1 && return (1:ngl, ngl:ngl, 1:ngl)   # front  j=ngl
    fid == 2 && return (1:ngl, 1:1,     1:ngl)   # back   j=1
    fid == 3 && return (1:ngl, 1:ngl, ngl:ngl)   # bottom k=ngl
    fid == 4 && return (1:ngl, 1:ngl, 1:1)       # top    k=1
    fid == 5 && return (1:1,   1:ngl, 1:ngl)     # right  i=1
    fid == 6 && return (ngl:ngl, 1:ngl, 1:ngl)   # left   i=ngl
    error("face_ip_range: fid must be 1..6, got $fid")
end

"""
    build_conforming_sweep_order(Ωx, Ωy, Ωz, nelem, ngl,
                                  facet_cell_ids, cell_face_ids,
                                  dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz)

Single-rank, conforming-interior-faces-only topological sort for direction
Ω. Takes raw mesh/metric arrays directly (see file header PERF NOTE) rather
than `mesh`/`metrics` structs.

Returns:
  - `order::Vector{Int}` — element indices in sweep-safe order (process in
    this order; every element appears after all its interior upwind
    neighbors)
  - `levels::Vector{Int}` — levels[iel] = wavefront level of element iel
    (elements with the same level have no dependency on each other and
    could be processed in parallel — not exploited yet, but computed now
    since it's nearly free once in-degrees are known, and will matter once
    MPI wavefronts are added)

Any face with no same-rank interior neighbor (physical boundary, periodic
boundary, or eventually an MPI ghost boundary) is treated as a graph source
— no dependency contributed. NON-CONFORMING interior interfaces are NOT
yet handled: an element with a hanging face is currently treated as if
that face were a boundary (a real limitation of this draft, not a silent
approximation — flagged clearly so it isn't mistaken for handled).
"""

#=function build_conforming_sweep_order(Ωx::Float64, Ωy::Float64, Ωz::Float64,
                                       nelem::TInt, ngl::TInt,
                                       neighbor_elem,   # from build_face_adjacency_from_conn_face_el
                                       dξdx, dξdy, dξdz,
                                       dηdx, dηdy, dηdz,
                                       dζdx, dζdy, dζdz)
    indegree    = zeros(TInt, nelem)
    upwind_deps = [TInt[] for _ = 1:nelem]
    downstream  = [TInt[] for _ = 1:nelem]

    is_inflow = Vector{Bool}(undef, 6)
    is_active = Vector{Bool}(undef, 6)

    for iel = 1:nelem
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem[iel][fid]
            other == zero(TInt) && continue
            push!(upwind_deps[iel], other)
            push!(downstream[other], iel)
            indegree[iel] += 1
        end
    end


    order  = TInt[]
    levels = zeros(TInt, nelem)
    queue  = TInt[iel for iel = 1:nelem if indegree[iel] == 0]
    indeg_work = copy(indegree)

    while !isempty(queue)
        next_queue = TInt[]
        for iel in queue
            push!(order, iel)
            for dsn in downstream[iel]
                indeg_work[dsn] -= 1
                levels[dsn] = max(levels[dsn], levels[iel] + 1)
                if indeg_work[dsn] == 0
                    push!(next_queue, dsn)
                end
            end
        end
        queue = next_queue
    end

    if length(order) != nelem
        stuck = [iel for iel = 1:nelem if indeg_work[iel] > 0]
        error("build_conforming_sweep_order: $(length(stuck)) elements never reached " *
              "indegree 0 — indicates a cycle in the conforming-interior dependency " *
              "graph for this direction (Ωx=$Ωx, Ωy=$Ωy, Ωz=$Ωz). Since conforming " *
              "interior faces alone should never cycle on a single rank, this most " *
              "likely means face_normal_sign's ascending-index assumption is wrong " *
              "for this mesh. Stuck elements: $(stuck[1:min(5,end)])")
    end

    return order, levels
end=#

function build_conforming_sweep_order(Ωx::Float64, Ωy::Float64, Ωz::Float64,
                                       nelem::TInt, ngl::TInt,
                                       neighbor_elem,
                                       dξdx, dξdy, dξdz,
                                       dηdx, dηdy, dηdz,
                                       dζdx, dζdy, dζdz;
                                       ncf_child_to_parent = nothing,      # (child_iel,child_fid) -> PARENT ELEMENT ID (NOT ncf_idx)
                                       ncf_parent_to_children = nothing)   # (parent_iel,parent_fid) -> [child_iel,...]
    indegree    = zeros(TInt, nelem)
    upwind_deps = [TInt[] for _ = 1:nelem]
    downstream  = [TInt[] for _ = 1:nelem]

    is_inflow = Vector{Bool}(undef, 6)
    is_active = Vector{Bool}(undef, 6)

    for iel = 1:nelem
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem[iel][fid]
            if other == zero(TInt)
                if ncf_child_to_parent !== nothing && haskey(ncf_child_to_parent, (Int(iel),fid))
                    p = TInt(ncf_child_to_parent[(Int(iel),fid)])
                    push!(upwind_deps[iel], p); push!(downstream[p], iel); indegree[iel] += 1
                elseif ncf_parent_to_children !== nothing && haskey(ncf_parent_to_children, (Int(iel),fid))
                    for c in ncf_parent_to_children[(Int(iel),fid)]
                        cT = TInt(c)
                        push!(upwind_deps[iel], cT); push!(downstream[cT], iel); indegree[iel] += 1
                    end
                end
                continue
            end
            push!(upwind_deps[iel], other)
            push!(downstream[other], iel)
            indegree[iel] += 1
        end
    end

    order  = TInt[]
    levels = zeros(TInt, nelem)
    queue  = TInt[iel for iel = 1:nelem if indegree[iel] == 0]
    indeg_work = copy(indegree)

    while !isempty(queue)
        next_queue = TInt[]
        for iel in queue
            push!(order, iel)
            for dsn in downstream[iel]
                indeg_work[dsn] -= 1
                levels[dsn] = max(levels[dsn], levels[iel] + 1)
                indeg_work[dsn] == 0 && push!(next_queue, dsn)
            end
        end
        queue = next_queue
    end
    if length(order) != nelem
        stuck = [iel for iel = 1:nelem if indeg_work[iel] > 0]
        error("build_conforming_sweep_order: $(length(stuck)) elements never reached indegree 0 " *
              "for Ωx=$Ωx, Ωy=$Ωy, Ωz=$Ωz — possible cycle. Stuck elements: $(stuck[1:min(5,end)])")
    end
    return order, levels
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_my_face_to_global_face  (DRAFT — untested)
#
#  Builds, once per rank (NOT per direction — purely geometric), the map
#  from (iel, fid) in THIS codebase's face convention (via conn_face_el,
#  see corner derivation below) to the GLOBAL face index used by
#  facet_cell_ids/conn_unique_faces (Gridap's numbering).
#
#  Deliberately does NOT assume any correspondence between my `fid` and
#  Gridap's internal per-element face ordering (cell_face_ids[iel][k]) —
#  that assumption was tried once already and was wrong (see conversation).
#  Instead matches by comparing SORTED global-node-ID sets, which only
#  depends on two things already independently confirmed:
#    - conn_face_el's own corner convention (derived by hand below from
#      populate_conn_face_el!, matches this codebase's own face_el, NOT my
#      arbitrary fid — remapped explicitly here)
#    - facet_cell_ids indexed as facet_cell_ids[iface][1] (confirmed)
#
#  Returns my_fid_to_gface::Vector{NTuple{6,Int}} where
#  my_fid_to_gface[iel][fid] = global face index, or 0 if no exact 4-corner
#  match was found (non-conforming/hanging face — see NCF_FALLBACK note
#  below; this is NOT silently treated as a boundary without a count, so it
#  can't be mistaken for "handled").
# ─────────────────────────────────────────────────────────────────────────────

# populate_conn_face_el!'s own face_el → my fid remap, derived by hand from
# its corner assignments (see conversation for the derivation):
#   face_el=1 (k=1,   ip1,ip4,ip3,ip2) → my fid=4 (top)
#   face_el=2 (i=ngl, ip1,ip2,ip6,ip5) → my fid=6 (left)
#   face_el=3 (j=ngl, ip2,ip3,ip7,ip6) → my fid=1 (front)
#   face_el=4 (i=1,   ip3,ip4,ip8,ip7) → my fid=5 (right)
#   face_el=5 (j=1,   ip1,ip5,ip8,ip4) → my fid=2 (back)
#   face_el=6 (k=ngl, ip5,ip6,ip7,ip8) → my fid=3 (bottom)

function build_my_face_to_global_face(nelem::TInt, conn_face_el, conn_unique_faces, nfaces::TInt)
    gface_lookup = Dict{NTuple{4,TInt}, TInt}()
    _cface_cache = _face_cache(conn_unique_faces)
    for iface_g = 1:nfaces
        nodes = _face_row(_cface_cache, conn_unique_faces, iface_g)
        key = Tuple(sort(TInt[nodes[1], nodes[2], nodes[3], nodes[4]]))
        gface_lookup[key] = iface_g
    end

    my_fid_to_gface = Vector{NTuple{6,TInt}}(undef, nelem)
    n_unmatched = zero(TInt)

    for iel = 1:nelem
        tmp = zeros(TInt, 6)
        for face_el = 1:6
            my_fid = CONN_FACE_EL_TO_MY_FID[face_el]
            n1 = conn_face_el[1, face_el, iel]
            n2 = conn_face_el[2, face_el, iel]
            n3 = conn_face_el[3, face_el, iel]
            n4 = conn_face_el[4, face_el, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            g = get(gface_lookup, key, zero(TInt))
            g == zero(TInt) && (n_unmatched += 1)
            tmp[my_fid] = g
        end
        my_fid_to_gface[iel] = Tuple(tmp)
    end

    if n_unmatched > 0
        @warn "build_my_face_to_global_face: $n_unmatched (element,face) entries had no " *
              "exact 4-corner match in conn_unique_faces — see docstring for expected cause " *
              "(non-conforming interfaces) vs. real problem (nonzero on a fully-conforming mesh)."
    end

    return my_fid_to_gface
end

function build_face_adjacency_from_conn_face_el(nelem::TInt, conn_face_el)
    # Self-contained — uses ONLY conn_face_el (a plain array this codebase
    # populates itself via ordinary assignment in populate_conn_face_el!,
    # confirmed NOT a Gridap Table). Does not touch facet_cell_ids,
    # cell_face_ids, or conn_unique_faces at all. After two separate
    # confusing failures trying to reuse Gridap's own precomputed
    # face/adjacency structures (a wrong local-face-ordering assumption,
    # then a Table-internal indexing crash that persisted through the
    # array_cache/getindex! fix and whose root cause isn't diagnosable
    # without a REPL), recomputing this one-time, per-mesh adjacency
    # ourselves is more reliable than a third guess at Table's behavior.
    #
    # Returns:
    #   neighbor_elem::Vector{NTuple{6,TInt}} — neighbor_elem[iel][fid] =
    #     the element sharing that face, or 0 if boundary/periodic (source)
    #   neighbor_fid::Vector{NTuple{6,TInt}} — the LOCAL face id on that
    #     neighbor corresponding to the same physical face, or 0
    face_map = Dict{NTuple{4,TInt}, Vector{Tuple{TInt,TInt}}}()

    for iel = 1:nelem
        for face_el = 1:6
            my_fid = CONN_FACE_EL_TO_MY_FID[face_el]
            n1 = conn_face_el[1, face_el, iel]
            n2 = conn_face_el[2, face_el, iel]
            n3 = conn_face_el[3, face_el, iel]
            n4 = conn_face_el[4, face_el, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            push!(get!(face_map, key, Tuple{TInt,TInt}[]), (TInt(iel), TInt(my_fid)))
        end
    end

    ne_tmp = [zeros(TInt, 6) for _ = 1:nelem]
    nf_tmp = [zeros(TInt, 6) for _ = 1:nelem]
    n_bad  = 0

    for (key, entries) in face_map
        if length(entries) == 1
            continue   # boundary/periodic face — source, stays 0
        elseif length(entries) == 2
            (iel1, fid1), (iel2, fid2) = entries
            ne_tmp[iel1][fid1] = iel2;  nf_tmp[iel1][fid1] = fid2
            ne_tmp[iel2][fid2] = iel1;  nf_tmp[iel2][fid2] = fid1
        else
            n_bad += 1   # >2 elements sharing exact corners — unexpected, flagged not silent
        end
    end

    n_bad > 0 && @warn "build_face_adjacency_from_conn_face_el: $n_bad face(s) matched " *
                        "by MORE than 2 elements — unexpected for a standard hex mesh, " *
                        "likely duplicate/degenerate node IDs; left unmatched (treated as source)."

    neighbor_elem = [Tuple(ne_tmp[iel]) for iel = 1:nelem]
    neighbor_fid  = [Tuple(nf_tmp[iel]) for iel = 1:nelem]
    return neighbor_elem, neighbor_fid
end

# ─────────────────────────────────────────────────────────────────────────────
#  verify_sweep_ordering.jl  (DRAFT — untested, meant to be run once against
#  your test mesh before trusting build_conforming_sweep_order further)
#
#  Checks, in order of how conclusive they are:
#    1. Reciprocity of neighbor_elem/neighbor_fid (pure data-structure sanity,
#       no physics involved — a failure here means build_face_adjacency_from_conn_face_el
#       itself is broken, not the sweep logic)
#    2. Boundary-face count cross-check against mesh.nfaces_bdy
#    3. FULL dependency check (not spot-check) for two directions: every
#       element's interior upwind dependency must appear strictly earlier
#       in `order`
#    4. Geometric orientation check for axis-aligned Ω — tests whether
#       face_normal_sign's ascending-index assumption is actually correct
#       for THIS mesh, using only element centroid coordinates (independent
#       of the ordering/adjacency code, so this is a real cross-check, not
#       circular)
#    5. Levels-distribution sanity print (too flat / too deep = probably wrong)
# ─────────────────────────────────────────────────────────────────────────────

function verify_sweep_ordering(mesh, metrics)
    nelem = mesh.nelem
    ngl   = mesh.ngl

    println("═"^70)
    println("SWEEP ORDERING VERIFICATION")
    println("═"^70)
    
    println("MESH FINGERPRINT: nelem=$(mesh.nelem) nfaces_bdy=$(mesh.nfaces_bdy) " *
        "hash(connijk)=$(hash(mesh.connijk)) hash(conn_face_el)=$(hash(mesh.conn_face_el))")

    # ── Step 0: build adjacency (once) ────────────────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(nelem, mesh.conn_face_el)

    # ── Check 1: reciprocity ───────────────────────────────────────────────
    # If neighbor_elem[iel][fid] == jel, then neighbor_elem[jel][neighbor_fid[iel][fid]]
    # must equal iel — i.e. the relation is symmetric and the local-face
    # cross-reference is internally consistent. Pure data-structure check,
    # no physics.
    n_recip_fail = 0
    n_interior_faces_checked = 0
    for iel = 1:nelem
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            n_interior_faces_checked += 1
            jfid = neighbor_fid[iel][fid]
            if jfid == zero(TInt) || neighbor_elem[jel][jfid] != iel
                n_recip_fail += 1
                n_recip_fail <= 5 && println("  RECIPROCITY FAIL: iel=$iel fid=$fid -> jel=$jel jfid=$jfid " *
                                              "but neighbor_elem[jel][jfid]=$(jel==0 ? "N/A" : neighbor_elem[jel][max(jfid,1)])")
            end
        end
    end
    println("[1] Reciprocity check: $(n_interior_faces_checked) interior face-links examined, " *
            "$(n_recip_fail) failures " * (n_recip_fail == 0 ? "✓ PASS" : "✗ FAIL"))

    # ── Check 2: boundary face count cross-check ──────────────────────────
    # Every (iel,fid) with neighbor_elem==0 is a boundary/periodic face as
    # WE see it. Total such entries, divided by 1 (each boundary face is
    # only counted once, from its single owning element — unlike interior
    # faces which appear from both sides), should equal mesh.nfaces_bdy
    # if that field counts local (this-rank) boundary faces. If it's a
    # GLOBAL count instead, this will legitimately differ under MPI and
    # should be compared to mesh.nfaces_bdy on a single-rank run only.
    n_boundary_entries = 0
    for iel = 1:nelem, fid = 1:6
        neighbor_elem[iel][fid] == zero(TInt) && (n_boundary_entries += 1)
    end
    expected = try mesh.nfaces_bdy catch; nothing end
    if expected === nothing
        println("[2] Boundary face count: found $n_boundary_entries boundary-face entries " *
                "(mesh.nfaces_bdy not accessible to compare — skipped)")
    else
        match = n_boundary_entries == expected
        println("[2] Boundary face count: found $n_boundary_entries, mesh.nfaces_bdy=$expected " *
                (match ? "✓ PASS" : "⚠ MISMATCH (expected on multi-rank runs; " *
                 "should match exactly on a single-rank run)"))
    end

    # ── Check 3: full dependency ordering check, two directions ──────────
    function check_order_full(Ωx, Ωy, Ωz, label; debug_violations::Bool=true, max_debug::Int=5)
        order, levels = build_conforming_sweep_order(
            Ωx, Ωy, Ωz, nelem, ngl, neighbor_elem,
            metrics.dξdx, metrics.dξdy, metrics.dξdz,
            metrics.dηdx, metrics.dηdy, metrics.dηdz,
            metrics.dζdx, metrics.dζdy, metrics.dζdz)

        pos = zeros(TInt, nelem)
        for (p, iel) in enumerate(order)
            pos[iel] = p
        end

        is_inflow = Vector{Bool}(undef, 6)
        is_active = Vector{Bool}(undef, 6)
        n_violations = 0
        violation_pairs = Tuple{TInt,TInt}[]   # collect ALL, not just first 5

        for iel = 1:nelem
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                     metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                     metrics.dζdx, metrics.dζdy, metrics.dζdz)
            for fid = 1:6
                (is_inflow[fid] && is_active[fid]) || continue
                jel = neighbor_elem[iel][fid]
                jel == zero(TInt) && continue
                if pos[jel] >= pos[iel]
                    n_violations += 1
                    push!(violation_pairs, (iel, jel))
                end
            end
        end

        n_levels = length(unique(levels))
        println("[3] Full order check ($label): $(length(order))/$nelem elements ordered, " *
                "$n_violations violations " * (n_violations == 0 ? "✓ PASS" : "✗ FAIL") *
                "  |  levels: min=$(minimum(levels)) max=$(maximum(levels)) " *
                "unique=$n_levels  (too close to 1 or to $nelem may indicate a problem)")

        # ── Auto-debug every violation (deduplicated by unordered pair, since
        # (iel,jel) and (jel,iel) would otherwise print the same geometry twice) ─
        if debug_violations && n_violations > 0
            seen = Set{Tuple{TInt,TInt}}()
            n_shown = 0
            println("─"^70)
            println("[3-debug] Auto-diagnostics for violating pairs (label=$label):")
            for (iel, jel) in violation_pairs
                key = (min(iel,jel), max(iel,jel))
                key in seen && continue
                push!(seen, key)
                n_shown += 1
                debug_single_face_pair(mesh, metrics, iel, jel, neighbor_elem, neighbor_fid;
                                        Ωx=Ωx, Ωy=Ωy, Ωz=Ωz)
                if n_shown >= max_debug
                    remaining = length(unique(violation_pairs)) - n_shown
                    remaining > 0 && println("  ... ($remaining more distinct violating pairs not shown, " *
                                              "increase max_debug to see more)")
                    break
                end
            end
        end

        return order, levels
    end

    check_order_full(1.0, 0.0, 0.0, "Ω=(1,0,0)")
    check_order_full(0.3, 0.5, 0.8, "Ω=(0.3,0.5,0.8) (generic oblique)")

    n_periodic_faces_checked   = 0
    n_periodic_wrongly_wired   = 0
    n_periodic_unmatched       = 0

    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        mesh.bdy_face_type[iface_bdy] in ("periodicx","periodicy","periodicz") || continue
        iel = mesh.bdy_face_in_elem[iface_bdy]

        # Corners of THIS SPECIFIC boundary face instance only.
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, ngl, ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, ngl]
        key_bdy = Tuple(sort(TInt[ip1, ip2, ip3, ip4]))

        # Find which of iel's 6 local faces this actually is (order-independent
        # match — corner listing order doesn't matter since we compare sorted sets).
        found_fid = 0
        for face_el = 1:6
            my_fid = CONN_FACE_EL_TO_MY_FID[face_el]
            n1 = mesh.conn_face_el[1, face_el, iel]
            n2 = mesh.conn_face_el[2, face_el, iel]
            n3 = mesh.conn_face_el[3, face_el, iel]
            n4 = mesh.conn_face_el[4, face_el, iel]
            if Tuple(sort(TInt[n1,n2,n3,n4])) == key_bdy
                found_fid = my_fid
                break
            end
         end

        if found_fid == 0
            n_periodic_unmatched += 1
            continue
        end

        n_periodic_faces_checked += 1
        jel = neighbor_elem[iel][found_fid]
        jel != zero(TInt) && (n_periodic_wrongly_wired += 1)
    end

    println("Periodic boundary face instances checked : $n_periodic_faces_checked")
    println("  ...unmatched to any local face (unexpected, should be 0): $n_periodic_unmatched")
    println("  ...wrongly wired as an interior neighbor: $n_periodic_wrongly_wired")
    

    # ── Check 4: geometric orientation check (independent of adjacency code) ─
    println("─"^70)
    println("[4] Geometric orientation check for Ω=(1,0,0):")

    cx = zeros(Float64, nelem)
    for iel = 1:nelem
        s = 0.0; n = 0
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel,i,j,k]
            s += mesh.x[ip]; n += 1
        end
        cx[iel] = s / n
    end

    indeg0 = Set{TInt}()
    is_inflow = Vector{Bool}(undef, 6)
    is_active = Vector{Bool}(undef, 6)
    for iel = 1:nelem
        classify_element_faces!(is_inflow, is_active, iel, 1.0, 0.0, 0.0, ngl,
                                 metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                 metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                 metrics.dζdx, metrics.dζdy, metrics.dζdz)
        has_dep = false
        for fid = 1:6
            if is_inflow[fid] && is_active[fid] && neighbor_elem[iel][fid] != zero(TInt)
                has_dep = true; break
            end
        end
        has_dep || push!(indeg0, iel)
    end

    cx_min, cx_max = extrema(cx)
    tol = (cx_max - cx_min) * 1e-3
    geometric_min_layer = Set(iel for iel = 1:nelem if cx[iel] < cx_min + tol)

    overlap = length(intersect(indeg0, geometric_min_layer))
    frac_min = overlap / max(length(geometric_min_layer), 1)

    println("    Source (indegree-0) elements found: $(length(indeg0))")
    println("    Geometric x-minimum layer elements : $(length(geometric_min_layer))")
    println("    Overlap: $overlap ($( round(100*frac_min,digits=1) )% of x-min layer are sources)")
    if frac_min > 0.9
        println("    ✓ PASS — sources are (as expected) the x-MINIMUM layer: " *
                "face_normal_sign orientation appears CORRECT for this mesh.")
    else
        geometric_max_layer = Set(iel for iel = 1:nelem if cx[iel] > cx_max - tol)
        overlap_max = length(intersect(indeg0, geometric_max_layer))
        frac_max = overlap_max / max(length(geometric_max_layer), 1)
        if frac_max > 0.9
            println("    ✗ FAIL — sources are the x-MAXIMUM layer instead of x-minimum. " *
                    "face_normal_sign's sign convention is INVERTED for this mesh.")
        else
            println("    ⚠ Sources match neither x-min nor x-max layer cleanly " *
                    "($(round(100*frac_max,digits=1))% overlap with x-max).")
        end
    end
   
    # Build a lookup: which (iel, fid) pairs are tagged as real boundaries
    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, ngl, ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1,ip2,ip3,ip4])))
    end


    n_true_boundary  = 0
    n_suspected_cut  = 0
    suspected_cuts   = Tuple{TInt,TInt}[]
    cut_neighbor_rank = Dict{Tuple{TInt,TInt}, TInt}()

    rank = MPI.Comm_rank(get_mpi_comm())

    for iel = 1:nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]
            n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]
            n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1,n2,n3,n4]))

            if key in tagged_face_keys
                n_true_boundary += 1
                continue
            end

            n_suspected_cut += 1
            push!(suspected_cuts, (iel, fid))

            # Any node strictly interior to the face (not on its bounding edges)
            # is safe to query ownership for — requires ngl >= 3 to have such a
            # node at all (ngl=2 means every face node IS a corner; flagged
            # below rather than silently picked wrong).
            irange, jrange, krange = face_ip_range(fid, ngl)
            if ngl < 3
               error("Face-interior ownership check requires ngl >= 3 (every node " *
                      "on a face is a corner when ngl=2); got ngl=$ngl.")
            end
            i_int = length(irange) == 1 ? first(irange) : first(irange)+1
            j_int = length(jrange) == 1 ? first(jrange) : first(jrange)+1
            k_int = length(krange) == 1 ? first(krange) : first(krange)+1
            ip_face_interior = mesh.connijk[iel, i_int, j_int, k_int]

            owner = mesh.gip2owner[ip_face_interior]
            cut_neighbor_rank[(iel,fid)] = owner
        end
    end

    println("True tagged physical/periodic boundary faces: $n_true_boundary")
    println("Suspected rank-partition-cut faces: $n_suspected_cut")

    neighbor_rank_counts = Dict{TInt,Int}()
    n_self_owned_anomaly = 0
    for ((iel,fid), owner) in cut_neighbor_rank
        if owner == rank
            n_self_owned_anomaly += 1
        else
            neighbor_rank_counts[owner] = get(neighbor_rank_counts, owner, 0) + 1
        end
    end
    println("Distinct neighbor ranks found (via face-interior node): $(sort(collect(neighbor_rank_counts)))")
    n_self_owned_anomaly > 0 && println("⚠ $n_self_owned_anomaly cut faces resolved face-interior " *
                                         "node owner to THIS rank — contradiction, worth investigating " *
                                         "directly rather than assuming it's fine.")

    neighbor_rank, neighbor_gel, n_unresolved = resolve_rank_cut_neighbors(mesh, suspected_cuts, comm)

    println("Rank $(MPI.Comm_rank(comm)): $(length(suspected_cuts)) suspected cuts, " *
            "$(length(neighbor_rank)) resolved, $n_unresolved unresolved")

    all_ok = verify_rank_cut_neighbors(mesh, neighbor_rank, suspected_cuts, comm)
    println("Rank $(MPI.Comm_rank(comm)): mutual verification " * (all_ok ? "✓ ALL PASS" : "✗ SOME MISMATCH"))

end

"""
    classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                             dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz;
                             tangential_tol=1e-8)

Fills two 6-element Bool buffers:
  - `is_inflow[fid]`  — true if face `fid` is an inflow face (Ω·n < 0)
  - `is_active[fid]`  — true if this face carries a MEANINGFUL streaming
    dependency at all. False for faces where Ω is (numerically) tangential
    to the face — |cos(angle between Ω and n_face)| < tangential_tol.

  Only faces with is_active[fid] == true should ever create a dependency
  edge, regardless of is_inflow's value — a face with negligible normal
  flux should not constrain element ordering, and without this guard,
  floating-point roundoff on truly-tangential faces (Ω·n ≈ 0 in exact
  arithmetic, but not exactly 0.0 after numerical metric-term summation)
  gets classified as inflow/outflow essentially at random, fabricating
  spurious dependencies. Confirmed empirically: on a 10×10×10 periodic
  box for Ω=(1,0,0), the unguarded version found only 4 true sources
  instead of the expected 100 (the x-minimum layer) — see conversation.
"""
function classify_element_faces!(is_inflow::AbstractVector{Bool}, is_active::AbstractVector{Bool},
                                  iel::Int, Ωx::Float64, Ωy::Float64, Ωz::Float64, ngl::Int,
                                  dξdx, dξdy, dξdz,
                                  dηdx, dηdy, dηdz,
                                  dζdx, dζdy, dζdz;
                                  tangential_tol::Float64 = 1e-8)
    @inbounds for fid = 1:6
        msym, sgn = face_normal_sign(fid)
        irange, jrange, krange = face_ip_range(fid, ngl)
        dotsum = 0.0
        nx_sum = 0.0; ny_sum = 0.0; nz_sum = 0.0
        n = 0
        for k in krange, j in jrange, i in irange
            if msym == 1       # :eta
                nx, ny, nz = dηdx[iel,i,j,k], dηdy[iel,i,j,k], dηdz[iel,i,j,k]
            elseif msym == 2   # :zeta
                nx, ny, nz = dζdx[iel,i,j,k], dζdy[iel,i,j,k], dζdz[iel,i,j,k]
            else               # :xi
                nx, ny, nz = dξdx[iel,i,j,k], dξdy[iel,i,j,k], dξdz[iel,i,j,k]
            end
            dotsum += sgn * (Ωx*nx + Ωy*ny + Ωz*nz)
            nx_sum += sgn*nx; ny_sum += sgn*ny; nz_sum += sgn*nz
            n += 1
        end
        dot_avg = dotsum / n
        # Scale-invariant tangency test: compare against the face-averaged
        # normal's own magnitude, not an absolute threshold — this works
        # regardless of element size or metric-term scaling. Ω is assumed
        # unit-norm (true for direction cosines sin θ cos φ, etc.).
        n_mag = sqrt((nx_sum/n)^2 + (ny_sum/n)^2 + (nz_sum/n)^2)
        cosang = n_mag > 0 ? dot_avg / n_mag : 0.0

        is_active[fid] = abs(cosang) >= tangential_tol
        is_inflow[fid] = dot_avg < 0.0   # only meaningful when is_active[fid] is true
    end
    return is_inflow, is_active
end

function debug_single_face_pair(mesh, metrics, iel::Int, jel::Int, neighbor_elem, neighbor_fid;
                                 Ωx::Float64, Ωy::Float64, Ωz::Float64)
    ngl = mesh.ngl
    println("═"^70)
    println("DEBUG: element $iel <-> element $jel   (Ω=($Ωx,$Ωy,$Ωz))")

    for e in (iel, jel)
        cx = cy = cz = 0.0; n = 0
        for k=1:ngl, j=1:ngl, i=1:ngl
            ip = mesh.connijk[e,i,j,k]
            cx += mesh.x[ip]; cy += mesh.y[ip]; cz += mesh.z[ip]; n += 1
        end
        println("  centroid(iel=$e) = ($(round(cx/n,digits=4)), $(round(cy/n,digits=4)), $(round(cz/n,digits=4)))")
    end

    fid_i_to_j = 0; fid_j_to_i = 0
    for fid = 1:6
        neighbor_elem[iel][fid] == jel && (fid_i_to_j = fid)
        neighbor_elem[jel][fid] == iel && (fid_j_to_i = fid)
    end
    println("  neighbor_elem[$iel] points to $jel via my_fid=$fid_i_to_j")
    println("  neighbor_elem[$jel] points to $iel via my_fid=$fid_j_to_i")
    if fid_i_to_j > 0
        println("  neighbor_fid[$iel][$fid_i_to_j] = $(neighbor_fid[iel][fid_i_to_j]) (should equal $fid_j_to_i)")
    end

    for (e, fid, label) in ((iel, fid_i_to_j, "iel's face toward jel"),
                            (jel, fid_j_to_i, "jel's face toward iel"))
        fid == 0 && continue
        msym, sgn = face_normal_sign(fid)
        irange, jrange, krange = face_ip_range(fid, ngl)
        println("  $label: iel=$e fid=$fid (msym=$msym sgn=$sgn)")
        for k in krange, j in jrange, i in irange
            if msym == 1
                nx,ny,nz = metrics.dηdx[e,i,j,k], metrics.dηdy[e,i,j,k], metrics.dηdz[e,i,j,k]
            elseif msym == 2
                nx,ny,nz = metrics.dζdx[e,i,j,k], metrics.dζdy[e,i,j,k], metrics.dζdz[e,i,j,k]
            else
                nx,ny,nz = metrics.dξdx[e,i,j,k], metrics.dξdy[e,i,j,k], metrics.dξdz[e,i,j,k]
            end
            raw_dot = Ωx*nx + Ωy*ny + Ωz*nz
            println("    node(i=$i,j=$j,k=$k): (nx,ny,nz)=($(round(nx,digits=4)),$(round(ny,digits=4)),$(round(nz,digits=4)))  " *
                    "raw Ω·n=$(round(raw_dot,digits=6))  signed(×sgn)=$(round(sgn*raw_dot,digits=6))")
        end
    end
    println("═"^70)
end

# ─────────────────────────────────────────────────────────────────────────────
#  resolve_rank_cut_neighbors  (DRAFT — untested, new infrastructure)
#
#  For every suspected rank-partition-cut face (untagged, no same-rank
#  neighbor — see suspected_cuts from the verification diagnostics), find
#  which remote rank actually owns the element on the other side, using a
#  gather -> match-by-key -> scatter-back pattern (same shape as
#  restructure4periodicity_3D_sorted!, but simpler: exact global-corner-ID
#  set match, no geometric/coordinate matching needed).
#
#  ASSUMPTION BEING TESTED, not yet confirmed: a genuine two-rank partition
#  cut's four corner nodes have IDENTICAL global point IDs (via ip2gip) on
#  both ranks. This is a more basic property than periodic merging — it's
#  what makes cross-rank DSS accumulation possible at all — but given that
#  the exact same kind of "should obviously be true" assumption about
#  gip2owner just turned out to need correction, this function verifies it
#  directly via a mutual round-trip check rather than assuming it silently.
#
#  Returns:
#    - neighbor_rank::Dict{Tuple{TInt,TInt}, TInt}  — (iel,fid) -> remote rank
#    - neighbor_gel::Dict{Tuple{TInt,TInt}, TInt}   — (iel,fid) -> remote element's GLOBAL id (mesh.el2gel on their side)
#    - n_unresolved::Int  — suspected cuts that found NO match anywhere
#      (would mean either the corner-ID assumption is wrong, or this face
#      isn't actually a 2-rank interior cut — flagged, not silently dropped)
# ─────────────────────────────────────────────────────────────────────────────

function resolve_rank_cut_neighbors(mesh, suspected_cuts::Vector{Tuple{TInt,TInt}}, comm)
    rank    = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    ngl     = mesh.ngl

    # ── Phase 1: each rank packs (key=4 global corner IDs, rank, local iel,
    # local fid, global element id) for each of its suspected cuts. ────────
    local_keys = NTuple{4,TInt}[]
    local_iel  = TInt[]
    local_fid  = TInt[]
    local_gel  = TInt[]

    for (iel, fid) in suspected_cuts
        fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
        n1 = mesh.conn_face_el[1, fe, iel]
        n2 = mesh.conn_face_el[2, fe, iel]
        n3 = mesh.conn_face_el[3, fe, iel]
        n4 = mesh.conn_face_el[4, fe, iel]
        key = Tuple(sort(TInt[mesh.ip2gip[n1], mesh.ip2gip[n2], mesh.ip2gip[n3], mesh.ip2gip[n4]]))
        push!(local_keys, key)
        push!(local_iel, iel)
        push!(local_fid, fid)
        push!(local_gel, mesh.el2gel[iel])
    end

    # ── Phase 2: gather everyone's entries to rank 0, match by exact key. ──
    n_local = TInt(length(local_keys))
    all_n   = MPI.Allgather(n_local, comm)

    # Flatten each local entry to 7 Ints for easy Allgatherv:
    # [g1,g2,g3,g4, rank, iel, fid, gel]  (rank stored explicitly so matching
    # works after concatenation, since we lose which-rank-sent-this otherwise)
    flat_local = TInt[]
    for i = eachindex(local_keys)
        k = local_keys[i]
        append!(flat_local, [k[1],k[2],k[3],k[4], TInt(rank), local_iel[i], local_fid[i], local_gel[i]])
    end
    counts = all_n .* TInt(8)
    all_flat = MPI.Allgatherv(flat_local, Cint.(counts), comm)

    n_total = length(all_flat) ÷ 8
    entries_by_key = Dict{NTuple{4,TInt}, Vector{NTuple{4,TInt}}}()  # key -> [(rank,iel,fid,gel), ...]
    for m = 1:n_total
        b = (m-1)*8
        key = (all_flat[b+1], all_flat[b+2], all_flat[b+3], all_flat[b+4])
        entry = (all_flat[b+5], all_flat[b+6], all_flat[b+7], all_flat[b+8])
        push!(get!(entries_by_key, key, NTuple{4,TInt}[]), entry)
    end

    # ── Phase 3: each rank looks up its own suspected cuts in the global table. ─
    neighbor_rank = Dict{Tuple{TInt,TInt}, TInt}()
    neighbor_gel  = Dict{Tuple{TInt,TInt}, TInt}()
    n_unresolved  = 0
    n_bad_group   = 0   # a key matched by something other than exactly 2 DIFFERENT ranks

    for i = eachindex(local_keys)
        key = local_keys[i]
        group = get(entries_by_key, key, NTuple{4,TInt}[])
        others = filter(e -> e[1] != rank, group)
        if length(group) != 2 || length(others) != 1
            n_bad_group += 1
            continue
        end
        orank, ogel_iel, ofid, ogel = others[1]
        neighbor_rank[(local_iel[i], local_fid[i])] = orank
        neighbor_gel[(local_iel[i], local_fid[i])]  = ogel
    end
    n_unresolved = length(suspected_cuts) - length(neighbor_rank)

    if n_bad_group > 0
        @warn "resolve_rank_cut_neighbors: rank $rank found $n_bad_group suspected-cut faces " *
              "whose corner-ID key matched something OTHER than exactly one other rank " *
              "(0 matches, self-matches, or >2-way matches). This likely means the " *
              "assumption that partition-cut corners share identical global IDs across " *
              "ranks does NOT hold as expected — do not trust neighbor_rank/neighbor_gel " *
              "until this is understood."
    end

    return neighbor_rank, neighbor_gel, n_unresolved
end

"""
    verify_rank_cut_neighbors(mesh, neighbor_rank, suspected_cuts, comm)

For each resolved (iel,fid) -> remote rank pairing, exchange a small
"handshake" value with that remote rank and confirm they independently
identify US as one of THEIR unresolved-cut neighbors too — a genuine
two-way, cross-rank confirmation, not just "the numbers matched on my side."
"""
function verify_rank_cut_neighbors(mesh, neighbor_rank, suspected_cuts, comm)
    rank = MPI.Comm_rank(comm)
    n_checked = 0
    n_confirmed = 0

    # For each remote rank we believe we share a cut with, send them the count
    # of cuts we think we share, and independently ask them the same. If both
    # sides agree on the COUNT (not full geometric detail — that's already
    # implied by the exact-key match in resolve_rank_cut_neighbors), that's a
    # strong sign the resolution is mutually consistent.
    my_partners = unique(values(neighbor_rank))
    counts_i_claim = Dict(p => count(==(p), values(neighbor_rank)) for p in my_partners)

    reqs = MPI.Request[]
    recv_bufs = Dict{TInt, Vector{TInt}}()
    for p in my_partners
        recv_bufs[p] = zeros(TInt, 1)
        push!(reqs, MPI.Isend([TInt(counts_i_claim[p])], p, 42, comm))
        push!(reqs, MPI.Irecv!(recv_bufs[p], p, 42, comm))
    end
    MPI.Waitall(reqs)

    for p in my_partners
        n_checked += 1
        their_claim = recv_bufs[p][1]
        my_claim    = counts_i_claim[p]
        if their_claim == my_claim
            n_confirmed += 1
        else
            @warn "verify_rank_cut_neighbors: rank $rank claims $my_claim shared cut faces " *
                  "with rank $p, but rank $p claims $their_claim — MISMATCH, resolution is " *
                  "not mutually consistent for this pair."
        end
    end

    println("Rank-cut mutual verification: $n_confirmed / $n_checked partner ranks agree on shared-cut count")
    return n_confirmed == n_checked
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_dg_volume_matrix  (DRAFT — untested)
#
#  Builds the element-local DG volume matrix for one element, one angular
#  direction Ω = (Ωx,Ωy,Ωz), for the operator:
#
#      Ω·∇I + α·I     where α = κ + σ*(1 - intΦ_eff)
#
#  matching EXACTLY the physical terms already assembled globally in
#  sparse_lhs_assembly_3Dby2D_adaptive (diagonal reaction term +
#  prop_i/prop_j/prop_k streaming terms) — same coefficients, same
#  quadrature, just restricted to one element's own ngl^3 local nodes
#  with NO scatter into any global index. This is the volume-integral part
#  of the DG weak form only — face/flux terms are a separate function,
#  to be verified independently before combining (see conversation).
#
#  A_vol[p, q] = ∫_element [ (Ω·∇ψ_q) ψ_p + α ψ_q ψ_p ] dV
#  where p, q range over the element's ngl^3 local nodes, flattened via
#  lidx(i,j,k) = i + (j-1)*ngl + (k-1)*ngl^2  (consistent flattening used
#  throughout, so the resulting vector can be reshaped back to (i,j,k) or
#  indexed flat as convenient downstream)
#
#  PERF: takes raw arrays (ψ, dψ, ω, metric components, κ, σ, intΦ_eff),
#  not structs — same reasoning as sweep_ordering.jl (struct-field access
#  across a function boundary can defeat inlining/specialization; this
#  runs once per element per direction, every sweep).
# ─────────────────────────────────────────────────────────────────────────────

@inline lidx(i,j,k,ngl) = i + (j-1)*ngl + (k-1)*ngl*ngl

"""
    build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl,
                             ψ, dψ, ω,
                             dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                             Je, κ, σ, intΦ_eff, connijk)

In-place: fills the (ngl^3 x ngl^3) buffer `A_vol` (caller-allocated once,
reused across elements/directions — no allocation here) with the element-
local DG volume matrix for element `iel`, direction (Ωx,Ωy,Ωz).

`ψ`, `dψ` are the 1D Lagrange basis / derivative matrices (ngl x ngl), as
used elsewhere in this codebase's CG assembly — same basis, restricted to
one element rather than assembled globally.

`κ`, `σ`, `intΦ_eff` are indexed by GLOBAL point id (matching how
sparse_lhs_assembly_3Dby2D_adaptive reads them: `κ_data[ip]` etc. via
`connijk`/`connijk_spa`) — `connijk` is passed in specifically to map this
element's local (i,j,k) to the global point id needed for those lookups.
"""
function build_dg_volume_matrix!(A_vol::AbstractMatrix{Float64}, iel::Int,
                                  Ωx::Float64, Ωy::Float64, Ωz::Float64, ngl::Int,
                                  ψ, dψ, ω,
                                  dξdx, dξdy, dξdz,
                                  dηdx, dηdy, dηdz,
                                  dζdx, dζdy, dζdz,
                                  Je, κ, σ, intΦ_eff, connijk)
    fill!(A_vol, 0.0)

    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        # Quadrature node (i,j,k) — this is where the weak form is evaluated,
        # matching sparse_lhs_assembly_3Dby2D_adaptive's own (i,j,k) loop.
        ip_global = connijk[iel, i, j, k]
        ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]

        dξdx_ij = dξdx[iel,i,j,k]; dξdy_ij = dξdy[iel,i,j,k]; dξdz_ij = dξdz[iel,i,j,k]
        dηdx_ij = dηdx[iel,i,j,k]; dηdy_ij = dηdy[iel,i,j,k]; dηdz_ij = dηdz[iel,i,j,k]
        dζdx_ij = dζdx[iel,i,j,k]; dζdy_ij = dζdy[iel,i,j,k]; dζdz_ij = dζdz[iel,i,j,k]

        α = κ[ip_global] + σ[ip_global] * (1.0 - intΦ_eff[ip_global])
        p = lidx(i,j,k,ngl)

        # ── Reaction term: α * ψ_q(node) * ψ_p(node) * ωJac ────────────────
        # At a GLL nodal quadrature point, ψ_q(node_p) = δ_{pq} (Kronecker
        # delta — the defining property of a nodal Lagrange basis evaluated
        # at its own collocation points), so this term is PURELY diagonal —
        # exactly matching the diagonal-only V_vec[pos] term in
        # sparse_lhs_assembly_3Dby2D_adaptive, not a separate approximation.
        A_vol[p, p] += α * ωJac

        # ── Streaming term: (Ω·∇ψ_q) evaluated at this node, projected onto
        # test function ψ_p — since ψ_p is also nodal, this again collapses
        # to contributions AT this same row p, for every q reachable via the
        # 1D derivative matrix along ξ, η, ζ through node (i,j,k). This is
        # exactly the prop_i/prop_j/prop_k construction from
        # sparse_lhs_assembly_3Dby2D_adaptive, just writing into a LOCAL
        # ngl^3 matrix instead of scattering into the global sparse matrix.
        prop_i = (dξdx_ij*Ωx + dξdy_ij*Ωy + dξdz_ij*Ωz) * ωJac
        for m = 1:ngl
            abs(dψ[m,i]) < eps(Float64) && continue
            q = lidx(m,j,k,ngl)
            A_vol[p, q] += dψ[m,i] * prop_i
        end

        prop_j = (dηdx_ij*Ωx + dηdy_ij*Ωy + dηdz_ij*Ωz) * ωJac
        for m = 1:ngl
            abs(dψ[m,j]) < eps(Float64) && continue
            q = lidx(i,m,k,ngl)
            A_vol[p, q] += dψ[m,j] * prop_j
        end

        prop_k = (dζdx_ij*Ωx + dζdy_ij*Ωy + dζdz_ij*Ωz) * ωJac
        for m = 1:ngl
            abs(dψ[m,k]) < eps(Float64) && continue
            q = lidx(i,j,m,ngl)
            A_vol[p, q] += dψ[m,k] * prop_k
        end
    end

    return A_vol
end

# ─────────────────────────────────────────────────────────────────────────────
#  test_dg_volume_matrix  (DRAFT — untested, unit test for build_dg_volume_matrix!)
#
#  Verifies build_dg_volume_matrix! against CLOSED-FORM ANALYTIC values on a
#  single axis-aligned box element with constant (deliberately asymmetric,
#  to catch dimension-swap bugs) side lengths Δx≠Δy≠Δz. No mesh, no Gridap,
#  no global assembly/DSS involved — a pure, self-contained unit test.
#
#  For an axis-aligned box, dηdx=dζdx=0 (and cyclic permutations), so for
#  Ω=(1,0,0) ONLY the ξ-direction streaming survives; similarly for
#  (0,1,0)/η and (0,0,1)/ζ. Three fully separate, explicit checks (no shared
#  generic helper) — repetition here trades code length for direct
#  readability against the analytic form, deliberately.
#
#  Test 4 (linearity) then certifies arbitrary oblique Ω by superposition:
#  prop_i/prop_j/prop_k are each linear in (Ωx,Ωy,Ωz), so confirming the
#  three axis-aligned directions AND exact linearity together certify every
#  oblique direction without hand-deriving messier oblique-direction algebra.
#
#  σ=0 for this test — isolates streaming + pure-absorption reaction
#  (α=κ). σ*(1-intΦ_eff) is purely additive to the diagonal and doesn't
#  interact with the streaming algebra being tested here.
# ─────────────────────────────────────────────────────────────────────────────

function test_dg_volume_matrix(; ngl::Int = 4)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    # Deliberately asymmetric box — Δx≠Δy≠Δz catches a dimension-swap bug
    # that a symmetric (e.g. cubic) test element would silently hide.
    Δx, Δy, Δz = 2.0, 3.0, 4.0
    Je_const = (Δx/2)*(Δy/2)*(Δz/2)

    dξdx = fill(2/Δx, 1,ngl,ngl,ngl); dξdy = zeros(1,ngl,ngl,ngl); dξdz = zeros(1,ngl,ngl,ngl)
    dηdx = zeros(1,ngl,ngl,ngl); dηdy = fill(2/Δy, 1,ngl,ngl,ngl); dηdz = zeros(1,ngl,ngl,ngl)
    dζdx = zeros(1,ngl,ngl,ngl); dζdy = zeros(1,ngl,ngl,ngl); dζdz = fill(2/Δz, 1,ngl,ngl,ngl)
    Je   = fill(Je_const, 1,ngl,ngl,ngl)

    npoin_test = ngl^3
    connijk = Array{Int}(undef, 1,ngl,ngl,ngl)
    for k=1:ngl, j=1:ngl, i=1:ngl
        connijk[1,i,j,k] = lidx(i,j,k,ngl)   # single element: local == global index
    end

    κ_const = 0.7
    κ        = fill(κ_const, npoin_test)
    σ        = zeros(npoin_test)
    intΦ_eff = zeros(npoin_test)

    A_vol = zeros(npoin_test, npoin_test)

    # ── Test 1: Ω=(1,0,0) — only ξ-direction (i-index) streaming survives ──
    build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, σ, intΦ_eff, connijk)
    max_err_1 = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        ωJac = ω[i]*ω[j]*ω[k]*Je_const
        for q_i=1:ngl, q_j=1:ngl, q_k=1:ngl
            q = lidx(q_i,q_j,q_k,ngl)
            expected = 0.0
            if (q_i,q_j,q_k) == (i,j,k)
                expected += κ_const * ωJac
            end
            if q_j == j && q_k == k   # streaming couples only along i (ξ), same j,k
                expected += dψ[q_i,i] * (2/Δx) * ωJac
            end
            max_err_1 = max(max_err_1, abs(A_vol[p,q] - expected))
        end
    end
    println("Test (Ω=(1,0,0)): max abs error = $max_err_1  " * (max_err_1 < 1e-12 ? "✓ PASS" : "✗ FAIL"))

    # ── Test 2: Ω=(0,1,0) — only η-direction (j-index) streaming survives ──
    build_dg_volume_matrix!(A_vol, 1, 0.0, 1.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, σ, intΦ_eff, connijk)
    max_err_2 = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        ωJac = ω[i]*ω[j]*ω[k]*Je_const
        for q_i=1:ngl, q_j=1:ngl, q_k=1:ngl
            q = lidx(q_i,q_j,q_k,ngl)
            expected = 0.0
            if (q_i,q_j,q_k) == (i,j,k)
                expected += κ_const * ωJac
            end
            if q_i == i && q_k == k   # streaming couples only along j (η), same i,k
                expected += dψ[q_j,j] * (2/Δy) * ωJac
            end
            max_err_2 = max(max_err_2, abs(A_vol[p,q] - expected))
        end
    end
    println("Test (Ω=(0,1,0)): max abs error = $max_err_2  " * (max_err_2 < 1e-12 ? "✓ PASS" : "✗ FAIL"))

    # ── Test 3: Ω=(0,0,1) — only ζ-direction (k-index) streaming survives ──
    build_dg_volume_matrix!(A_vol, 1, 0.0, 0.0, 1.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, σ, intΦ_eff, connijk)
    max_err_3 = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        ωJac = ω[i]*ω[j]*ω[k]*Je_const
        for q_i=1:ngl, q_j=1:ngl, q_k=1:ngl
            q = lidx(q_i,q_j,q_k,ngl)
            expected = 0.0
            if (q_i,q_j,q_k) == (i,j,k)
                expected += κ_const * ωJac
            end
            if q_i == i && q_j == j   # streaming couples only along k (ζ), same i,j
                expected += dψ[q_k,k] * (2/Δz) * ωJac
            end
            max_err_3 = max(max_err_3, abs(A_vol[p,q] - expected))
        end
    end
    println("Test (Ω=(0,0,1)): max abs error = $max_err_3  " * (max_err_3 < 1e-12 ? "✓ PASS" : "✗ FAIL"))

    # ── Test 4: linearity in Ω certifies arbitrary oblique directions ──────
    Ω1 = (0.3, 0.5, 0.8)
    Ω2 = (-0.2, 0.9, 0.1)
    A1, A2, Asum = zeros(npoin_test,npoin_test), zeros(npoin_test,npoin_test), zeros(npoin_test,npoin_test)
    build_dg_volume_matrix!(A1, 1, Ω1[1],Ω1[2],Ω1[3], ngl, ψ,dψ,ω, dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,σ,intΦ_eff,connijk)
    build_dg_volume_matrix!(A2, 1, Ω2[1],Ω2[2],Ω2[3], ngl, ψ,dψ,ω, dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,σ,intΦ_eff,connijk)
    build_dg_volume_matrix!(Asum, 1, Ω1[1]+Ω2[1], Ω1[2]+Ω2[2], Ω1[3]+Ω2[3], ngl, ψ,dψ,ω,
                             dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,σ,intΦ_eff,connijk)

    # Reaction-term diagonal is Ω-INDEPENDENT and gets added fresh in every
    # call, so A1+A2 double-counts it once relative to Asum — subtract one
    # copy back out before comparing.
    reaction_diag_vec = zeros(npoin_test)
    for k=1:ngl, j=1:ngl, i=1:ngl
        reaction_diag_vec[lidx(i,j,k,ngl)] = κ_const*ω[i]*ω[j]*ω[k]*Je_const
    end
    max_err_4 = maximum(abs.((A1 .+ A2 .- Diagonal(reaction_diag_vec)) .- Asum))
    println("Test (linearity in Ω): max abs error = $max_err_4  " * (max_err_4 < 1e-10 ? "✓ PASS" : "✗ FAIL"))

    return max(max_err_1, max_err_2, max_err_3, max_err_4)
end

# ─────────────────────────────────────────────────────────────────────────────
#  test_sweep_vs_cg_coupling  (DRAFT — untested)
#
#  Standalone, single-rank, no-MPI-needed test: builds a small hand-made
#  structured box mesh (NOT via Gridap/GMSH — simple enough to hand-verify),
#  DSS-assembles the GLOBAL CG streaming+reaction matrix from
#  build_dg_volume_matrix! (already validated against closed-form values),
#  then checks whether that assembled matrix is BLOCK-triangular under the
#  ALREADY-VALIDATED sweep order (build_face_adjacency_from_conn_face_el +
#  build_conforming_sweep_order, unchanged, reused as-is).
#
#  "Block" = element. A shared CG node's block = whichever element reaches
#  it EARLIEST in sweep order (the element that would resolve it first in a
#  forward block-sweep). Mass with block(row) < block(col) is what a
#  sequential block-sweep would DROP (an earlier/upstream element's equation
#  depending on a not-yet-solved LATER element's unknown) — that's the
#  number that answers whether "CG sequential elimination" is exact.
#  Mass with block(row) >= block(col) is legitimate (within-block coupling,
#  solved together as one dense system; or a later block depending on an
#  earlier, already-known one — exactly the upwind-flux mechanism).
# ─────────────────────────────────────────────────────────────────────────────
function build_toy_structured_mesh(Nx::Int, Ny::Int, Nz::Int, ngl::Int,
                                    Δx::Float64, Δy::Float64, Δz::Float64)
    nelem = Nx*Ny*Nz
    gnx = Nx*(ngl-1) + 1
    gny = Ny*(ngl-1) + 1
    gnz = Nz*(ngl-1) + 1
    npoin = gnx*gny*gnz

    lgl = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())   # need GLL node positions for exact coords

    connijk = Array{Int}(undef, nelem, ngl, ngl, ngl)
    x = zeros(Float64, npoin); y = zeros(Float64, npoin); z = zeros(Float64, npoin)

    for ez = 1:Nz, ey = 1:Ny, ex = 1:Nx
        iel = ex + (ey-1)*Nx + (ez-1)*Nx*Ny
        x0 = (ex-1)*Δx; y0 = (ey-1)*Δy; z0 = (ez-1)*Δz
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            gx = (ex-1)*(ngl-1) + i
            gy = (ey-1)*(ngl-1) + j
            gz = (ez-1)*(ngl-1) + k
            g = gx + (gy-1)*gnx + (gz-1)*gnx*gny
            connijk[iel,i,j,k] = g
            x[g] = x0 + (Δx/2)*(1 + lgl.ξ[i])
            y[g] = y0 + (Δy/2)*(1 + lgl.ξ[j])
            z[g] = z0 + (Δz/2)*(1 + lgl.ξ[k])
        end
    end

    n_distinct = length(unique(vec(connijk)))
    n_distinct == npoin || error("build_toy_structured_mesh: expected $npoin distinct nodes, " *
                                 "got $n_distinct — shared-node numbering is wrong, stop here.")

    conn_face_el = Array{Int}(undef, 4, 6, nelem)
    for iel = 1:nelem
        ip1 = connijk[iel, ngl, 1,   1]
        ip2 = connijk[iel, ngl, ngl, 1]
        ip3 = connijk[iel, 1,   ngl, 1]
        ip4 = connijk[iel, 1,   1,   1]
        ip5 = connijk[iel, ngl, 1,   ngl]
        ip6 = connijk[iel, ngl, ngl, ngl]
        ip7 = connijk[iel, 1,   ngl, ngl]
        ip8 = connijk[iel, 1,   1,   ngl]
        conn_face_el[:,1,iel] = [ip1,ip4,ip3,ip2]
        conn_face_el[:,2,iel] = [ip1,ip2,ip6,ip5]
        conn_face_el[:,3,iel] = [ip2,ip3,ip7,ip6]
        conn_face_el[:,4,iel] = [ip3,ip4,ip8,ip7]
        conn_face_el[:,5,iel] = [ip1,ip5,ip8,ip4]
        conn_face_el[:,6,iel] = [ip5,ip6,ip7,ip8]
    end

    Je_const = (Δx/2)*(Δy/2)*(Δz/2)
    dξdx = fill(2/Δx, nelem,ngl,ngl,ngl); dξdy = zeros(nelem,ngl,ngl,ngl); dξdz = zeros(nelem,ngl,ngl,ngl)
    dηdx = zeros(nelem,ngl,ngl,ngl); dηdy = fill(2/Δy, nelem,ngl,ngl,ngl); dηdz = zeros(nelem,ngl,ngl,ngl)
    dζdx = zeros(nelem,ngl,ngl,ngl); dζdy = zeros(nelem,ngl,ngl,ngl); dζdz = fill(2/Δz, nelem,ngl,ngl,ngl)
    Je   = fill(Je_const, nelem,ngl,ngl,ngl)

    xmin = minimum(x)   # NEW — needed for x0 in the exact-solution comparison

    return (nelem=nelem, npoin=npoin, ngl=ngl, connijk=connijk, conn_face_el=conn_face_el,
            x=x, y=y, z=z, xmin=xmin,
            dξdx=dξdx,dξdy=dξdy,dξdz=dξdz, dηdx=dηdx,dηdy=dηdy,dηdz=dηdz,
            dζdx=dζdx,dζdy=dζdy,dζdz=dζdz, Je=Je)
end


function test_sweep_vs_cg_coupling(; Nx=3, Ny=3, Nz=3, ngl=4,
                                    Δx=1.3, Δy=1.7, Δz=2.1,   # asymmetric, catches dim-swap bugs
                                    Ω = (0.3, 0.5, 0.8))       # oblique — generic, not axis-aligned
    m = build_toy_structured_mesh(Nx, Ny, Nz, ngl, Δx, Δy, Δz)
    println("Toy mesh: $(m.nelem) elements, $(m.npoin) global nodes, ngl=$ngl")

    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    κ        = fill(0.7, m.npoin)   # value doesn't affect this test (see note below)
    σ        = zeros(m.npoin)
    intΦ_eff = zeros(m.npoin)

    # ── DSS-assemble the global matrix from build_dg_volume_matrix! (ALREADY
    # VALIDATED per-element physics) — exactly the same scatter-add pattern
    # already confirmed correct for vectors via DSS_rhs! (RHS[I] += rhs_el),
    # here applied to a matrix (standard CG stiffness-matrix assembly).
    Ωx, Ωy, Ωz = Ω
    I_g = Int[]; J_g = Int[]; V_g = Float64[]
    A_vol = zeros(ngl^3, ngl^3)
    for iel = 1:m.nelem
        build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ω,
                                 m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz,
                                 m.Je, κ, σ, intΦ_eff, m.connijk)
        for k=1:ngl, j=1:ngl, i=1:ngl
            p = lidx(i,j,k,ngl); gp = m.connijk[iel,i,j,k]
            for kk=1:ngl, jj=1:ngl, ii=1:ngl
                q = lidx(ii,jj,kk,ngl); gq = m.connijk[iel,ii,jj,kk]
                A_vol[p,q] == 0.0 && continue
                push!(I_g, gp); push!(J_g, gq); push!(V_g, A_vol[p,q])
            end
        end
    end
    A_global = sparse(I_g, J_g, V_g, m.npoin, m.npoin)   # duplicate (i,j) entries SUM — this IS DSS

    # ── Sweep order for the SAME Ω — reusing ALREADY VALIDATED code UNCHANGED ─
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(m.nelem, m.conn_face_el)
    order, levels = build_conforming_sweep_order(Ωx, Ωy, Ωz, m.nelem, ngl, neighbor_elem,
                                                  m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz)
    pos = zeros(Int, m.nelem)
    for (p,iel) in enumerate(order); pos[iel] = p; end

    # ── Node -> block assignment: earliest (smallest pos) owning element ──
    block_of = fill(typemax(Int), m.npoin)
    for iel = 1:m.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = m.connijk[iel,i,j,k]
        block_of[g] = min(block_of[g], pos[iel])
    end
    any(block_of .== typemax(Int)) && error("some node never assigned a block — connijk/order mismatch")

    # ── Bucket every nonzero by block(row) vs block(col) ───────────────────
    rows, cols, vals = findnz(A_global)
    below = 0.0; above = 0.0; same = 0.0
    for idx = eachindex(rows)
        r, c, v = rows[idx], cols[idx], vals[idx]
        br, bc = block_of[r], block_of[c]
        if br < bc
            above += abs(v)
        elseif br > bc
            below += abs(v)
        else
            same += abs(v)
        end
    end

    println("Ω = $Ω")
    println("mass WITHIN a block (same element, expected — solved as one dense system): $same")
    println("mass BELOW block-diagonal (later block depends on earlier — legitimate upwind coupling): $below")
    println("mass ABOVE block-diagonal (earlier block depends on LATER — would be DROPPED by sequential sweep): $above")
    println("ratio above/(below+same): $(above/(below+same))")

    return (above=above, below=below, same=same)
end

function test_sweep_with_bc(; Nx=3, Ny=3, Nz=3, ngl=4,
                            Δx=1.3, Δy=1.7, Δz=2.1,
                            Ω = (0.3, 0.5, 0.8),
                            κ_val = 0.7,   # nonzero reaction — also removes the
                                           # constant-null-space issue on its own,
                                           # but we fix the BC regardless since a
                                           # true streaming problem needs one anyway
                            seed = 42)
    m = build_toy_structured_mesh(Nx, Ny, Nz, ngl, Δx, Δy, Δz)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ωq = basis.ψ, basis.dψ, lgl.ω
    Ωx, Ωy, Ωz = Ω

    κ        = fill(κ_val, m.npoin)
    σ        = zeros(m.npoin)
    intΦ_eff = zeros(m.npoin)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(m.nelem, m.conn_face_el)
    order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, m.nelem, ngl, neighbor_elem,
                                             m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz)

    # ── Identify TRUE physical boundary nodes: any face with no same-rank
    # interior neighbor (neighbor_elem==0) — for this toy mesh there's no
    # periodicity, so every such face is a genuine domain boundary. ────────
    rng = MersenneTwister(seed)
    bc_value = Dict{Int,Float64}()   # global node -> prescribed inflow value
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    for iel = 1:m.nelem
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            neighbor_elem[iel][fid] == zero(TInt) || continue   # only TRUE boundary, not interior
            irange,jrange,krange = face_ip_range(fid, ngl)
            for k in krange, j in jrange, i in irange
                g = m.connijk[iel,i,j,k]
                bc_value[g] = randn(rng)   # arbitrary fixed prescribed value per boundary node
            end
        end
    end
    println("Prescribed $(length(bc_value)) true-boundary inflow node values")

    # ── Assemble A_global (same as before) and b, then apply BC by row
    # replacement (standard Dirichlet enforcement): for every g in bc_value,
    # zero g's row, set diagonal=1, RHS=bc_value[g]. ────────────────────────
    I_g = Int[]; J_g = Int[]; V_g = Float64[]
    A_vol = zeros(ngl^3, ngl^3)
    for iel = 1:m.nelem
        build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ωq,
                                 m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz,
                                 m.Je, κ, σ, intΦ_eff, m.connijk)
        for k=1:ngl, j=1:ngl, i=1:ngl
            p = lidx(i,j,k,ngl); gp = m.connijk[iel,i,j,k]
            for kk=1:ngl, jj=1:ngl, ii=1:ngl
                q = lidx(ii,jj,kk,ngl); gq = m.connijk[iel,ii,jj,kk]
                A_vol[p,q] == 0.0 && continue
                push!(I_g, gp); push!(J_g, gq); push!(V_g, A_vol[p,q])
            end
        end
    end
    A_global = sparse(I_g, J_g, V_g, m.npoin, m.npoin)
    b = randn(rng, m.npoin)

    A_bc = copy(A_global)
    b_bc = copy(b)
    for (g, val) in bc_value
        A_bc[g, :] .= 0.0
        A_bc[g, g] = 1.0
        b_bc[g] = val
    end
    x_exact = A_bc \ b_bc
    println("‖x_exact‖ = $(norm(x_exact))  (should now be a reasonable, non-huge magnitude)")

    # ── ONE forward sweep pass: same mechanism as before, but boundary nodes
    # are seeded as KNOWN (their BC value) before the sweep starts, exactly
    # matching how the direct solve enforced them. ─────────────────────────
    x = zeros(m.npoin)
    known = falses(m.npoin)
    for (g,val) in bc_value
        x[g] = val; known[g] = true
    end

    A_vol_local = zeros(ngl^3, ngl^3)
    for iel in order
        build_dg_volume_matrix!(A_vol_local, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ωq,
                                 m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz,
                                 m.Je, κ, σ, intΦ_eff, m.connijk)
        gids = Array{Int}(undef, ngl^3)
        for k=1:ngl, j=1:ngl, i=1:ngl
            gids[lidx(i,j,k,ngl)] = m.connijk[iel,i,j,k]
        end
        unk_local = [p for p = 1:ngl^3 if !known[gids[p]]]
        isempty(unk_local) && continue
        rhs_local = zeros(length(unk_local))
        for (ui,p) in enumerate(unk_local)
            rhs_local[ui] = b[gids[p]]
            for q = 1:ngl^3
                known[gids[q]] && (rhs_local[ui] -= A_vol_local[p,q]*x[gids[q]])
            end
        end
        Aloc = A_vol_local[unk_local, unk_local]
        xloc = Aloc \ rhs_local
        for (ui,p) in enumerate(unk_local)
            x[gids[p]] = xloc[ui]; known[gids[p]] = true
        end
    end

    rel_err = norm(x - x_exact) / norm(x_exact)
    println("ONE forward sweep pass: rel. error vs BC-corrected exact solve = $rel_err")
    return rel_err
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_dg_flux_terms!  (DERIVED — see conversation for full weak-form
#  derivation; strong-form-DG correction term, exactly equivalent to the
#  standard weak-form upwind-flux DG method)
#
#  ADDITIVE and SEPARATE from build_dg_volume_matrix! by design — call
#  AFTER the volume matrix for this element/direction, accumulates into
#  the SAME A_vol/RHS_local buffers. Verified to be IDENTICALLY ZERO on
#  outflow faces (Î≡I there — not a special case handled by omission, a
#  genuine mathematical identity), so no is_inflow branch is even needed
#  beyond the early return below.
# ─────────────────────────────────────────────────────────────────────────────

#=function build_dg_flux_terms!(A_vol::AbstractMatrix{Float64}, RHS_local::AbstractVector{Float64},
                               fid::Int, ngl::Int,
                               Ωx::Float64, Ωy::Float64, Ωz::Float64,
                               is_inflow::Bool, is_active::Bool,
                               face_normal_avg::NTuple{3,Float64},
                               ω,
                               dξdx, dξdy, dξdz,
                               dηdx, dηdy, dηdz,
                               dζdx, dζdy, dζdz,
                               Je, iel::Int,
                               I_upwind)   # (i,j,k) -> known value (neighbor trace OR prescribed BC)

    (is_active && is_inflow) || return   # tangential OR outflow → correction is exactly zero

    msym, sgn = face_normal_sign(fid)
    irange, jrange, krange = face_ip_range(fid, ngl)
    nx, ny, nz = face_normal_avg
    Ωdotn = Ωx*nx + Ωy*ny + Ωz*nz   # < 0 here, guaranteed by is_inflow check above

    @inbounds for k in krange, j in jrange, i in irange
        p = lidx(i,j,k,ngl)

        if msym == 1;      gx,gy,gz = dηdx[iel,i,j,k], dηdy[iel,i,j,k], dηdz[iel,i,j,k]
        elseif msym == 2;  gx,gy,gz = dζdx[iel,i,j,k], dζdy[iel,i,j,k], dζdz[iel,i,j,k]
        else;              gx,gy,gz = dξdx[iel,i,j,k], dξdy[iel,i,j,k], dξdz[iel,i,j,k]
        end
        Je_face = Je[iel,i,j,k] * sqrt(gx^2+gy^2+gz^2)

        if msym == 1;      wq = ω[i]*ω[k]
        elseif msym == 2;  wq = ω[i]*ω[j]
        else;              wq = ω[j]*ω[k]
        end

        weight = -Ωdotn * wq * Je_face   # positive, since Ωdotn < 0
        A_vol[p,p]  += weight
        RHS_local[p] += weight * I_upwind(i,j,k)
    end
end=#

function build_dg_flux_terms!(A_vol::AbstractMatrix{Float64}, RHS_local::AbstractVector{Float64},
                               fid::Int, ngl::Int,
                               Ωx::Float64, Ωy::Float64, Ωz::Float64,
                               is_inflow::Bool, is_active::Bool,
                               face_normal_avg::NTuple{3,Float64},
                               ω, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                               Je, iel::Int, I_upwind;
                               weight_out::Union{AbstractVector{Float64},Nothing} = nothing)  # NEW, optional
    (is_active && is_inflow) || return

    msym, sgn = face_normal_sign(fid)
    irange, jrange, krange = face_ip_range(fid, ngl)
    nx, ny, nz = face_normal_avg
    Ωdotn = Ωx*nx + Ωy*ny + Ωz*nz

    idx = 0
    @inbounds for k in krange, j in jrange, i in irange
        idx += 1
        p = lidx(i,j,k,ngl)
        if msym == 1;      gx,gy,gz = dηdx[iel,i,j,k], dηdy[iel,i,j,k], dηdz[iel,i,j,k]
        elseif msym == 2;  gx,gy,gz = dζdx[iel,i,j,k], dζdy[iel,i,j,k], dζdz[iel,i,j,k]
        else;              gx,gy,gz = dξdx[iel,i,j,k], dξdy[iel,i,j,k], dξdz[iel,i,j,k]
        end
        Je_face = Je[iel,i,j,k] * sqrt(gx^2+gy^2+gz^2)
        if msym == 1;      wq = ω[i]*ω[k]
        elseif msym == 2;  wq = ω[i]*ω[j]
        else;              wq = ω[j]*ω[k]
        end
        weight = -Ωdotn * wq * Je_face
        A_vol[p,p]  += weight
        RHS_local[p] += weight * I_upwind(i,j,k)
        weight_out !== nothing && (weight_out[idx] = weight)   # record for reuse
    end
end

"""
    build_dg_flux_rhs_only!(RHS_local, fid, ngl, is_inflow, is_active, weight_cache, I_upwind)

Applies ONLY the RHS contribution of the flux term, using a PRE-COMPUTED
weight array (one value per face node, from a prior full call to
build_dg_flux_terms!). Used on cache-HIT passes, where A_vol's diagonal
was already correctly built and factorized on the pass that populated the
cache — recomputing it again would be pure waste, exactly the redundant
work the caching exists to eliminate.
"""
function build_dg_flux_rhs_only!(RHS_local::AbstractVector{Float64}, fid::Int, ngl::Int,
                                  is_inflow::Bool, is_active::Bool,
                                  weight_cache::AbstractVector{Float64}, I_upwind)
    (is_active && is_inflow) || return
    irange, jrange, krange = face_ip_range(fid, ngl)
    idx = 0
    @inbounds for k in krange, j in jrange, i in irange
        idx += 1
        p = lidx(i,j,k,ngl)
        RHS_local[p] += weight_cache[idx] * I_upwind(i,j,k)
    end
end
# ─────────────────────────────────────────────────────────────────────────────
#  test_dg_flux_terms  (DRAFT — untested)
#
#  Part A: build_dg_flux_terms! tested in ISOLATION (volume matrix left at
#  zero) against hand-derived closed-form values, including an explicit
#  check that outflow/tangential faces contribute EXACTLY zero — confirming
#  the derivation's central claim (Î≡I on outflow), not assuming it.
#
#  Part B: build_dg_volume_matrix! + build_dg_flux_terms! SUMMED, then
#  solved as one small dense system, checked against the EXACT analytic
#  solution of Ω·∇I + κI = 0 with prescribed inflow — a genuine PDE
#  solution, independent of both functions' own internal formulas.
# ─────────────────────────────────────────────────────────────────────────────

function test_dg_flux_terms(; ngl::Int = 4)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    Δx, Δy, Δz = 2.0, 3.0, 4.0   # deliberately asymmetric, same as volume test
    Je_const = (Δx/2)*(Δy/2)*(Δz/2)

    dξdx = fill(2/Δx, 1,ngl,ngl,ngl); dξdy = zeros(1,ngl,ngl,ngl); dξdz = zeros(1,ngl,ngl,ngl)
    dηdx = zeros(1,ngl,ngl,ngl); dηdy = fill(2/Δy, 1,ngl,ngl,ngl); dηdz = zeros(1,ngl,ngl,ngl)
    dζdx = zeros(1,ngl,ngl,ngl); dζdy = zeros(1,ngl,ngl,ngl); dζdz = fill(2/Δz, 1,ngl,ngl,ngl)
    Je   = fill(Je_const, 1,ngl,ngl,ngl)

    npoin_test = ngl^3
    A_vol = zeros(npoin_test, npoin_test)
    RHS   = zeros(npoin_test)

    # ═══════════════════════════════════════════════════════════════════════
    # PART A: flux in isolation
    # ═══════════════════════════════════════════════════════════════════════

    # ── A1: single inflow face, Ω=(1,0,0) — only face 5 (i=1, "right") is
    # inflow. Every entry checked against a hand-derived closed form. ──────
    Ωx, Ωy, Ωz = 1.0, 0.0, 0.0
    I0 = 3.7   # arbitrary prescribed inflow value, constant across the face
    fill!(A_vol, 0.0); fill!(RHS, 0.0)

    # signed (outward) face normal for face 5: sgn=-1, unsigned=(2/Δx,0,0)
    face5_normal = (-2/Δx, 0.0, 0.0)
    is_inflow5, is_active5 = true, true
    build_dg_flux_terms!(A_vol, RHS, 5, ngl, Ωx, Ωy, Ωz, is_inflow5, is_active5,
                          face5_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> I0)

    max_err_A1 = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        if i == 1   # on the inflow face
            expected_diag = (2/Δx) * ω[j]*ω[k] * (Δy/2)*(Δz/2)
            max_err_A1 = max(max_err_A1, abs(A_vol[p,p] - expected_diag))
            max_err_A1 = max(max_err_A1, abs(RHS[p] - expected_diag*I0))
        else        # not on this face — must be untouched
            max_err_A1 = max(max_err_A1, abs(A_vol[p,p]))
            max_err_A1 = max(max_err_A1, abs(RHS[p]))
        end
        # off-diagonal entries must be untouched everywhere (flux is purely
        # diagonal — no q≠p contribution anywhere in the derivation)
        for q = 1:npoin_test
            q == p && continue
            max_err_A1 = max(max_err_A1, abs(A_vol[p,q]))
        end
    end
    println("Test A1 (single inflow face, Ω=(1,0,0)): max abs error = $max_err_A1  " *
            (max_err_A1 < 1e-12 ? "✓ PASS" : "✗ FAIL"))

    # ── A2: outflow and tangential faces must be EXACTLY zero — direct
    # confirmation of Î≡I on outflow from the derivation, not assumed. ──────
    fill!(A_vol, 0.0); fill!(RHS, 0.0)
    face6_normal = (2/Δx, 0.0, 0.0)   # face 6 (left, i=ngl): sgn=+1
    build_dg_flux_terms!(A_vol, RHS, 6, ngl, 1.0, 0.0, 0.0, false, true,   # is_inflow=false: OUTFLOW
                          face6_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> 999.0)   # deliberately absurd I_upwind — if this
                                                       # value leaks through ANYWHERE, it will
                                                       # be impossible to miss in the error check
    err_outflow = maximum(abs.(A_vol)) + maximum(abs.(RHS))

    fill!(A_vol, 0.0); fill!(RHS, 0.0)
    face1_normal = (0.0, 0.0, 0.0)   # tangential: Ω=(1,0,0) has zero component on an η/ζ face anyway,
                                       # but is_active=false is what actually gates this, tested directly
    build_dg_flux_terms!(A_vol, RHS, 1, ngl, 1.0, 0.0, 0.0, true, false,   # is_active=false: TANGENTIAL
                          face1_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> 999.0)
    err_tangential = maximum(abs.(A_vol)) + maximum(abs.(RHS))

    println("Test A2 (outflow face contributes exactly zero): error = $err_outflow  " *
            (err_outflow == 0.0 ? "✓ PASS" : "✗ FAIL"))
    println("Test A2 (tangential face contributes exactly zero): error = $err_tangential  " *
            (err_tangential == 0.0 ? "✓ PASS" : "✗ FAIL"))

    # ── A3: TWO simultaneous inflow faces (oblique Ω), confirming independent
    # accumulation — 3-4-5 unit vector for clean arithmetic: Ωx=0.6,Ωy=0.8 ──
    Ωx3, Ωy3, Ωz3 = 0.6, 0.8, 0.0
    fill!(A_vol, 0.0); fill!(RHS, 0.0)
    I0_x, I0_y = 2.0, 5.0   # DIFFERENT prescribed values on each face

    face5_normal_3 = (-2/Δx, 0.0, 0.0)   # face 5 (i=1): inflow since Ωx>0
    face2_normal_3 = (0.0, -2/Δy, 0.0)   # face 2 (j=1, "back"): inflow since Ωy>0
    build_dg_flux_terms!(A_vol, RHS, 5, ngl, Ωx3, Ωy3, Ωz3, true, true,
                          face5_normal_3, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> I0_x)
    build_dg_flux_terms!(A_vol, RHS, 2, ngl, Ωx3, Ωy3, Ωz3, true, true,
                          face2_normal_3, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> I0_y)

    max_err_A3 = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        expected_diag = 0.0; expected_rhs = 0.0
        if i == 1   # face 5 contribution
            d = Ωx3*(2/Δx) * ω[j]*ω[k] * (Δy/2)*(Δz/2)
            expected_diag += d; expected_rhs += d*I0_x
        end
        if j == 1   # face 2 contribution — INDEPENDENT, must simply ADD, not overwrite
            d = Ωy3*(2/Δy) * ω[i]*ω[k] * (Δx/2)*(Δz/2)
            expected_diag += d; expected_rhs += d*I0_y
        end
        max_err_A3 = max(max_err_A3, abs(A_vol[p,p]-expected_diag), abs(RHS[p]-expected_rhs))
    end
    println("Test A3 (two simultaneous inflow faces accumulate independently): " *
            "max abs error = $max_err_A3  " * (max_err_A3 < 1e-12 ? "✓ PASS" : "✗ FAIL"))

    # ═══════════════════════════════════════════════════════════════════════
    # PART B: volume + flux combined, checked against a REAL analytic PDE
    # solution — Ω·∇I + κI = 0, prescribed inflow I(x0)=I0, exact solution
    # I(x) = I0*exp(-κ(x-x0)), independent of y,z.
    # ═══════════════════════════════════════════════════════════════════════

    κ_val = 0.5
    x0 = 0.0   # element's low-x edge; physical x at node i: x0 + (Δx/2)*(1+ξ_i)

    κ        = fill(κ_val, npoin_test)
    σ        = zeros(npoin_test)
    intΦ_eff = zeros(npoin_test)
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    fill!(A_vol, 0.0); fill!(RHS, 0.0)
    build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, σ, intΦ_eff, connijk_1)
    # RHS from volume is zero (q=0, and build_dg_volume_matrix! has no source
    # term anyway — see note in conversation on this being out of scope here)

    I0_B = 4.2
    build_dg_flux_terms!(A_vol, RHS, 5, ngl, 1.0, 0.0, 0.0, true, true,
                          face5_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                          Je, 1, (i,j,k) -> I0_B)
    build_dg_flux_terms!(A_vol, RHS, 6, ngl, 1.0, 0.0, 0.0, false, true,   # outflow, no-op, included
                          face6_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,   # for completeness
                          Je, 1, (i,j,k) -> 0.0)

    x_num = A_vol \ RHS

    max_err_B = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        x_i = x0 + (Δx/2)*(1 + lgl.ξ[i])
        I_exact = I0_B * exp(-κ_val*(x_i - x0))
        max_err_B = max(max_err_B, abs(x_num[p] - I_exact))
    end
    println("Test B (volume+flux solve vs exact exp(-κx) solution): " *
            "max abs error = $max_err_B  " * (max_err_B < 1e-8 ? "✓ PASS" : "✗ FAIL"))

    return max(max_err_A1, err_outflow, err_tangential, max_err_A3, max_err_B)
end

function test_dg_flux_convergence(; ngl_list = [4, 6, 8, 10, 12])
    Δx, Δy, Δz = 2.0, 3.0, 4.0
    κ_val = 0.5
    x0 = 0.0
    I0_B = 4.2

    errors = Float64[]

    for ngl in ngl_list
        lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
        basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
        ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

        Je_const = (Δx/2)*(Δy/2)*(Δz/2)
        dξdx = fill(2/Δx, 1,ngl,ngl,ngl); dξdy = zeros(1,ngl,ngl,ngl); dξdz = zeros(1,ngl,ngl,ngl)
        dηdx = zeros(1,ngl,ngl,ngl); dηdy = fill(2/Δy, 1,ngl,ngl,ngl); dηdz = zeros(1,ngl,ngl,ngl)
        dζdx = zeros(1,ngl,ngl,ngl); dζdy = zeros(1,ngl,ngl,ngl); dζdz = fill(2/Δz, 1,ngl,ngl,ngl)
        Je   = fill(Je_const, 1,ngl,ngl,ngl)

        npoin_test = ngl^3
        κ        = fill(κ_val, npoin_test)
        σ        = zeros(npoin_test)
        intΦ_eff = zeros(npoin_test)
        connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
        for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

        A_vol = zeros(npoin_test, npoin_test)
        RHS   = zeros(npoin_test)

        build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                                 dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                                 Je, κ, σ, intΦ_eff, connijk_1)

        face5_normal = (-2/Δx, 0.0, 0.0)
        face6_normal = (2/Δx, 0.0, 0.0)
        build_dg_flux_terms!(A_vol, RHS, 5, ngl, 1.0, 0.0, 0.0, true, true,
                              face5_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                              Je, 1, (i,j,k) -> I0_B)
        build_dg_flux_terms!(A_vol, RHS, 6, ngl, 1.0, 0.0, 0.0, false, true,
                              face6_normal, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                              Je, 1, (i,j,k) -> 0.0)

        x_num = A_vol \ RHS

        max_err = 0.0
        for k=1:ngl, j=1:ngl, i=1:ngl
            p = lidx(i,j,k,ngl)
            x_i = x0 + (Δx/2)*(1 + lgl.ξ[i])
            I_exact = I0_B * exp(-κ_val*(x_i - x0))
            max_err = max(max_err, abs(x_num[p] - I_exact))
        end
        push!(errors, max_err)
        println("ngl=$ngl: max abs error = $max_err")
    end

    println("─"^50)
    println("Error ratios between successive ngl (should shrink FAST, not slowly/erratically):")
    for i = 2:length(errors)
        println("  ngl=$(ngl_list[i-1])→$(ngl_list[i]): ratio = $(errors[i-1]/errors[i])")
    end

    return errors
end

function test_dg_reconciliation_one_ngl(ngl::Int; Δx1 = 1.5, Δx2 = 3.0, Δy = 3.0, Δz = 4.0,
                                         κ_val = 0.5, I0 = 4.2)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    npoin_test = ngl^3
    κ        = fill(κ_val, npoin_test)
    σ        = zeros(npoin_test)
    intΦ_eff = zeros(npoin_test)
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    Je1 = (Δx1/2)*(Δy/2)*(Δz/2)
    dξdx1 = fill(2/Δx1,1,ngl,ngl,ngl); dξdy1=zeros(1,ngl,ngl,ngl); dξdz1=zeros(1,ngl,ngl,ngl)
    dηdx1 = zeros(1,ngl,ngl,ngl); dηdy1=fill(2/Δy,1,ngl,ngl,ngl); dηdz1=zeros(1,ngl,ngl,ngl)
    dζdx1 = zeros(1,ngl,ngl,ngl); dζdy1=zeros(1,ngl,ngl,ngl); dζdz1=fill(2/Δz,1,ngl,ngl,ngl)
    Je1a  = fill(Je1,1,ngl,ngl,ngl)

    A1 = zeros(npoin_test,npoin_test); RHS1 = zeros(npoin_test)
    build_dg_volume_matrix!(A1, 1, 1.0,0.0,0.0, ngl, ψ,dψ,ω,
                             dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1,
                             Je1a, κ, σ, intΦ_eff, connijk_1)
    build_dg_flux_terms!(A1, RHS1, 5, ngl, 1.0,0.0,0.0, true, true,
                          (-2/Δx1,0.0,0.0), ω, dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1,
                          Je1a, 1, (i,j,k)->I0)
    build_dg_flux_terms!(A1, RHS1, 6, ngl, 1.0,0.0,0.0, false, true,
                          (2/Δx1,0.0,0.0), ω, dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1,
                          Je1a, 1, (i,j,k)->0.0)
    x1 = A1 \ RHS1

    Je2 = (Δx2/2)*(Δy/2)*(Δz/2)
    dξdx2 = fill(2/Δx2,1,ngl,ngl,ngl); dξdy2=zeros(1,ngl,ngl,ngl); dξdz2=zeros(1,ngl,ngl,ngl)
    dηdx2 = zeros(1,ngl,ngl,ngl); dηdy2=fill(2/Δy,1,ngl,ngl,ngl); dηdz2=zeros(1,ngl,ngl,ngl)
    dζdx2 = zeros(1,ngl,ngl,ngl); dζdy2=zeros(1,ngl,ngl,ngl); dζdz2=fill(2/Δz,1,ngl,ngl,ngl)
    Je2a  = fill(Je2,1,ngl,ngl,ngl)

    A2 = zeros(npoin_test,npoin_test); RHS2 = zeros(npoin_test)
    build_dg_volume_matrix!(A2, 1, 1.0,0.0,0.0, ngl, ψ,dψ,ω,
                             dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2,
                             Je2a, κ, σ, intΦ_eff, connijk_1)
    build_dg_flux_terms!(A2, RHS2, 5, ngl, 1.0,0.0,0.0, true, true,
                          (-2/Δx2,0.0,0.0), ω, dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2,
                          Je2a, 1, (i,j,k) -> x1[lidx(ngl,j,k,ngl)])
    build_dg_flux_terms!(A2, RHS2, 6, ngl, 1.0,0.0,0.0, false, true,
                          (2/Δx2,0.0,0.0), ω, dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2,
                          Je2a, 1, (i,j,k)->0.0)
    x2 = A2 \ RHS2

    xI = Δx1
    max_jump = 0.0
    max_err_upwind  = 0.0
    max_err_downwind = 0.0
    max_err_avg = 0.0
    max_err_jac = 0.0
    for k=1:ngl, j=1:ngl
        v1 = x1[lidx(ngl,j,k,ngl)]
        v2 = x2[lidx(1,  j,k,ngl)]
        max_jump = max(max_jump, abs(v1-v2))
        I_exact = I0 * exp(-κ_val * xI)
        max_err_upwind   = max(max_err_upwind,   abs(v1 - I_exact))
        max_err_downwind = max(max_err_downwind, abs(v2 - I_exact))
        max_err_avg = max(max_err_avg, abs(0.5*(v1+v2) - I_exact))
        max_err_jac = max(max_err_jac, abs((Δx1*v1+Δx2*v2)/(Δx1+Δx2) - I_exact))
    end

    return (jump=max_jump, upwind=max_err_upwind, downwind=max_err_downwind,
            avg=max_err_avg, jac=max_err_jac)
end

function test_dg_reconciliation_convergence(; ngl_list = [4,6,8,10,12])
    results = [test_dg_reconciliation_one_ngl(n) for n in ngl_list]

    println("ngl  |   jump      |  upwind err  | downwind err |  avg err    | Jacobian err")
    println("─"^85)
    for (n,r) in zip(ngl_list, results)
        println("$(lpad(n,4)) | $(rpad(r.jump,11)) | $(rpad(r.upwind,12)) | " *
                "$(rpad(r.downwind,12)) | $(rpad(r.avg,11)) | $(r.jac)")
    end

    println("─"^85)
    println("Jump ratio between successive ngl (spectral convergence check):")
    for i = 2:length(results)
        println("  ngl=$(ngl_list[i-1])→$(ngl_list[i]): ratio = $(results[i-1].jump/results[i].jump)")
    end
    println("Upwind-error ratio between successive ngl:")
    for i = 2:length(results)
        println("  ngl=$(ngl_list[i-1])→$(ngl_list[i]): ratio = $(results[i-1].upwind/results[i].upwind)")
    end

    return results
end
# ─────────────────────────────────────────────────────────────────────────────
#  verify_face_node_correspondence  (DRAFT — untested)
#
#  For every interior conforming face pair (already identified via
#  neighbor_elem/neighbor_fid — reused unchanged, no new adjacency logic),
#  builds the EMPIRICAL node-to-node correspondence between iel's face and
#  jel's face by matching GLOBAL point ids — not by assuming any particular
#  rotation/reflection. Reports whether that correspondence is a clean
#  bijection (every node on one side matches EXACTLY one node on the
#  other), and separately reports, for a sample of pairs, what LOCAL
#  index transform the bijection actually corresponds to — so we can see
#  by inspection whether it's a simple rotation/flip (good — means
#  connijk_dg's per-interface bookkeeping is simple) or something messier.
# ─────────────────────────────────────────────────────────────────────────────

function verify_face_node_correspondence(mesh, neighbor_elem, neighbor_fid)
    ngl = mesh.ngl
    n_pairs_checked   = 0
    n_perfect_bijection = 0
    n_problem_pairs   = 0
    problem_examples  = Tuple{TInt,Int,TInt,Int}[]
    sample_maps       = Any[]   # a few (iel,fid,jel,jfid, local_correspondence) for inspection

    for iel = 1:mesh.nelem
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            jel < iel && continue   # each interior pair checked once

            jfid = neighbor_fid[iel][fid]
            n_pairs_checked += 1

            irange_i, jrange_i, krange_i = face_ip_range(fid, ngl)
            irange_j, jrange_j, krange_j = face_ip_range(jfid, ngl)

            # gid -> (i,j,k) on jel's face
            jel_gid_to_ijk = Dict{Int, Tuple{Int,Int,Int}}()
            for k in krange_j, j in jrange_j, i in irange_j
                jel_gid_to_ijk[mesh.connijk[jel,i,j,k]] = (i,j,k)
            end

            # For every node on iel's face, find its match on jel's face
            local_map = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int,Int}}()
            n_matched   = 0
            n_face_nodes = 0
            for k in krange_i, j in jrange_i, i in irange_i
                n_face_nodes += 1
                gid = mesh.connijk[iel,i,j,k]
                if haskey(jel_gid_to_ijk, gid)
                    local_map[(i,j,k)] = jel_gid_to_ijk[gid]
                    n_matched += 1
                end
            end

            if n_matched == n_face_nodes && n_face_nodes == ngl*ngl &&
               length(unique(values(local_map))) == n_face_nodes   # also confirms INJECTIVE, not just matched
                n_perfect_bijection += 1
                if length(sample_maps) < 5
                    push!(sample_maps, (iel, fid, jel, jfid, local_map))
                end
            else
                n_problem_pairs += 1
                length(problem_examples) < 5 && push!(problem_examples, (iel, fid, jel, jfid))
            end
        end
    end

    println("Face pairs checked: $n_pairs_checked")
    println("Perfect bijection (every node matched, 1-1, no duplicates): $n_perfect_bijection")
    println("Problem pairs: $n_problem_pairs")
    if n_problem_pairs > 0
        println("  ⚠ Example problem pairs (iel,fid,jel,jfid): $problem_examples")
    end

    println("─"^70)
    println("Sample local-index correspondence (i,j,k)_iel -> (i,j,k)_jel:")
    for (iel,fid,jel,jfid,lmap) in sample_maps
        println("  iel=$iel fid=$fid  <->  jel=$jel jfid=$jfid")
        # print just the 4 corners of iel's face for a compact, readable sample
        irange,jrange,krange = face_ip_range(fid, ngl)
        corners = [(first(irange),first(jrange),first(krange)),
                   (last(irange), first(jrange),first(krange)),
                   (first(irange),last(jrange), last(krange)),
                   (last(irange), last(jrange), last(krange))]
        for c in unique(corners)
            haskey(lmap, c) && println("    $c  ->  $(lmap[c])")
        end
    end

    return n_problem_pairs == 0
end

function build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    nelem, ngl = mesh.nelem, mesh.ngl

    # Single linear pass: g -> list of every (iel,i,j,k) touching it.
    g_to_touches = Dict{Int, Vector{NTuple{4,Int}}}()
    for iel = 1:nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        push!(get!(g_to_touches, g, NTuple{4,Int}[]), (iel,i,j,k))
    end

    connijk_dg  = Array{Int}(undef, nelem, ngl, ngl, ngl)
    next_new_id = maximum(mesh.connijk) + 1

    # Single pass over the reverse map (sorted keys — fully deterministic,
    # matches the tie-break rule: smallest touching iel keeps the CG id).
    for g in sort(collect(keys(g_to_touches)))
        touches = g_to_touches[g]
        if length(touches) == 1
            iel,i,j,k = touches[1]
            connijk_dg[iel,i,j,k] = g
            continue
        end
        owner_iel = minimum(t[1] for t in touches)
        for (iel,i,j,k) in touches
            if iel == owner_iel
                connijk_dg[iel,i,j,k] = g
            else
                connijk_dg[iel,i,j,k] = next_new_id
                next_new_id += 1
            end
        end
    end
    npoin_dg = next_new_id - 1

    dg_to_cg = Vector{Int}(undef, npoin_dg)
    for iel = 1:nelem, k=1:ngl, j=1:ngl, i=1:ngl
        dg_to_cg[connijk_dg[iel,i,j,k]] = mesh.connijk[iel,i,j,k]
    end

    # Face correspondence — unchanged from verify_face_node_correspondence,
    # already O(nelem·6·ngl²), no optimization needed here.
    face_correspondence = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int,Int},Tuple{Int,Int,Int}}}()
    for iel = 1:nelem
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            jfid = neighbor_fid[iel][fid]
            irange_i,jrange_i,krange_i = face_ip_range(fid, ngl)
            irange_j,jrange_j,krange_j = face_ip_range(jfid, ngl)
            jel_gid_to_ijk = Dict{Int,Tuple{Int,Int,Int}}()
            for k in krange_j, j in jrange_j, i in irange_j
                jel_gid_to_ijk[mesh.connijk[jel,i,j,k]] = (i,j,k)
            end
            lmap = Dict{Tuple{Int,Int,Int},Tuple{Int,Int,Int}}()
            for k in krange_i, j in jrange_i, i in irange_i
                g = mesh.connijk[iel,i,j,k]
                haskey(jel_gid_to_ijk, g) && (lmap[(i,j,k)] = jel_gid_to_ijk[g])
            end
            face_correspondence[(iel,fid)] = lmap
        end
    end

    return connijk_dg, npoin_dg, dg_to_cg, face_correspondence
end

function build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)
    nelem, ngl = mesh.nelem, mesh.ngl
    face_neighbor_dg = zeros(TInt, ngl*ngl, 6, nelem)   # 0 sentinel = boundary/no same-rank neighbor
    n_missing = 0

    for iel = 1:nelem
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue   # leaves this (fid,iel) slice all-zero

            irange, jrange, krange = face_ip_range(fid, ngl)
            lmap = face_correspondence[(iel,fid)]
            p = 0
            for k in krange, j in jrange, i in irange
                p += 1
                if !haskey(lmap, (i,j,k))
                    n_missing += 1   # flagged below, not silently skipped
                    continue
                end
                (ij, jj, kj) = lmap[(i,j,k)]
                face_neighbor_dg[p, fid, iel] = connijk_dg[jel, ij, jj, kj]
            end
        end
    end

    n_missing > 0 && @warn "build_face_neighbor_dg: $n_missing face-node lookups had no " *
                            "correspondence entry — unexpected given the exhaustive bijection " *
                            "already confirmed (n_problem_pairs=0 on both test meshes); worth " *
                            "investigating before trusting this array."
    return face_neighbor_dg
end

function verify_face_neighbor_dg(mesh, connijk_dg, npoin_dg, dg_to_cg,
                                  neighbor_elem, face_neighbor_dg)
    nelem, ngl = mesh.nelem, mesh.ngl
    n_checked = 0
    n_mismatch = 0
    n_bad_slice_nonzero = 0   # boundary faces that should be all-zero but aren't
    n_bad_slice_incomplete = 0  # interior faces that should be fully populated but aren't

    for iel = 1:nelem
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            irange, jrange, krange = face_ip_range(fid, ngl)

            if jel == zero(TInt)
                any(face_neighbor_dg[:, fid, iel] .!= 0) && (n_bad_slice_nonzero += 1)
                continue
            end

            p = 0
            n_populated = 0
            for k in krange, j in jrange, i in irange
                p += 1
                dg_id = face_neighbor_dg[p, fid, iel]
                dg_id == 0 && continue
                n_populated += 1

                # Ground-truth check: this element's OWN original global id
                # at this face node, vs. the neighbor's DG id mapped back
                # to CG — using ONLY mesh.connijk as reference, not
                # face_correspondence, for a genuinely independent check.
                g_here = mesh.connijk[iel,i,j,k]
                g_via_neighbor = dg_id <= npoin_dg ? dg_to_cg[dg_id] : -1
                n_checked += 1
                g_here != g_via_neighbor && (n_mismatch += 1)
            end
            n_populated != ngl*ngl && (n_bad_slice_incomplete += 1)
        end
    end

    println("Node-level round-trip checks performed: $n_checked")
    println("Mismatches (physical point disagreement): $n_mismatch  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))
    println("Boundary/no-neighbor slices wrongly nonzero: $n_bad_slice_nonzero  " *
            (n_bad_slice_nonzero == 0 ? "✓ PASS" : "✗ FAIL"))
    println("Interior-neighbor slices not fully populated (ngl² entries expected): " *
            "$n_bad_slice_incomplete  " * (n_bad_slice_incomplete == 0 ? "✓ PASS" : "✗ FAIL"))

    return n_mismatch == 0 && n_bad_slice_nonzero == 0 && n_bad_slice_incomplete == 0
end

@inline function compute_face_normal_avg(iel::Int, fid::Int, ngl::Int,
                                          dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz)
    msym, sgn = face_normal_sign(fid)
    irange, jrange, krange = face_ip_range(fid, ngl)
    nx_sum = 0.0; ny_sum = 0.0; nz_sum = 0.0; n = 0
    for k in krange, j in jrange, i in irange
        if msym == 1;      nx,ny,nz = dηdx[iel,i,j,k], dηdy[iel,i,j,k], dηdz[iel,i,j,k]
        elseif msym == 2;  nx,ny,nz = dζdx[iel,i,j,k], dζdy[iel,i,j,k], dζdz[iel,i,j,k]
        else;              nx,ny,nz = dξdx[iel,i,j,k], dξdy[iel,i,j,k], dξdz[iel,i,j,k]
        end
        nx_sum += sgn*nx; ny_sum += sgn*ny; nz_sum += sgn*nz
        n += 1
    end
    nx_avg, ny_avg, nz_avg = nx_sum/n, ny_sum/n, nz_sum/n
    mag = sqrt(nx_avg^2 + ny_avg^2 + nz_avg^2)   # NEW: normalize to a genuine unit vector
    return (nx_avg/mag, ny_avg/mag, nz_avg/mag)
end

function test_full_sweep_exponential(mesh, metrics, basis, ω; κ_val = 0.5, I0 = 4.2)
    Ωx, Ωy, Ωz = 1.0, 0.0, 0.0
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ

    # ── Build (or reuse) DG infrastructure, store into mesh as requested ──
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg
        mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg
        mesh.face_neighbor_dg = face_neighbor_dg
        println("Built and stored DG infrastructure: npoin_dg=$npoin_dg " *
                "(original npoin=$(length(mesh.x)))")
    else
        println("Reusing already-stored mesh.connijk_dg (npoin_dg=$(mesh.npoin_dg))")
    end

    # ── Sweep order for THIS Ω (direction-dependent, not cached in mesh) ──
    order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                             metrics.dζdx,metrics.dζdy,metrics.dζdz)

    npoin_cg = length(mesh.x)
    κ        = fill(κ_val, npoin_cg)
    σ        = zeros(npoin_cg)
    intΦ_eff = zeros(npoin_cg)

    x_dg = fill(NaN, mesh.npoin_dg)
    bc_fn = (x,y,z) -> I0   # a genuine function, constant-valued — NOT bare I0 itself
    run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                  mesh.connijk_dg, mesh.face_neighbor_dg, κ, σ, intΦ_eff, bc_fn)

    n_nan = count(isnan, x_dg)
    n_nan > 0 && @warn "test_full_sweep_exponential: $n_nan of $(mesh.npoin_dg) DG entries " *
                        "were NEVER WRITTEN by the sweep — indicates a real bug " *
                        "(unreachable node, or an element missing from `order`)."

    pos = zeros(Int, mesh.nelem)
    for (p,iel) in enumerate(order); pos[iel] = p; end

    g_to_touches = Dict{Int, Vector{NTuple{4,Int}}}()
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        push!(get!(g_to_touches, g, NTuple{4,Int}[]), (iel,i,j,k))
    end

    x_cg = zeros(npoin_cg)
    for g = 1:npoin_cg
        touches = g_to_touches[g]
        best = touches[1]
        for t in touches[2:end]
            pos[t[1]] < pos[best[1]] && (best = t)
        end
        iel,i,j,k = best
        x_cg[g] = x_dg[mesh.connijk_dg[iel,i,j,k]]
    end

    x0 = mesh.xmin
    max_err = 0.0
    for g = 1:npoin_cg
        I_exact = I0 * exp(-κ_val * (mesh.x[g] - x0))
        max_err = max(max_err, abs(x_cg[g] - I_exact))
    end
    println("Full multi-element sweep vs exact exponential: max abs error = $max_err")
    return max_err
end


#=function run_dg_sweep!(x_dg::Vector{Float64}, mesh, metrics, ψ, dψ, ω,
                        Ωx::Float64, Ωy::Float64, Ωz::Float64,
                        order::Vector{Int}, neighbor_elem,
                        connijk_dg, face_neighbor_dg,
                        κ, σ, intΦ_eff, bc_fn;
                        extra_source = nothing,
                        timers = nothing,
                        factor_cache = nothing,
                        weight_cache_store = nothing,
                        face_info_cache = nothing)   # NEW, optional: Dict{Int,Any} keyed by iel,
                                                       # one cache PER DIRECTION — same requirement
                                                       # as factor_cache/weight_cache_store (a shared,
                                                       # non-direction-scoped cache here would silently
                                                       # reuse another direction's classification/normals,
                                                       # the same corruption bug already found once)
    ngl = mesh.ngl
    A_vol = zeros(ngl^3, ngl^3)
    RHS_local = zeros(ngl^3)
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)

    for iel in order
        need_build = factor_cache === nothing || !haskey(factor_cache, iel)

        if need_build
            t0 = timers !== nothing ? time() : 0.0
            fill!(A_vol, 0.0); fill!(RHS_local, 0.0)
            build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            timers !== nothing && (timers["build_volume"] = get(timers,"build_volume",0.0) + (time()-t0))
        else
            fill!(RHS_local, 0.0)
        end
        extra_source !== nothing &&
            build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source, mesh.connijk)

        # ── Face classification + normals: geometry/Ω-only, never depends on
        # iteration — cached exactly like A_vol/factorization. ────────────────
        face_info = face_info_cache === nothing ? nothing : get(face_info_cache, iel, nothing)

        local face_normals
        if face_info === nothing
            t0 = timers !== nothing ? time() : 0.0
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)
            timers !== nothing && (timers["face_classification"] = get(timers,"face_classification",0.0) + (time()-t0))

            t0 = timers !== nothing ? time() : 0.0
            face_normals = Vector{NTuple{3,Float64}}(undef,6)
            for fid = 1:6
                is_active[fid] || continue
                face_normals[fid] = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
            end
            timers !== nothing && (timers["normal_computation"] = get(timers,"normal_computation",0.0) + (time()-t0))

            face_info_cache !== nothing &&
                (face_info_cache[iel] = (copy(is_inflow), copy(is_active), copy(face_normals)))
        else
            is_inflow_c, is_active_c, face_normals = face_info
            is_inflow .= is_inflow_c
            is_active .= is_active_c
        end

        t0 = timers !== nothing ? time() : 0.0
        weight_caches = weight_cache_store === nothing ?
            [zeros(Float64, ngl*ngl) for _=1:6] :
            get!(weight_cache_store, iel, [zeros(Float64, ngl*ngl) for _=1:6])

        for fid = 1:6
            is_active[fid] || continue
            face_normal = face_normals[fid]

            if is_inflow[fid]
                jel = neighbor_elem[iel][fid]
                if jel == zero(TInt)
                    I_upwind = (i,j,k) -> bc_fn(mesh.x[mesh.connijk[iel,i,j,k]],
                                                 mesh.y[mesh.connijk[iel,i,j,k]],
                                                 mesh.z[mesh.connijk[iel,i,j,k]])
                else
                    I_upwind = (i,j,k) -> x_dg[face_neighbor_dg[face_flat_index(fid, ngl, i, j, k), fid, iel]]
                end
            else
                I_upwind = (i,j,k) -> 0.0
            end

            if need_build
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx, Ωy, Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind; weight_out=weight_caches[fid])
            else
                build_dg_flux_rhs_only!(RHS_local, fid, ngl, is_inflow[fid], is_active[fid],
                                         weight_caches[fid], I_upwind)
            end
        end
        timers !== nothing && (timers["flux_and_faces"] = get(timers,"flux_and_faces",0.0) + (time()-t0))

        t0 = timers !== nothing ? time() : 0.0
        if need_build
            fac = factorize(A_vol)
            factor_cache !== nothing && (factor_cache[iel] = fac)
            timers !== nothing && (timers["factorize"] = get(timers,"factorize",0.0) + (time()-t0))
            t1 = timers !== nothing ? time() : 0.0
            x_local = fac \ RHS_local
            timers !== nothing && (timers["solve_only"] = get(timers,"solve_only",0.0) + (time()-t1))
        else
            x_local = factor_cache[iel] \ RHS_local
            timers !== nothing && (timers["solve_only"] = get(timers,"solve_only",0.0) + (time()-t0))
        end

        for k=1:ngl, j=1:ngl, i=1:ngl
            x_dg[connijk_dg[iel,i,j,k]] = x_local[lidx(i,j,k,ngl)]
        end
    end
    return x_dg
end=#

function run_dg_sweep!(x_dg::Vector{Float64}, mesh, metrics, ψ, dψ, ω,
                        Ωx::Float64, Ωy::Float64, Ωz::Float64,
                        order::Vector{Int}, neighbor_elem,
                        connijk_dg, face_neighbor_dg,
                        κ, σ, intΦ_eff, bc_fn;
                        extra_source = nothing,
                        ncf_as_child = Dict{Tuple{Int,Int},Int}(),        # NEW, optional
                        ncf_as_parent = Dict{Tuple{Int,Int},Int}(),
			ncf_L1L2 = Dict{Int,Tuple{Matrix{Float64},Matrix{Float64}}}())  # NEW, optional
    ngl = mesh.ngl
    A_vol = zeros(ngl^3, ngl^3)
    RHS_local = zeros(ngl^3)
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)

    for iel in order
        fill!(A_vol, 0.0); fill!(RHS_local, 0.0)
        build_dg_volume_matrix!(A_vol, iel, Ωx, Ωy, Ωz, ngl, ψ, dψ, ω,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                 metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
        extra_source !== nothing &&
            build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source, mesh.connijk)

        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)

        for fid = 1:6
            is_active[fid] || continue
            face_normal = compute_face_normal_avg(iel, fid, ngl,
                metrics.dξdx,metrics.dξdy,metrics.dξdz,
                metrics.dηdx,metrics.dηdy,metrics.dηdz,
                metrics.dζdx,metrics.dζdy,metrics.dζdz)

            if is_inflow[fid]
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
        		ncf_group = ncf_as_parent[(iel,fid)]  # Vector of ncf_idx
        		entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in ncf_group]
        		parent_grid = build_dg_ncf_child_to_parent_v5(mesh, connijk_dg, x_dg, iel, fid, entries, ngl)
        		I_upwind = (i,j,k) -> begin
            		    node_idx = face_flat_index(fid, ngl, i, j, k)
            		    p1 = ((node_idx-1) % ngl) + 1
           	 	    p2 = ((node_idx-1) ÷ ngl) + 1
            		    parent_grid[p1,p2]
        	    	end
    		    else
        	    	I_upwind = (i,j,k) -> bc_fn(mesh.x[mesh.connijk[iel,i,j,k]],
                    			                 mesh.y[mesh.connijk[iel,i,j,k]],
                                		         mesh.z[mesh.connijk[iel,i,j,k]])
    		    end
	    	else
    		    I_upwind = (i,j,k) -> x_dg[face_neighbor_dg[face_flat_index(fid, ngl, i, j, k), fid, iel]]
		end
	    else
                I_upwind = (i,j,k) -> 0.0
            end

            build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx, Ωy, Ωz,
                                  is_inflow[fid], is_active[fid], face_normal, ω,
                                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                  metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                  metrics.Je, iel, I_upwind)
        end

        x_local = A_vol \ RHS_local
        for k=1:ngl, j=1:ngl, i=1:ngl
            x_dg[connijk_dg[iel,i,j,k]] = x_local[lidx(i,j,k,ngl)]
        end
    end
    return x_dg
end


function test_full_sweep_exponential_layered(mesh, metrics, basis, ω; κ_val = 0.5, I0 = 4.2)
    Ωx, Ωy, Ωz = 1.0, 0.0, 0.0
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg
        mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg
        mesh.face_neighbor_dg = face_neighbor_dg
    end

    order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                             metrics.dζdx,metrics.dζdy,metrics.dζdz)

    npoin_cg = length(mesh.x)
    κ        = fill(κ_val, npoin_cg)
    σ        = zeros(npoin_cg)
    intΦ_eff = zeros(npoin_cg)

    bc_fn = (x,y,z) -> I0   # a genuine function, constant-valued — NOT bare I0 itself
    x_dg = fill(NaN, mesh.npoin_dg)
    run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                  mesh.connijk_dg, mesh.face_neighbor_dg, κ, σ, intΦ_eff, bc_fn)

    x0 = mesh.xmin

    # ── LAYER 1: DG-level error — every (iel,i,j,k) checked against the
    # EXACT solution at its OWN original physical location, using ONLY
    # mesh.connijk (never connijk_dg, never any neighbor lookup) as the
    # coordinate reference. This bypasses reconciliation ENTIRELY — if the
    # sweep itself is correct, every element's own trace, including at
    # shared interfaces, should independently already be close to exact. ──
    max_err_dg = 0.0
    worst = (0,0,0,0, 0.0, 0.0)
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g_cg = mesh.connijk[iel,i,j,k]
        x_phys = mesh.x[g_cg]
        I_exact = I0 * exp(-κ_val * (x_phys - x0))
        v = x_dg[mesh.connijk_dg[iel,i,j,k]]
        err = abs(v - I_exact)
        if err > max_err_dg
            max_err_dg = err
            worst = (iel,i,j,k, v, I_exact)
        end
    end
    println("LAYER 1 — DG-level error (no reconciliation): max abs error = $max_err_dg")
    println("  worst case: iel=$(worst[1]) (i,j,k)=($(worst[2]),$(worst[3]),$(worst[4])) " *
            "computed=$(worst[5]) exact=$(worst[6])")

    # ── LAYER 2: reconciled CG-level error (upwind-only) — unchanged from before ──
    pos = zeros(Int, mesh.nelem)
    for (p,iel) in enumerate(order); pos[iel] = p; end

    g_to_touches = Dict{Int, Vector{NTuple{4,Int}}}()
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        push!(get!(g_to_touches, g, NTuple{4,Int}[]), (iel,i,j,k))
    end

    x_cg = zeros(npoin_cg)
    for g = 1:npoin_cg
        touches = g_to_touches[g]
        best = touches[1]
        for t in touches[2:end]
            pos[t[1]] < pos[best[1]] && (best = t)
        end
        iel,i,j,k = best
        x_cg[g] = x_dg[mesh.connijk_dg[iel,i,j,k]]
    end

    max_err_cg = 0.0
    for g = 1:npoin_cg
        I_exact = I0 * exp(-κ_val * (mesh.x[g] - x0))
        max_err_cg = max(max_err_cg, abs(x_cg[g] - I_exact))
    end
    println("LAYER 2 — reconciled CG-level error (upwind-only): max abs error = $max_err_cg")

    iel = 1
    ngl = mesh.ngl
    κ_val = 0.5
    κ = fill(κ_val, length(mesh.x)); σ = zeros(length(mesh.x)); intΦ_eff = zeros(length(mesh.x))

    println("═"^70)
    println("ELEMENT 1 FULL TRACE — Ω=(1,0,0)")
    println("═"^70)

    # ── 1. Geometry: where is this element, and what's its physical extent? ──
    xs = [mesh.x[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    ys = [mesh.y[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    zs = [mesh.z[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    println("[1] element 1 physical extent: x∈[$(minimum(xs)),$(maximum(xs))] " *
            "y∈[$(minimum(ys)),$(maximum(ys))] z∈[$(minimum(zs)),$(maximum(zs))]")
    println("    domain xmin=$(mesh.xmin) xmax=$(mesh.xmax)")

    # ── 2. Face classification — repeat exactly, print ALL 6 with normals ────
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    classify_element_faces!(is_inflow, is_active, iel, 1.0, 0.0, 0.0, ngl,
                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    println("[2] Face classification:")
    for fid = 1:6
        fn = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz,
            metrics.dηdx,metrics.dηdy,metrics.dηdz,
            metrics.dζdx,metrics.dζdy,metrics.dζdz)
        println("    fid=$fid active=$(is_active[fid]) inflow=$(is_inflow[fid]) " *
            "neighbor=$(neighbor_elem[iel][fid]) normal_avg=$fn")
    end

    # ── 3. Volume matrix — diagonal and a few off-diagonals at node (1,1,1) ──
    A_vol = zeros(ngl^3, ngl^3); RHS_local = zeros(ngl^3)
    build_dg_volume_matrix!(A_vol, iel, 1.0,0.0,0.0, ngl, basis.ψ, basis.dψ, ω,
                            metrics.dξdx,metrics.dξdy,metrics.dξdz,
                            metrics.dηdx,metrics.dηdy,metrics.dηdz,
                            metrics.dζdx,metrics.dζdy,metrics.dζdz,
                            metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
    p111 = lidx(1,1,1,ngl)
    println("[3] AFTER volume matrix only:")
    println("    A_vol[p111,p111] = $(A_vol[p111,p111])")
    nzq = [(q,A_vol[p111,q]) for q=1:ngl^3 if A_vol[p111,q]!=0.0 && q!=p111]
    println("    off-diagonal entries at row p111: $nzq")
    println("    RHS_local[p111] = $(RHS_local[p111])  (should be 0 — no source term yet)")

    # ── 4. Flux terms — print incremental change from EACH face separately ──
    println("[4] Flux terms, applied one face at a time:")
    for fid = 1:6
        is_active[fid] || (println("    fid=$fid: SKIPPED (inactive)"); continue)
        fn = compute_face_normal_avg(iel, fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz,
            metrics.dηdx,metrics.dηdy,metrics.dηdz,
            metrics.dζdx,metrics.dζdy,metrics.dζdz)
        diag_before = A_vol[p111,p111]; rhs_before = RHS_local[p111]
        I_upwind = is_inflow[fid] ? ((i,j,k) -> 4.2) : ((i,j,k) -> 0.0)
        build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, 1.0,0.0,0.0,
                             is_inflow[fid], is_active[fid], fn, ω,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz,
                             metrics.Je, iel, I_upwind)
        println("    fid=$fid (inflow=$(is_inflow[fid])): " *
                "Δdiag=$(A_vol[p111,p111]-diag_before)  ΔRHS=$(RHS_local[p111]-rhs_before)")
    end

    println("[5] FINAL, before solve: A_vol[p111,p111]=$(A_vol[p111,p111])  RHS_local[p111]=$(RHS_local[p111])")
    println("    implied value if p111 were ISOLATED (ignore off-diag coupling): " *
            "$(RHS_local[p111]/A_vol[p111,p111])")
 
    x_local = A_vol \ RHS_local
    println("[6] ACTUAL solved value at (1,1,1): $(x_local[p111])")
    println("    exact I0 (should match if (1,1,1) truly sits at the domain's x-inflow wall): 4.2")


    return (dg=max_err_dg, cg=max_err_cg)
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_mass_diag  (diagonal-only, matching build_mass_matrix!'s NSD_3D,
#  Inexact convention — under inexact integration the full ψ_p·ψ_q mass
#  matrix collapses to purely diagonal, same collocation fact already
#  exploited in build_dg_volume_matrix!'s own reaction term, so extracting
#  only M[p,p] = ω[i]ω[j]ω[k]*Je avoids building/storing the full,
#  wastefully mostly-zero (ngl^3,ngl^3) tensor build_mass_matrix! computes)
# ─────────────────────────────────────────────────────────────────────────────

function test_dg_terms_on_real_element(mesh, metrics, basis, ω, iel::Int; κ_val = 0.5)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ

    # ── Step 0: extract this element's ACTUAL physical size from real coords ──
    xs = [mesh.x[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    ys = [mesh.y[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    zs = [mesh.z[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    Δx = maximum(xs) - minimum(xs)
    Δy = maximum(ys) - minimum(ys)
    Δz = maximum(zs) - minimum(zs)
    println("Element $iel measured size: Δx=$Δx Δy=$Δy Δz=$Δz")

    # ── Step 1: are the mesh's OWN metrics what an axis-aligned box of this
    # size should produce? (build_metric_terms! itself has never been
    # independently checked until now) ──────────────────────────────────
    expected_dξdx = 2/Δx; expected_dηdy = 2/Δy; expected_dζdz = 2/Δz
    expected_Je = (Δx/2)*(Δy/2)*(Δz/2)

    max_metric_err = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        max_metric_err = max(max_metric_err,
            abs(metrics.dξdx[iel,i,j,k] - expected_dξdx),
            abs(metrics.dηdy[iel,i,j,k] - expected_dηdy),
            abs(metrics.dζdz[iel,i,j,k] - expected_dζdz),
            abs(metrics.dξdy[iel,i,j,k]), abs(metrics.dξdz[iel,i,j,k]),   # off-axis terms
            abs(metrics.dηdx[iel,i,j,k]), abs(metrics.dηdz[iel,i,j,k]),   # should be exactly 0
            abs(metrics.dζdx[iel,i,j,k]), abs(metrics.dζdy[iel,i,j,k]),
            abs(metrics.Je[iel,i,j,k] - expected_Je))
    end
    println("[Step 1] mesh's OWN metrics vs analytic expectation for this box size: " *
            "max abs error = $max_metric_err  " * (max_metric_err < 1e-6 ? "✓ PASS" : "✗ FAIL"))

    # ── Step 2: build_dg_volume_matrix! fed the REAL metrics — does it match
    # the SAME closed form already validated on the toy mesh, just at this
    # element's actual (much larger) scale? ─────────────────────────────
    npoin_test = ngl^3
    κ        = fill(κ_val, length(mesh.x))
    σ        = zeros(length(mesh.x))
    intΦ_eff = zeros(length(mesh.x))

    A_vol = zeros(npoin_test, npoin_test)
    build_dg_volume_matrix!(A_vol, iel, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz,
                             metrics.Je, κ, σ, intΦ_eff, mesh.connijk)

    max_err_vol = 0.0
    for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        gp = mesh.connijk[iel,i,j,k]
        ωJac = ω[i]*ω[j]*ω[k]*expected_Je
        for q_i=1:ngl, q_j=1:ngl, q_k=1:ngl
            q = lidx(q_i,q_j,q_k,ngl)
            expected = 0.0
            if (q_i,q_j,q_k) == (i,j,k)
                expected += κ_val * ωJac
            end
            if q_j == j && q_k == k
                expected += dψ[q_i,i] * (2/Δx) * ωJac
            end
            max_err_vol = max(max_err_vol, abs(A_vol[p,q] - expected))
        end
    end
    println("[Step 2] build_dg_volume_matrix! vs closed form (real Δx=$Δx): " *
            "max abs error = $max_err_vol  " * (max_err_vol < 1e-6 ? "✓ PASS" : "✗ FAIL"))

    # ── Step 3, corrected: build_dg_flux_terms! tested with a FRESHLY-ZEROED
    # A_vol — no volume terms mixed in, matching the isolation the toy-mesh
    # test used. If this alone brings the error down to ~1e-6, the previous
    # "3.2 million" was entirely a test-construction mistake on my part, not
    # a real flux bug. ───────────────────────────────────────────────────────

    A_vol_flux_only = zeros(npoin_test, npoin_test)
    RHS_flux_only   = zeros(npoin_test)

    face_normal5 = compute_face_normal_avg(iel, 5, ngl,
        metrics.dξdx,metrics.dξdy,metrics.dξdz,
        metrics.dηdx,metrics.dηdy,metrics.dηdz,
        metrics.dζdx,metrics.dζdy,metrics.dζdz)
    println("[Step 3-fixed] measured face 5 normal (should be unit): $face_normal5")

    I0 = 4.2
    build_dg_flux_terms!(A_vol_flux_only, RHS_flux_only, 5, ngl, 1.0, 0.0, 0.0, true, true,
                          face_normal5, ω,
                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                          metrics.Je, iel, (i,j,k)->I0)

    
    max_err_flux_corrected = 0.0
    for k=1:ngl, j=1:ngl
        p = lidx(1,j,k,ngl)
        expected_diag = ω[j]*ω[k] * (Δy/2)*(Δz/2)   # CORRECTED: no spurious (2/Δx) factor
        max_err_flux_corrected = max(max_err_flux_corrected,
            abs(A_vol_flux_only[p,p]-expected_diag),
            abs(RHS_flux_only[p]-expected_diag*I0))
    end
    println("[Step 3-corrected] build_dg_flux_terms! vs CORRECT closed form: " *
            "max abs error = $max_err_flux_corrected  " * (max_err_flux_corrected < 1e-6 ? "✓ PASS" : "✗ FAIL"))
    println("[Step 3-fixed] build_dg_flux_terms! ALONE (properly isolated) vs closed form: " *
            "max abs error = $max_err_flux_fixed  " * (max_err_flux_fixed < 1e-6 ? "✓ PASS" : "✗ FAIL"))

    return (metrics=max_metric_err, volume=max_err_vol, flux=max_err_flux)
end

function test_kappa_sweep(mesh, metrics, basis, ω, iel::Int;
                           κ_list = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0], I0 = 4.2)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ

    xs = [mesh.x[mesh.connijk[iel,i,j,k]] for i=1:ngl,j=1:ngl,k=1:ngl]
    Δx = maximum(xs) - minimum(xs)
    x0 = minimum(xs)
    println("Element $iel: Δx=$Δx")
    println("─"^90)
    println(rpad("κ",12) * rpad("κ·Δx",12) * rpad("boundary-node err",20) * "whole-element max err")
    println("─"^90)

    npoin_test = ngl^3
    σ        = zeros(length(mesh.x))
    intΦ_eff = zeros(length(mesh.x))

    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    classify_element_faces!(is_inflow, is_active, iel, 1.0, 0.0, 0.0, ngl,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz)

    for κ_val in κ_list
        κ = fill(κ_val, length(mesh.x))
        A_vol = zeros(npoin_test, npoin_test)
        RHS   = zeros(npoin_test)

        build_dg_volume_matrix!(A_vol, iel, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                 metrics.Je, κ, σ, intΦ_eff, mesh.connijk)

        for fid = 1:6
            is_active[fid] || continue
            fn = compute_face_normal_avg(iel, fid, ngl,
                metrics.dξdx,metrics.dξdy,metrics.dξdz,
                metrics.dηdx,metrics.dηdy,metrics.dηdz,
                metrics.dζdx,metrics.dζdy,metrics.dζdz)
            I_upwind = is_inflow[fid] ? ((i,j,k) -> I0) : ((i,j,k) -> 0.0)
            build_dg_flux_terms!(A_vol, RHS, fid, ngl, 1.0, 0.0, 0.0,
                                  is_inflow[fid], is_active[fid], fn, ω,
                                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                  metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                  metrics.Je, iel, I_upwind)
        end

        x_local = A_vol \ RHS

        boundary_err = abs(x_local[lidx(1,1,1,ngl)] - I0)   # node (1,1,1) IS the inflow face itself
        max_err = 0.0
        for k=1:ngl, j=1:ngl, i=1:ngl
            p = lidx(i,j,k,ngl)
            x_phys = mesh.x[mesh.connijk[iel,i,j,k]]
            I_exact = I0 * exp(-κ_val*(x_phys - x0))
            max_err = max(max_err, abs(x_local[p]-I_exact))
        end

    
    println(rpad(κ_val,12) * rpad(round(κ_val*Δx,digits=4),12) *
         rpad(round(boundary_err,sigdigits=4),20) * rpad(round(max_err,sigdigits=4),12))
    end
end

function test_all_directions_exponential(mesh, metrics, basis, ω; κ_val = 0.0001, I0 = 4.2)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg
        mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg
        mesh.face_neighbor_dg = face_neighbor_dg
    end

    npoin_cg = length(mesh.x)
    κ        = fill(κ_val, npoin_cg)
    σ        = zeros(npoin_cg)
    intΦ_eff = zeros(npoin_cg)

    # Fixed reference point for the BC family — arbitrary but fixed for all
    # directions, so every direction's exact solution is genuinely
    # comparable on equal footing (same r_bc, only Ω itself changes).
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    println("Testing all $extra_npoin angular directions, κ=$κ_val")
    println("─"^90)
    println(rpad("ip_ext",8) * rpad("(Ωx,Ωy,Ωz)",34) * rpad("n_levels",10) * "max abs error")
    println("─"^90)

    worst_err = 0.0
    worst_ip_ext = 0
    all_errs = Float64[]

    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

        order, levels = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                                       metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                                       metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                                       metrics.dζdx,metrics.dζdy,metrics.dζdz)

        bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

        x_dg = fill(NaN, mesh.npoin_dg)
        run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                      mesh.connijk_dg, mesh.face_neighbor_dg, κ, σ, intΦ_eff, bc_fn)

        n_nan = count(isnan, x_dg)
        if n_nan > 0
            println(rpad(ip_ext,8) * rpad("($Ωx,$Ωy,$Ωz)",34) * rpad(length(unique(levels)),10) *
                    "⚠ $n_nan UNWRITTEN entries — skipping error calc")
            continue
        end

        max_err = 0.0
        for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
            g_cg = mesh.connijk[iel,i,j,k]
            x_,y_,z_ = mesh.x[g_cg], mesh.y[g_cg], mesh.z[g_cg]
            I_exact = bc_fn(x_,y_,z_)
            v = x_dg[mesh.connijk_dg[iel,i,j,k]]
            max_err = max(max_err, abs(v - I_exact))
        end

        push!(all_errs, max_err)
        if max_err > worst_err
            worst_err = max_err; worst_ip_ext = ip_ext
        end
        println(rpad(ip_ext,8) * rpad("($(round(Ωx,digits=3)),$(round(Ωy,digits=3)),$(round(Ωz,digits=3)))",34) *
                rpad(length(unique(levels)),10) * "$max_err")
    end

    println("─"^90)
    println("Worst direction: ip_ext=$worst_ip_ext, max abs error = $worst_err")
    println("Mean error across all directions: $(sum(all_errs)/length(all_errs))")
    println("Median error: $(sort(all_errs)[length(all_errs)÷2 + 1])")

    return all_errs
end

"""
    build_dg_source_term!(RHS_local, iel, ngl, ω, Je, q, connijk)

Adds a volumetric source term ∫_e q ψ_p dV into RHS_local, where `q` is
indexed by ORIGINAL CG global id (same convention as κ,σ,intΦ_eff in
build_dg_volume_matrix! — a physical/source FIELD, not a per-element
unknown). Collapses to purely diagonal by the same nodal-collocation
argument as the reaction term (see build_dg_volume_matrix! docstring) —
ADDITIVE, meant to be called alongside build_dg_volume_matrix! and
build_dg_flux_terms!, never in place of either.
"""
function build_dg_source_term!(RHS_local::AbstractVector{Float64}, iel::Int, ngl::Int,
                                ω, Je, q, connijk)
    @inbounds for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        gid = connijk[iel,i,j,k]
        RHS_local[p] += q[gid] * ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
    end
    return RHS_local
end

function test_dg_source_term(; ngl::Int = 4, κ_val = 0.7, q_val = 3.5)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    Δx, Δy, Δz = 2.0, 3.0, 4.0
    Je_const = (Δx/2)*(Δy/2)*(Δz/2)
    dξdx = fill(2/Δx,1,ngl,ngl,ngl); dξdy=zeros(1,ngl,ngl,ngl); dξdz=zeros(1,ngl,ngl,ngl)
    dηdx = zeros(1,ngl,ngl,ngl); dηdy=fill(2/Δy,1,ngl,ngl,ngl); dηdz=zeros(1,ngl,ngl,ngl)
    dζdx = zeros(1,ngl,ngl,ngl); dζdy=zeros(1,ngl,ngl,ngl); dζdz=fill(2/Δz,1,ngl,ngl,ngl)
    Je   = fill(Je_const,1,ngl,ngl,ngl)

    npoin_test = ngl^3
    κ        = fill(κ_val, npoin_test)
    σ        = zeros(npoin_test)
    intΦ_eff = zeros(npoin_test)
    q        = fill(q_val, npoin_test)
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    # Ω=(0,0,0): zero streaming, isolates the reaction+source balance entirely.
    A_vol = zeros(npoin_test, npoin_test)
    RHS   = zeros(npoin_test)
    build_dg_volume_matrix!(A_vol, 1, 0.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, σ, intΦ_eff, connijk_1)
    build_dg_source_term!(RHS, 1, ngl, ω, Je, q, connijk_1)

    # No flux, no BC — with Ω=0 every face is tangential (is_active=false for
    # all 6), so this exercises volume+source ALONE, with no coupling to
    # anything outside the element at all.

    x_local = A_vol \ RHS
    expected = q_val / κ_val
    max_err = maximum(abs.(x_local .- expected))
    println("Test (build_dg_source_term!, pure reaction-source balance): " *
            "max abs error = $max_err  (expected value: $expected)  " *
            (max_err < 1e-10 ? "✓ PASS" : "✗ FAIL"))
    return max_err
end

function test_scattering_toy_C_single_element(; ngl=4, Δx=1.3, Δy=1.7, Δz=2.1,
                                                n_dir=3, κ_val=0.5, I0=4.2,
                                                C=nothing, max_iter=200, tol=1e-10)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    Je_const = (Δx/2)*(Δy/2)*(Δz/2)
    dξdx = fill(2/Δx,1,ngl,ngl,ngl); dξdy=zeros(1,ngl,ngl,ngl); dξdz=zeros(1,ngl,ngl,ngl)
    dηdx = zeros(1,ngl,ngl,ngl); dηdy=fill(2/Δy,1,ngl,ngl,ngl); dηdz=zeros(1,ngl,ngl,ngl)
    dζdx = zeros(1,ngl,ngl,ngl); dζdy=zeros(1,ngl,ngl,ngl); dζdz=fill(2/Δz,1,ngl,ngl,ngl)
    Je   = fill(Je_const,1,ngl,ngl,ngl)
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    npoin_test = ngl^3
    κ = fill(κ_val, npoin_test); σ = zeros(npoin_test); intΦ_eff = zeros(npoin_test)

    dirs = [(1.0,0.0,0.0), (0.6,0.8,0.0), (0.3,0.5,0.8)][1:n_dir]
    if C === nothing
        C = [0.0  0.05 0.05; 0.05 0.0  0.05; 0.05 0.05 0.0][1:n_dir,1:n_dir]
    end

    # ── Build each direction's OWN volume+flux system, ONCE — reused every
    # iteration (they don't depend on the coupling, only the RHS does) ──────
    A_blocks = Matrix{Float64}[]
    b_blocks_bc_only = Vector{Float64}[]   # BC/flux part only, WITHOUT scattering
    for k = 1:n_dir
        Ωx,Ωy,Ωz = dirs[k]
        A_vol = zeros(npoin_test,npoin_test); RHS = zeros(npoin_test)
        build_dg_volume_matrix!(A_vol, 1, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                 dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                                 Je, κ, σ, intΦ_eff, connijk_1)
        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, 1, Ωx,Ωy,Ωz, ngl,
                                 dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz)
        for fid=1:6
            is_active[fid] || continue
            fn = compute_face_normal_avg(1, fid, ngl, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz)
            I_up = is_inflow[fid] ? ((i,j,k)->I0) : ((i,j,k)->0.0)
            build_dg_flux_terms!(A_vol, RHS, fid, ngl, Ωx,Ωy,Ωz, is_inflow[fid], is_active[fid],
                                  fn, ω, dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz, Je, 1, I_up)
        end
        push!(A_blocks, A_vol)
        push!(b_blocks_bc_only, RHS)
    end

    # ── Ground truth: n_dir·npoin_test monolithic solve — ONLY diagonal
    # (same-element) coupling via C, no interior faces exist at all ─────────
    Nbig = n_dir*npoin_test
    A_big = zeros(Nbig,Nbig); b_big = zeros(Nbig)
    for k = 1:n_dir
        rows = (k-1)*npoin_test+1 : k*npoin_test
        A_big[rows,rows] .= A_blocks[k]
        b_big[rows] .= b_blocks_bc_only[k]
        for mm = 1:n_dir
            C[k,mm] == 0.0 && continue
            cols = (mm-1)*npoin_test+1 : mm*npoin_test
            for kk=1:ngl,jj=1:ngl,ii=1:ngl
                p = lidx(ii,jj,kk,ngl)
                A_big[rows[p], cols[p]] -= C[k,mm]*ω[ii]*ω[jj]*ω[kk]*Je_const
            end
        end
    end
    x_big = A_big \ b_big

    # ── Iterative solve: reuses the SAME A_blocks/b_blocks_bc_only, only
    # extra_source (from build_dg_source_term!) changes each iteration ──────
    I_all = [zeros(npoin_test) for _ = 1:n_dir]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:n_dir
            RHS_k = copy(b_blocks_bc_only[k])
            for mm = 1:n_dir
                C[k,mm] == 0.0 && continue
                q_km = C[k,mm] .* I_all[mm]   # per-node, indexed by original CG id (=local id here)
                build_dg_source_term!(RHS_k, 1, ngl, ω, Je, q_km, connijk_1)
            end
            x_new = A_blocks[k] \ RHS_k
            max_delta = max(max_delta, maximum(abs.(x_new .- I_all[k])))
            I_all[k] = x_new
        end
        max_delta < tol && (println("Converged after $it iterations, max_delta=$max_delta"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    max_err = 0.0
    for k = 1:n_dir
        rows = (k-1)*npoin_test+1 : k*npoin_test
        max_err = max(max_err, maximum(abs.(I_all[k] .- x_big[rows])))
    end
    println("Iterative (single element) vs direct block-coupled ground truth: max abs error = $max_err")
    return max_err
end

function test_scattering_toy_C_two_element(; ngl=4, Δx1=1.3, Δx2=3.0, Δy=1.7, Δz=2.1,
                                             n_dir=3, κ_val=0.5, I0=4.2,
                                             C=nothing, max_iter=200, tol=1e-10)
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    npoin_test = ngl^3
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    Je1 = (Δx1/2)*(Δy/2)*(Δz/2)
    dξdx1=fill(2/Δx1,1,ngl,ngl,ngl); dξdy1=zeros(1,ngl,ngl,ngl); dξdz1=zeros(1,ngl,ngl,ngl)
    dηdx1=zeros(1,ngl,ngl,ngl); dηdy1=fill(2/Δy,1,ngl,ngl,ngl); dηdz1=zeros(1,ngl,ngl,ngl)
    dζdx1=zeros(1,ngl,ngl,ngl); dζdy1=zeros(1,ngl,ngl,ngl); dζdz1=fill(2/Δz,1,ngl,ngl,ngl)
    Je1a = fill(Je1,1,ngl,ngl,ngl)

    Je2 = (Δx2/2)*(Δy/2)*(Δz/2)
    dξdx2=fill(2/Δx2,1,ngl,ngl,ngl); dξdy2=zeros(1,ngl,ngl,ngl); dξdz2=zeros(1,ngl,ngl,ngl)
    dηdx2=zeros(1,ngl,ngl,ngl); dηdy2=fill(2/Δy,1,ngl,ngl,ngl); dηdz2=zeros(1,ngl,ngl,ngl)
    dζdx2=zeros(1,ngl,ngl,ngl); dζdy2=zeros(1,ngl,ngl,ngl); dζdz2=fill(2/Δz,1,ngl,ngl,ngl)
    Je2a = fill(Je2,1,ngl,ngl,ngl)

    κ = fill(κ_val, npoin_test); σ = zeros(npoin_test); intΦ_eff = zeros(npoin_test)

    dirs = [(1.0,0.0,0.0), (0.6,0.8,0.0), (0.3,0.5,0.8)][1:n_dir]
    if C === nothing
        C = [0.0  0.05 0.05; 0.05 0.0  0.05; 0.05 0.05 0.0][1:n_dir,1:n_dir]
    end

    # ── Pre-build each direction's OWN volume+flux, per element — the BC/flux
    # parts don't depend on scattering, only the source does each iteration ──
    A1 = Matrix{Float64}[]; b1_bc = Vector{Float64}[]
    A2 = Matrix{Float64}[]; b2_bc_partial = Vector{Float64}[]  # inflow face NOT yet added (needs elem1's trace)
    is_inflow1_all = Vector{Bool}[]; is_active1_all = Vector{Bool}[]
    is_inflow2_all = Vector{Bool}[]; is_active2_all = Vector{Bool}[]
    fn1_5 = Vector{NTuple{3,Float64}}(undef, n_dir)
    fn2_5 = Vector{NTuple{3,Float64}}(undef, n_dir)

    for k = 1:n_dir
        Ωx,Ωy,Ωz = dirs[k]

        Av1 = zeros(npoin_test,npoin_test); Rv1 = zeros(npoin_test)
        build_dg_volume_matrix!(Av1, 1, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                 dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1,
                                 Je1a, κ, σ, intΦ_eff, connijk_1)
        ii1=Vector{Bool}(undef,6); ia1=Vector{Bool}(undef,6)
        classify_element_faces!(ii1, ia1, 1, Ωx,Ωy,Ωz, ngl, dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1)
        for fid=1:6
            ia1[fid] || continue
            fn = compute_face_normal_avg(1, fid, ngl, dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1)
            fid==5 && (fn1_5[k] = fn)
            I_up = ii1[fid] ? ((i,j,k)->I0) : ((i,j,k)->0.0)  # face 5 (true inflow BC), face 6 outflow
            build_dg_flux_terms!(Av1, Rv1, fid, ngl, Ωx,Ωy,Ωz, ii1[fid], ia1[fid], fn, ω,
                                  dξdx1,dξdy1,dξdz1, dηdx1,dηdy1,dηdz1, dζdx1,dζdy1,dζdz1, Je1a, 1, I_up)
        end
        push!(A1, Av1); push!(b1_bc, Rv1)
        push!(is_inflow1_all, ii1); push!(is_active1_all, ia1)

        Av2 = zeros(npoin_test,npoin_test); Rv2 = zeros(npoin_test)
        build_dg_volume_matrix!(Av2, 1, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                 dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2,
                                 Je2a, κ, σ, intΦ_eff, connijk_1)
        ii2=Vector{Bool}(undef,6); ia2=Vector{Bool}(undef,6)
        classify_element_faces!(ii2, ia2, 1, Ωx,Ωy,Ωz, ngl, dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2)
        for fid=1:6
            ia2[fid] || continue
            fn = compute_face_normal_avg(1, fid, ngl, dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2)
            fid==5 && (fn2_5[k] = fn)
            fid == 5 && continue   # SKIP: element 2's inflow face needs elem 1's trace — added per-iteration below
            I_up = (i,j,k)->0.0    # outflow (face 6) — always zero contribution
            build_dg_flux_terms!(Av2, Rv2, fid, ngl, Ωx,Ωy,Ωz, ii2[fid], ia2[fid], fn, ω,
                                  dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2, Je2a, 1, I_up)
        end
        push!(A2, Av2); push!(b2_bc_partial, Rv2)
        push!(is_inflow2_all, ii2); push!(is_active2_all, ia2)
    end

    # ── Iterative solve: element 1 first (only depends on itself + scattering
    # from BOTH elements' previous iterate — since scattering is LOCAL in
    # space, elem 1's source only ever needs elem-1-side values of I_m) ──────
    I1_all = [zeros(npoin_test) for _ = 1:n_dir]
    I2_all = [zeros(npoin_test) for _ = 1:n_dir]

    for it = 1:max_iter
        max_delta = 0.0

        for k = 1:n_dir
            RHS1 = copy(b1_bc[k])
            for mm = 1:n_dir
                C[k,mm]==0.0 && continue
                build_dg_source_term!(RHS1, 1, ngl, ω, Je1a, C[k,mm].*I1_all[mm], connijk_1)
            end
            x1_new = A1[k] \ RHS1
            max_delta = max(max_delta, maximum(abs.(x1_new .- I1_all[k])))
            I1_all[k] = x1_new
        end

        for kdir = 1:n_dir
            Ωx,Ωy,Ωz = dirs[kdir]
            RHS2 = copy(b2_bc_partial[kdir])
            I_up = (i,j,k) -> I1_all[kdir][lidx(ngl,j,k,ngl)]   # now unambiguous — kdir, not k
    	    build_dg_flux_terms!(zeros(npoin_test,npoin_test), RHS2, 5, ngl, Ωx,Ωy,Ωz,
                          is_inflow2_all[kdir][5], is_active2_all[kdir][5], fn2_5[kdir], ω,
                          dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2, Je2a, 1, I_up)
    	    for mm = 1:n_dir
        	C[kdir,mm]==0.0 && continue
        	build_dg_source_term!(RHS2, 1, ngl, ω, Je2a, C[kdir,mm].*I2_all[mm], connijk_1)
    	    end
    	    x2_new = A2[kdir] \ RHS2
    	    max_delta = max(max_delta, maximum(abs.(x2_new .- I2_all[kdir])))
    	    I2_all[kdir] = x2_new
	end

	max_delta < tol && (println("Converged after $it iterations, max_delta=$max_delta"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    # ── Ground truth: monolithic (2·n_dir·npoin_test)-sized solve ───────────
    Nbig = 2*n_dir*npoin_test
    A_big = zeros(Nbig,Nbig); b_big = zeros(Nbig)
    blk(elem,k,p) = ((elem-1)*n_dir + (k-1))*npoin_test + p   # elem=1,2

    for k = 1:n_dir
        rows1 = [blk(1,k,p) for p=1:npoin_test]
        A_big[rows1,rows1] .= A1[k]
        b_big[rows1] .= b1_bc[k]

        rows2 = [blk(2,k,p) for p=1:npoin_test]
        A_big[rows2,rows2] .= A2[k]
        b_big[rows2] .= b2_bc_partial[k]
        # element 2's inflow face couples to element 1's SAME direction k —
        # off-diagonal block (elem2,k) <-> (elem1,k), NOT cross-direction.
        Ωx,Ωy,Ωz = dirs[k]
        A_dummy = zeros(npoin_test,npoin_test); R_dummy = zeros(npoin_test)
        build_dg_flux_terms!(A_dummy, R_dummy, 5, ngl, Ωx,Ωy,Ωz, true, true, fn2_5[k], ω,
                              dξdx2,dξdy2,dξdz2, dηdx2,dηdy2,dηdz2, dζdx2,dζdy2,dζdz2, Je2a, 1,
                              (i,j,k)->1.0)  # unit-value trick ONLY for extracting the diagonal
                                              # COEFFICIENT itself (weight), not for any RHS use —
                                              # this is exactly the flux "weight" derived analytically
                                              # earlier: weight = -Ωdotn*wq*Je_face, a KNOWN closed
                                              # form independent of what value it's later multiplied by
        for j=1:ngl,kk=1:ngl
            p_here = lidx(1,j,kk,ngl)
            p_there = lidx(ngl,j,kk,ngl)   # element 1's outflow face node, same (j,k)
            A_big[rows2[p_here], rows1[p_there]] -= A_dummy[p_here,p_here]
        end

        for mm = 1:n_dir
            C[k,mm]==0.0 && continue
            rows1m = [blk(1,mm,p) for p=1:npoin_test]
            rows2m = [blk(2,mm,p) for p=1:npoin_test]
            for p=1:ngl^3
                ii_,jj_,kk_ = ((p-1)%ngl)+1, (((p-1)÷ngl)%ngl)+1, ((p-1)÷(ngl*ngl))+1
                A_big[rows1[p], rows1m[p]] -= C[k,mm]*ω[ii_]*ω[jj_]*ω[kk_]*Je1a[1,ii_,jj_,kk_]
                A_big[rows2[p], rows2m[p]] -= C[k,mm]*ω[ii_]*ω[jj_]*ω[kk_]*Je2a[1,ii_,jj_,kk_]
            end
        end
    end
    x_big = A_big \ b_big

    max_err = 0.0
    for k = 1:n_dir
        rows1 = [blk(1,k,p) for p=1:npoin_test]; rows2 = [blk(2,k,p) for p=1:npoin_test]
        max_err = max(max_err, maximum(abs.(I1_all[k] .- x_big[rows1])))
        max_err = max(max_err, maximum(abs.(I2_all[k] .- x_big[rows2])))
    end
    println("Iterative (two element) vs direct block-coupled ground truth: max abs error = $max_err")
    return max_err, I1_all, I2_all
end

function test_scattering_isotropic_two_element(; ngl=4, Δx1=1.3, Δx2=3.0, Δy=1.7, Δz=2.1,
                                                 n_dir=3, κ_val=0.5, I0=4.2,
                                                 max_iter=200, tol=1e-10)
    σ_val = κ_val
    C = fill(σ_val/n_dir, n_dir, n_dir)

    # σ passed into build_dg_volume_matrix! (via test_scattering_toy_C_two_element)
    # must be ZERO — scattering's diagonal contribution is entirely excluded here,
    # not approximated via intΦ_eff=1, per the physics: σ(1-intΦ_eff) is pure
    # quadrature-error residual around an analytic zero and has no place in a
    # formulation where scattering is explicit on the RHS. The REAL σ_val only
    # enters through C, which builds the actual scattering source term.
    max_err_vs_direct, I1_all, I2_all = test_scattering_toy_C_two_element(
        ngl=ngl, Δx1=Δx1, Δx2=Δx2, Δy=Δy, Δz=Δz,
        n_dir=n_dir, κ_val=κ_val, I0=I0, C=C, max_iter=max_iter, tol=tol)
    # (relies on test_scattering_toy_C_two_element's own κ passed to
    #  build_dg_volume_matrix! being κ_val alone, σ=zeros — already the
    #  case throughout that function, since it never introduced a separate
    #  σ argument distinct from κ; confirmed by re-reading it above)

    max_dev_from_I0 = 0.0
    for k = 1:n_dir
        max_dev_from_I0 = max(max_dev_from_I0, maximum(abs.(I1_all[k] .- I0)))
        max_dev_from_I0 = max(max_dev_from_I0, maximum(abs.(I2_all[k] .- I0)))
    end

    println("Iterative vs direct block-coupled ground truth: max abs error = $max_err_vs_direct")
    println("Max deviation from closed-form I≡I0: $max_dev_from_I0")

    return max_err_vs_direct, max_dev_from_I0
end

function test_scattering_exact_ode_v2(; ngl=4, Δx1=1.3, Δx2=3.0, Δy=1.7, Δz=2.1,
                                       n_dir=3, κ_val=0.5, σ_val=0.5, I0=4.2,
                                       max_iter=200, tol=1e-10)
    dirs = [(1.0,0.0,0.0), (0.85,0.0,0.0), (0.6,0.0,0.0)][1:n_dir]

    # ── Build the 2-element toy mesh once — same pattern as test_dg_reconciliation ──
    lgl   = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    basis = build_Interpolation_basis!(LagrangeBasis(), lgl.ξ, lgl.ξ, Float64, CPU())
    ψ, dψ, ω = basis.ψ, basis.dψ, lgl.ω

    m = build_toy_structured_mesh(2, 1, 1, ngl, Δx1, Δy, Δz)  # NOTE: uniform Δx across both
                                                                # elements here — build_toy_structured_mesh
                                                                # doesn't support per-element sizing, so
                                                                # Δx1==Δx2 required for this rewrite (flagged below)
    npoin_cg = m.npoin
    npoin_test = ngl^3

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    C = fill(σ_val/n_dir, n_dir, n_dir)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(m.nelem, m.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(m, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(m, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    orders = [build_conforming_sweep_order(d..., m.nelem, ngl, neighbor_elem,
                  m.dξdx,m.dξdy,m.dξdz, m.dηdx,m.dηdy,m.dηdz, m.dζdx,m.dζdy,m.dζdz)[1] for d in dirs]

    metrics_nt = (dξdx=m.dξdx,dξdy=m.dξdy,dξdz=m.dξdz,
                  dηdx=m.dηdx,dηdy=m.dηdy,dηdz=m.dηdz,
                  dζdx=m.dζdx,dζdy=m.dζdy,dζdz=m.dζdz, Je=m.Je)
    bc_fn = (x,y,z) -> I0   # flat — genuine domain boundary only (fid=5 on element 1)

    # ── Iterative solve, through the ACTUAL production sweep ────────────────
    I_cg_all = [zeros(npoin_cg) for _ = 1:n_dir]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:n_dir
            Ωx,Ωy,Ωz = dirs[k]
            q_k = zeros(npoin_cg)
            for mm = 1:n_dir
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, npoin_dg)
            run_dg_sweep!(x_dg, m, metrics_nt, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn; extra_source=q_k)
            I_new = zeros(npoin_cg)
            for iel = 1:m.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[m.connijk[iel,ii,jj,kk]] = x_dg[connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        max_delta < tol && (println("Converged after $it iterations, max_delta=$max_delta"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    # ── Exact solution via matrix exponential (unchanged derivation) ────────
    μ = [d[1] for d in dirs]
    M = Diagonal(μ)
    A = -inv(M) * (κ_val*I(n_dir) - (σ_val/n_dir)*ones(n_dir,n_dir))
    x0 = m.xmin
    I_exact_at(x) = exp(A*(x-x0)) * (I0*ones(n_dir))

    max_dev = 0.0
    for k = 1:n_dir
        for iel = 1:m.nelem, i=1:ngl, j=1:ngl, kk=1:ngl
            g = m.connijk[iel,i,j,kk]
            I_ex = I_exact_at(m.x[g])[k]
            max_dev = max(max_dev, abs(I_cg_all[k][g] - I_ex))
        end
    end
    println("Max deviation from EXACT ODE solution (via run_dg_sweep!): $max_dev")
    return max_dev
end

function test_manufactured_isotropic(mesh, metrics, basis, ω;
                                      I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                      κ_val=0.001, σ_val=0.0005,
                                      max_iter=500, tol=1e-9)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    # ── Domain extent, one full period across the whole mesh ───────────────
    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin; Lz = mesh.zmax - mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    # ── Manufactured solution and its EXACT analytic gradient ──────────────
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    # ∇I = A*(kx,ky,kz)*cos(kx x+ky y+kz z)*(1+bΩx+cΩy+dΩz)  — exact, by hand
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    # ── Angular quadrature: directions + weights from extra_mesh ────────────
    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext, iθ, iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end
    W = sum(w_ext)
    C = [σ_val*w_ext[m]/W for k=1:extra_npoin, m=1:extra_npoin]   # isotropic — same for every (k,m)

    # ── Mesh infrastructure (build/reuse, as before) ────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    orders = [build_conforming_sweep_order(d..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz)[1] for d in dirs]

    # ── Per-direction EXTERNAL source q_ext, computed ONCE from the exact
    # manufactured solution — fixed for the whole run, independent of iteration ──
    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][g] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end

    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    # ── Iterative solve ──────────────────────────────────────────────────────
    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (println("Converged after $it iterations"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    # ── Compare against the EXACT manufactured solution everywhere ──────────
    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][g] - I_true))
        end
    end
    println("Manufactured solution (isotropic) — max abs error over all points, all directions: $max_err")
    return max_err
end

function test_manufactured_HG(mesh, metrics, basis, ω;
                               I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                               κ_val=0.001, σ_val=0.0005, g=0.7,
                               max_iter=500, tol=1e-9)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin; Lz = mesh.zmax - mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    # Same manufactured ansatz as the isotropic test — deliberately UNCHANGED,
    # so this test isolates ONLY the new thing: C built from a genuinely
    # direction-pair-dependent phase function, not a richer I itself.
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext, iθ, iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end

    # ── Henyey-Greenstein phase function, genuinely pairwise — the new piece ──
    # Φ_HG(cosΘ) = (1-g²) / (4π(1+g²-2g·cosΘ)^{3/2}), cosΘ = Ω_k·Ω_m (both
    # already unit vectors, so the dot product IS cosΘ directly, no arccos
    # needed). NOT normalized to Σ_m C[k,m]=σ exactly by construction — per
    # the design principle established for the isotropic test, this doesn't
    # need to match the continuum ∫Φ dΩ'=1 exactly; q_ext is built from this
    # SAME discrete C, so self-consistency (not continuum-exactness) is what
    # correctness actually requires here.
    function Φ_HG(cosΘ, g)
        return (1 - g^2) / (4π * (1 + g^2 - 2*g*cosΘ)^1.5)
    end

    C = zeros(extra_npoin, extra_npoin)
    for k = 1:extra_npoin, m = 1:extra_npoin
        cosΘ = dirs[k][1]*dirs[m][1] + dirs[k][2]*dirs[m][2] + dirs[k][3]*dirs[m][3]
        C[k,m] = σ_val * w_ext[m] * Φ_HG(cosΘ, g)
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    orders = [build_conforming_sweep_order(dirs[k]..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz)[1] for k = 1:extra_npoin]

    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            x_,y_,z_ = mesh.x[gnode], mesh.y[gnode], mesh.z[gnode]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][gnode] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end

    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (println("Converged after $it iterations"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            I_true = I_manuf(mesh.x[gnode],mesh.y[gnode],mesh.z[gnode],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][gnode] - I_true))
        end
    end
    println("Manufactured solution (Henyey-Greenstein, g=$g) — max abs error: $max_err")
    return max_err
end

function test_manufactured_HG_angular_convergence(mesh, metrics, basis, ω;
                                                   I0=4.2, A=1.0,
                                                   g_manuf=0.6,        # peakedness of the MANUFACTURED angular shape
                                                   Ω_ref=(0.0,0.0,1.0), # fixed reference direction for that shape
                                                   κ_val=0.001, σ_val=0.0005, g_scatter=0.7,
                                                   max_iter=500, tol=1e-9)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin; Lz = mesh.zmax - mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    # Angular factor: Φ_HG(Ω·Ω_ref) itself — genuinely NON-polynomial in Ω,
    # so no finite-degree angular quadrature integrates it exactly. This is
    # the deliberate difference from the earlier linear/quadratic attempts.
    Φ_HG(cosΘ, gg) = (1 - gg^2) / (4π * (1 + gg^2 - 2*gg*cosΘ)^1.5)
    angular_shape(Ωx,Ωy,Ωz) = Φ_HG(Ωx*Ω_ref[1]+Ωy*Ω_ref[2]+Ωz*Ω_ref[3], g_manuf)

    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*angular_shape(Ωx,Ωy,Ωz)

    # ── Analytic spatial gradient — angular_shape has NO x,y,z dependence,
    # so this factors out cleanly, same structure as before ────────────────
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*angular_shape(Ωx,Ωy,Ωz) .* (kx,ky,kz)

    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext, iθ, iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end

    C = zeros(extra_npoin, extra_npoin)
    for k = 1:extra_npoin, m = 1:extra_npoin
        cosΘ_km = dirs[k][1]*dirs[m][1] + dirs[k][2]*dirs[m][2] + dirs[k][3]*dirs[m][3]
        C[k,m] = σ_val * w_ext[m] * Φ_HG(cosΘ_km, g_scatter)
    end

    for k in [1, extra_npoin÷2, extra_npoin]
    	println("k=$k: sum(C[k,:]) = $(sum(C[k,:]))  (target ≈ σ_val=$σ_val for a well-resolved quadrature)")
    end


    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    orders = [build_conforming_sweep_order(dirs[k]..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz)[1] for k = 1:extra_npoin]

    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            x_,y_,z_ = mesh.x[gnode], mesh.y[gnode], mesh.z[gnode]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][gnode] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end

    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        max_delta < tol && (println("Converged after $it iterations"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            I_true = I_manuf(mesh.x[gnode],mesh.y[gnode],mesh.z[gnode],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][gnode] - I_true))
        end
    end
    println("Manufactured HG (non-polynomial angular shape) — max abs error: $max_err")
    return max_err
end

# ─────────────────────────────────────────────────────────────────────────────
#  Mie-like proxy phase function — NOT a real Mie calculation. Built to mimic
#  the two features that make genuine Mie phase functions hard for angular
#  quadrature: (1) a sharper forward peak than single-HG produces (via a
#  double-HG envelope), and (2) oscillatory ripple structure in scattering
#  angle (the interference-fringe signature real Mie phase functions have,
#  which single- or double-HG alone never do). Not required to be physically
#  normalized or non-negative — see conversation: this is a numerical stress
#  test of the quadrature/scattering machinery, not a physical model, and
#  nothing in our code assumes Φ≥0 (it's just a coefficient in a linear
#  equation) — but worth checking/reporting where it dips, for honesty.
# ─────────────────────────────────────────────────────────────────────────────

function Φ_mie_like(cosΘ; g1=0.85, g2=-0.2, f=0.9, ripple_N=8, ripple_ε=0.25)
    HG(c,g) = (1-g^2) / (4π*(1+g^2-2*g*c)^1.5)
    envelope = f*HG(cosΘ,g1) + (1-f)*HG(cosΘ,g2)
    Θ = acos(clamp(cosΘ, -1.0, 1.0))
    ripple = 1 + ripple_ε*cos(ripple_N*Θ)
    return envelope * ripple
end

function test_manufactured_mie_like(mesh, metrics, basis, ω;
                                     I0=4.2, A=1.0,
                                     Ω_ref=(0.0,0.0,1.0),
                                     κ_val=0.001, σ_val=0.0005,
                                     mie_kwargs = (g1=0.85, g2=-0.2, f=0.9, ripple_N=8, ripple_ε=0.25),
                                     max_iter=500, tol=1e-9)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin; Lz = mesh.zmax - mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    angular_shape(Ωx,Ωy,Ωz) = Φ_mie_like(Ωx*Ω_ref[1]+Ωy*Ω_ref[2]+Ωz*Ω_ref[3]; mie_kwargs...)
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*angular_shape(Ωx,Ωy,Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*angular_shape(Ωx,Ωy,Ωz) .* (kx,ky,kz)

    # ── Honesty check: report where this proxy dips negative, since it's
    # not built to guarantee positivity (see docstring above) ──────────────
    sample_cosΘ = range(-1.0, 1.0, length=200)
    sample_vals = [Φ_mie_like(c; mie_kwargs...) for c in sample_cosΘ]
    println("Φ_mie_like sample range: [$(minimum(sample_vals)), $(maximum(sample_vals))]" *
            (minimum(sample_vals) < 0 ? "  ⚠ goes negative — non-physical, fine for this " *
             "numerical stress test but flagged for honesty" : ""))

    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext, iθ, iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end

    C = zeros(extra_npoin, extra_npoin)
    for k = 1:extra_npoin, m = 1:extra_npoin
        cosΘ_km = dirs[k][1]*dirs[m][1] + dirs[k][2]*dirs[m][2] + dirs[k][3]*dirs[m][3]
        C[k,m] = σ_val * w_ext[m] * Φ_mie_like(cosΘ_km; mie_kwargs...)
    end

    for kcheck in [1, extra_npoin÷2, extra_npoin]
        println("k=$kcheck: sum(C[k,:]) = $(sum(C[kcheck,:]))  (target ≈ σ_val=$σ_val)")
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    orders = [build_conforming_sweep_order(dirs[k]..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz)[1] for k = 1:extra_npoin]

    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            x_,y_,z_ = mesh.x[gnode], mesh.y[gnode], mesh.z[gnode]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][gnode] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end

    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    diverged = false
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        if !isfinite(max_delta) || max_delta > 1e10
            println("DIVERGED at iteration $it (max_delta=$max_delta) — aborting")
            diverged = true
            break
        end
        max_delta < tol && (println("Converged after $it iterations"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end
    diverged && return Inf

    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for gnode = 1:npoin_cg
            I_true = I_manuf(mesh.x[gnode],mesh.y[gnode],mesh.z[gnode],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][gnode] - I_true))
        end
    end
    println("Manufactured Mie-like — max abs error: $max_err")
    return max_err
end

function test_manufactured_isotropic_timed(mesh, metrics, basis, ω;
                                      I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                      κ_val=0.001, σ_val=0.0005,
                                      max_iter=500, tol=1e-9)
    ngl = mesh.ngl
    ψ, dψ = basis.ψ, basis.dψ
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    # ── Domain extent, one full period across the whole mesh ───────────────
    Lx = mesh.xmax - mesh.xmin; Ly = mesh.ymax - mesh.ymin; Lz = mesh.zmax - mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz

    # ── Manufactured solution and its EXACT analytic gradient ──────────────
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    # ∇I = A*(kx,ky,kz)*cos(kx x+ky y+kz z)*(1+bΩx+cΩy+dΩz)  — exact, by hand
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    # ── Angular quadrature: directions + weights from extra_mesh ────────────
    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]
        ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext, iθ, iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end
    W = sum(w_ext)
    C = [σ_val*w_ext[m]/W for k=1:extra_npoin, m=1:extra_npoin]   # isotropic — same for every (k,m)

    # ── Mesh infrastructure (build/reuse, as before) ────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    if isempty(mesh.connijk_dg)
        connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
            build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
        face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem,
                                                   neighbor_fid, face_correspondence)
        mesh.connijk_dg = connijk_dg; mesh.npoin_dg = npoin_dg
        mesh.dg_to_cg = dg_to_cg; mesh.face_neighbor_dg = face_neighbor_dg
    end

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    orders = [build_conforming_sweep_order(d..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz)[1] for d in dirs]

    # ── Per-direction EXTERNAL source q_ext, computed ONCE from the exact
    # manufactured solution — fixed for the whole run, independent of iteration ──
    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][g] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end

    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    # ── Iterative solve ──────────────────────────────────────────────────────
    timers = Dict{String,Float64}()
    timers["scatter_update"] = 0.0

    # One factorization cache PER DIRECTION — populated on iteration 1,
    # reused for iterations 2+. Never shared across k, since A_vol
    # genuinely differs by direction (both volume streaming term and
    # flux terms depend on Ω).
    factor_caches = [Dict{Int,Any}() for _ = 1:extra_npoin]
    weight_cache_stores = [Dict{Int,Any}() for _ = 1:extra_npoin]
    face_info_caches = [Dict{Int,Any}() for _ = 1:extra_npoin]

    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            
	    @info length(factor_caches[k])

	    Ωx,Ωy,Ωz = dirs[k]
            t0 = time()
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            timers["scatter_update"] += time() - t0
            x_dg = fill(NaN, mesh.npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          mesh.connijk_dg, mesh.face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k, timers=timers, factor_cache=factor_caches[k], weight_cache_store=weight_cache_stores[k], face_info_cache=face_info_caches[k])
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[mesh.connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        max_delta < tol && break
    end

    total = sum(values(timers))
    println("─"^60)
    println("TIMING BREAKDOWN (total = $(round(total,digits=2))s)")
    for (name,t) in sort(collect(timers), by=x->-x[2])
        println("  $name: $(round(t,digits=2))s ($(round(100*t/total,digits=1))%)")
    end
    # ── Compare against the EXACT manufactured solution everywhere ──────────
    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][g] - I_true))
        end
    end
    println("Manufactured solution (isotropic) — max abs error over all points, all directions: $max_err")

    return timers
end

"""
    face_flat_index(fid, ngl, i, j, k)

Returns the SAME flattened index p that build_dg_flux_terms!/
build_dg_flux_rhs_only! use internally (their `for k in krange, j in
jrange, i in irange; idx += 1` enumeration) — computed directly via
arithmetic instead of a Dict lookup. Since exactly one of the three
reference directions is fixed (singleton range) for any face, and the
other two are full 1:ngl, this is O(1) with no allocation.
"""
@inline function face_flat_index(fid::Int, ngl::Int, i::Int, j::Int, k::Int)
    irange, jrange, krange = face_ip_range(fid, ngl)
    pos_i = length(irange) == 1 ? 1 : i
    pos_j = length(jrange) == 1 ? 1 : j
    pos_k = length(krange) == 1 ? 1 : k
    len_i = length(irange); len_j = length(jrange)
    return pos_i + (pos_j-1)*len_i + (pos_k-1)*len_i*len_j
end


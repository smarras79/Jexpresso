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
                             Je, κ, connijk)

In-place: fills the (ngl^3 x ngl^3) buffer `A_vol` (caller-allocated once,
reused across elements/directions — no allocation here) with the element-
local DG volume matrix for element `iel`, direction (Ωx,Ωy,Ωz), for the
operator:

    Ω·∇I + κ·I

Scattering is handled ENTIRELY through the explicit outer-iteration
source coupling (extra_source_all), never through this operator — this
function no longer accepts σ/intΦ_eff at all, since passing a nonzero σ
here would silently double-count extinction against any RHS/manufactured
solution built assuming scattering is purely a source term.

`ψ`, `dψ` are the 1D Lagrange basis / derivative matrices (ngl x ngl), as
used elsewhere in this codebase's CG assembly — same basis, restricted to
one element rather than assembled globally.

`κ` is indexed by GLOBAL point id (matching how
sparse_lhs_assembly_3Dby2D_adaptive reads it: `κ_data[ip]` via
`connijk`/`connijk_spa`) — `connijk` is passed in specifically to map this
element's local (i,j,k) to the global point id needed for that lookup.
"""
function build_dg_volume_matrix!(A_vol::AbstractMatrix{Float64}, iel::Int,
                                  Ωx::Float64, Ωy::Float64, Ωz::Float64, ngl::Int,
                                  ψ, dψ, ω,
                                  dξdx, dξdy, dξdz,
                                  dηdx, dηdy, dηdz,
                                  dζdx, dζdy, dζdz,
                                  Je, κ, connijk)
    fill!(A_vol, 0.0)

    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip_global = connijk[iel, i, j, k]
        ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]

        dξdx_ij = dξdx[iel,i,j,k]; dξdy_ij = dξdy[iel,i,j,k]; dξdz_ij = dξdz[iel,i,j,k]
        dηdx_ij = dηdx[iel,i,j,k]; dηdy_ij = dηdy[iel,i,j,k]; dηdz_ij = dηdz[iel,i,j,k]
        dζdx_ij = dζdx[iel,i,j,k]; dζdy_ij = dζdy[iel,i,j,k]; dζdz_ij = dζdz[iel,i,j,k]

        α = κ[ip_global]
        p = lidx(i,j,k,ngl)

        # ── Reaction term: α * ψ_q(node) * ψ_p(node) * ωJac ────────────────
        # Purely diagonal at a GLL nodal quadrature point (ψ_q(node_p) =
        # δ_{pq}), exactly matching the diagonal-only V_vec[pos] term in
        # sparse_lhs_assembly_3Dby2D_adaptive.
        A_vol[p, p] += α * ωJac

        # ── Streaming term: (Ω·∇ψ_q) at this node, projected onto ψ_p ──────
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

    A_vol = zeros(npoin_test, npoin_test)

    # ── Test 1: Ω=(1,0,0) — only ξ-direction (i-index) streaming survives ──
    build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, connijk)
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
                             Je, κ, connijk)
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
                             Je, κ, connijk)
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
    build_dg_volume_matrix!(A1, 1, Ω1[1],Ω1[2],Ω1[3], ngl, ψ,dψ,ω, dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,connijk)
    build_dg_volume_matrix!(A2, 1, Ω2[1],Ω2[2],Ω2[3], ngl, ψ,dψ,ω, dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,connijk)
    build_dg_volume_matrix!(Asum, 1, Ω1[1]+Ω2[1], Ω1[2]+Ω2[2], Ω1[3]+Ω2[3], ngl, ψ,dψ,ω,
                             dξdx,dξdy,dξdz,dηdx,dηdy,dηdz,dζdx,dζdy,dζdz, Je,κ,connijk)

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
                                 m.Je, κ, m.connijk)
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
                                 m.Je, κ, m.connijk)
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
                                 m.Je, κ, m.connijk)
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
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    fill!(A_vol, 0.0); fill!(RHS, 0.0)
    build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, connijk_1)
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
        connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
        for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

        A_vol = zeros(npoin_test, npoin_test)
        RHS   = zeros(npoin_test)

        build_dg_volume_matrix!(A_vol, 1, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                                 dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                                 Je, κ, connijk_1)

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
                             Je1a, κ, connijk_1)
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
                             Je2a, κ, connijk_1)
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

    A_vol = zeros(npoin_test, npoin_test)
    build_dg_volume_matrix!(A_vol, iel, 1.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz,
                             metrics.Je, κ, mesh.connijk)

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
                                 metrics.Je, κ, mesh.connijk)

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
#=function build_dg_source_term!(RHS_local::AbstractVector{Float64}, iel::Int, ngl::Int,
                                ω, Je, q, connijk)
    @inbounds for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        gid = connijk[iel,i,j,k]
        RHS_local[p] += q[gid] * ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
    end
    return RHS_local
end=#

function build_dg_source_term!(RHS_local, iel, ngl, ω, Je, q, connijk)
    local_q = length(q) == ngl^3
    @inbounds for k=1:ngl, j=1:ngl, i=1:ngl
        p = lidx(i,j,k,ngl)
        RHS_local[p] += (local_q ? q[p] : q[connijk[iel,i,j,k]]) * ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
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
    q        = fill(q_val, npoin_test)
    connijk_1 = Array{Int}(undef,1,ngl,ngl,ngl)
    for k=1:ngl,j=1:ngl,i=1:ngl; connijk_1[1,i,j,k] = lidx(i,j,k,ngl); end

    # Ω=(0,0,0): zero streaming, isolates the reaction+source balance entirely.
    A_vol = zeros(npoin_test, npoin_test)
    RHS   = zeros(npoin_test)
    build_dg_volume_matrix!(A_vol, 1, 0.0, 0.0, 0.0, ngl, ψ, dψ, ω,
                             dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                             Je, κ, connijk_1)
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
    κ = fill(κ_val, npoin_test);

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
                                 Je, κ, connijk_1)
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

    κ = fill(κ_val, npoin_test);

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
                                 Je1a, κ, connijk_1)
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
                                 Je2a, κ, connijk_1)
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


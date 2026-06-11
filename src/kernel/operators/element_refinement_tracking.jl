# Element-level refinement state tracking for the spatial-angular AMR solver.
#
# Invariant enforced here:
#   A pair of adjacent elements may be non-conforming in AT MOST ONE dimension
#   (spatial or angular, never both at the same interface).
#
# Consequence: angular refinement of element `iel` is only safe when EVERY
# spatial neighbor that shares a FACE or EDGE with `iel` (local and cross-rank)
# is at the same spatial refinement level.  Corner-only contacts are excluded
# because the angular non-conformity handler does not need to communicate across
# corners.  Face and edge contacts both require matching spatial nodes, so a
# spatially non-conforming edge is just as problematic as a non-conforming face.

"""
Refinement state snapshot for one spatial element, including the refinement
levels of every known spatial neighbor (face- and edge-adjacent only; corners
excluded).

Fields
------
- `spatial_level`   : `mesh.ad_lvl[iel]`
- `max_angular_level`: `maximum(extra_meshes_ref_level[iel])`
- `local_neighbor_elems`          : local (same-rank) face+edge neighbour element indices
- `local_neighbor_spatial_levels` : `mesh.ad_lvl` for each local neighbour
- `local_neighbor_max_ang_levels` : `maximum(ref_level)` for each local neighbour
- `cross_rank_neighbor_spatial_levels` : spatial levels of face+edge cross-rank neighbours,
  discovered via AllGather bbox exchange
- `cross_rank_neighbor_max_ang_levels` : angular levels of those same neighbours
- `can_refine_angular` : true iff the element is not at a NCF face AND no face/edge
  neighbour (local or cross-rank) has a different spatial refinement level
"""
struct ElementRefinementRecord
    spatial_level::Int
    max_angular_level::Int

    local_neighbor_elems         ::Vector{Int}
    local_neighbor_spatial_levels::Vector{Int}
    local_neighbor_max_ang_levels::Vector{Int}

    cross_rank_neighbor_spatial_levels::Vector{Int}
    cross_rank_neighbor_max_ang_levels::Vector{Int}

    can_refine_angular::Bool
end

# ─── internal helpers ─────────────────────────────────────────────────────────

function _element_bbox(iel, connijk, x, y, z, ngl)
    xlo = xhi = x[connijk[iel,1,1,1]]
    ylo = yhi = y[connijk[iel,1,1,1]]
    zlo = zhi = z[connijk[iel,1,1,1]]
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = connijk[iel,i,j,k]
        xlo = min(xlo, x[ip]); xhi = max(xhi, x[ip])
        ylo = min(ylo, y[ip]); yhi = max(yhi, y[ip])
        zlo = min(zlo, z[ip]); zhi = max(zhi, z[ip])
    end
    return xlo, xhi, ylo, yhi, zlo, zhi
end

# Contact classification for one coordinate dimension.
#   2 → positive-length overlap (interior shared region)
#   1 → touching (degenerate, zero-length overlap = shared face/edge/corner)
#   0 → no contact (gap)
function _contact_dim(lo_a, hi_a, lo_b, hi_b)
    tol = 1e-5
    lo = max(lo_a, lo_b)
    hi = min(hi_a, hi_b)
    if hi > lo + tol
        return 2
    elseif hi >= lo - tol
        return 1
    else
        return 0
    end
end

# Classify the shared topology between two axis-aligned bounding boxes.
# Returns :face, :edge, :corner, or :none.
#
# This handles non-conforming (different-size) bounding boxes correctly because
# it uses interval overlap rather than exact bound matching.
#
# Face  : 2 dims overlap, 1 dim touches
# Edge  : 1 dim overlaps, 2 dims touch
# Corner: 0 dims overlap, 3 dims touch
function _bbox_topology(xl,xh,yl,yh,zl,zh, xil,xih,yil,yih,zil,zih)
    cx = _contact_dim(xl, xh, xil, xih)
    cy = _contact_dim(yl, yh, yil, yih)
    cz = _contact_dim(zl, zh, zil, zih)
    (cx == 0 || cy == 0 || cz == 0) && return :none
    n_overlap = (cx == 2) + (cy == 2) + (cz == 2)
    n_touch   = (cx == 1) + (cy == 1) + (cz == 1)
    if n_touch == 1 && n_overlap == 2; return :face
    elseif n_touch == 2 && n_overlap == 1; return :edge
    elseif n_touch == 3; return :corner
    else; return :none
    end
end

# ─── public API ───────────────────────────────────────────────────────────────

"""
    build_element_refinement_records(mesh, extra_meshes_ref_level, nelem, ngl, comm) -> Vector{ElementRefinementRecord}

`extra_meshes_ref_level` may be `nothing` when the angular mesh is uniform (no
angular AMR); all angular refinement levels are then treated as 0.

Build a per-element refinement record for every local spatial element.

Algorithm
---------
1. Compute local spatial levels (`mesh.ad_lvl`) and max angular levels.
2. Flag elements that appear in any `mesh.non_conforming_facets*` list as
   spatially non-conforming (authoritative from p4est); those immediately set
   `can_refine_angular = false`.
3. Compute bounding boxes for all local elements.
4. AllGather bounding boxes + levels from all ranks (8 Float64 per element).
5. For each local element, scan all global elements via `_bbox_topology`:
   - :face and :edge contacts → tracked as neighbours (corners excluded)
   - Same-rank matches → `local_neighbor_elems`
   - Cross-rank matches → `cross_rank_neighbor_*_levels`
6. `can_refine_angular` = `!spa_nc[iel]` AND no face/edge neighbour
   (local or cross-rank) has a different spatial level.
   (Differs from old code: edge neighbours with different spatial level now also
   block angular refinement, since the angular solver requires matching nodes on
   every shared edge.)
"""
function build_element_refinement_records(
    mesh, extra_meshes_ref_level, nelem::Int, ngl::Int, comm
)::Vector{ElementRefinementRecord}
    _ang_ref = extra_meshes_ref_level

    # ── 1. Self levels ────────────────────────────────────────────────────────
    spatial_levels = Vector{Int}(undef, nelem)
    max_ang_levels = Vector{Int}(undef, nelem)
    for iel = 1:nelem
        spatial_levels[iel] = Int(mesh.ad_lvl[iel])
        max_ang_levels[iel] = _ang_ref === nothing ? 0 :
                              (isempty(_ang_ref[iel]) ? 0 : maximum(Int, _ang_ref[iel]))
    end

    # ── 2. Spatially non-conforming element sets (authoritative from p4est) ───
    spa_nc = falses(nelem)
    if !ismissing(mesh.num_ncf) && mesh.num_ncf > 0
        for ncf in mesh.non_conforming_facets
            cid = Int(ncf[1]); pid = Int(ncf[2])
            1 <= cid <= nelem && (spa_nc[cid] = true)
            1 <= pid <= nelem && (spa_nc[pid] = true)
        end
    end
    if !ismissing(mesh.num_ncf_pg) && mesh.num_ncf_pg > 0
        for ncf in mesh.non_conforming_facets_parents_ghost
            cid = Int(ncf[1])
            1 <= cid <= nelem && (spa_nc[cid] = true)
        end
    end
    if !ismissing(mesh.num_ncf_cg) && mesh.num_ncf_cg > 0
        for ncf in mesh.non_conforming_facets_children_ghost
            pid = Int(ncf[1])
            1 <= pid <= nelem && (spa_nc[pid] = true)
        end
    end

    # ── 3. Bounding boxes for all local elements ──────────────────────────────
    bboxes = Vector{NTuple{6,Float64}}(undef, nelem)
    for iel = 1:nelem
        bboxes[iel] = _element_bbox(iel, mesh.connijk, mesh.x, mesh.y, mesh.z, ngl)
    end

    # ── 4. AllGather bounding boxes + levels from all ranks ───────────────────
    # Each element contributes 8 Float64: (xl,xh,yl,yh,zl,zh,spatial_lvl,ang_lvl)
    n_fields = 8
    local_flat = Vector{Float64}(undef, nelem * n_fields)
    for iel = 1:nelem
        xl,xh,yl,yh,zl,zh = bboxes[iel]
        s = n_fields*(iel-1)
        local_flat[s+1] = xl;  local_flat[s+2] = xh
        local_flat[s+3] = yl;  local_flat[s+4] = yh
        local_flat[s+5] = zl;  local_flat[s+6] = zh
        local_flat[s+7] = Float64(spatial_levels[iel])
        local_flat[s+8] = Float64(max_ang_levels[iel])
    end
    n_loc  = Int32(nelem)
    n_all  = MPI.Allgather([n_loc], comm)
    all_flat = MPI.Allgatherv(local_flat, Int32.(n_all .* n_fields), comm)

    local_start = Int(sum(n_all[1:MPI.Comm_rank(comm)])) + 1
    local_end   = local_start + nelem - 1

    # ── 5. Find face+edge neighbours; track edge neighbours separately ─────────
    n_global = length(all_flat) ÷ n_fields
    local_nbrs      = [Int[] for _ = 1:nelem]
    local_edge_nbrs = [Int[] for _ = 1:nelem]   # subset: edge contacts only
    cr_spa_nbr      = [Int[] for _ = 1:nelem]
    cr_ang_nbr      = [Int[] for _ = 1:nelem]
    cr_edge_spa_nbr = [Int[] for _ = 1:nelem]   # cross-rank edge-contact spatial levels

    for iel = 1:nelem
        xl,xh,yl,yh,zl,zh = bboxes[iel]
        for k = 1:n_global
            s = n_fields*(k-1)
            xil = all_flat[s+1]; xih = all_flat[s+2]
            yil = all_flat[s+3]; yih = all_flat[s+4]
            zil = all_flat[s+5]; zih = all_flat[s+6]
            slvl = Int(round(all_flat[s+7]))
            alvl = Int(round(all_flat[s+8]))

            topo = _bbox_topology(xl,xh,yl,yh,zl,zh, xil,xih,yil,yih,zil,zih)
            # Corners excluded: angular solver does not communicate across corners
            (topo == :none || topo == :corner) && continue

            if local_start <= k <= local_end
                jel = k - local_start + 1
                jel == iel && continue
                push!(local_nbrs[iel], jel)
                topo == :edge && push!(local_edge_nbrs[iel], jel)
            else
                push!(cr_spa_nbr[iel], slvl)
                push!(cr_ang_nbr[iel], alvl)
                topo == :edge && push!(cr_edge_spa_nbr[iel], slvl)
            end
        end
    end

    # ── 6. Assemble records ───────────────────────────────────────────────────
    records = Vector{ElementRefinementRecord}(undef, nelem)
    for iel = 1:nelem
        lnbrs  = local_nbrs[iel]
        ln_spa = Int[spatial_levels[j] for j in lnbrs]
        ln_ang = Int[max_ang_levels[j] for j in lnbrs]

        # can_refine_angular: blocked if element is at a NCF face OR any
        # face/edge neighbour has a different spatial refinement level.
        # The edge-neighbour check is new: a non-conforming edge requires the
        # same matching treatment as a non-conforming face, so we must block it.
        my_slvl = spatial_levels[iel]
        edge_ncf = any(spatial_levels[j] != my_slvl for j in local_edge_nbrs[iel]) ||
                   any(slvl             != my_slvl for slvl in cr_edge_spa_nbr[iel])

        records[iel] = ElementRefinementRecord(
            my_slvl,
            max_ang_levels[iel],
            lnbrs,
            ln_spa,
            ln_ang,
            cr_spa_nbr[iel],
            cr_ang_nbr[iel],
            !spa_nc[iel] && !edge_ncf,
        )
    end

    return records
end

"""
    build_angular_refinement_mask(records) -> BitVector

Return a `BitVector` of length `nelem` where entry `iel` is `true` iff
element `iel` may be angularly refined without violating the spatial-angular
non-conformity constraint.
"""
function build_angular_refinement_mask(records::Vector{ElementRefinementRecord})::BitVector
    mask = BitVector(undef, length(records))
    for (iel, rec) in enumerate(records)
        mask[iel] = rec.can_refine_angular
    end
    return mask
end

"""
    verify_element_refinement_records(records, nelem, comm, rank;
                                      uniform_spatial=false)

Run consistency checks on refinement records and print a diagnostic report.

Checks (always)
---------------
- Local neighbour symmetry: if `j` appears in `records[i].local_neighbor_elems`
  then `i` must appear in `records[j].local_neighbor_elems`.
- Total neighbour count (local + cross-rank) in [1, 26]: rank-independent
  because cross-rank neighbours are found via AllGather.

Checks (uniform_spatial=true only)
------------------------------------
- All elements must have `can_refine_angular = true`.
- No cross-rank neighbour may have a different spatial level.

Returns `true` if all checks pass.
"""
function verify_element_refinement_records(
    records::Vector{ElementRefinementRecord},
    nelem::Int, comm, rank::Int;
    uniform_spatial::Bool = false
)::Bool
    ok = true

    # ── 1. Uniform-spatial assertions ─────────────────────────────────────────
    if uniform_spatial
        n_blocked = count(!r.can_refine_angular for r in records)
        if n_blocked > 0
            @warn "[$rank] verify_refinement: $n_blocked/$nelem elements blocked " *
                  "on a spatially uniform mesh — expected 0"
            ok = false
        end
        n_cross_nc = count(records) do r
            any(lvl != r.spatial_level for lvl in r.cross_rank_neighbor_spatial_levels)
        end
        if n_cross_nc > 0
            @warn "[$rank] verify_refinement: $n_cross_nc elements have cross-rank " *
                  "neighbours at a different spatial level on a uniform mesh — expected 0"
            ok = false
        end
    end

    # ── 2. Local neighbour symmetry ───────────────────────────────────────────
    nbr_sets = [Set(records[iel].local_neighbor_elems) for iel = 1:nelem]
    n_asym = 0
    for iel = 1:nelem
        for jel in nbr_sets[iel]
            iel ∉ nbr_sets[jel] && (n_asym += 1)
        end
    end
    if n_asym > 0
        @warn "[$rank] verify_refinement: $n_asym asymmetric local neighbour " *
              "pairs (j in nbrs[i] but i ∉ nbrs[j])"
        ok = false
    end

    # ── 3. Total neighbour count (local + cross-rank, face+edge only) ─────────
    counts = [length(records[iel].local_neighbor_elems) +
              length(records[iel].cross_rank_neighbor_spatial_levels)
              for iel = 1:nelem]
    n_zero = count(==(0), counts)
    n_zero > 0 && @warn "[$rank] verify_refinement: $n_zero elements have 0 total neighbours"

    local_min  = nelem > 0 ? minimum(counts) : 0
    local_max  = nelem > 0 ? maximum(counts) : 0
    local_sum  = sum(counts)
    global_min = MPI.Allreduce(local_min, MPI.MIN, comm)
    global_max = MPI.Allreduce(local_max, MPI.MAX, comm)
    global_sum = MPI.Allreduce(local_sum, MPI.SUM, comm)
    global_n   = MPI.Allreduce(nelem,     MPI.SUM, comm)

    if rank == 0
        @info "verify_refinement: neighbour count (face+edge, corners excluded) " *
              "min=$global_min max=$global_max " *
              "mean=$(round(global_sum/global_n, digits=2)) " *
              "(should match serial)"
        if global_max > 18   # 6 face + 12 edge neighbours max in 3D
            @warn "verify_refinement: max total neighbour count $global_max > 18 " *
                  "(face+edge max in 3D is 18: 6 face + 12 edge)"
            ok = false
        end
    end

    # ── 4. can_refine_angular stats ───────────────────────────────────────────
    n_blocked_local  = count(!r.can_refine_angular for r in records)
    n_blocked_global = MPI.Allreduce(n_blocked_local, MPI.SUM, comm)
    if rank == 0
        @info "verify_refinement: " *
              "$(global_n - n_blocked_global)/$global_n elements eligible " *
              "for angular refinement across all ranks"
    end

    global_ok = Bool(MPI.Allreduce(Int32(ok), MPI.MIN, comm))
    if rank == 0
        global_ok ? @info("verify_refinement: all checks PASSED") :
                    @warn("verify_refinement: some checks FAILED — see warnings above")
    end
    return global_ok
end

export ElementRefinementRecord
export build_element_refinement_records
export build_angular_refinement_mask
export verify_element_refinement_records

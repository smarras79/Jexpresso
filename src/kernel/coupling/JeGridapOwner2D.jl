using MPI
using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads

# --------------------------
# Robust geometry utilities
# --------------------------

@inline function _dist2_point_seg(px, py, ax, ay, bx, by)
    vx = bx - ax;  vy = by - ay
    wx = px - ax;  wy = py - ay
    seg2 = vx*vx + vy*vy
    if seg2 == 0.0
        dx = px - ax; dy = py - ay
        return dx*dx + dy*dy
    else
        t = (wx*vx + wy*vy) / seg2
        t = ifelse(t < 0.0, 0.0, ifelse(t > 1.0, 1.0, t))
        qx = ax + t*vx; qy = ay + t*vy
        dx = px - qx; dy = py - qy
        return dx*dx + dy*dy
    end
end

# Edge-aware polygon test (ray-casting + on-edge check), scale-aware tol
function _point_in_polygon_with_edges(px::Float64, py::Float64, poly)::Bool
    xs = (p->p[1]).(poly); ys = (p->p[2]).(poly)
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    span = max(xmax - xmin, ymax - ymin, 1.0)
    tol  = 1e-12 * span
    if px < xmin - tol || px > xmax + tol || py < ymin - tol || py > ymax + tol
        return false
    end
    tol2 = tol*tol
    @inbounds for i in 1:length(poly)
        p = poly[i]; q = poly[i == length(poly) ? 1 : i+1]
        if _dist2_point_seg(px, py, p[1], p[2], q[1], q[2]) <= tol2
            return true
        end
    end
    inside = false
    @inbounds for i in 1:length(poly)
        p = poly[i]; q = poly[i == length(poly) ? 1 : i+1]
        if ((p[2] > py) != (q[2] > py))
            xint = p[1] + (py - p[2]) * (q[1] - p[1]) / (q[2] - p[2])
            if xint >= px - tol
                inside = !inside
            end
        end
    end
    return inside
end

# --------------------------
# Fast per-rank owner locator
# --------------------------

"""
OwnerLocator2D

Per-rank spatial index to locate owners of arbitrary (x,y) points quickly.

Fields:
  • cell_coords :: Vector{Vector{Point{2,Float64}}}   # local polygons (owned+ghost)
  • bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax       # per-cell AABBs
  • l2owner_part :: Vector{Int32}                     # local cell -> owner (PartitionedArrays part id)
  • bins :: Dict{Tuple{Int,Int}, Vector{Int}}         # (ix,iy) -> candidate cell indices
  • xmin, xmax, ymin, ymax :: Float64                 # local patch AABB
  • nx, ny :: Int                                     # number of bins in x/y
"""
struct OwnerLocator2D
    cell_coords::Vector{Vector{Point{2,Float64}}}
    bbox_xmin::Vector{Float64}
    bbox_xmax::Vector{Float64}
    bbox_ymin::Vector{Float64}
    bbox_ymax::Vector{Float64}
    l2owner_part::Vector{Int32}
    bins::Dict{Tuple{Int,Int}, Vector{Int}}
    xmin::Float64; xmax::Float64; ymin::Float64; ymax::Float64
    nx::Int; ny::Int
end

@inline _clampi(x::Int, a::Int, b::Int) = x < a ? a : (x > b ? b : x)

# Heuristic: nbins from cell count and aspect ratio
function _choose_bins(ncells::Int, xmin, xmax, ymin, ymax)
    if ncells ≤ 0
        return 1, 1
    end
    side = max(1, floor(Int, sqrt(ncells)))             # ~√N bins total
    dx = max(xmax - xmin, 1e-9); dy = max(ymax - ymin, 1e-9)
    ar = dx / dy
    nx = max(1, floor(Int, side * sqrt(ar)))
    ny = max(1, floor(Int, side / max(sqrt(ar), 1e-9)))
    return nx, ny
end

"""
    build_owner_locator2d(partitioned_model; nbins=:auto)

Build `OwnerLocator2D` on the **current rank** from a GridapDistributed model.

- Uses `Triangulation(model)` → `local_views` to get the **serial, overlapped** (owned+ghost) triangulation per rank.  
- Extracts geometry via `get_cell_coordinates(::Triangulation)` and ownership via
  `get_cell_gids(model)` + `PartitionedArrays.local_to_owner`. [2](https://p4est.github.io/api/p4est-2.8.6/p8est__ghost_8h.html)[1](https://github.com/cburstedde/p4est/blob/master/src/p8est_search.h)

If `nbins=:auto`, picks a good bin count from the number of local cells and aspect ratio.
"""
function build_owner_locator2d(partitioned_model; nbins=:auto)
    Ωd        = Triangulation(partitioned_model)
    cell_gids = get_cell_gids(partitioned_model)

    loc_ref = Ref{OwnerLocator2D}()

    PartitionedArrays.map(local_views(Ωd), PartitionedArrays.partition(cell_gids)) do Ωloc, part
        # Serial geometry on overlapped triangulation  [2](https://p4est.github.io/api/p4est-2.8.6/p8est__ghost_8h.html)
        cell_coords = get_cell_coordinates(Ωloc)
        nloc = length(cell_coords)

        # Ownership maps (local idx -> owner part id)  [1](https://github.com/cburstedde/p4est/blob/master/src/p8est_search.h)
        l2owner_part = PartitionedArrays.local_to_owner(part)

        # Per-cell AABBs + local patch AABB
        bbox_xmin = Vector{Float64}(undef, nloc)
        bbox_xmax = Vector{Float64}(undef, nloc)
        bbox_ymin = Vector{Float64}(undef, nloc)
        bbox_ymax = Vector{Float64}(undef, nloc)
        xmin = +Inf; xmax = -Inf; ymin = +Inf; ymax = -Inf

        @inbounds for iloc in 1:nloc
            poly = cell_coords[iloc]
            xs = (p->p[1]).(poly); ys = (p->p[2]).(poly)
            xmn, xmx = minimum(xs), maximum(xs)
            ymn, ymx = minimum(ys), maximum(ys)
            bbox_xmin[iloc] = xmn; bbox_xmax[iloc] = xmx
            bbox_ymin[iloc] = ymn; bbox_ymax[iloc] = ymx
            xmin = min(xmin, xmn); xmax = max(xmax, xmx)
            ymin = min(ymin, ymn); ymax = max(ymax, ymx)
        end

        # Bin grid
        nx, ny = nbins === :auto ? _choose_bins(nloc, xmin, xmax, ymin, ymax) : nbins
        nx = max(nx,1); ny = max(ny,1)
        bins = Dict{Tuple{Int,Int}, Vector{Int}}()

        # Helpers
        dx = max(xmax - xmin, 1e-12);  dy = max(ymax - ymin, 1e-12)
        _ix(x) = _clampi(1 + floor(Int, (x - xmin) / dx * nx), 1, nx)
        _iy(y) = _clampi(1 + floor(Int, (y - ymin) / dy * ny), 1, ny)

        # Insert each cell AABB into overlapping bins
        @inbounds for iloc in 1:nloc
            ia = _ix(bbox_xmin[iloc]); ib = _ix(bbox_xmax[iloc])
            ja = _iy(bbox_ymin[iloc]); jb = _iy(bbox_ymax[iloc])
            for i in ia:ib, j in ja:jb
                push!(get!(bins, (i,j), Int[]), iloc)
            end
        end

        loc_ref[] = OwnerLocator2D(cell_coords, bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax,
                                   l2owner_part, bins, xmin, xmax, ymin, ymax, nx, ny)
        nothing
    end

    return loc_ref[]
end

# --------------------------
# Fast batched query
# --------------------------

"""
    owns_points_batch!(locator::OwnerLocator2D, comm, X::AbstractMatrix{<:Real};
                       unique_to_owner=false, threaded=true)

For a batch of M points stored as a matrix `X` with size (M,2) (columns: x,y), return 3 vectors:

  • found::Vector{Bool}         – whether the point lies in some local cell (owned+ghost)
  • owner_mpi::Vector{Int32}    – MPI rank that **owns** the containing cell (0-based)
  • is_mine::Vector{Bool}       – owner_mpi == my rank

If `unique_to_owner=true`, we report `found[i]=false` for points whose owner is not my rank
(this way **only the owner** returns true, neighbors with ghosts return false).

Notes:
  • Uses bin lookup + AABB + robust polygon test per candidate.
  • Thread-safe; set `threaded=false` to disable threading.
"""
function owns_points_batch!(locator::OwnerLocator2D, comm, X::AbstractMatrix{<:Real};
                            unique_to_owner::Bool=false, threaded::Bool=true)

    myrank = MPI.Comm_rank(comm)
    M = size(X, 1)
    found      = falses(M)
    owner_mpi  = fill(Int32(-1), M)
    is_mine    = falses(M)

    # Local helpers
    @inline function _bin_of(x::Float64, y::Float64)
        dx = max(locator.xmax - locator.xmin, 1e-12)
        dy = max(locator.ymax - locator.ymin, 1e-12)
        i = _clampi(1 + floor(Int, (x - locator.xmin) / dx * locator.nx), 1, locator.nx)
        j = _clampi(1 + floor(Int, (y - locator.ymin) / dy * locator.ny), 1, locator.ny)
        return i, j
    end

    # Per-point worker (capturing locator by closure is fine; it's read-only)
    function _process_range(rng)
        @inbounds for k in rng
            px = Float64(X[k,1]); py = Float64(X[k,2])

            # Primary bin
            i, j = _bin_of(px, py)
            candidates = get(locator.bins, (i,j), Int[])

            # If empty, peek neighbors (3x3)
            if isempty(candidates)
                for ii in max(1,i-1):min(locator.nx,i+1), jj in max(1,j-1):min(locator.ny,j+1)
                    candidates = get(locator.bins, (ii,jj), Int[])
                    !isempty(candidates) && break
                end
            end

            # Search candidates
            local_found = false
            local_owner_mpi = Int32(-1)

            for iloc in candidates
                # AABB prune
                if px < locator.bbox_xmin[iloc] || px > locator.bbox_xmax[iloc] ||
                   py < locator.bbox_ymin[iloc] || py > locator.bbox_ymax[iloc]
                    continue
                end
                # Exact polygon test (edge-aware)
                if _point_in_polygon_with_edges(px, py, locator.cell_coords[iloc])
                    local_found = true
                    owner_part = locator.l2owner_part[iloc]      # 1-based part id  [1](https://github.com/cburstedde/p4est/blob/master/src/p8est_search.h)
                    local_owner_mpi = Int32(owner_part - 1)      # 0-based MPI
                    break
                end
            end

            if local_found
                if unique_to_owner && local_owner_mpi != myrank
                    # suppress ghost hits so only owner returns true
                    found[k] = false; owner_mpi[k] = Int32(-1); is_mine[k] = false
                else
                    found[k]     = true
                    owner_mpi[k] = local_owner_mpi
                    is_mine[k]   = (local_owner_mpi == myrank)
                end
            end
        end
    end

    if threaded && nthreads() > 1
        # Chunk the batch across threads
        chunk = cld(M, nthreads())
        @threads for t in 1:nthreads()
            lo = (t-1)*chunk + 1
            hi = min(t*chunk, M)
            lo <= hi && _process_range(lo:hi)
        end
    else
        _process_range(1:M)
    end

    return found, owner_mpi, is_mine
end

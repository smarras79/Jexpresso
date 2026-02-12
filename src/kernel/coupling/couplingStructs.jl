using Gridap
using GridapDistributed
using PartitionedArrays
using Base.Threads

function couplingAlloc(nrank1, nrank2, T;)

    couple = zeros(T, nrank1, nrank2)

    return couple
end

"""
    St_coupling

Holds all inter-code coupling state needed during the time loop:
  - world communicator for inter-code MPI
  - local communicator for Jexpresso-internal MPI
  - rank/size info on both communicators
  - pre-allocated send/receive buffers for field exchange with Alya
  - coupling metadata received during initialization (ndime, remote grid info, etc.)
"""
mutable struct St_coupling
    comm_world::MPI.Comm
    comm_local::MPI.Comm
    wrank::Int
    wsize::Int
    lrank::Int
    lsize::Int
    nranks_alya::Int

    # Remote grid metadata from Alya
    ndime::Int32
    rem_min::Vector{Float32}
    rem_max::Vector{Float32}
    rem_nx::Vector{Int32}
    alya2world::Vector{Int32}

    # Pre-allocated exchange buffers (sized after mesh is known)
    send_buf::Vector{Float64}
    recv_buf::Vector{Float64}

    # Coupling allocation matrix
    couple::Matrix{Int64}
end

"""
    je_couplingSetup(inputs) -> Union{St_coupling, Nothing}

Initialize coupling data structures using the pre-stored coupling data
from `Jexpresso-mini-coupled.jl`.

Returns `St_coupling` if coupling is active, `nothing` otherwise.
"""
function je_couplingSetup(inputs)

    if !inputs[:lcoupling]
        return nothing
    end

    comm_local = get_mpi_comm()
    lrank = MPI.Comm_rank(comm_local)
    lsize = MPI.Comm_size(comm_local)

    println("[Jexpresso rank $lrank] Coupling setup..."); flush(stdout)

    # Retrieve coupling data stored by Jexpresso-mini-coupled.jl
    cdata = get_coupling_data()
    if cdata === nothing
        if lrank == 0
            @warn "[Jexpresso] lcoupling=true but no coupling data found. " *
                  "Are you running via Jexpresso-mini-coupled.jl?"
        end
        return nothing
    end

    comm_world  = get_mpi_comm_world()
    wrank       = MPI.Comm_rank(comm_world)
    wsize       = MPI.Comm_size(comm_world)
    nranks_alya = wsize - lsize

    ndime      = cdata[:ndime]
    rem_min    = cdata[:rem_min]
    rem_max    = cdata[:rem_max]
    rem_nx     = cdata[:rem_nx]
    alya2world = cdata[:alya2world]

    couple = couplingAlloc(wsize, wsize, Int64)

    # Initial empty exchange buffers (will be resized when mesh coupling map is built)
    send_buf = Float64[]
    recv_buf = Float64[]

    coupling = St_coupling(
        comm_world, comm_local,
        wrank, wsize,
        lrank, lsize,
        nranks_alya,
        ndime, rem_min, rem_max, rem_nx, alya2world,
        send_buf, recv_buf,
        couple
    )

    println("[Jexpresso rank $lrank] Coupling setup DONE " *
            "(world_rank=$wrank, world_size=$wsize, nranks_alya=$nranks_alya, ndime=$ndime)"); flush(stdout)

    return coupling
end

#------------------------------------------------------------------------------------
# Coupling exchange: send data to Alya and receive data from Alya
# Called within the time loop at every coupling step.
#
# This uses MPI_Sendrecv on COMM_WORLD between the Jexpresso local root
# and the Alya root (world rank 0).
# Non-root Jexpresso ranks participate via Bcast on the local communicator.
#------------------------------------------------------------------------------------

"""
    coupling_send_recv!(coupling::St_coupling, send_data::AbstractVector{Float64},
                        recv_data::AbstractVector{Float64};
                        alya_root::Int=0, tag::Int=100)

Exchange field data with Alya:
  1. Jexpresso local root (lrank==0) does MPI.Sendrecv! with the Alya root on comm_world
  2. The received data is then broadcast to all Jexpresso local ranks via comm_local

`send_data` and `recv_data` must be the same length on the root.
On non-root ranks, `recv_data` is filled by the broadcast.
"""
function coupling_send_recv!(coupling::St_coupling,
                             send_data::AbstractVector{Float64},
                             recv_data::AbstractVector{Float64};
                             alya_root::Int=0, tag::Int=100)

    if coupling.lrank == 0
        # Root exchanges with Alya root via world communicator
        MPI.Sendrecv!(send_data, recv_data, coupling.comm_world;
                      dest=alya_root, sendtag=tag,
                      source=alya_root, recvtag=tag)
    end

    # Broadcast received data to all local Jexpresso ranks
    MPI.Bcast!(recv_data, 0, coupling.comm_local)
end

"""
    coupling_send!(coupling::St_coupling, send_data::AbstractVector{Float64};
                   alya_root::Int=0, tag::Int=200)

Send field data to Alya (non-blocking from Jexpresso root).
Only the Jexpresso local root sends.
"""
function coupling_send!(coupling::St_coupling,
                        send_data::AbstractVector{Float64};
                        alya_root::Int=0, tag::Int=200)

    if coupling.lrank == 0
        MPI.Send(send_data, coupling.comm_world; dest=alya_root, tag=tag)
    end
end

"""
    coupling_recv!(coupling::St_coupling, recv_data::AbstractVector{Float64};
                   alya_root::Int=0, tag::Int=200)

Receive field data from Alya (blocking on Jexpresso root),
then broadcast to all local ranks.
"""
function coupling_recv!(coupling::St_coupling,
                        recv_data::AbstractVector{Float64};
                        alya_root::Int=0, tag::Int=200)

    if coupling.lrank == 0
        MPI.Recv!(recv_data, coupling.comm_world; source=alya_root, tag=tag)
    end

    MPI.Bcast!(recv_data, 0, coupling.comm_local)
end

#------------------------------------------------------------------------------------
# Receive the grid coordinate from Alya (structured and regular only)
#------------------------------------------------------------------------------------
function distribute_and_count!(
    rem_nx,
    rem_min,
    rem_max,
    ndime,
    nranks2,
    in_my_rank,
    a,
    wrank,
    alya2world,
    mesh)

    #
    # is_in_my_rank must be
    # found from Jexpresso
    #

    nx, ny, nz = rem_nx
    nxy        = nx * ny
    nmax       = nxy * nz
    rem_dx     = similar(rem_max)

    @assert ndime in (2, 3) "ndime is typically 2 or 3"
    @assert length(rem_min) ≥ ndime && length(rem_nx) ≥ ndime

    r     = mod(nmax, nranks2 - 1)
    npoin = nmax ÷ (nranks2 - 1)

    ri = zeros(Int32, ndime)
    x  = zeros(Float64, ndime)

    rem_dx[1:ndime] = (rem_max[1:ndime] .- rem_min[1:ndime])./(rem_nx[1:ndime] .- 1)

    #
    # Local max min on this rank:
    #
    lxmin = min(mesh.x); lxmax = max(mesh.x)
    lymin = min(mesh.y); lymax = max(mesh.y)
    lzmin = min(mesh.z); lzmax = max(mesh.z)

    my_boundingBox = St_myBoundingBox(lxmin, lxmax, lymin, lymax, lzmin, lzmax)

    @inbounds for ipoin in 1:nmax
        i0   = ipoin - 1
        iz   = i0 ÷ nxy
        remz = i0 - iz * nxy
        iy   = remz ÷ nx
        ix   = mod(remz - iy * nx, nx)

        if ndime ≥ 1; ri[1] = ix; end
        if ndime ≥ 2; ri[2] = iy; end
        if ndime ≥ 3; ri[3] = iz; end

        x[1:ndime] = rem_min[1:ndime] .+ Float32.(ri[1:ndime]) .* rem_dx[1:ndime]

        ###
        ### in, yin, zin = is_on_my_rank(local::St_myBoundingBox, global::St_myBoundingBox, x, y, z; atol=default_atol(local))

        if is_in_my_rank
            alya_rank = if ipoin ≤ r * (npoin + 1)
                (ipoin - 1) ÷ (npoin + 1) + 1
            else
                r + (ipoin - r * (npoin + 1) - 1) ÷ npoin + 1
            end

            world_rank = alya2world[alya_rank]
            a[wrank, world_rank] += 1
        end
    end
    println("$a")
    return a
end

struct St_myBoundingBox{T}
    xmin::T; xmax::T
    ymin::T; ymax::T
    zmin::T; zmax::T
end

# Scale-aware tolerance (important for large coordinate magnitudes)
@inline function default_atol(xlocal::St_myBoundingBox{T}) where {T<:Real}
    Lx = xlocal.xmax - xlocal.xmin
    Ly = xlocal.ymax - xlocal.ymin
    Lz = xlocal.zmax - xlocal.zmin
    scale = max(Lx, Ly, Lz, one(T))
    return sqrt(eps(T)) * scale   # ~1e-8 * scale for Float64
end

@inline function is_on_my_rank(xlocal::St_myBoundingBox, xglobal::St_myBoundingBox, x, y, z; atol=default_atol(xlocal))
    # Identify whether this subdomain touches the global max on each axis.
    # Those ranks get inclusive upper bound to avoid "nobody owns global max boundary".
    lastx = abs(xlocal.xmax - xglobal.xmax) <= atol
    lasty = abs(xlocal.ymax - xglobal.ymax) <= atol
    lastz = abs(xlocal.zmax - xglobal.zmax) <= atol

    # Half-open: [min, max) for interior ranks; [min, max] for ranks on global max boundary.
    xin = (x >= xlocal.xmin - atol) & (lastx ? (x <= xlocal.xmax + atol) : (x < xlocal.xmax - atol))
    yin = (y >= xlocal.ymin - atol) & (lasty ? (y <= xlocal.ymax + atol) : (y < xlocal.ymax - atol))
    zin = (z >= xlocal.zmin - atol) & (lastz ? (z <= xlocal.zmax + atol) : (z < xlocal.zmax - atol))

    return xin & yin & zin
end

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
  * cell_coords :: Vector{Vector{Point{2,Float64}}}   # local polygons (owned+ghost)
  * bbox_xmin, bbox_xmax, bbox_ymin, bbox_ymax       # per-cell AABBs
  * l2owner_part :: Vector{Int32}                     # local cell -> owner (PartitionedArrays part id)
  * bins :: Dict{Tuple{Int,Int}, Vector{Int}}         # (ix,iy) -> candidate cell indices
  * xmin, xmax, ymin, ymax :: Float64                 # local patch AABB
  * nx, ny :: Int                                     # number of bins in x/y
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
"""
function build_owner_locator2d(partitioned_model; nbins=:auto)
    Ωd        = Triangulation(partitioned_model)
    cell_gids = get_cell_gids(partitioned_model)

    loc_ref = Ref{OwnerLocator2D}()

    PartitionedArrays.map(local_views(Ωd), PartitionedArrays.partition(cell_gids)) do Ωloc, part
        cell_coords = get_cell_coordinates(Ωloc)
        nloc = length(cell_coords)

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

  * found::Vector{Bool}         - whether the point lies in some local cell (owned+ghost)
  * owner_mpi::Vector{Int32}    - MPI rank that **owns** the containing cell (0-based)
  * is_mine::Vector{Bool}       - owner_mpi == my rank
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
                    owner_part = locator.l2owner_part[iloc]      # 1-based part id
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

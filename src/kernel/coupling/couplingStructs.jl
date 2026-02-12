using .Jexpresso
using Gridap
using GridapDistributed
using PartitionedArrays
using PartitionedArrays
using Base.Threads

function couplingAlloc(nrank1, nrank2, T;)
    
    couple = zeros(T, nrank1, nrank2)
    
    return couple
end

function je_couplingSetup(je_mesh, lcouple)

    if !lcoupling
        return nothing
    end
    
    println("Coupling setup..."); flush(stdout)
    
    world, wrank, wsize = je_mpi_init()
    println("[Jexpresso rank $wrank] World size: $wsize"); flush(stdout)

    # Read APPID (0 or 1) from environment
    appid = try parse(Int, get(ENV, "APPID", "2")) catch; 2 end
    println("[appid: $appid"); flush(stdout)
    if appid < 0
        if wrank == 0
            println("[Jexpresso] ERROR: APPID not set. Launch with -x APPID=0 (Fortran) and -x APPID=1 (Julia).")
        end
        MPI.Abort(world, 1)
    end

    # Split WORLD into per-app local comms
    println("[Split before $wrank"); flush(stdout)
    local_comm = MPI.Comm_split(world, appid, wrank)
    println("[Split after $wrank"); flush(stdout)

    lrank        = MPI.Comm_rank(local_comm)
    lsize        = MPI.Comm_size(local_comm)
    nranks1      = lsize                      # Jexpresso
    nranks2      = wsize - lsize              # Other code
    local_chars  = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
    recv_buffer  = nothing
    is_jexpresso = (appid == 2)

    MPI.Gather!(local_chars, recv_buffer, 0, world)

    # Set coupling mode to prevent auto-execution on module load
    ENV["JEXPRESSO_COUPLING_MODE"] = "false"

    # Set command line arguments for Jexpresso
    #push!(empty!(ARGS), "CompEuler", "wave1d")
    
    # Load Jexpresso module (setup doesn't run yet because of JEXPRESSO_COUPLING_MODE)
    #println("[Jexpresso rank $wrank] Loading Jexpresso module (JIT compilation may take minutes)..."); flush(stdout)
    #include("./src/Jexpresso.jl")
    #println("[Jexpresso rank $wrank] Jexpresso module loaded."); flush(stdout)

    # Set the custom local communicator
    Jexpresso.set_mpi_comm(local_comm)
    println("[AAAAAAA Jexpresso rank $wrank] MPI communicator set (local_comm, size=$lsize)."); flush(stdout)

    #--------------------------------------------------------------------------------------------
    # Receive ndime from Alya via Bcast on COMM_WORLD (all ranks must participate)
    # Alya: call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    # Use Int32 to match Fortran's MPI_INTEGER (4 bytes)
    #--------------------------------------------------------------------------------------------
    couple = couplingAlloc(wsize, wsize, Int64;)

    ndime_buf = Vector{Int32}(undef, 1)
    MPI.Bcast!(ndime_buf, 0, world)
    
    ndime = ndime_buf[1]
    println("[Jexpresso rank $wrank] Received ndime = $ndime from Alya"); flush(stdout)

    rem_min  = Vector{Float32}(undef, 3)
    rem_max  = Vector{Float32}(undef, 3)
    rem_nx   = Vector{Int32}(undef, 3)
    for idime in 1:3
        MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
        MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
    end
    println("[Jexpresso rank $wrank] Received ndime      = $ndime      from Alya"); flush(stdout)
    println("[Jexpresso rank $wrank] Received rem_min    = $rem_min    from Alya"); flush(stdout)
    println("[Jexpresso rank $wrank] Received rem_max    = $rem_max    from Alya"); flush(stdout)
    println("[Jexpresso rank $wrank] Received rem_nx     = $rem_nx     from Alya"); flush(stdout)

    alya2world_l = zeros(Int32, nranks2)
    alya2world   = MPI.Allreduce(alya2world_l,MPI.SUM,world)

    println("[Jexpresso rank $wrank] Received alya2world = $alya2world from Alya"); flush(stdout)

    a_l = zeros(Int32, wsize,wsize)
    a   = MPI.Allreduce(a_l,MPI.SUM,world)
#=
    distribute_and_count!(
        rem_nx,
        rem_min,
        rem_max,
        ndime,
        nranks2,
        is_jexpresso,
        a,
        wrank,
        alya2world,
        mesh)
=#
    #--------------------------------------------------------------------------------------------
    # END Receive ndime from Alya
    #--------------------------------------------------------------------------------------------
    
   # MPI.Finalize()
    
    println("Coupling setup..."); flush(stdout)
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

    @info " ASSAASSAAS"
    @info my_boundingBox.xmin, my_boundingBox.ymin, my_boundingBox.zmin
    @info my_boundingBox.xmax, my_boundingBox.ymax, my_boundingBox.zmax
    @info "EWWEEWWEWEWEW"
    
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

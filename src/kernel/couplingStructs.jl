using .Jexpresso

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

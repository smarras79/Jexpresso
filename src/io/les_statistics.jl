Base.@kwdef mutable struct LESStatCache{T <: AbstractFloat, nz, dims_buf, dims_stress, backend, VT1, VT2}
    z_levels        ::Vector{Float64}       = zeros(Float64, nz)
    z_groups        ::Vector{Vector{Int64}} = [Int64[] for _ in 1:nz]
    npts_per_z      ::Vector{Int64}         = zeros(Int64,   nz)
    local_sum       ::VT2                   = KernelAbstractions.zeros(backend, T, dims_buf...)
    global_sum      ::VT2                   = KernelAbstractions.zeros(backend, T, dims_buf...)
    prof            ::VT1                   = KernelAbstractions.zeros(backend, T, dims_buf[2])
    local_stress    ::VT2                   = KernelAbstractions.zeros(backend, T, dims_stress...)
    global_stress   ::VT2                   = KernelAbstractions.zeros(backend, T, dims_stress...)
    stress_prof     ::VT1                   = KernelAbstractions.zeros(backend, T, dims_stress[2])
    # temporal accumulation (running sums over statistics output steps)
    n_samples       ::Int                   = 0
    sum_global_sum    ::VT2                 = KernelAbstractions.zeros(backend, T, dims_buf...)
    sum_global_stress ::VT2                 = KernelAbstractions.zeros(backend, T, dims_stress...)
    # online accumulation (Approach 2): local running sums per time step, Allreduce only at end
    n_online_samples     ::Int  = 0
    local_online_sum     ::VT2  = KernelAbstractions.zeros(backend, T, dims_buf...)
    local_online_stress  ::VT2  = KernelAbstractions.zeros(backend, T, dims_stress...)
    global_online_sum    ::VT2  = KernelAbstractions.zeros(backend, T, dims_buf...)
    global_online_stress ::VT2  = KernelAbstractions.zeros(backend, T, dims_stress...)
end


# ---- Spectral cache: element-wise Lagrange interpolation GLL → uniform y-grid ----
mutable struct LESSpectralCache
    interp_mat ::Matrix{Float64}  # (ngl × ngl) Lagrange matrix: GLL → uniform
    nel_y      ::Int              # number of y-elements per strip (global)
    N_unif     ::Int              # uniform y-points per strip = nel_y*(ngl-1)+1
    nk         ::Int              # one-sided spectral bins = N_unif÷2+1
end

Base.@kwdef mutable struct LESCrossSection{T <: AbstractFloat}
    xz_coords    ::Vector{Tuple{Float64,Float64}}  # sorted unique (x,z) pairs
    xz_groups    ::Vector{Vector{Int64}}            # local point indices per (x,z)
    npts_per_xz  ::Vector{Int64}                   # global y-count per (x,z)
    local_mean   ::Matrix{T}                        # nxz × nprofiles
    global_mean  ::Matrix{T}                        # nxz × nprofiles
    local_stress ::Matrix{T}                        # nxz × nstress
    global_stress::Matrix{T}                        # nxz × nstress
    prof_buf     ::Vector{T}                        # scratch, length nprofiles
    stress_buf   ::Vector{T}                        # scratch, length nstress
    # temporal accumulation (running sums over statistics output steps)
    n_samples    ::Int                              # number of snapshots accumulated
    sum_mean     ::Matrix{T}                        # nxz × nprofiles: running sum of global_mean
    sum_stress   ::Matrix{T}                        # nxz × nstress:   running sum of global_stress
    # VTK cell connectivity — works for both flat and terrain-following meshes
    xz_cells     ::Vector{NTuple{4,Int64}}          # quad corner indices (into xz_coords)
    # spectral cache (built once; nothing if spectra disabled)
    spec_cache   ::Union{Nothing, LESSpectralCache}
    # online accumulation (Approach 2): local running sums per time step, Allreduce only at end
    n_online_samples     ::Int
    local_online_mean    ::Matrix{T}
    local_online_stress  ::Matrix{T}
    global_online_mean   ::Matrix{T}
    global_online_stress ::Matrix{T}
end

function build_les_stat_cache(mesh, nprofiles::Int, nstress::Int, T, backend)
    nprofiles == 0 && return nothing

    comm = MPI.COMM_WORLD

    z         = Array(mesh.z)  # ensure CPU array
    gip2owner = Array(mesh.gip2owner)

    # Gather all unique z values across all ranks
    local_z_unique = sort(unique(round.(z; digits=8)))
    rank = MPI.Comm_rank(comm)

    # Gather variable-length arrays to rank 0, merge, then broadcast back
    all_z_arrays = MPI.gather(local_z_unique, comm)
    if rank == 0
        z_levels = sort(unique(round.(vcat(all_z_arrays...); digits=8)))
    else
        z_levels = Float64[]
    end
    nz_global = MPI.bcast(length(z_levels), 0, comm)
    if rank != 0
        resize!(z_levels, nz_global)
    end
    MPI.Bcast!(z_levels, 0, comm)

    nz = length(z_levels)

    # Build local index groups for each global z-level — owned points only
    # so that Allreduce counts each shared boundary node exactly once.
    z_groups = Vector{Vector{Int64}}(undef, nz)
    local_counts = zeros(Int64, nz)
    for iz in 1:nz
        ids = findall(ip -> abs(z[ip] - z_levels[iz]) < 1e-6 && gip2owner[ip] == rank,
                      1:length(z))
        z_groups[iz] = ids
        local_counts[iz] = length(ids)
    end

    # Get global point counts per z-level
    npts_per_z = zeros(Int64, nz)
    MPI.Allreduce!(local_counts, npts_per_z, MPI.SUM, comm)

    dims_buf    = (Int64(nz), Int64(nprofiles))
    dims_stress = (Int64(nz), Int64(nstress))
    VT1   = typeof(KernelAbstractions.zeros(backend, T, dims_buf[2]))
    VT2   = typeof(KernelAbstractions.zeros(backend, T, dims_buf...))
    cache = LESStatCache{T, Int64(nz), dims_buf, dims_stress, backend, VT1, VT2}()
    copyto!(cache.z_levels,   z_levels)
    copyto!(cache.z_groups,   z_groups)
    copyto!(cache.npts_per_z, npts_per_z)
    return cache
end

"""
    horizontal_mean!(cache, uaux, qe, ET, comm)

For each z-level, call `user_les_profiles!` at every local point to accumulate
all profile sums into `cache.local_sum` (nz × nprofiles), then perform a single
`MPI.Allreduce!` and divide by `npts_per_z`. Results are stored in `cache.global_sum`.
"""
function horizontal_mean!(cache, uaux, qe, ET, comm)
    nz        = length(cache.z_levels)
    nprofiles = size(cache.local_sum,    2)
    nstress   = size(cache.local_stress, 2)
    fill!(cache.local_sum,    0.0)
    fill!(cache.local_stress, 0.0)

    means_buf = cache.prof        # scratch: means (length nprofiles)
    prof_buf  = cache.stress_prof # scratch: raw products (length nstress)

    for iz in 1:nz
        for ip in cache.z_groups[iz]
            user_les_profiles!(means_buf, prof_buf, @view(uaux[ip,:]), @view(qe[ip,:]), ET)
            for k in 1:nprofiles
                cache.local_sum[iz, k] += means_buf[k]
            end
            for k in 1:nstress
                cache.local_stress[iz, k] += prof_buf[k]
            end
        end
    end

    MPI.Allreduce!(cache.local_sum,    cache.global_sum,    MPI.SUM, comm)
    MPI.Allreduce!(cache.local_stress, cache.global_stress, MPI.SUM, comm)
    for iz in 1:nz
        npts = cache.npts_per_z[iz]
        for k in 1:nprofiles
            cache.global_sum[iz, k] /= npts
        end
        for k in 1:nstress
            cache.global_stress[iz, k] /= npts
        end
    end
end

"""
    horizontal_stress!(cache, uaux, qe, ET, comm)

Second pass: uses `cache.global_sum` (means from `horizontal_mean!`) to compute
perturbations at each point, calls `user_les_stress!`, and accumulates into
`cache.local_stress`. A single `MPI.Allreduce!` gives `cache.global_stress`.
"""
function horizontal_stress!(cache)
    nz      = length(cache.z_levels)
    nstress = size(cache.global_stress, 2)
    profp   = cache.stress_prof  # scratch (length nstress)

    for iz in 1:nz
        means = @view cache.global_sum[iz, :]
        prof  = @view cache.global_stress[iz, :]
        user_les_stress!(profp, prof, means)
        for k in 1:nstress
            cache.global_stress[iz, k] = profp[k]
        end
    end
end

function les_statistics(u, params, ::Any)

    isnothing(params.les_stat_cache) && return

    mesh  = params.mesh
    npoin = mesh.npoin
    neqs  = params.neqs
    cache = params.les_stat_cache
    ET    = params.inputs[:SOL_VARS_TYPE]
    comm  = MPI.COMM_WORLD

    uaux = params.uaux
    qe   = params.qp.qe
    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    # Single pass: accumulate means + raw products, two Allreduces
    horizontal_mean!(cache, uaux, qe, ET, comm)
    # Post-reduction: convert raw products → fluctuation statistics (no MPI)
    horizontal_stress!(cache)

    # Temporal accumulation for 1D profiles
    cache.sum_global_sum    .+= cache.global_sum
    cache.sum_global_stress .+= cache.global_stress
    cache.n_samples         += 1

    # XZ cross section (y-averaged 2D fields)
    cs = params.les_cross_section
    if !isnothing(cs)
        compute_xz_cross_section!(cs, uaux, qe, ET, comm)
        # Temporal accumulation for xz cross-sections
        cs.sum_mean   .+= cs.global_mean
        cs.sum_stress .+= cs.global_stress
        cs.n_samples  += 1
    end
end

# ================================================================================
# XZ Cross Section (y-averaged 2D fields)
# ================================================================================

"""
    build_les_spectral_cache(cs, mesh, basis)

Build a `LESSpectralCache`:
  - Uses global point count from `cs.npts_per_xz` (set by Allreduce in
    `build_les_cross_section`) to compute `nel_y` and `N_unif` correctly
    for any MPI decomposition in y.
  - Precomputes the Lagrange interpolation matrix from GLL points to
    equally-spaced points (dy = dy_el/(ngl-1)) within each element.
"""
function build_les_spectral_cache(cs::LESCrossSection, mesh, basis)
    ngl = mesh.ngl

    # Use global count — npts_per_xz is already Allreduced across all ranks
    Ny_global = maximum(cs.npts_per_xz)
    nel_y  = (Ny_global - 1) ÷ (ngl - 1)
    N_unif = nel_y * (ngl - 1) + 1
    nk     = N_unif ÷ 2 + 1

    # GLL reference points ξ ∈ [-1,1] from the SEM basis
    ξ_gll  = Array(basis.ξ)  # length ngl

    # Uniformly spaced reference points: η_i = -1 + 2(i-1)/(ngl-1)
    η_unif = [-1.0 + 2.0*(i-1)/(ngl-1) for i in 1:ngl]

    # Lagrange interpolation matrix: I[i,j] = L_j(η_i)
    #   L_j(ξ) = Π_{k≠j} (ξ - ξ_k) / (ξ_j - ξ_k)
    interp_mat = zeros(ngl, ngl)
    for i in 1:ngl
        for j in 1:ngl
            val = 1.0
            for k in 1:ngl
                k == j && continue
                val *= (η_unif[i] - ξ_gll[k]) / (ξ_gll[j] - ξ_gll[k])
            end
            interp_mat[i, j] = val
        end
    end

    return LESSpectralCache(interp_mat, nel_y, N_unif, nk)
end

"""
    compute_les_spectra!(cs, uaux, qe, mesh, ET, comm, params, t, iout)

MPI-aware turbulent power spectra computation in the spanwise (y) direction.

Algorithm:
  1. Each rank packs its local (ixz, y, u, v, w, θ) data as a flat Float64 array
     (6 floats per point).
  2. `MPI.gather` collects all local arrays to rank 0.
  3. Rank 0 decodes, sorts each (x,z) strip by y, applies element-wise Lagrange
     interpolation (GLL → uniform grid), calls `user_les_spectral!`, and writes
     one output file per statistics output step.
"""
function compute_les_spectra!(cs, uaux, qe, mesh, ET, comm, params)
    isnothing(cs.spec_cache) && return
    sc   = cs.spec_cache
    rank = MPI.Comm_rank(comm)
    nxz  = length(cs.xz_coords)
    y    = Array(mesh.y)

    # --- Pack local data: 6 floats per point = (ixz_float, y, u, v, w, θ) ---
    n_local   = sum(length(cs.xz_groups[ixz]) for ixz in 1:nxz)
    local_buf = zeros(Float64, 6 * n_local)
    ptr = 0
    for ixz in 1:nxz
        for ip in cs.xz_groups[ixz]
            if ET == PERT()
                ρ   = uaux[ip,1] + qe[ip,1]
                u_v = uaux[ip,2] / ρ
                v_v = uaux[ip,3] / ρ
                w_v = uaux[ip,4] / ρ
                θ_v = (uaux[ip,5]+qe[ip,5])/ρ - qe[ip,5]/qe[ip,1]
            else
                ρ   = uaux[ip,1]
                u_v = uaux[ip,2] / ρ
                v_v = uaux[ip,3] / ρ
                w_v = uaux[ip,4] / ρ
                θ_v = uaux[ip,5] / ρ
            end
            local_buf[ptr+1] = Float64(ixz)
            local_buf[ptr+2] = y[ip]
            local_buf[ptr+3] = u_v
            local_buf[ptr+4] = v_v
            local_buf[ptr+5] = w_v
            local_buf[ptr+6] = θ_v
            ptr += 6
        end
    end

    # --- Gather all local buffers to rank 0 ---
    all_bufs = MPI.gather(local_buf, comm)
    rank != 0 && return

    # --- On rank 0: decode into per-strip arrays ---
    combined = vcat(all_bufs...)
    n_global = length(combined) ÷ 6

    strip_ys = [Float64[] for _ in 1:nxz]
    strip_us = [Float64[] for _ in 1:nxz]
    strip_vs = [Float64[] for _ in 1:nxz]
    strip_ws = [Float64[] for _ in 1:nxz]
    strip_θs = [Float64[] for _ in 1:nxz]

    for i in 1:n_global
        base  = (i-1)*6
        ixz   = round(Int, combined[base+1])
        1 <= ixz <= nxz || continue
        push!(strip_ys[ixz], combined[base+2])
        push!(strip_us[ixz], combined[base+3])
        push!(strip_vs[ixz], combined[base+4])
        push!(strip_ws[ixz], combined[base+5])
        push!(strip_θs[ixz], combined[base+6])
    end

    # --- Interpolate each strip and compute spectra ---
    lesspectra_vars = get(params.inputs, :lesspectra_vars, String[])
    if isempty(lesspectra_vars)
        return nothing, nothing
    end

    ngl          = size(sc.interp_mat, 1)
    u_unif       = zeros(Float64, sc.N_unif, 4)
    kappa        = zeros(Float64, sc.nk)
    spectra      = zeros(Float64, sc.nk, 4)
    u_gll        = zeros(Float64, ngl)
    u_tmp        = zeros(Float64, ngl)
    all_spectra  = zeros(Float64, nxz, sc.nk, 4)
    kappa_global = zeros(Float64, sc.nk)

    for ixz in 1:nxz
        isempty(strip_ys[ixz]) && continue
        perm      = sortperm(strip_ys[ixz])
        ys_sorted = strip_ys[ixz][perm]
        vals      = hcat(strip_us[ixz][perm],
                         strip_vs[ixz][perm],
                         strip_ws[ixz][perm],
                         strip_θs[ixz][perm])  # (Ny_global × 4)

        Ly = ys_sorted[end] - ys_sorted[1]
        Ly <= 0.0 && continue

        # Element-wise Lagrange interpolation: GLL → uniform
        fill!(u_unif, 0.0)
        nel_y = sc.nel_y
        for ivar in 1:4
            for e in 1:nel_y
                for k in 1:ngl
                    idx      = (e-1)*(ngl-1) + k
                    u_gll[k] = vals[idx, ivar]
                end
                fill!(u_tmp, 0.0)
                for i in 1:ngl, j in 1:ngl
                    u_tmp[i] += sc.interp_mat[i,j] * u_gll[j]
                end
                i_start = (e-1)*(ngl-1) + 1
                n_store = (e < nel_y) ? ngl-1 : ngl
                for k in 1:n_store
                    u_unif[i_start+k-1, ivar] = u_tmp[k]
                end
            end
        end

        fill!(spectra, 0.0)
        fill!(kappa,   0.0)
        user_les_spectral!(spectra, kappa, u_unif, Ly)

        all_spectra[ixz, :, :] = spectra
        ixz == 1 && (kappa_global .= kappa)
    end

    return all_spectra, kappa_global
end


function build_les_cross_section(mesh, basis, nprofiles::Int, nstress::Int, T)
    nprofiles == 0 && return nothing

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    x         = Array(mesh.x)
    z         = Array(mesh.z)
    gip2owner = Array(mesh.gip2owner)
    npoin_local = length(x)

    xr = round.(x; digits=6)
    zr = round.(z; digits=6)

    # Build local (x,z) → point index groups — owned points only
    # so that Allreduce counts each shared boundary node exactly once.
    xz_to_local = Dict{Tuple{Float64,Float64}, Vector{Int64}}()
    for ip in 1:npoin_local
        gip2owner[ip] == rank || continue
        key = (xr[ip], zr[ip])
        push!(get!(xz_to_local, key, Int64[]), ip)
    end

    # Flatten local unique pairs for MPI gather
    local_pairs = sort(collect(keys(xz_to_local)))
    local_flat  = Float64[v for p in local_pairs for v in p]

    # Gather all pairs to rank 0, merge unique, broadcast back
    all_flats = MPI.gather(local_flat, comm)
    if rank == 0
        merged       = vcat(all_flats...)
        all_pairs    = [(merged[2i-1], merged[2i]) for i in 1:length(merged)÷2]
        unique_pairs = sort(unique(all_pairs))
        xz_flat      = Float64[v for p in unique_pairs for v in p]
    else
        xz_flat = Float64[]
    end

    nxz_global = MPI.bcast(rank == 0 ? length(xz_flat) ÷ 2 : 0, 0, comm)
    if rank != 0
        resize!(xz_flat, 2 * nxz_global)
    end
    MPI.Bcast!(xz_flat, 0, comm)

    xz_coords = [(xz_flat[2i-1], xz_flat[2i]) for i in 1:nxz_global]
    nxz       = nxz_global

    # Build xz_groups from the local dict
    xz_groups    = Vector{Vector{Int64}}(undef, nxz)
    local_counts = zeros(Int64, nxz)
    for ixz in 1:nxz
        key               = xz_coords[ixz]
        ids               = get(xz_to_local, key, Int64[])
        xz_groups[ixz]    = ids
        local_counts[ixz] = length(ids)
    end

    npts_per_xz = zeros(Int64, nxz)
    MPI.Allreduce!(local_counts, npts_per_xz, MPI.SUM, comm)

    # Build VTK quad connectivity from element structure.
    # Using connijk[iel, ix, 1, iz] (fixing iy=1) extracts the xz face of each
    # element; works for both flat and terrain-following meshes.  Quads from
    # all y-elements at the same (ie,ke) map to identical ixz corner indices
    # and are deduplicated via a Set.
    coord_to_ixz = Dict{Tuple{Float64,Float64}, Int64}(
        xz_coords[ixz] => ixz for ixz in 1:nxz)

    connijk_arr = Array(mesh.connijk)  # (nelem, ngl, ngl, ngl)
    ngl_c       = mesh.ngl
    x_arr       = Array(mesh.x)
    z_arr       = Array(mesh.z)

    seen_quads = Set{NTuple{4,Int64}}()
    xz_cells   = NTuple{4,Int64}[]

    for iel in 1:mesh.nelem
        for ix in 1:ngl_c-1, iz in 1:ngl_c-1
            ip1 = connijk_arr[iel, ix,   1, iz  ]
            ip2 = connijk_arr[iel, ix+1, 1, iz  ]
            ip3 = connijk_arr[iel, ix+1, 1, iz+1]
            ip4 = connijk_arr[iel, ix,   1, iz+1]
            k1 = (round(x_arr[ip1]; digits=6), round(z_arr[ip1]; digits=6))
            k2 = (round(x_arr[ip2]; digits=6), round(z_arr[ip2]; digits=6))
            k3 = (round(x_arr[ip3]; digits=6), round(z_arr[ip3]; digits=6))
            k4 = (round(x_arr[ip4]; digits=6), round(z_arr[ip4]; digits=6))
            i1 = get(coord_to_ixz, k1, 0)
            i2 = get(coord_to_ixz, k2, 0)
            i3 = get(coord_to_ixz, k3, 0)
            i4 = get(coord_to_ixz, k4, 0)
            (i1 > 0 && i2 > 0 && i3 > 0 && i4 > 0) || continue
            dedup_key = tuple(sort!(Int64[i1, i2, i3, i4])...)
            dedup_key in seen_quads && continue
            push!(seen_quads, dedup_key)
            push!(xz_cells, (i1, i2, i3, i4))
        end
    end

    # Gather cells from all ranks to rank 0 — VTK is written by rank 0 only, but
    # each rank only has its local elements, so rank 0 would only see its own slice.
    local_cells_flat = Int64[v for t in xz_cells for v in t]
    all_cells_flat   = MPI.gather(local_cells_flat, comm)
    if rank == 0
        merged_flat = vcat(all_cells_flat...)
        ncells_all  = length(merged_flat) ÷ 4
        seen_global = Set{NTuple{4,Int64}}()
        xz_cells    = NTuple{4,Int64}[]
        for i in 1:ncells_all
            cell = (merged_flat[4i-3], merged_flat[4i-2], merged_flat[4i-1], merged_flat[4i])
            key  = tuple(sort!(Int64[cell...])...)
            key in seen_global && continue
            push!(seen_global, key)
            push!(xz_cells, cell)
        end
    else
        xz_cells = NTuple{4,Int64}[]
    end

    cs = LESCrossSection{T}(
        xz_coords, xz_groups, npts_per_xz,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nprofiles),
        zeros(T, nxz, nstress),   zeros(T, nxz, nstress),
        zeros(T, nprofiles),      zeros(T, nstress),
        0,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress),
        xz_cells,
        nothing,
        0,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress),
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress)
    )
    # cs.spec_cache = build_les_spectral_cache(cs, mesh, basis)
    return cs
end

"""
    compute_xz_cross_section!(cs, uaux, qe, ET, comm)

Two-pass y-averaging on the xz plane:
  Pass 1 — y-average of profiles via `user_les_profiles!` → `cs.global_mean`
  Pass 2 — y-average of stress products via `user_les_stress!` → `cs.global_stress`
"""
function compute_xz_cross_section!(cs, uaux, qe, ET, comm)
    nxz       = length(cs.xz_coords)
    nprofiles = size(cs.local_mean,   2)
    nstress   = size(cs.local_stress, 2)

    # Single pass: accumulate means + raw products
    fill!(cs.local_mean,   0.0)
    fill!(cs.local_stress, 0.0)
    for ixz in 1:nxz
        for ip in cs.xz_groups[ixz]
            user_les_profiles!(cs.prof_buf, cs.stress_buf, @view(uaux[ip,:]), @view(qe[ip,:]), ET)
            for k in 1:nprofiles
                cs.local_mean[ixz, k] += cs.prof_buf[k]
            end
            for k in 1:nstress
                cs.local_stress[ixz, k] += cs.stress_buf[k]
            end
        end
    end
    MPI.Allreduce!(cs.local_mean,   cs.global_mean,   MPI.SUM, comm)
    MPI.Allreduce!(cs.local_stress, cs.global_stress, MPI.SUM, comm)
    for ixz in 1:nxz
        npts = cs.npts_per_xz[ixz]
        for k in 1:nprofiles
            cs.global_mean[ixz, k] /= npts
        end
        for k in 1:nstress
            cs.global_stress[ixz, k] /= npts
        end
    end

    # Post-process: convert averaged raw products → fluctuation statistics in-place
    profp = cs.stress_buf  # scratch (length nstress)
    for ixz in 1:nxz
        means = @view cs.global_mean[ixz, :]
        prof  = @view cs.global_stress[ixz, :]
        user_les_stress!(profp, prof, means)
        for k in 1:nstress
            cs.global_stress[ixz, k] = profp[k]
        end
    end
end

function write_xz_cross_section_vtk(cs, params, time, iout, all_spectra=nothing, kappa_global=nothing)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank != 0 && return

    nprofiles       = size(cs.global_mean,   2)
    nstress         = size(cs.global_stress, 2)
    lesprofile_vars = params.inputs[:lesprofile_vars]
    lesstress_vars  = params.inputs[:lesstress_vars]
    nxz             = length(cs.xz_coords)

    # VTK_QUAD cells — connectivity built from element structure at init time,
    # so this works for both flat and terrain-following meshes.
    xz_cells_vtk = [MeshCell(VTKCellTypes.VTK_QUAD, [i1, i2, i3, i4])
                    for (i1, i2, i3, i4) in cs.xz_cells]

    # Point coordinates — xz plane, y = 0
    x_pts = [p[1] for p in cs.xz_coords]
    y_pts = zeros(Float64, nxz)
    z_pts = [p[2] for p in cs.xz_coords]

    # ---- instantaneous xz snapshot ----
    # fout_name = joinpath(params.inputs[:output_dir], @sprintf("les_xz_%06d", iout))
    # vtk = vtk_grid(fout_name, x_pts, y_pts, z_pts, xz_cells_vtk)
    # for k in 1:nprofiles
    #     vtk[lesprofile_vars[k], VTKPointData()] = cs.global_mean[:, k]
    # end
    # for k in 1:nstress
    #     vtk[lesstress_vars[k], VTKPointData()] = cs.global_stress[:, k]
    # end
    # outfiles = vtk_save(vtk)

    # pvd_path = joinpath(params.inputs[:output_dir], "les_xz.pvd")
    # if !isfile(pvd_path); init_pvd_file(pvd_path); end
    # append_pvd_entry(pvd_path, time, basename(outfiles[1]))

    # ---- time-averaged xz (running average, overwritten each call) ----
    ns = cs.n_samples
    fout_tavg = joinpath(params.inputs[:output_dir], "les_xz_tavg")
    vtk_tavg  = vtk_grid(fout_tavg, x_pts, y_pts, z_pts, xz_cells_vtk)
    for k in 1:nprofiles
        vtk_tavg[lesprofile_vars[k], VTKPointData()] = cs.sum_mean[:, k] ./ ns
    end
    for k in 1:nstress
        vtk_tavg[lesstress_vars[k], VTKPointData()] = cs.sum_stress[:, k] ./ ns
    end
    vtk_save(vtk_tavg)

    # ---- turbulent spectra (3-D grid: x × κ × z) ----
    isnothing(all_spectra) && return
    lesspectra_vars = get(params.inputs, :lesspectra_vars, String[])
    isempty(lesspectra_vars) && return

    nk    = length(kappa_global)
    n_pts = nxz * nk
    xs_pts = zeros(Float64, n_pts)
    yk_pts = zeros(Float64, n_pts)   # κ in the y slot
    zs_pts = zeros(Float64, n_pts)
    for ik in 1:nk, ixz in 1:nxz
        idx         = (ik-1)*nxz + ixz
        xs_pts[idx] = cs.xz_coords[ixz][1]
        yk_pts[idx] = kappa_global[ik]
        zs_pts[idx] = cs.xz_coords[ixz][2]
    end

    # VTK_HEXAHEDRON cells: each xz-quad extruded between κ-slices ik and ik+1
    spec_cells = MeshCell[]
    for ik in 1:nk-1
        for (i1, i2, i3, i4) in cs.xz_cells
            b1=(ik-1)*nxz+i1; b2=(ik-1)*nxz+i2; b3=(ik-1)*nxz+i3; b4=(ik-1)*nxz+i4
            t1= ik   *nxz+i1; t2= ik   *nxz+i2; t3= ik   *nxz+i3; t4= ik   *nxz+i4
            push!(spec_cells, MeshCell(VTKCellTypes.VTK_HEXAHEDRON,
                                       [b1, b2, b3, b4, t1, t2, t3, t4]))
        end
    end

    # fout_spec = joinpath(params.inputs[:output_dir], @sprintf("les_spectra_%06d", iout))
    # vtk_spec  = vtk_grid(fout_spec, xs_pts, yk_pts, zs_pts, spec_cells)
    # for (ivar, vname) in enumerate(lesspectra_vars)
    #     data = [all_spectra[ixz, ik, ivar] for ik in 1:nk for ixz in 1:nxz]
    #     vtk_spec[vname, VTKPointData()] = data
    # end
    # spec_outfiles = vtk_save(vtk_spec)

    # pvd_spec = joinpath(params.inputs[:output_dir], "les_spectra.pvd")
    # if !isfile(pvd_spec); init_pvd_file(pvd_spec); end
    # append_pvd_entry(pvd_spec, time, basename(spec_outfiles[1]))
end

# ================================================================================
# Finalize (Approach 1): write final time-and-space averages once at the end
# ================================================================================

"""
    les_finalize!(params, t)

Write the final time-and-space averaged 1D profiles and xz cross-section VTK
using all snapshots accumulated by `les_statistics`. Call once at the end of
the statistics window. No additional MPI needed: all sums were already
Allreduced during accumulation in `horizontal_mean!`.
"""
function les_finalize!(params, t)
    isnothing(params.les_stat_cache) && return

    cache = params.les_stat_cache
    ns    = cache.n_samples
    ns == 0 && return

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank != 0 && return

    nz        = length(cache.z_levels)
    nprofiles = size(cache.sum_global_sum, 2)
    nstress   = size(cache.sum_global_stress, 2)
    lesprofile_vars = params.inputs[:lesprofile_vars]
    lesstress_vars  = params.inputs[:lesstress_vars]

    outfile_les_tavg = joinpath(params.inputs[:output_dir], "les_statistics_tavg.dat")
    open(outfile_les_tavg, "w") do io
        print(io, "# time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
        for k in 1:nprofiles
            print(io, "  ", lesprofile_vars[k])
        end
        println(io)
        for iz in 1:nz
            @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
            for k in 1:nprofiles
                @printf(io, "  %.6e", cache.sum_global_sum[iz, k] / ns)
            end
            println(io)
        end
    end

    outfile_stress_tavg = joinpath(params.inputs[:output_dir], "les_stress_tavg.dat")
    open(outfile_stress_tavg, "w") do io
        print(io, "# time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
        for k in 1:nstress
            print(io, "  ", lesstress_vars[k])
        end
        println(io)
        for iz in 1:nz
            @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
            for k in 1:nstress
                @printf(io, "  %.6e", cache.sum_global_stress[iz, k] / ns)
            end
            println(io)
        end
    end

    # XZ cross-section final VTK
    cs = params.les_cross_section
    isnothing(cs) && return
    ns_cs = cs.n_samples
    ns_cs == 0 && return

    nprofiles_cs = size(cs.sum_mean, 2)
    nstress_cs   = size(cs.sum_stress, 2)
    nxz          = length(cs.xz_coords)
    xz_cells_vtk = [MeshCell(VTKCellTypes.VTK_QUAD, [i1, i2, i3, i4])
                    for (i1, i2, i3, i4) in cs.xz_cells]
    x_pts = [p[1] for p in cs.xz_coords]
    y_pts = zeros(Float64, nxz)
    z_pts = [p[2] for p in cs.xz_coords]

    fout_tavg = joinpath(params.inputs[:output_dir], "les_xz_tavg")
    vtk_tavg  = vtk_grid(fout_tavg, x_pts, y_pts, z_pts, xz_cells_vtk)
    for k in 1:nprofiles_cs
        vtk_tavg[lesprofile_vars[k], VTKPointData()] = cs.sum_mean[:, k] ./ ns_cs
    end
    for k in 1:nstress_cs
        vtk_tavg[lesstress_vars[k], VTKPointData()] = cs.sum_stress[:, k] ./ ns_cs
    end
    vtk_save(vtk_tavg)
end

# ================================================================================
# Online accumulation (Approach 2): no MPI per step, single Allreduce at end
# ================================================================================

"""
    les_accumulate_online!(u, params)

Accumulate local running sums for online recursive averaging (Approach 2).
Called every time step in the statistics window; performs NO MPI communication.
Call `les_finalize_online!` once at the end to Allreduce and write output.
"""
function les_accumulate_online!(u, params)
    isnothing(params.les_stat_cache) && return

    cache  = params.les_stat_cache
    cs     = params.les_cross_section
    mesh   = params.mesh
    npoin  = mesh.npoin
    neqs   = params.neqs
    ET     = params.SOL_VARS_TYPE
    uaux   = params.uaux
    qe     = params.qp.qe

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    local_online_sum    = cache.local_online_sum
    local_online_stress = cache.local_online_stress
    # 1D profile accumulation: sum over local owned points (no Allreduce)
    nz        = length(cache.z_levels)
    nprofiles = size(local_online_sum, 2)
    nstress   = size(local_online_stress, 2)
    means_buf = cache.prof
    prof_buf  = cache.stress_prof
    z_groups  = cache.z_groups

    for iz in 1:nz
        for ip in z_groups[iz]
            user_les_profiles!(means_buf, prof_buf, @view(uaux[ip,:]), @view(qe[ip,:]), ET)
            for k in 1:nprofiles
                local_online_sum[iz, k] += means_buf[k]
            end
            for k in 1:nstress
                local_online_stress[iz, k] += prof_buf[k]
            end
        end
    end
    cache.n_online_samples += 1

    # XZ cross-section accumulation: sum over local points (no Allreduce)
    isnothing(cs) && return

    cslocal_online_mean   = cs.local_online_mean
    cslocal_online_stress = cs.local_online_stress
    nxz          = length(cs.xz_coords)
    nprofiles_cs = size(cslocal_online_mean, 2)
    nstress_cs   = size(cslocal_online_stress, 2)
    csprof_buf     = cs.prof_buf
    csstress_buf   = cs.stress_buf
    csxz_groups    = cs.xz_groups

    for ixz in 1:nxz
        for ip in csxz_groups[ixz]
            user_les_profiles!(csprof_buf, csstress_buf, @view(uaux[ip,:]), @view(qe[ip,:]), ET)
            for k in 1:nprofiles_cs
                cslocal_online_mean[ixz, k] += csprof_buf[k]
            end
            for k in 1:nstress_cs
                cslocal_online_stress[ixz, k] += csstress_buf[k]
            end
        end
    end
    cs.n_online_samples += 1
end

"""
    les_finalize_online!(params, t)

Finalize online statistics: perform a single MPI.Allreduce across all ranks,
divide by (npts × n_online_samples) to get time-and-y averages, apply the
fluctuation decomposition `<a'b'> = <ab> - <a><b>` via `user_les_stress!`,
then write output files. Safe to call even if no steps were accumulated.
"""
function les_finalize_online!(params, t)
    isnothing(params.les_stat_cache) && return

    cache = params.les_stat_cache
    ns    = cache.n_online_samples
    ns == 0 && return

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    nz        = length(cache.z_levels)
    nprofiles = size(cache.local_online_sum, 2)
    nstress   = size(cache.local_online_stress, 2)

    # Single Allreduce for 1D profiles
    MPI.Allreduce!(cache.local_online_sum,    cache.global_online_sum,    MPI.SUM, comm)
    MPI.Allreduce!(cache.local_online_stress, cache.global_online_stress, MPI.SUM, comm)

    # Divide by (npts_per_z × n_online_samples) → time-and-y averaged means and raw products
    for iz in 1:nz
        denom = Float64(cache.npts_per_z[iz]) * ns
        for k in 1:nprofiles
            cache.global_online_sum[iz, k] /= denom
        end
        for k in 1:nstress
            cache.global_online_stress[iz, k] /= denom
        end
    end

    # Convert raw products → fluctuation statistics using time-and-y means
    profp = cache.stress_prof
    for iz in 1:nz
        means = @view cache.global_online_sum[iz, :]
        prof  = @view cache.global_online_stress[iz, :]
        user_les_stress!(profp, prof, means)
        for k in 1:nstress
            cache.global_online_stress[iz, k] = profp[k]
        end
    end

    if rank == 0
        lesprofile_vars = params.inputs[:lesprofile_vars]
        lesstress_vars  = params.inputs[:lesstress_vars]

        outfile = joinpath(params.inputs[:output_dir], "les_statistics_online.dat")
        open(outfile, "w") do io
            print(io, "# online_tavg  time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
            for k in 1:nprofiles; print(io, "  ", lesprofile_vars[k]); end
            println(io)
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nprofiles; @printf(io, "  %.6e", cache.global_online_sum[iz, k]); end
                println(io)
            end
        end

        outfile_stress = joinpath(params.inputs[:output_dir], "les_stress_online.dat")
        open(outfile_stress, "w") do io
            print(io, "# online_tavg  time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
            for k in 1:nstress; print(io, "  ", lesstress_vars[k]); end
            println(io)
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nstress; @printf(io, "  %.6e", cache.global_online_stress[iz, k]); end
                println(io)
            end
        end
    end

    # XZ cross-section finalization
    cs = params.les_cross_section
    isnothing(cs) && return

    ns_cs = cs.n_online_samples
    ns_cs == 0 && return

    nxz          = length(cs.xz_coords)
    nprofiles_cs = size(cs.local_online_mean, 2)
    nstress_cs   = size(cs.local_online_stress, 2)

    MPI.Allreduce!(cs.local_online_mean,   cs.global_online_mean,   MPI.SUM, comm)
    MPI.Allreduce!(cs.local_online_stress, cs.global_online_stress, MPI.SUM, comm)

    for ixz in 1:nxz
        denom = Float64(cs.npts_per_xz[ixz]) * ns_cs
        for k in 1:nprofiles_cs
            cs.global_online_mean[ixz, k] /= denom
        end
        for k in 1:nstress_cs
            cs.global_online_stress[ixz, k] /= denom
        end
    end

    profp_cs = cs.stress_buf
    for ixz in 1:nxz
        means = @view cs.global_online_mean[ixz, :]
        prof  = @view cs.global_online_stress[ixz, :]
        user_les_stress!(profp_cs, prof, means)
        for k in 1:nstress_cs
            cs.global_online_stress[ixz, k] = profp_cs[k]
        end
    end

    rank != 0 && return

    lesprofile_vars = params.inputs[:lesprofile_vars]
    lesstress_vars  = params.inputs[:lesstress_vars]
    xz_cells_vtk = [MeshCell(VTKCellTypes.VTK_QUAD, [i1, i2, i3, i4])
                    for (i1, i2, i3, i4) in cs.xz_cells]
    x_pts = [p[1] for p in cs.xz_coords]
    y_pts = zeros(Float64, nxz)
    z_pts = [p[2] for p in cs.xz_coords]

    fout = joinpath(params.inputs[:output_dir], "les_xz_online")
    vtk  = vtk_grid(fout, x_pts, y_pts, z_pts, xz_cells_vtk)
    for k in 1:nprofiles_cs
        vtk[lesprofile_vars[k], VTKPointData()] = cs.global_online_mean[:, k]
    end
    for k in 1:nstress_cs
        vtk[lesstress_vars[k], VTKPointData()] = cs.global_online_stress[:, k]
    end
    vtk_save(vtk)
end

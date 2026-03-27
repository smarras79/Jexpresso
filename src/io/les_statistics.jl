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
    # spectral cache (built once; nothing if spectra disabled)
    spec_cache   ::Union{Nothing, LESSpectralCache}
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

function les_statistics(u, params, time)

    isnothing(params.les_stat_cache) && return
    t = time
    iout  = findfirst(x -> AlmostEqual(x, time), params.inputs[:statistics_time])

    mesh      = params.mesh
    npoin     = mesh.npoin
    neqs      = params.neqs
    cache     = params.les_stat_cache
    nz        = length(cache.z_levels)
    nprofiles = size(cache.global_sum, 2)
    ngl       = mesh.ngl
    ET        = params.inputs[:SOL_VARS_TYPE]

    wθ     = params.WM.wθ
    wqv    = params.WM.wqv

    nfaces_bdy       = mesh.nfaces_bdy
    bdy_face_type    = mesh.bdy_face_type
    poin_in_bdy_face = mesh.poin_in_bdy_face

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

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
        # all_spectra, kappa_global = compute_les_spectra!(cs, uaux, qe, mesh, ET, comm, params)
        # all_spectra, kappa_global = nothing
        write_xz_cross_section_vtk(cs, params, t, iout)
        # write_xz_cross_section_vtk(cs, params, t, iout, all_spectra, kappa_global)
    end

    # MOST surface fluxes
    # PhysConst   = PhysicalConst{Float64}()
    # wθ_sum      = 0.0
    # wqv_sum     = 0.0
    # cnt         = 0.0
    # for iface = 1:nfaces_bdy
    #     if bdy_face_type[iface] == "MOST"
    #         for i = 1:ngl, j = 1:ngl
    #             ip       = poin_in_bdy_face[iface,i,j]
    #             ρ        = uaux[ip, 1]
    #             wθ_sum  += ρ * PhysConst.cp * wθ[iface,i,j,1]
    #             wqv_sum += ρ * PhysConst.Lc * wqv[iface,i,j,1]
    #             cnt     += 1.0
    #         end
    #     end
    # end
    # most_sum_local  = [wθ_sum, wqv_sum, cnt]
    # most_sum_global = zeros(3)
    # MPI.Allreduce!(most_sum_local, most_sum_global, MPI.SUM, comm)

    if rank == 0
        lesprofile_vars = params.inputs[:lesprofile_vars]

        outfile_les = joinpath(params.inputs[:output_dir], "les_statistics.dat")
        if !isfile(outfile_les)
            open(outfile_les, "w") do io
                print(io, "# time  z")
                for k in 1:nprofiles
                    print(io, "  ", lesprofile_vars[k])
                end
                println(io)
            end
        end
        open(outfile_les, "a") do io
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nprofiles
                    @printf(io, "  %.6e", cache.global_sum[iz, k])
                end
                println(io)
            end
            println(io)
        end

        nstress          = size(cache.global_stress, 2)
        lesstress_vars   = params.inputs[:lesstress_vars]
        outfile_stress   = joinpath(params.inputs[:output_dir], "les_stress.dat")
        if !isfile(outfile_stress)
            open(outfile_stress, "w") do io
                print(io, "# time  z")
                for k in 1:nstress
                    print(io, "  ", lesstress_vars[k])
                end
                println(io)
            end
        end
        open(outfile_stress, "a") do io
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nstress
                    @printf(io, "  %.6e", cache.global_stress[iz, k])
                end
                println(io)
            end
            println(io)
        end

        # outfile_most = joinpath(params.inputs[:output_dir], "MOST_statistics.dat")
        # if !isfile(outfile_most)
        #     open(outfile_most, "w") do io
        #         println(io, "# time  LHF  SHF")
        #     end
        # end
        # open(outfile_most, "a") do io
        #     @printf(io, "%.6e  %.6e  %.6e\n", t, most_sum_global[2]/most_sum_global[3], most_sum_global[1]/most_sum_global[3])
        # end

        # Time-averaged 1D profiles (running average over all accumulated snapshots)
        ns = cache.n_samples
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

        nstress_tavg   = size(cache.global_stress, 2)
        lesstress_vars = params.inputs[:lesstress_vars]
        outfile_stress_tavg = joinpath(params.inputs[:output_dir], "les_stress_tavg.dat")
        open(outfile_stress_tavg, "w") do io
            print(io, "# time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
            for k in 1:nstress_tavg
                print(io, "  ", lesstress_vars[k])
            end
            println(io)
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nstress_tavg
                    @printf(io, "  %.6e", cache.sum_global_stress[iz, k] / ns)
                end
                println(io)
            end
        end
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

    cs = LESCrossSection{T}(
        xz_coords, xz_groups, npts_per_xz,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nprofiles),
        zeros(T, nxz, nstress),   zeros(T, nxz, nstress),
        zeros(T, nprofiles),      zeros(T, nstress),
        0,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress),
        nothing
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

    # Unique sorted x and z values — shared by xz and spectra grids
    xs = sort(unique(p[1] for p in cs.xz_coords))
    zs = sort(unique(p[2] for p in cs.xz_coords))
    nx = length(xs)
    nz = length(zs)

    coord_to_idx = Dict{Tuple{Float64,Float64}, Int}(p => i for (i, p) in enumerate(cs.xz_coords))

    # VTK_QUAD cells on the (x,z) plane
    xz_cells = MeshCell[]
    for ix in 1:nx-1, iz in 1:nz-1
        i1 = get(coord_to_idx, (xs[ix],   zs[iz]  ), 0)
        i2 = get(coord_to_idx, (xs[ix+1], zs[iz]  ), 0)
        i3 = get(coord_to_idx, (xs[ix+1], zs[iz+1]), 0)
        i4 = get(coord_to_idx, (xs[ix],   zs[iz+1]), 0)
        i1 > 0 && i2 > 0 && i3 > 0 && i4 > 0 &&
            push!(xz_cells, MeshCell(VTKCellTypes.VTK_QUAD, [i1, i2, i3, i4]))
    end

    # Point coordinates — xz plane, y = 0
    x_pts = [p[1] for p in cs.xz_coords]
    y_pts = zeros(Float64, nxz)
    z_pts = [p[2] for p in cs.xz_coords]

    # ---- instantaneous xz snapshot ----
    fout_name = joinpath(params.inputs[:output_dir], @sprintf("les_xz_%06d", iout))
    vtk = vtk_grid(fout_name, x_pts, y_pts, z_pts, xz_cells)
    for k in 1:nprofiles
        vtk[lesprofile_vars[k], VTKPointData()] = cs.global_mean[:, k]
    end
    for k in 1:nstress
        vtk[lesstress_vars[k], VTKPointData()] = cs.global_stress[:, k]
    end
    outfiles = vtk_save(vtk)

    pvd_path = joinpath(params.inputs[:output_dir], "les_xz.pvd")
    if !isfile(pvd_path); init_pvd_file(pvd_path); end
    append_pvd_entry(pvd_path, time, basename(outfiles[1]))

    # ---- time-averaged xz (running average, overwritten each call) ----
    ns = cs.n_samples
    fout_tavg = joinpath(params.inputs[:output_dir], "les_xz_tavg")
    vtk_tavg  = vtk_grid(fout_tavg, x_pts, y_pts, z_pts, xz_cells)
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
        for ix in 1:nx-1, iz in 1:nz-1
            i1 = get(coord_to_idx, (xs[ix],   zs[iz]  ), 0)
            i2 = get(coord_to_idx, (xs[ix+1], zs[iz]  ), 0)
            i3 = get(coord_to_idx, (xs[ix+1], zs[iz+1]), 0)
            i4 = get(coord_to_idx, (xs[ix],   zs[iz+1]), 0)
            (i1>0 && i2>0 && i3>0 && i4>0) || continue
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

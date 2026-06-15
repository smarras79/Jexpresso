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
    # scratch buffers for fill_sgs_cache! — pre-allocated to avoid per-call heap allocation
    sgs_uprim  ::Array{Float64,4} = zeros(Float64, 1, 1, 1, 1)
    sgs_ncount ::Vector{Int32}    = Int32[]
end


# ---- Bottom-wall u* cache ----
mutable struct LESBottomCache
    x_coords           ::Vector{Float64}       # nx sorted x-positions at bottom wall
    wall_groups        ::Vector{Vector{Int64}} # local owned ip1 (first-wall-node) indices per x
    z_inside           ::Dict{Int64,Float64}   # ip1 → normal distance from surface (|Δcoords·n̂|)
    npts_per_x         ::Vector{Int64}         # global y-count per x
    z0_m               ::Float64               # aerodynamic roughness length [m]
    κ_v                ::Float64               # von Kármán constant
    # explicit-path accumulation (Allreduce per step, consistent with horizontal_mean!)
    sum_ustar          ::Vector{Float64}       # time-accumulated y-mean u*(x), rank-0 valid
    n_samples          ::Int
    # online-path accumulation (local only; Allreduce at finalization)
    local_ustar_online ::Vector{Float64}       # local accumulated sum over all steps and y-nodes
    n_online_samples   ::Int
    # scratch buffers — pre-allocated to avoid per-call heap allocation
    local_step         ::Vector{Float64}
    global_step        ::Vector{Float64}
end

"""
    build_les_bottom_cache(mesh, inputs)

Build `LESBottomCache` by scanning boundary faces tagged "MOST" to identify
first-wall-layer nodes above the surface, exactly as `BCs.jl` does.

`z_inside` per node mirrors the `z_inside` in `BCs.jl` (line 761):
    z_inside = |Δcoords · n̂|   (normal distance from surface node to first-wall node)

Returns `nothing` if `:lwall_model` is false or no MOST faces are found.
"""
function build_les_bottom_cache(mesh, metrics, inputs)
    get(inputs, :lwall_model, false) || return nothing

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)

    ifirst_wall_node = get(inputs, :ifirst_wall_node_index, 2)::Int
    z0_m = Float64(get(inputs, :z0_m, 0.1))
    PhysConst_b = PhysicalConst{Float64}()
    κ_v = PhysConst_b.karman

    coords       = Array(mesh.coords)      # (npoin, 3)
    gip2owner    = Array(mesh.gip2owner)
    connijk_arr  = Array(mesh.connijk)
    poin_bdy     = Array(mesh.poin_in_bdy_face)
    face_in_elem = Array(mesh.bdy_face_in_elem)
    face_type    = mesh.bdy_face_type
    # For Laguerre meshes sem.metrics is a tuple; take first element
    met          = metrics isa Tuple ? metrics[1] : metrics
    metrics_nx   = Array(met.nx)           # (nfaces_bdy, ngl, ngl)
    metrics_ny   = Array(met.ny)
    metrics_nz   = Array(met.nz)
    nfaces_bdy   = mesh.nfaces_bdy
    ngl          = mesh.ngl

    xr = round.(coords[:, 1]; digits=6)

    # Scan MOST faces; mirror BCs.jl exactly:
    #   ip_sfc = bottom surface node (connijk[e,i,j,1])
    #   ip1    = first-wall-layer node (connijk[e,i,j,ifirst_wall_node])
    #   z_inside[ip1] = |Δcoords · n̂|  (same formula as BCs.jl line 761)
    x_to_ip1   = Dict{Float64, Vector{Int64}}()
    ip1_zin    = Dict{Int64,   Float64}()      # ip1 → z_inside

    for iface in 1:nfaces_bdy
        face_type[iface] == "MOST" || continue
        e = face_in_elem[iface]
        for i in 1:ngl, j in 1:ngl
            ip_sfc = poin_bdy[iface, i, j]
            gip2owner[ip_sfc] == rank || continue
            ip1 = connijk_arr[e, i, j, ifirst_wall_node]
            nx_f = metrics_nx[iface, i, j]
            ny_f = metrics_ny[iface, i, j]
            nz_f = metrics_nz[iface, i, j]
            Δx = coords[ip1, 1] - coords[ip_sfc, 1]
            Δy = coords[ip1, 2] - coords[ip_sfc, 2]
            Δz = coords[ip1, 3] - coords[ip_sfc, 3]
            z_in = abs(Δx*nx_f + Δy*ny_f + Δz*nz_f)
            ip1_zin[ip1] = z_in
            xk = xr[ip_sfc]
            push!(get!(x_to_ip1, xk, Int64[]), ip1)
        end
    end

    # Gather unique x-coords globally
    local_xs  = sort(collect(keys(x_to_ip1)))
    all_xs    = MPI.gather(local_xs, comm)
    if rank == 0
        x_flat = sort(unique(vcat(all_xs...)))
    else
        x_flat = Float64[]
    end
    nx_global = MPI.bcast(rank == 0 ? length(x_flat) : 0, 0, comm)
    rank != 0 && resize!(x_flat, nx_global)
    MPI.Bcast!(x_flat, 0, comm)

    nx       = nx_global
    x_coords = x_flat

    wall_groups  = Vector{Vector{Int64}}(undef, nx)
    local_counts = zeros(Int64, nx)
    for ix in 1:nx
        ids = get(x_to_ip1, x_coords[ix], Int64[])
        wall_groups[ix]  = ids
        local_counts[ix] = length(ids)
    end
    npts_per_x = zeros(Int64, nx)
    MPI.Allreduce!(local_counts, npts_per_x, MPI.SUM, comm)

    return LESBottomCache(
        x_coords, wall_groups, ip1_zin, npts_per_x,
        z0_m, κ_v,
        zeros(Float64, nx), 0,
        zeros(Float64, nx), 0,
        zeros(Float64, nx), zeros(Float64, nx)
    )
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
    # temporal accumulation of spectra (rank-0 only; nothing until spec_cache is built)
    n_spec_samples::Int
    sum_spectra  ::Union{Nothing, Array{Float64,3}}  # (nxz, nk, 4), rank-0 only
    kappa_global ::Union{Nothing, Vector{Float64}}   # (nk,), filled on first call
    # online accumulation (Approach 2): local running sums per time step, Allreduce only at end
    n_online_samples     ::Int
    local_online_mean    ::Matrix{T}
    local_online_stress  ::Matrix{T}
    global_online_mean   ::Matrix{T}
    global_online_stress ::Matrix{T}
end

function build_les_stat_cache(mesh, nprofiles::Int, nstress::Int, T, backend)
    nprofiles == 0 && return nothing

    comm = get_mpi_comm()

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
    cache.sgs_uprim  = zeros(Float64, mesh.ngl, mesh.ngl, mesh.ngl, 5)
    cache.sgs_ncount = zeros(Int32, mesh.npoin)
    return cache
end

"""
    fill_sgs_cache!(params)

Loop over all elements, compute velocity and temperature gradients at each GLL
point, call `compute_sij_and_mu_turb` (dispatches on SMAG/VREM), and fill
`params.sgs_stress[ip, 1:12]`:
  [1:6]  → kinematic SGS momentum stress components −2ν_t Sij  (m²/s²)
  [7]    → 0 (tptp_sfs, not modeled in Smagorinsky/Vreman)
  [8:10] → SGS heat-flux components −κ_t ∂θ/∂xᵢ  (K·m/s)
  [11]   → TKE dissipation rate ε = 2ν_eff Sij Sij  (m²/s³)
  [12]   → scalar-variance dissipation rate εθ = κ_eff |∇θ|²  (K²/s)
Called once per statistics output step (Approach 1), not every time step.
"""
function fill_sgs_cache!(params)
    isnothing(params.les_stat_cache) && return
    VT = params.VT
    (VT isa SMAG || VT isa VREM) || return

    sgs_stress = params.sgs_stress
    uaux       = params.uaux
    qe         = params.qp.qe
    mesh       = params.mesh
    connijk    = mesh.connijk
    ngl        = mesh.ngl
    nelem      = mesh.nelem
    npoin      = mesh.npoin
    dψ         = params.basis.dψ
    dξdx = params.metrics.dξdx;  dξdy = params.metrics.dξdy;  dξdz = params.metrics.dξdz
    dηdx = params.metrics.dηdx;  dηdy = params.metrics.dηdy;  dηdz = params.metrics.dηdz
    dζdx = params.metrics.dζdx;  dζdy = params.metrics.dζdy;  dζdz = params.metrics.dζdz
    Δ        = Float64(mesh.Δeffective_l)
    ad_lvl   = mesh.ad_lvl
    PhysConst = PhysicalConst{Float64}()
    Pr_t = PhysConst.Pr_t;  μ_mol = PhysConst.μ_mol;  κ_mol = PhysConst.κ_mol
    ET   = params.SOL_VARS_TYPE

    uprim  = params.les_stat_cache.sgs_uprim
    # Track how many elements contribute to each node's gradient estimate.
    # At shared boundary nodes, adjacent elements give different gradient values
    # (CG continuity guarantees continuous solution but not continuous derivative).
    # We average all element contributions so that boundary nodes get the
    # mean gradient from both sides rather than the last-writer's value.
    ncount = params.les_stat_cache.sgs_ncount

    _fill_sgs_inner!(sgs_stress, ncount, uaux, qe, connijk, uprim,
                     dψ, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                     ad_lvl, Δ, PhysConst, Pr_t, μ_mol, κ_mol,
                     ngl, nelem, npoin, VT, ET)
end

function _fill_sgs_inner!(sgs_stress, ncount, uaux, qe, connijk, uprim,
                           dψ, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                           ad_lvl, Δ, PhysConst, Pr_t, μ_mol, κ_mol,
                           ngl, nelem, npoin, VT, ET)
    fill!(sgs_stress, 0.0)
    fill!(ncount, Int32(0))

    for iel in 1:nelem
        Δ2 = (ldexp(Δ, -ad_lvl[iel]))^2

        for k in 1:ngl, j in 1:ngl, i in 1:ngl
            ip = connijk[iel,i,j,k]
            if ET == PERT()
                ρ = uaux[ip,1] + qe[ip,1]
                uprim[i,j,k,1] = ρ
                uprim[i,j,k,2] = uaux[ip,2] / ρ
                uprim[i,j,k,3] = uaux[ip,3] / ρ
                uprim[i,j,k,4] = uaux[ip,4] / ρ
                uprim[i,j,k,5] = (uaux[ip,5]+qe[ip,5])/ρ - qe[ip,5]/qe[ip,1]
            else
                ρ = uaux[ip,1]
                uprim[i,j,k,1] = ρ
                uprim[i,j,k,2] = uaux[ip,2] / ρ
                uprim[i,j,k,3] = uaux[ip,3] / ρ
                uprim[i,j,k,4] = uaux[ip,4] / ρ
                uprim[i,j,k,5] = uaux[ip,5] / ρ
            end
        end

        for k in 1:ngl, j in 1:ngl, i in 1:ngl
            ip = connijk[iel,i,j,k]
            ρ  = uprim[i,j,k,1]

            dudξ=0.0; dudη=0.0; dudζ=0.0
            dvdξ=0.0; dvdη=0.0; dvdζ=0.0
            dwdξ=0.0; dwdη=0.0; dwdζ=0.0
            dθdξ=0.0; dθdη=0.0; dθdζ=0.0
            for ii in 1:ngl
                dudξ += dψ[ii,i]*uprim[ii,j,k,2];  dudη += dψ[ii,j]*uprim[i,ii,k,2];  dudζ += dψ[ii,k]*uprim[i,j,ii,2]
                dvdξ += dψ[ii,i]*uprim[ii,j,k,3];  dvdη += dψ[ii,j]*uprim[i,ii,k,3];  dvdζ += dψ[ii,k]*uprim[i,j,ii,3]
                dwdξ += dψ[ii,i]*uprim[ii,j,k,4];  dwdη += dψ[ii,j]*uprim[i,ii,k,4];  dwdζ += dψ[ii,k]*uprim[i,j,ii,4]
                dθdξ += dψ[ii,i]*uprim[ii,j,k,5];  dθdη += dψ[ii,j]*uprim[i,ii,k,5];  dθdζ += dψ[ii,k]*uprim[i,j,ii,5]
            end

            Jdξdx=dξdx[iel,i,j,k]; Jdξdy=dξdy[iel,i,j,k]; Jdξdz=dξdz[iel,i,j,k]
            Jdηdx=dηdx[iel,i,j,k]; Jdηdy=dηdy[iel,i,j,k]; Jdηdz=dηdz[iel,i,j,k]
            Jdζdx=dζdx[iel,i,j,k]; Jdζdy=dζdy[iel,i,j,k]; Jdζdz=dζdz[iel,i,j,k]

            dudx = dudξ*Jdξdx + dudη*Jdηdx + dudζ*Jdζdx
            dudy = dudξ*Jdξdy + dudη*Jdηdy + dudζ*Jdζdy
            dudz = dudξ*Jdξdz + dudη*Jdηdz + dudζ*Jdζdz
            dvdx = dvdξ*Jdξdx + dvdη*Jdηdx + dvdζ*Jdζdx
            dvdy = dvdξ*Jdξdy + dvdη*Jdηdy + dvdζ*Jdζdy
            dvdz = dvdξ*Jdξdz + dvdη*Jdηdz + dvdζ*Jdζdz
            dwdx = dwdξ*Jdξdx + dwdη*Jdηdx + dwdζ*Jdζdx
            dwdy = dwdξ*Jdξdy + dwdη*Jdηdy + dwdζ*Jdζdy
            dwdz = dwdξ*Jdξdz + dwdη*Jdηdz + dwdζ*Jdζdz
            dθdx = dθdξ*Jdξdx + dθdη*Jdηdx + dθdζ*Jdζdx
            dθdy = dθdξ*Jdξdy + dθdη*Jdηdy + dθdζ*Jdζdy
            dθdz = dθdξ*Jdξdz + dθdη*Jdηdz + dθdζ*Jdζdz

            μ_t, S11, S22, S33, S12, S13, S23, SijSij =
                compute_sij_and_mu_turb(ρ, dudx, dudy, dudz,
                                        dvdx, dvdy, dvdz,
                                        dwdx, dwdy, dwdz,
                                        PhysConst, Δ2, VT)

            ν_t   = μ_t / ρ
            κ_t   = μ_t / (ρ * Pr_t)
            ν_eff = μ_mol/ρ + ν_t
            κ_eff = κ_mol + κ_t

            sgs_stress[ip, 1]  += -2 * ν_t * S11
            sgs_stress[ip, 2]  += -2 * ν_t * S12
            sgs_stress[ip, 3]  += -2 * ν_t * S13
            sgs_stress[ip, 4]  += -2 * ν_t * S22
            sgs_stress[ip, 5]  += -2 * ν_t * S23
            sgs_stress[ip, 6]  += -2 * ν_t * S33
            # sgs_stress[ip, 7] stays 0 (θθ sfs not modeled in Smagorinsky/Vreman)
            sgs_stress[ip, 8]  += -κ_t * dθdx
            sgs_stress[ip, 9]  += -κ_t * dθdy
            sgs_stress[ip, 10] += -κ_t * dθdz
            sgs_stress[ip, 11] += 2 * ν_eff * SijSij
            sgs_stress[ip, 12] += κ_eff * (dθdx*dθdx + dθdy*dθdy + dθdz*dθdz)
            ncount[ip] += one(Int32)
        end
    end

    # Average contributions at shared element-boundary nodes so that
    # gradient-derived SGS quantities are the mean of all adjacent elements
    # rather than the last writer's value.
    @inbounds for ip in 1:npoin
        n = ncount[ip]
        n < 2 && continue
        inv_n = 1.0 / n
        for c in 1:12
            sgs_stress[ip, c] *= inv_n
        end
    end
end

"""
    horizontal_mean!(cache, uaux, qe, sgs_stress, ET, comm)

For each z-level, call `user_les_profiles!` at every local point to accumulate
all profile sums into `cache.local_sum` (nz × nprofiles), then perform a single
`MPI.Allreduce!` and divide by `npts_per_z`. Results are stored in `cache.global_sum`.
SFS stress components (prof[11:20]) and dissipation rates (prof[24:25]) are
overridden from the pre-computed `sgs_stress` array after the user callback.
"""
function horizontal_mean!(cache, uaux, qe, sgs_stress, ET, comm)
    nz        = length(cache.z_levels)
    nprofiles = size(cache.local_sum,    2)
    nstress   = size(cache.local_stress, 2)
    fill!(cache.local_sum,    0.0)
    fill!(cache.local_stress, 0.0)

    means_buf = cache.prof        # scratch: means (length nprofiles)
    prof_buf  = cache.stress_prof # scratch: raw products (length nstress)

    for iz in 1:nz
        for ip in cache.z_groups[iz]
            user_les_profiles!(means_buf, prof_buf, @view(uaux[ip,:]), @view(qe[ip,:]), @view(sgs_stress[ip,:]), ET)
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
    ET    = params.SOL_VARS_TYPE
    comm  = get_mpi_comm()

    uaux = params.uaux
    qe   = params.qp.qe
    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    # Fill SGS stress/dissipation cache (gradient computation at statistics time)
    fill_sgs_cache!(params)
    sgs_stress = params.sgs_stress

    # Single pass: accumulate means + raw products, two Allreduces
    horizontal_mean!(cache, uaux, qe, sgs_stress, ET, comm)
    # Accumulate raw products (NOT yet fluctuations) so that Reynolds decomposition
    # can be applied once at finalization using the long-time mean as reference.
    # Applying decomposition per-snapshot and then time-averaging gives
    # <u²>_{t,y} - <<u>_y²>_t  ≠  <u²>_{t,y} - <u>_{t,y}²  (wrong for coarse meshes).
    cache.sum_global_sum    .+= cache.global_sum
    cache.sum_global_stress .+= cache.global_stress
    cache.n_samples         += 1

    # XZ cross section (y-averaged 2D fields)
    # Bottom-wall u*(x): replicate BCs.jl/CM_MOST! — wall-normal velocity is projected
    # out before computing speed, then u* = κ_v * |u_tangential| / ln(z_inside/z0_m).
    # For flat terrain (n = (0,0,1)) the tangential speed reduces to sqrt(u²+v²).
    bc = params.les_bottom_cache
    if !isnothing(bc)
        nx_bc       = length(bc.x_coords)
        local_step  = bc.local_step
        global_step = bc.global_step
        fill!(local_step,  0.0)
        fill!(global_step, 0.0)
        for ix in 1:nx_bc
            for ip1 in bc.wall_groups[ix]
                z_in = get(bc.z_inside, ip1, bc.z0_m * 10.0)
                z_in <= bc.z0_m && continue
                if ET == PERT()
                    ρ = uaux[ip1,1] + qe[ip1,1]
                    u = (uaux[ip1,2] + qe[ip1,2]) / ρ
                    v = (uaux[ip1,3] + qe[ip1,3]) / ρ
                else
                    ρ = uaux[ip1,1]
                    u = uaux[ip1,2] / ρ
                    v = uaux[ip1,3] / ρ
                end
                local_step[ix] += bc.κ_v * sqrt(u*u + v*v) / log(z_in / bc.z0_m)
            end
        end
        MPI.Allreduce!(local_step, global_step, MPI.SUM, comm)
        if MPI.Comm_rank(comm) == 0
            for ix in 1:nx_bc
                npts = bc.npts_per_x[ix]
                npts > 0 && (bc.sum_ustar[ix] += global_step[ix] / npts)
            end
        end
        bc.n_samples += 1
    end

    cs = params.les_cross_section
    if !isnothing(cs)
        compute_xz_cross_section!(cs, uaux, qe, sgs_stress, ET, comm)
        # Temporal accumulation for xz cross-sections
        cs.sum_mean   .+= cs.global_mean
        cs.sum_stress .+= cs.global_stress
        cs.n_samples  += 1

        # Spectra disabled — uncomment when FFTW is loaded
        # spec_result = compute_les_spectra!(cs, uaux, qe, mesh, ET, comm, params)
        # if spec_result isa Tuple && !isnothing(spec_result[1]) && !isnothing(cs.sum_spectra)
        #     all_spec, kappa = spec_result
        #     cs.sum_spectra .+= all_spec
        #     cs.n_spec_samples == 0 && (cs.kappa_global .= kappa)
        #     cs.n_spec_samples += 1
        # end
    end
end

# ================================================================================
# XZ Cross Section (y-averaged 2D fields)
# ================================================================================

"""
    build_les_spectral_cache(cs, mesh)

Build a `LESSpectralCache`:
  - Uses global point count from `cs.npts_per_xz` (set by Allreduce in
    `build_les_cross_section`) to compute `nel_y` and `N_unif` correctly
    for any MPI decomposition in y.
  - Precomputes the Lagrange interpolation matrix from GLL points to
    equally-spaced points (dy = dy_el/(ngl-1)) within each element.
"""
function build_les_spectral_cache(cs::LESCrossSection, mesh)
    ngl = mesh.ngl

    # Use global count — npts_per_xz is already Allreduced across all ranks
    Ny_global = maximum(cs.npts_per_xz)
    nel_y  = (Ny_global - 1) ÷ (ngl - 1)
    N_unif = nel_y * (ngl - 1) + 1
    nk     = N_unif ÷ 2 + 1

    # GLL reference points ξ ∈ [-1,1]: St_Lagrange only stores ψ/dψ, not ξ.
    # Recompute the LGL nodes from ngl — same nodes used to build the basis.
    ξ_gll = Array(basis_structs_ξ_ω!(LGL(), ngl - 1, CPU()).ξ)

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


function build_les_cross_section(mesh, nprofiles::Int, nstress::Int, T)
    nprofiles == 0 && return nothing

    comm = get_mpi_comm()
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
        0, nothing, nothing,
        0,
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress),
        zeros(T, nxz, nprofiles), zeros(T, nxz, nstress)
    )
    cs.spec_cache = build_les_spectral_cache(cs, mesh)
    if !isnothing(cs.spec_cache) && MPI.Comm_rank(get_mpi_comm()) == 0
        nk = cs.spec_cache.nk
        cs.sum_spectra  = zeros(Float64, nxz, nk, 4)
        cs.kappa_global = zeros(Float64, nk)
    end
    return cs
end

"""
    compute_xz_cross_section!(cs, uaux, qe, sgs_stress, ET, comm)

Two-pass y-averaging on the xz plane:
  Pass 1 — y-average of profiles via `user_les_profiles!` → `cs.global_mean`
  Pass 2 — y-average of stress products via `user_les_stress!` → `cs.global_stress`
SFS stress components (stress_buf[11:20]) and dissipation rates ([24:25]) are
overridden from `sgs_stress` after the user callback.
"""
function compute_xz_cross_section!(cs, uaux, qe, sgs_stress, ET, comm)
    nxz       = length(cs.xz_coords)
    nprofiles = size(cs.local_mean,   2)
    nstress   = size(cs.local_stress, 2)

    # Single pass: accumulate means + raw products
    fill!(cs.local_mean,   0.0)
    fill!(cs.local_stress, 0.0)
    for ixz in 1:nxz
        for ip in cs.xz_groups[ixz]
            user_les_profiles!(cs.prof_buf, cs.stress_buf, @view(uaux[ip,:]), @view(qe[ip,:]), @view(sgs_stress[ip,:]), ET)
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

    # Do NOT apply user_les_stress! here: raw products are left in cs.global_stress
    # so that les_finalize! can apply Reynolds decomposition once using the
    # long-time mean (sum_mean/ns) rather than per-snapshot y-means.
end

function write_xz_cross_section_vtk(cs, params, time, iout, all_spectra=nothing, kappa_global=nothing)
    comm = get_mpi_comm()
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
    # Apply Reynolds decomposition using long-time y-mean to avoid per-snapshot
    # y-mean noise that creates element-boundary artifacts.
    ns = cs.n_samples
    profp = zeros(Float64, nstress)
    stress_xz = zeros(Float64, nxz, nstress)
    for ixz in 1:nxz
        means_tavg = cs.sum_mean[ixz, :] ./ ns
        prof_raw   = cs.sum_stress[ixz, :] ./ ns
        user_les_stress!(profp, prof_raw, means_tavg)
        stress_xz[ixz, :] .= profp
    end
    fout_tavg = joinpath(params.inputs[:output_dir], "les_xz_tavg")
    vtk_tavg  = vtk_grid(fout_tavg, x_pts, y_pts, z_pts, xz_cells_vtk)
    for k in 1:nprofiles
        vtk_tavg[lesprofile_vars[k], VTKPointData()] = cs.sum_mean[:, k] ./ ns
    end
    for k in 1:nstress
        vtk_tavg[lesstress_vars[k], VTKPointData()] = stress_xz[:, k]
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

    comm = get_mpi_comm()
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

    # Apply Reynolds decomposition using time-averaged means as reference.
    # sum_global_stress / ns = <u²>_{t,x,y}  (time-mean of raw products)
    # sum_global_sum   / ns = <u>_{t,x,y}    (time-mean)
    # Correct: <u'u'> = <u²>_{t,x,y} - <u>_{t,x,y}²
    let profp = zeros(Float64, nstress)
        stress_out = zeros(Float64, nz, nstress)
        for iz in 1:nz
            means_tavg = cache.sum_global_sum[iz, :] ./ ns
            prof_raw   = cache.sum_global_stress[iz, :] ./ ns
            user_les_stress!(profp, prof_raw, means_tavg)
            stress_out[iz, :] .= profp
        end
        outfile_stress_tavg = joinpath(params.inputs[:output_dir], "les_stress_tavg.dat")
        open(outfile_stress_tavg, "w") do io
            print(io, "# time_end=", @sprintf("%.6e", t), "  n_samples=", ns, "  z")
            for k in 1:nstress; print(io, "  ", lesstress_vars[k]); end
            println(io)
            for iz in 1:nz
                @printf(io, "%.6e  %.6e", t, cache.z_levels[iz])
                for k in 1:nstress; @printf(io, "  %.6e", stress_out[iz, k]); end
                println(io)
            end
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

    # Apply Reynolds decomposition once using long-time y-mean as reference.
    # Correct: <u'u'>(x,z) = <u²>_{t,y}(x,z) - <u>_{t,y}(x,z)²
    let profp_cs = zeros(Float64, nstress_cs), stress_xz = zeros(Float64, nxz, nstress_cs)
        for ixz in 1:nxz
            means_tavg = cs.sum_mean[ixz, :] ./ ns_cs
            prof_raw   = cs.sum_stress[ixz, :] ./ ns_cs
            user_les_stress!(profp_cs, prof_raw, means_tavg)
            stress_xz[ixz, :] .= profp_cs
        end
        fout_tavg = joinpath(params.inputs[:output_dir], "les_xz_tavg")
        vtk_tavg  = vtk_grid(fout_tavg, x_pts, y_pts, z_pts, xz_cells_vtk)
        for k in 1:nprofiles_cs
            vtk_tavg[lesprofile_vars[k], VTKPointData()] = cs.sum_mean[:, k] ./ ns_cs
        end
        for k in 1:nstress_cs
            vtk_tavg[lesstress_vars[k], VTKPointData()] = stress_xz[:, k]
        end
        vtk_save(vtk_tavg)
    end

    # ---- time-averaged spanwise spectra: les_xz_spectra_tavg.vtu ----
    ns_spec = cs.n_spec_samples
    if ns_spec > 0 && !isnothing(cs.sum_spectra) && !isnothing(cs.kappa_global)
        lesspectra_vars = get(params.inputs, :lesspectra_vars, String[])
        if !isempty(lesspectra_vars)
            nk       = length(cs.kappa_global)
            n_pts    = nxz * nk
            xs_pts   = zeros(Float64, n_pts)
            yk_pts   = zeros(Float64, n_pts)
            zs_pts   = zeros(Float64, n_pts)
            for ik in 1:nk, ixz in 1:nxz
                idx         = (ik-1)*nxz + ixz
                xs_pts[idx] = cs.xz_coords[ixz][1]
                yk_pts[idx] = cs.kappa_global[ik]
                zs_pts[idx] = cs.xz_coords[ixz][2]
            end

            spec_cells = MeshCell[]
            for ik in 1:nk-1
                for (i1, i2, i3, i4) in cs.xz_cells
                    b1=(ik-1)*nxz+i1; b2=(ik-1)*nxz+i2; b3=(ik-1)*nxz+i3; b4=(ik-1)*nxz+i4
                    t1= ik   *nxz+i1; t2= ik   *nxz+i2; t3= ik   *nxz+i3; t4= ik   *nxz+i4
                    push!(spec_cells, MeshCell(VTKCellTypes.VTK_HEXAHEDRON,
                                               [b1, b2, b3, b4, t1, t2, t3, t4]))
                end
            end

            fout_spec = joinpath(params.inputs[:output_dir], "les_xz_spectra_tavg")
            vtk_spec  = vtk_grid(fout_spec, xs_pts, yk_pts, zs_pts, spec_cells)
            for (ivar, vname) in enumerate(lesspectra_vars)
                data = [cs.sum_spectra[ixz, ik, ivar] / ns_spec
                        for ik in 1:nk for ixz in 1:nxz]
                vtk_spec[vname, VTKPointData()] = data
            end
            vtk_save(vtk_spec)
        end
    end

    # ---- time-averaged surface friction velocity: les_ustar_tavg.dat ----
    bc = params.les_bottom_cache
    if !isnothing(bc) && bc.n_samples > 0
        ns_bc = bc.n_samples
        outfile_ustar = joinpath(params.inputs[:output_dir], "les_ustar_tavg.dat")
        open(outfile_ustar, "w") do io
            println(io, "# time_end=", @sprintf("%.6e", t),
                    "  n_samples=", ns_bc, "  x  ustar")
            for ix in 1:length(bc.x_coords)
                @printf(io, "%.6e  %.6e\n", bc.x_coords[ix], bc.sum_ustar[ix] / ns_bc)
            end
        end
    end

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
    uaux       = params.uaux
    qe         = params.qp.qe
    sgs_stress = params.sgs_stress

    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)
    fill_sgs_cache!(params)

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
            user_les_profiles!(means_buf, prof_buf, @view(uaux[ip,:]), @view(qe[ip,:]), @view(sgs_stress[ip,:]), ET)
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
            user_les_profiles!(csprof_buf, csstress_buf, @view(uaux[ip,:]), @view(qe[ip,:]), @view(sgs_stress[ip,:]), ET)
            for k in 1:nprofiles_cs
                cslocal_online_mean[ixz, k] += csprof_buf[k]
            end
            for k in 1:nstress_cs
                cslocal_online_stress[ixz, k] += csstress_buf[k]
            end
        end
    end
    cs.n_online_samples += 1

    # Bottom-wall u* online accumulation (local sum; no Allreduce per step)
    bc = params.les_bottom_cache
    if !isnothing(bc)
        for ix in 1:length(bc.x_coords)
            for ip1 in bc.wall_groups[ix]
                z_in = get(bc.z_inside, ip1, bc.z0_m * 10.0)
                z_in <= bc.z0_m && continue
                if ET == PERT()
                    ρ = uaux[ip1,1] + qe[ip1,1]
                    u = (uaux[ip1,2] + qe[ip1,2]) / ρ
                    v = (uaux[ip1,3] + qe[ip1,3]) / ρ
                else
                    ρ = uaux[ip1,1]
                    u = uaux[ip1,2] / ρ
                    v = uaux[ip1,3] / ρ
                end
                bc.local_ustar_online[ix] += bc.κ_v * sqrt(u*u + v*v) / log(z_in / bc.z0_m)
            end
        end
        bc.n_online_samples += 1
    end
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

    comm = get_mpi_comm()
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

    # Bottom-wall u* online finalization
    bc = params.les_bottom_cache
    if !isnothing(bc) && bc.n_online_samples > 0
        nx_bc   = length(bc.x_coords)
        ns_bc   = bc.n_online_samples
        global_ustar = zeros(Float64, nx_bc)
        MPI.Allreduce!(bc.local_ustar_online, global_ustar, MPI.SUM, comm)
        rank != 0 && return
        outfile_ustar = joinpath(params.inputs[:output_dir], "les_ustar_online.dat")
        open(outfile_ustar, "w") do io
            println(io, "# online_tavg  time_end=", @sprintf("%.6e", t),
                    "  n_samples=", ns_bc, "  x  ustar")
            for ix in 1:nx_bc
                denom = Float64(bc.npts_per_x[ix]) * ns_bc
                @printf(io, "%.6e  %.6e\n", bc.x_coords[ix], global_ustar[ix] / denom)
            end
        end
    end
end

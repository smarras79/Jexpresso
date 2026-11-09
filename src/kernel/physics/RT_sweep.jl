Base.@kwdef struct RadiativeCase
    kappa_fn::Function                    # (x,y,z,ip) -> κ
    sigma_fn::Function                    # (x,y,z,ip) -> σ
    phase_fn::Function = (θ,θ1,ϕ,ϕ1,g) -> 1/(4π)
    g_const::Float64 = 0.0
    g_field::Union{Function,Nothing} = nothing   # (x,y,z,ip) -> g
    bc_fn::Function                       # (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_bc
    source_fn::Function = (x,y,z,ip,θ,ϕ) -> 0.0
end

Base.@kwdef struct SweepOptions
    async::Bool = true
    upwind_order::Bool = true
    incremental_pack::Bool = true
    face_plans::Bool = true
    local_plans::Bool = true
    flat_faces::Bool = true
    angular_plans::Bool = true
    gs_mode::Symbol = :cached_upwind
    gs_cyclic::Bool = true
    gs_cycle_cache::Bool = true
    gs_skip::Bool = true
    gs_inplace::Bool = true
    gs_gc_off::Bool = true
    gs_tol::Float64 = 1e-9
    lu_cache_mb::Float64 = -1.0        # -1 = from available memory (auto_cache_budget_mb)
    gc_mode::Symbol = :sync
    gc_every::Int = 0                  # 0 = auto: 16 for <=32 ranks, 8 above (measured optima)
    async_chunk::Int = 256
    async_interval_ms::Float64 = 0.0   # 0 = auto: 20 ms for <=32 ranks, 5 ms above
    round_cap::Int = -64               # used only by the synchronous sweep (async=false)
end

function apply_sweep_options!(o::SweepOptions, nranks::Int)
    USE_ASYNC_SWEEP[] = o.async;           USE_UPWIND_ORDER[] = o.upwind_order
    USE_INCREMENTAL_PACK[] = o.incremental_pack
    USE_FACE_PLANS[] = o.face_plans;       USE_LOCAL_NCF_PLANS[] = o.local_plans
    USE_FLAT_FACES[] = o.flat_faces;       USE_ANGULAR_PLANS[] = o.angular_plans
    GS_MODE[] = o.gs_mode;                 USE_GS_CYCLIC[] = o.gs_cyclic
    USE_GS_CYCLE_CACHE[] = o.gs_cycle_cache; USE_GS_SKIP[] = o.gs_skip
    USE_GS_INPLACE[] = o.gs_inplace;       USE_GS_GC_OFF[] = o.gs_gc_off
    ASYNC_CHUNK[] = o.async_chunk
    ASYNC_INTERVAL[] = (o.async_interval_ms > 0 ? o.async_interval_ms : (nranks <= 32 ? 20.0 : 5.0)) / 1000
    gc_every = o.gc_every > 0 ? o.gc_every : (nranks <= 32 ? 16 : 8)
    return (gc_mode = o.gc_mode, gc_every = gc_every, gs_tol = o.gs_tol,
            lu_cache_mb = o.lu_cache_mb, round_cap = o.round_cap)
end

function build_case_fields(mesh, case::RadiativeCase)
    npoin_cg = length(mesh.x)
    κ = zeros(npoin_cg); σ = zeros(npoin_cg)
    for ip = 1:npoin_cg
        κ[ip] = case.kappa_fn(mesh.x[ip], mesh.y[ip], mesh.z[ip], ip)
        σ[ip] = case.sigma_fn(mesh.x[ip], mesh.y[ip], mesh.z[ip], ip)
    end

    ngl = mesh.ngl
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            bc_fn_all[(iel,d)] = (x,y,z,ip,nx,ny,nz) -> case.bc_fn(x,y,z,ip,nx,ny,nz,θ,ϕ)
            #=q_vec = zeros(npoin_cg)
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                q_vec[gip] = case.source_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip, θ, ϕ)
            end=#
	    q_vec = zeros(ngl^3)
	    for k=1:ngl, j=1:ngl, i=1:ngl
    		gip = mesh.connijk[iel,i,j,k]
    		q_vec[lidx(i,j,k,ngl)] = case.source_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip, θ, ϕ)
	    end

            q_ext_all[(iel,d)] = q_vec
        end
    end
    return κ, σ, bc_fn_all, q_ext_all
end

function build_phase_coupling_matrices_elementwise(mesh, case::RadiativeCase)
    ngl = mesh.ngl
    C_iel_all = Dict{Int, Matrix{Float64}}()

    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        npd = em.extra_npoin
        w_ext = zeros(Float64, npd)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        thetas = em.extra_coords[1,:]; phis = em.extra_coords[2,:]

        g_iel = if case.g_field !== nothing
            gsum = 0.0
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                gsum += case.g_field(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip)
            end
            gsum / ngl^3
        else
            case.g_const
        end

        C = zeros(npd, npd)
        for dk = 1:npd, m = 1:npd
            C[dk,m] = case.phase_fn(thetas[dk], thetas[m], phis[dk], phis[m], g_iel) * w_ext[m]
        end
        C_iel_all[iel] = C
    end
    return C_iel_all
end

function build_longwave_case(atmos_data, κ_lw::Vector{Float64}, σ_lw::Vector{Float64}, bdy;
                              ε_surface=0.97, T_space=0.0)
    lw = LWParams(ε_surface, T_space)
    kappa_fn = (x,y,z,ip) -> κ_lw[ip] + σ_lw[ip]
    sigma_fn = (x,y,z,ip) -> σ_lw[ip]
    bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> user_rad_bc_longwave(x,y,z,θ,ϕ,nx,ny,nz,ip,atmos_data,bdy,lw)
    source_fn = (x,y,z,ip,θ,ϕ) -> user_rhs_longwave(x,y,z,θ,ϕ,ip,κ_lw,atmos_data)
    # LW scattering treated as isotropic per atmos_to_rad_longwave's own docs
    # (thermal-IR scattering is weak; ω₀_cloud=0.05 already folds the
    # MAGNITUDE into σ, only the ANGULAR redistribution is assumed isotropic).
    return RadiativeCase(kappa_fn=kappa_fn, sigma_fn=sigma_fn,
                          phase_fn=(θ,θ1,ϕ,ϕ1,g) -> 1/(4π), g_const=0.0,
                          bc_fn=bc_fn, source_fn=source_fn)
end

"""
Beer-Lambert direct beam, evaluated at every spatial node from the
precomputed vertical optical-depth profile. F_dir is a FLUX (W/m²);
G_dir = F_dir/μ₀ is the corresponding beam RADIANCE contribution to G,
matching radiative_heating.jl's convention exactly.
"""
function compute_direct_beam(mesh, κ_sw::Vector{Float64}, σ_sw::Vector{Float64}, sw::SWParams)
    npoin = length(mesh.x)
    κ_ext_sw = κ_sw .+ σ_sw
    F_dir = zeros(npoin)
    for ip = 1:npoin
        τ_ip = interp_optical_depth(mesh.z[ip], sw.z_prof, sw.τ_from_TOA)
        F_dir[ip] = sw.S₀_flux * sw.μ₀ * exp(-τ_ip / sw.μ₀)
    end
    G_dir = F_dir ./ sw.μ₀
    Q_dir = κ_sw .* G_dir
    return F_dir, G_dir, Q_dir
end

function build_shortwave_diffuse_case(κ_sw::Vector{Float64}, σ_sw::Vector{Float64}, g_eff::Vector{Float64},
                                       sw::SWParams, F_dir::Vector{Float64}, bdy; surface_albedo=0.15)
    κ_ext_sw = κ_sw .+ σ_sw
    ω0_lateral = sum(σ_sw) / sum(max.(κ_ext_sw, 1e-30))   # domain-mean SSA, per user_rad_bc_shortwave_diffuse's docstring

    kappa_fn = (x,y,z,ip) -> κ_ext_sw[ip]
    sigma_fn = (x,y,z,ip) -> σ_sw[ip]
    g_field  = (x,y,z,ip) -> g_eff[ip]
    bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> begin
        τ_z = interp_optical_depth(z, sw.z_prof, sw.τ_from_TOA)
        user_rad_bc_shortwave_diffuse(x,y,z,θ,ϕ,bdy,sw,F_dir[ip],τ_z,ω0_lateral,g_eff[ip]; surface_albedo=surface_albedo)
    end
    source_fn = (x,y,z,ip,θ,ϕ) -> user_rhs_shortwave_diffuse(x,y,z,θ,ϕ,ip,F_dir,σ_sw,sw,g_eff[ip])

    return RadiativeCase(kappa_fn=kappa_fn, sigma_fn=sigma_fn,
                          phase_fn=phase_function_hg, g_field=g_field,
                          bc_fn=bc_fn, source_fn=source_fn)
end

function compute_rt_radiative_heating_unified(mesh, I_cg_all::Dict{Tuple{Int,Int},Vector{Float64}},
                                               atmos_data, κ::Vector{Float64}, σ::Vector{Float64},
                                               swlw::Bool, comm;
                                               G_dir::Union{Vector{Float64},Nothing}=nothing,
                                               Q_dir::Union{Vector{Float64},Nothing}=nothing,
                                               sw_μ₀::Float64=0.0, sw_φ₀::Float64=0.0)
    PhysConst = PhysicalConst{Float64}()
    npoin = mesh.npoin
    σ_SB = 5.670374419e-8
    cₚ = PhysConst.cp
    ngl = mesh.ngl

    w_ext_all = build_angular_quadrature_weights(mesh)

    # ── GLOBAL node multiplicity, EXACTLY the CG solver's convention:
    # counted ONCE per (iel,i,j,k) triple (spatial only, no direction
    # index), reduced via ip2gip so inter-rank shares sum correctly. ────────
    node_div_local = zeros(Int, npoin)
    for iel = 1:mesh.nelem
        for k=1:ngl, j=1:ngl, i=1:ngl
            node_div_local[mesh.connijk[iel,i,j,k]] += 1
        end
    end
    g_count = zeros(Int, mesh.gnpoin)
    for ip = 1:npoin
        g_count[mesh.ip2gip[ip]] += node_div_local[ip]
    end
    MPI.Allreduce!(g_count, MPI.SUM, comm)
    node_div = zeros(Int, npoin)
    for ip = 1:npoin
        node_div[ip] = g_count[mesh.ip2gip[ip]]
    end

    # ── RAW (undivided) local sums — divide ONCE, after reduction. ──────────
    G_accum = zeros(npoin); Fx_accum = zeros(npoin); Fy_accum = zeros(npoin); Fz_accum = zeros(npoin)
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            haskey(I_cg_all, (iel,d)) || continue
            Ivec = I_cg_all[(iel,d)]
            dΩ = w_ext_all[iel][d]
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                Ival = Ivec[lidx(i,j,k,ngl)]
                G_accum[gip]  += Ival * dΩ
                Fx_accum[gip] += Ival * Ωx * dΩ
                Fy_accum[gip] += Ival * Ωy * dΩ
                Fz_accum[gip] += Ival * Ωz * dΩ
            end
        end
    end

    if MPI.Comm_size(comm) > 1
        gnpoin = mesh.gnpoin
        g_G = zeros(gnpoin); g_Fx = zeros(gnpoin); g_Fy = zeros(gnpoin); g_Fz = zeros(gnpoin)
        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            g_G[gip] += G_accum[ip]; g_Fx[gip] += Fx_accum[ip]
            g_Fy[gip] += Fy_accum[ip]; g_Fz[gip] += Fz_accum[ip]
        end
        MPI.Allreduce!(g_G, MPI.SUM, comm); MPI.Allreduce!(g_Fx, MPI.SUM, comm)
        MPI.Allreduce!(g_Fy, MPI.SUM, comm); MPI.Allreduce!(g_Fz, MPI.SUM, comm)
        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            G_accum[ip] = g_G[gip]; Fx_accum[ip] = g_Fx[gip]
            Fy_accum[ip] = g_Fy[gip]; Fz_accum[ip] = g_Fz[gip]
        end
    end

    for ip = 1:npoin
        G_accum[ip] /= node_div[ip]; Fx_accum[ip] /= node_div[ip]
        Fy_accum[ip] /= node_div[ip]; Fz_accum[ip] /= node_div[ip]
    end

    splitting = G_dir !== nothing && Q_dir !== nothing
    if splitting
        θ_sun = π - acos(clamp(sw_μ₀, 0.0, 1.0))
        Ωx_sun = sin(θ_sun)*cos(sw_φ₀); Ωy_sun = sin(θ_sun)*sin(sw_φ₀); Ωz_sun = cos(θ_sun)
        for ip = 1:npoin
            G_accum[ip] += G_dir[ip]
            F_dir_ip = G_dir[ip] * sw_μ₀
            Fx_accum[ip] += F_dir_ip*Ωx_sun; Fy_accum[ip] += F_dir_ip*Ωy_sun; Fz_accum[ip] += F_dir_ip*Ωz_sun
        end
    end

    Q_rad = zeros(npoin); dTdt_rad = zeros(npoin); F_net = zeros(npoin)
    for ip = 1:npoin
        T = atmos_data.t_current[ip]; ρ = atmos_data.rho[ip]
        κ_a = κ[ip]-σ[ip]; G = G_accum[ip]
        if swlw
            Q_rad[ip] = κ_a * (G - 4.0*σ_SB*T^4)
        else
            Q_rad[ip] = splitting ? Q_dir[ip] + κ_a*(G_accum[ip]-G_dir[ip]) : κ_a*G
        end
        F_net[ip] = sqrt(Fx_accum[ip]^2 + Fy_accum[ip]^2 + Fz_accum[ip]^2)
        dTdt_rad[ip] = Q_rad[ip] / (ρ*cₚ)
    end
    return Q_rad, dTdt_rad, F_net, G_accum
end

function build_angular_quadrature_weights(mesh)
    w_ext_all = Dict{Int, Vector{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        w_ext_all[iel] = w_ext
    end
    return w_ext_all
end

"""
    build_unified_ncf_infrastructure(mesh, sem, comm)

Builds every neighbor/correspondence/exchange table run_dg_sweep_multirank_unified!
needs, from the mesh AS GIVEN — no refinement is applied here. Call this
once per mesh topology (at initialization, and again whenever an external
AMR step changes the mesh), not once per RT solve.
"""
function build_unified_ncf_infrastructure(mesh, sem, comm)
    clear_face_plans!()
    rank = MPI.Comm_rank(comm)
    ngl = mesh.ngl

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_L1L2 = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
    end
    ncf_as_child_keys = Set(keys(ncf_as_child))
    ncf_as_parent_keys = Set(keys(ncf_as_parent))

    rigid_face_keys, periodic_face_keys = build_tagged_face_keys_periodic_aware(mesh)
    suspected_faces = Tuple{TInt,TInt}[]
    for iel = 1:mesh.nelem
        for fid = 1:6
            neighbor_elem[iel][fid] == zero(TInt) || continue
            (iel,fid) in ncf_as_child_keys && continue
            (iel,fid) in ncf_as_parent_keys && continue
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1 = mesh.conn_face_el[1, fe, iel]; n2 = mesh.conn_face_el[2, fe, iel]
            n3 = mesh.conn_face_el[3, fe, iel]; n4 = mesh.conn_face_el[4, fe, iel]
            key = Tuple(sort(TInt[n1, n2, n3, n4]))
            key in rigid_face_keys && continue
            push!(suspected_faces, (iel, fid))
        end
    end

    neighbor_rank, neighbor_kind, neighbor_gel_group, neighbor_is_periodic,
        same_rank_kind, same_rank_partner, same_rank_is_periodic, n_unresolved =
        resolve_periodic_or_rankcut_neighbors_ncf(mesh, suspected_faces, periodic_face_keys, comm)
    println("Rank $rank: unresolved=$n_unresolved  " * (n_unresolved==0 ? "✓" : "✗"))
    n_unresolved > 0 && error("Rank $rank: unresolved faces — aborting")

    neighbor_gel = Dict{Tuple{TInt,TInt}, TInt}()
    for ((iel,fid), kind) in neighbor_kind
        kind == :conforming || continue
        neighbor_gel[(TInt(iel),TInt(fid))] = neighbor_gel_group[(iel,fid)][1]
    end

    periodic_face_correspondence = build_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)

    remote_gid_coords = exchange_ncf_face_coordinates(mesh, neighbor_rank, comm)
    remote_angular_mesh = exchange_angular_mesh_across_rank_cut(mesh, neighbor_rank, comm)
    angular_deps_cr = build_angular_dependency_graph_crossrank(mesh, neighbor_rank, neighbor_gel, remote_angular_mesh)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    #gid_exch = build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm)
    gid_exch = augment_gid_exchange(build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm), ngl)
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    #periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)
    periodic_gid_exch = augment_gid_exchange(build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm), ngl)

    return (mesh=mesh, sem=sem,
            neighbor_elem=neighbor_elem, neighbor_fid=neighbor_fid, face_correspondence=face_correspondence,
            periodic_face_correspondence=periodic_face_correspondence,
            neighbor_rank=neighbor_rank, neighbor_kind=neighbor_kind, neighbor_gel=neighbor_gel,
            neighbor_gel_group=neighbor_gel_group, neighbor_is_periodic=neighbor_is_periodic,
            same_rank_kind=same_rank_kind, same_rank_partner=same_rank_partner, same_rank_is_periodic=same_rank_is_periodic,
            ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2, ncf_as_parent=ncf_as_parent,
            remote_gid_coords=remote_gid_coords,
            angular_deps=angular_deps, angular_deps_cr=angular_deps_cr, remote_angular_mesh=remote_angular_mesh,
            gid_exch=gid_exch, periodic_gid_exch=periodic_gid_exch, n_dirs=n_dirs)
end

function run_unified_case_solve(setup, metrics, comm, case::RadiativeCase;
                                 max_iter::Int=200, tol::Float64=1e-5,
                                 I_cg_all_init::Union{Dict{Tuple{Int,Int},Vector{Float64}},Nothing}=nothing,
                                 solver::Symbol=:jacobi,
                                 gmres_restart::Int=20, gmres_max_restarts::Int=10,
                                 opts::SweepOptions=SweepOptions(),
                                 report::Union{Nothing,Dict{Symbol,Any}}=nothing)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    mesh = setup.mesh
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs

    ro = apply_sweep_options!(opts, np)          # sets the switches; returns the keyword values

    κ, σ, bc_fn_all, q_ext_all = build_case_fields(mesh, case)
    C_iel_all = build_phase_coupling_matrices_elementwise(mesh, case)

    I_cg_init = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
        I_cg_init[(iel,dd)] = I_cg_all_init !== nothing && haskey(I_cg_all_init,(iel,dd)) ?
                               copy(I_cg_all_init[(iel,dd)]) : zeros(ngl^3)
    end
    order = sort(collect(keys(I_cg_init)))
    len_each = ngl^3

    # ---- LU cache budget: explicit MB, or from available memory when < 0 ----
        krylov = (solver == :gmres || solver == :gmres_precond) ? gmres_restart : 0
    cache_mb = ro.lu_cache_mb; budget_lines = String[]
    if cache_mb < 0
        ni_ = node_info(comm)
        AUTO_BUDGET_FN[] = () -> begin
            mb, lines = auto_cache_budget_mb(setup, ni_, comm; krylov_vectors=krylov)
            MPI.Comm_rank(comm) == 0 && (foreach(println, lines); flush(stdout))
            mb
        end
        cache_mb = 0.0                              # the sweep sets the real budget after discovery
    end
    cache = make_elem_cache(mesh, metrics, setup.n_dirs, ngl; budget_mb=cache_mb)
    memprobe("after make_elem_cache", comm)

    t_src = Float64[]; t_swp = Float64[]; t_tot = Float64[]
    gc_prev = ro.gc_mode == :sync ? GC.enable(false) : true      # sync mode: the sweep collects on a schedule

    function apply_G(flat::Vector{Float64})
        t_a = time_ns()
        I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        unflatten_I!(I_cg_all, flat, order, len_each)
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]; nd = em.extra_npoin; C_iel = C_iel_all[iel]
            gips  = [mesh.connijk[iel,i,j,k] for i=1:ngl, j=1:ngl, k=1:ngl][:]      # i fastest = lidx order
            σ_loc = σ[gips]
            Imat = Matrix{Float64}(undef, ngl^3, nd)
            for m = 1:nd; Imat[:,m] .= I_cg_all[(iel,m)]; end
            S = Imat * transpose(C_iel)
            for dk = 1:nd
                extra_source_all[(iel,dk)] = q_ext_all[(iel,dk)] .+ σ_loc .* view(S,:,dk)
            end
        end
        t_b = time_ns()
        x_dg, _ = run_dg_sweep_multirank_unified!(mesh, metrics, ψ, dψ, ω,
            setup.neighbor_elem, setup.neighbor_fid, setup.face_correspondence, setup.periodic_face_correspondence,
            setup.neighbor_rank, setup.neighbor_kind, setup.neighbor_gel, setup.neighbor_gel_group, setup.neighbor_is_periodic,
            setup.same_rank_kind, setup.same_rank_partner, setup.same_rank_is_periodic,
            setup.ncf_as_child, setup.ncf_L1L2, setup.ncf_as_parent,
            setup.remote_gid_coords,
            setup.angular_deps, setup.angular_deps_cr, setup.remote_angular_mesh,
            setup.gid_exch, setup.periodic_gid_exch,
            κ, bc_fn_all, comm;
            extra_source_all=extra_source_all, I_cg_all_prev=I_cg_all, elem_cache=cache,
            gc_mode=ro.gc_mode, gc_every=ro.gc_every, round_cap=ro.round_cap, gs_tol=ro.gs_tol)
        t_c = time_ns()
        out = flatten_I(x_dg, order)
        #if ro.gc_mode == :sync; GC.enable(true); GC.gc(false); GC.enable(false); end
        if ro.gc_mode == :sync; GC.enable(true); GC.gc(true); GC.enable(false); end
	t_d = time_ns()
        push!(t_src, (t_b - t_a)*1e-9); push!(t_swp, (t_c - t_b)*1e-9); push!(t_tot, (t_d - t_a)*1e-9)
        # memory per matvec (collective: every rank must reach these Allreduces)
        mv = Float64[Sys.maxrss()/1e9, Base.gc_live_bytes()/1e9]
        mmax = MPI.Allreduce(mv, MPI.MAX, comm)
        if rank == 0
            println("memory after matvec $(length(t_tot)): peak RSS max $(round(mmax[1]; digits=2)) GB, " *
                    "live heap max $(round(mmax[2]; digits=2)) GB")
        end
	if length(t_tot) == 1
            cm = MPI.Allreduce(Float64[cache.bytes/1e6, length(cache.factors)], MPI.MAX, comm)
            rank == 0 && (println("elem cache after matvec 1: max $(round(cm[1]; digits=0)) MB, $(Int(cm[2])) factors on the fullest rank"); flush(stdout))
        end
	return out
    end

    local result
    t_solve0 = time_ns()
    try
        if solver == :gmres
            x_final = gmres_fixed_point_solve(apply_G, flatten_I(I_cg_init, order), comm;
                                               restart=gmres_restart, max_restarts=gmres_max_restarts, tol=tol)
            I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(I_cg_all, x_final, order, len_each)
            result = (I_cg_all, κ, σ)
        elseif solver == :gmres_precond
            κ_mean = sum(κ)/max(length(κ),1)
            periodic_axes = detect_periodic_axes(mesh)
            isempty(periodic_axes) && error("gmres_precond requested but mesh has no periodic-tagged faces at all")
            scales_flat = build_scales_flat(compute_direction_preconditioner_scales(mesh, κ_mean, periodic_axes), order, len_each)
            x_final = gmres_fixed_point_solve_precond(apply_G, flatten_I(I_cg_init, order), scales_flat, comm;
                                                   restart=gmres_restart, max_restarts=gmres_max_restarts, tol=tol)
            I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(I_cg_all, x_final, order, len_each)
            result = (I_cg_all, κ, σ)
        elseif solver == :jacobi_extrap || solver == :jacobi_extrap_block
            mode = solver == :jacobi_extrap ? :global : :block
            x_final = jacobi_with_extrapolation(apply_G, flatten_I(I_cg_init, order), comm;
                                             len_each=len_each, max_iter=max_iter, tol=tol, mode=mode)
            I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(I_cg_all, x_final, order, len_each)
            result = (I_cg_all, κ, σ)
        else                                                                   # plain Jacobi
            I_cg_all = I_cg_init; result = (I_cg_all, κ, σ)
            for it = 1:max_iter
                flat_new = apply_G(flatten_I(I_cg_all, order))
                x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(x_dg, flat_new, order, len_each)
                max_delta_local = 0.0
                for iel = 1:mesh.nelem, dd = 1:n_dirs[iel]
                    max_delta_local = max(max_delta_local, maximum(abs.(x_dg[(iel,dd)] .- I_cg_all[(iel,dd)])))
                    I_cg_all[(iel,dd)] = copy(x_dg[(iel,dd)])
                end
                max_delta_global = MPI.Allreduce(max_delta_local, MPI.MAX, comm)
                rank == 0 && println("  RT iter $it: max_delta=$max_delta_global")
                if max_delta_global < tol
                    rank == 0 && println("RT converged after $it iterations")
                    break
                end
            end
        end
    finally
        GC.enable(gc_prev)
    end
    t_solve = (time_ns() - t_solve0) * 1e-9

    # ---- report (collective) ----
    nmv = length(t_tot)
    T = MPI.Allreduce(vcat(t_src, t_swp, t_tot), MPI.MAX, comm)              # max over ranks, per matvec
    rss = MPI.Allreduce(Sys.maxrss()/1e9, MPI.MAX, comm)
    if rank == 0
        src = T[1:nmv]; swp = T[nmv+1:2nmv]; tot = T[2nmv+1:3nmv]
        st = nmv >= 3 ? (3:nmv) : (1:nmv)
        foreach(println, budget_lines)
        println("coupled solve [$(solver)]: $nmv matvecs, solve time $(round(t_solve; digits=2)) s, peak RSS max $(round(rss; digits=2)) GB; " *
                "matvec 1 $(round(tot[1]; digits=2)) s (includes discovery, cache build); " *
                "steady (matvecs $(first(st))-$(nmv)) mean $(round(sum(tot[st])/length(st); digits=2)) s " *
                "(source $(round(sum(src[st])/length(st); digits=2)), sweep $(round(sum(swp[st])/length(st); digits=2))), " *
                "min $(round(minimum(tot[st]); digits=2)) max $(round(maximum(tot[st]); digits=2)); " *
        	"options: async=$(opts.async) cache_mb=$(round(cache.budget_bytes/1e6; digits=0)) cached_mb=$(round(cache.bytes/1e6; digits=0)) gc=$(ro.gc_mode)/$(ro.gc_every)")
	println("matvec times [s]: " * join([string(round(x; digits=2)) for x in tot], " "))
    end
    report === nothing || (report[:matvec_times] = tot_dummy(T, nmv); report[:solve_time] = t_solve)
    return result
end
tot_dummy(T, nmv) = T[2nmv+1:3nmv]


"""
Standard Type-I Anderson acceleration (Walker & Ni 2011; Fang & Saad 2009),
with explicit mixing parameter β and a conditioning safeguard.

x_{k+1} = (1-β)·(Σ α_i x_i) + β·(Σ α_i g_i)

where g_i = G(x_i) (the RAW fixed-point map output, i.e. x_new_raw at
step i), x_i are past ITERATES (not increments), and α minimizes
‖Σ α_i f_i‖ subject to Σα_i=1, with f_i = g_i - x_i.
"""
mutable struct AndersonAccelerator
    depth::Int
    β::Float64
    X::Vector{Vector{Float64}}   # past iterates x_i
    G::Vector{Vector{Float64}}   # past raw map outputs g_i = G(x_i)
end
AndersonAccelerator(depth::Int; β::Float64=0.5) = AndersonAccelerator(depth, β, Vector{Float64}[], Vector{Float64}[])

function anderson_step!(acc::AndersonAccelerator, x_prev::Vector{Float64}, g_new::Vector{Float64})
    push!(acc.X, copy(x_prev)); push!(acc.G, copy(g_new))
    if length(acc.X) > acc.depth
        popfirst!(acc.X); popfirst!(acc.G)
    end
    m = length(acc.X)
    if m == 1
        return (1-acc.β).*x_prev .+ acc.β.*g_new   # plain damped step, no history yet
    end

    F = [acc.G[i] .- acc.X[i] for i = 1:m]          # residuals f_i
    ΔF = hcat([F[i] .- F[m] for i = 1:m-1]...)       # (n, m-1)

    # ── Safeguard: if ΔF is (near-)rank-deficient, fall back to the plain
    # damped step rather than risk an ill-conditioned solve blowing up. ──────
    local γ
    try
        γ = ΔF \ (-F[m])
        any(!isfinite, γ) && error("non-finite γ")
    catch
        return (1-acc.β).*x_prev .+ acc.β.*g_new
    end

    α = zeros(m); α[1:m-1] .= γ; α[m] = 1 - sum(γ)

    x_mix = sum(α[i] .* acc.X[i] for i = 1:m)
    g_mix = sum(α[i] .* acc.G[i] for i = 1:m)
    return (1-acc.β).*x_mix .+ acc.β.*g_mix
end

flatten_I(I_cg_all, order) = vcat([I_cg_all[k] for k in order]...)
function unflatten_I!(I_cg_all, flat::Vector{Float64}, order, len_each::Int)
    for (idx,k) in enumerate(order)
        I_cg_all[k] = flat[(idx-1)*len_each+1 : idx*len_each]
    end
end

function memprobe(tag::String, comm)
    mv = Float64[Sys.maxrss()/1e9, Base.gc_live_bytes()/1e9]
    mx = MPI.Allreduce(mv, MPI.MAX, comm)
    if MPI.Comm_rank(comm) == 0
        println("memprobe [$tag]: peak RSS max $(round(mx[1]; digits=2)) GB, live heap max $(round(mx[2]; digits=2)) GB")
        flush(stdout)
    end
    return nothing
end

function invalidate_angular_caches!()
    clear_face_plans!()                 # plans, local/flat/angular plans, ASYNC_LAYOUTS, ASYNC_DEPS, GS_CYCLE_CACHE
    return nothing
end

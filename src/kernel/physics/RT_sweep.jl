Base.@kwdef struct RadiativeCase
    kappa_fn::Function                    # (x,y,z,ip) -> κ
    sigma_fn::Function                    # (x,y,z,ip) -> σ
    phase_fn::Function = (θ,θ1,ϕ,ϕ1,g) -> 1/(4π)
    g_const::Float64 = 0.0
    g_field::Union{Function,Nothing} = nothing   # (x,y,z,ip) -> g
    bc_fn::Function                       # (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_bc
    source_fn::Function = (x,y,z,ip,θ,ϕ) -> 0.0
    # Normalisation of the discrete scattering operator C[k,m] = Φ(Ω_k,Ω_m) w_m: :none, :row or :symmetric
    phase_norm::Symbol = :none
    # Direct-beam scattering source, added per node in build_case_fields:
    #   q_k(x) = coeff(x) Φ(Ω_k, Ω_beam; g(x)) / N(x),  N(x) = Σ_k w_k Φ(Ω_k, Ω_beam; g(x)) if normalize, else 1
    # fields: θ, ϕ (beam propagation direction), coeff::Vector (σ·I_dir per node), g::Vector, normalize::Bool
    beam_source::Union{Nothing,NamedTuple} = nothing
    # Fixed values that replace entries of the lagged iterate the sweep reads across PERIODIC faces:
    # (iel,d) => (lidx, values). Used to turn periodic faces into walls with a prescribed inflow
    # (lw_lateral_prev_override). nothing = ordinary periodic coupling.
    prev_override::Union{Nothing,Dict{Tuple{Int,Int},Tuple{Vector{Int},Vector{Float64}}}} = nothing
    # Lambertian reflection of the DIFFUSE downward flux at the bottom boundary: (albedo = A, is_bottom = z -> Bool).
    # Adds (A/π) F↓(x,y) to the inflow of every upward direction at bottom nodes, with F↓ = Σ_{Ωz<0} w |Ωz| I taken
    # from the iterate apply_G is called with (element-local), so G stays affine. Any fixed part of the reflected
    # radiance (e.g. A F_dir / π for the direct beam) stays in bc_fn. nothing = no diffuse reflection.
    surface_reflection::Union{Nothing,NamedTuple} = nothing
    # Asymmetry used for the scattering matrix: :element (one σ-weighted g per element) or :nodal (each node's own g;
    # nodes sharing a g share one matrix, so only cloud-edge nodes add matrices)
    phase_g_mode::Symbol = :element
    # source_fn independent of direction (LW emission): build_case_fields evaluates it once per element and every
    # direction shares that vector (only without a beam source, whose per-direction term is added in place)
    source_isotropic::Bool = false
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
    # gc_mode :sync only: also run a young-generation collection whenever this many GB have been allocated since the
    # last one (inside the sweep and after the source build). 0 = off (the schedule above only). The scheduled GC alone
    # let ~31 GB/rank of per-pair garbage build up within one sweep at 1.2 M pairs/rank (BOMEX 64x64x60 6x6O4, OOM).
    gc_alloc_gb::Float64 = 4.0
    # element solver on the workspace path: :dense (LU) or :fastdiag (sweep_fastdiag.jl; affine elements, falls back to
    # :dense per pair). ENV JEXPRESSO_RT_FASTDIAG=1/0 overrides.
    elem_solver::Symbol = :dense
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
    gb = something(tryparse(Float64, get(ENV, "JEXPRESSO_RT_GC_ALLOC_GB", "")), o.gc_alloc_gb)   # env override (tuning/tests)
    GC_ALLOC_LIMIT[] = o.gc_mode == :sync ? gb * 1e9 : 0.0
    GC_SYNC_MIN[] = something(tryparse(Float64, get(ENV, "JEXPRESSO_RT_GC_SYNC_MIN_GB", "")), 2.0) * 1e9
    GC_ALLOC_MARK[] = Base.gc_bytes()
    fdenv = get(ENV, "JEXPRESSO_RT_FASTDIAG", "")
    USE_FASTDIAG[] = fdenv == "" ? o.elem_solver == :fastdiag : fdenv == "1"
    o.elem_solver in (:dense, :fastdiag) || error("SweepOptions.elem_solver must be :dense or :fastdiag")
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
    share_q = case.source_isotropic && case.beam_source === nothing
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        q_shared = nothing
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            bc_fn_all[(iel,d)] = (x,y,z,ip,nx,ny,nz) -> case.bc_fn(x,y,z,ip,nx,ny,nz,θ,ϕ)
            if share_q && q_shared !== nothing          # same values for every direction: share one vector (read-only)
                q_ext_all[(iel,d)] = q_shared
                continue
            end
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
            share_q && (q_shared = q_vec)
        end
    end

    # Direct-beam scattering source (single beam direction): per node, Φ(Ω_k, Ω_beam; g) over the element's angular
    # nodes, optionally divided by its discrete integral Σ_k w_k Φ so the scattered beam energy is exactly σ·I_dir.
    bs = case.beam_source
    if bs !== nothing
        w_all = build_angular_quadrature_weights(mesh)
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]; nd = em.extra_npoin; w = w_all[iel]
            θd = em.extra_coords[1,1:nd]; ϕd = em.extra_coords[2,1:nd]
            Φd = zeros(nd)
            gtab = case.phase_norm == :symmetric_g ? phase_g_beam_table(θd, ϕd, w, Float64(bs.θ), Float64(bs.ϕ), phase_function_hg) : nothing
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                c = (haskey(bs, :coeff_elem) && bs.coeff_elem !== nothing) ? bs.coeff_elem[lidx(i,j,k,ngl), iel] : bs.coeff[gip]
                c == 0.0 && continue
                gb = gtab === nothing ? bs.g[gip] : phase_g_lookup(gtab, Float64(bs.g[gip]))
                for d = 1:nd; Φd[d] = phase_function_hg(θd[d], bs.θ, ϕd[d], bs.ϕ, gb); end
                N = bs.normalize ? sum(w[d]*Φd[d] for d = 1:nd) : 1.0
                l = lidx(i,j,k,ngl)
                for d = 1:nd; q_ext_all[(iel,d)][l] += c * Φd[d] / N; end
            end
        end
    end
    return κ, σ, bc_fn_all, q_ext_all
end

"""
    surface_reflection_setup!(bc_fn_all, mesh, case) -> state or nothing

For case.surface_reflection = (albedo A, is_bottom): finds, per element, the nodes on the bottom boundary, the downward
directions with weights w_d |cos θ_d| (Ωz = cos θ) and wraps bc_fn_all[(iel,d)] of every upward direction so its bottom
inflow gains (A/π) F↓[iel][node]. F↓ is filled by surface_reflection_update! from the iterate (element-local values,
so no communication). Upward-direction inflow at a bottom node is only read when that node is on an inflow face.
"""
function surface_reflection_setup!(bc_fn_all, mesh, case::RadiativeCase)
    sr = case.surface_reflection
    sr === nothing && return nothing
    A = sr.albedo; isb = sr.is_bottom; ngl = mesh.ngl
    w_all = build_angular_quadrature_weights(mesh)
    elems = Dict{Int, NamedTuple}()
    for iel = 1:mesh.nelem
        lids = Int[]; gips = Int[]
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            gip = mesh.connijk[iel,i,j,k]
            isb(mesh.z[gip]) && (push!(lids, lidx(i,j,k,ngl)); push!(gips, gip))
        end
        isempty(lids) && continue
        em = mesh.extra_mesh[iel]; nd = em.extra_npoin; w = w_all[iel]
        dn = Int[]; wdn = Float64[]
        for d = 1:nd
            μ = cos(em.extra_coords[1,d])
            μ < 0.0 && (push!(dn, d); push!(wdn, w[d] * (-μ)))
        end
        Fdn = Dict{Int,Float64}(g => 0.0 for g in gips)
        elems[iel] = (lids = lids, gips = gips, dn = dn, wdn = wdn, Fdn = Fdn)
        for d = 1:nd
            cos(em.extra_coords[1,d]) > 0.0 || continue
            base = bc_fn_all[(iel,d)]
            bc_fn_all[(iel,d)] = (x,y,z,ip,nx,ny,nz) -> begin
                v = base(x,y,z,ip,nx,ny,nz)
                isb(z) ? v + A / π * get(Fdn, ip, 0.0) : v
            end
        end
    end
    return (albedo = A, elems = elems)
end

function surface_reflection_update!(refl, I_cg_all)
    for (iel, e) in refl.elems
        for q in eachindex(e.lids)
            l = e.lids[q]; F = 0.0
            for m in eachindex(e.dn)
                F += e.wdn[m] * I_cg_all[(iel, e.dn[m])][l]
            end
            e.Fdn[e.gips[q]] = F
        end
    end
    return nothing
end

"""
    normalize_phase_matrix(C, w, mode) -> C̃

C[k,m] = Φ(Ω_k,Ω_m) w_m is the discrete scattering operator; the exact one conserves scattered energy,
Σ_m C[k,m] = 1. mode = :row divides each row by its sum. mode = :symmetric writes C = K diag(w) with K symmetric and
finds d > 0 such that C̃ = diag(d) K diag(d) diag(w) has unit row sums (Sinkhorn: d ← sqrt(d ./ (K (d .* w)))), so the
kernel stays symmetric (reciprocity Φ(Ω,Ω′) = Φ(Ω′,Ω)).
"""
function normalize_phase_matrix(C::Matrix{Float64}, w::Vector{Float64}, mode::Symbol)
    mode == :none && return C
    if mode == :row
        return C ./ sum(C; dims=2)
    elseif mode == :symmetric || mode == :symmetric_g
        K = C ./ transpose(w)
        d = ones(length(w))
        for it = 1:500
            d = sqrt.(d ./ (K * (d .* w)))
            maximum(abs.(d .* (K * (d .* w)) .- 1.0)) < 1e-13 && break
        end
        return (d .* K .* transpose(d)) .* transpose(w)
    else
        error("normalize_phase_matrix: unknown mode $(mode)")
    end
end

# :symmetric_g (asymmetry-preserving normalisation). A pointwise-sampled, Sinkhorn-normalised HG kernel on the angular
# grid conserves energy but scatters too far forward when the peak is narrow compared with the node spacing: its discrete
# asymmetry Σ_m C_km Ω_k·Ω_m exceeds g (HG g' = 0.840 -> 0.864 on 4x4O4, +15 % of 1-g'; 6x6O4 +5 %; 6x6O6 +1.2 %), and
# σ(1-g) controls diffusion through the cloud. For each distinct angular grid this tabulates the discrete asymmetry of
# the normalised kernel as a function of the HG parameter, then builds each kernel with the parameter whose discrete
# asymmetry equals the target g, so the 0th and 1st moments are both exact (Hunter & Guo-type normalisation).
const PHASE_GMAP_CACHE = Dict{UInt, Tuple{Vector{Float64}, Vector{Float64}}}()
function phase_g_table(thetas, phis, w, phase_fn)
    key = hash((thetas, phis, w))
    return get!(PHASE_GMAP_CACHE, key) do
        n = length(w)
        X = [sin(thetas[k])*cos(phis[k]) for k = 1:n]; Y = [sin(thetas[k])*sin(phis[k]) for k = 1:n]; Z = cos.(thetas)
        cosm = [clamp(X[k]*X[m] + Y[k]*Y[m] + Z[k]*Z[m], -1.0, 1.0) for k = 1:n, m = 1:n]
        gs = collect(0.0:0.025:0.975); χ = similar(gs); W = sum(w)
        for (i, gin) in enumerate(gs)
            C = [phase_fn(thetas[k], thetas[m], phis[k], phis[m], gin) * w[m] for k = 1:n, m = 1:n]
            Cn = normalize_phase_matrix(C, w, :symmetric)
            χ[i] = sum(w[k] * Cn[k,m] * cosm[k,m] for k = 1:n, m = 1:n) / W
        end
        all(diff(χ) .> 0) || @warn "phase_g_table: discrete asymmetry not monotone in g on this grid" χ
        (χ, gs)
    end
end
# Beam-source analogue of phase_g_table: the direct-beam scattering source Φ(Ω_k, Ω_beam; g)/N over the grid nodes has
# its own discrete asymmetry Σ_k w_k Φ Ω_k·Ω_beam / N, which the matrix table does not fix: the beam direction need not be
# a grid node, and when it falls on one (μ0 = 0.5 is a θ element boundary of the 6x6 grids) the sampled peak is biased
# (6x6O4 with the matrix map: <cos> 0.836 for g' 0.847, -7.5 % of 1-g'; median-column in-cloud absorption -0.8 %).
# Tabulated per (grid, beam direction) and inverted with phase_g_lookup.
const PHASE_BEAM_GMAP_CACHE = Dict{UInt, Tuple{Vector{Float64}, Vector{Float64}}}()
function phase_g_beam_table(thetas, phis, w, θb::Float64, ϕb::Float64, phase_fn)
    key = hash((thetas, phis, w, θb, ϕb))
    return get!(PHASE_BEAM_GMAP_CACHE, key) do
        n = length(w)
        c = [clamp(sin(thetas[k])*sin(θb)*cos(phis[k] - ϕb) + cos(thetas[k])*cos(θb), -1.0, 1.0) for k = 1:n]
        gs = collect(0.0:0.025:0.975); χ = similar(gs); Φ = zeros(n)
        for (i, gin) in enumerate(gs)
            for k = 1:n; Φ[k] = phase_fn(thetas[k], θb, phis[k], ϕb, gin); end
            χ[i] = sum(w .* Φ .* c) / sum(w .* Φ)
        end
        all(diff(χ) .> 0) || @warn "phase_g_beam_table: discrete beam asymmetry not monotone in g" χ
        (χ, gs)
    end
end

# HG parameter whose normalised discrete kernel has asymmetry gt (linear interpolation; clamped at the table ends).
function phase_g_lookup(tab, gt::Float64)
    χ, gs = tab
    gt <= χ[1] && return max(gt - χ[1] + gs[1], 0.0)
    gt >= χ[end] && return gs[end]
    i = searchsortedlast(χ, gt)
    return gs[i] + (gs[i+1] - gs[i]) * (gt - χ[i]) / (χ[i+1] - χ[i])
end

# Phase matrices shared between elements: C depends only on the element's angular grid, its g, the phase function and
# the normalisation, so elements with the same key share ONE matrix (read-only everywhere). LW (isotropic): one matrix for
# all elements instead of one per element (6x6O6: 1332^2 x 8 B = 14 MB each, 13.6 GB/rank at 960 elements); SW: all
# clear-air elements (Rayleigh, g = 0) share one, cloud-interior elements another. The previous call's matrices are kept
# and reused when the key repeats; only the current call's set is retained, so memory stays bounded.
const PHASE_MATRIX_CACHE = Ref(Dict{Tuple{UInt,Float64,Symbol,UInt},Matrix{Float64}}())

# One cached matrix per (grid, g, normalisation, phase function); see PHASE_MATRIX_CACHE.
# g is rounded to PHASE_G_QUANTUM before keying and building: nodal and σ-weighted g vary continuously (cloud/Rayleigh mix:
# 0.8449-0.8471 in the BOMEX median cloud), which would otherwise give one matrix per cloud node or element. 1e-3 changes
# 1-g by <= 0.3 % in cloud (1-g ~ 0.15).
const PHASE_G_QUANTUM = 1e-3
function _phase_matrix!(cur, prev, case::RadiativeCase, thetas, phis, w_ext, gridkey::UInt, g::Float64, cnt)
    g = round(g / PHASE_G_QUANTUM) * PHASE_G_QUANTUM
    key = (gridkey, g, case.phase_norm, objectid(case.phase_fn))
    C = get(cur, key, nothing)
    C === nothing || return C
    C = get(prev, key, nothing)
    if C === nothing
        npd = length(w_ext)
        g_use = case.phase_norm == :symmetric_g ?
                phase_g_lookup(phase_g_table(thetas, phis, w_ext, case.phase_fn), g) : g
        Cr = zeros(npd, npd)
        for dk = 1:npd, m = 1:npd
            Cr[dk,m] = case.phase_fn(thetas[dk], thetas[m], phis[dk], phis[m], g_use) * w_ext[m]
        end
        C = normalize_phase_matrix(Cr, w_ext, case.phase_norm); cnt[1] += 1
    else
        cnt[2] += 1
    end
    cur[key] = C
    return C
end

function _angular_grid(em)
    npd = em.extra_npoin
    w_ext = zeros(Float64, npd)
    for e_ext = 1:em.extra_nelem
        nop = em.extra_nop[e_ext]
        for iθ=1:nop+1, iϕ=1:nop+1
            ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
            w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
        end
    end
    thetas = em.extra_coords[1,1:npd]; phis = em.extra_coords[2,1:npd]
    return thetas, phis, w_ext, hash((thetas, phis, w_ext))
end

"""
    build_phase_coupling_matrices_elementwise(mesh, case) -> (C_iel_all, node_groups)

C_iel_all[iel]: the element's scattering matrix C[k,m] = Φ(Ω_k·Ω_m; g) w_m (normalised per case.phase_norm), with g the
σ-weighted mean of the nodal g. node_groups (case.phase_g_mode == :nodal, else nothing): node_groups[iel] = list of
(local node indices, matrix) with each node's own g; nodes with σ = 0 join the first group (their source is zero).
"""
function build_phase_coupling_matrices_elementwise(mesh, case::RadiativeCase)
    ngl = mesh.ngl
    C_iel_all = Dict{Int, Matrix{Float64}}()
    nodal = case.phase_g_mode == :nodal && case.g_field !== nothing
    case.phase_g_mode in (:element, :nodal) || error("phase_g_mode must be :element or :nodal")
    groups = nodal ? Dict{Int, Vector{Tuple{Vector{Int},Matrix{Float64}}}}() : nothing
    prev = PHASE_MATRIX_CACHE[]; cur = Dict{Tuple{UInt,Float64,Symbol,UInt},Matrix{Float64}}()
    cnt = [0, 0]; ngroups_max = 0
    gn = zeros(ngl^3); sn = zeros(ngl^3)

    for iel = 1:mesh.nelem
        thetas, phis, w_ext, gridkey = _angular_grid(mesh.extra_mesh[iel])

        # Element asymmetry: σ-weighted mean of the nodal g (the scattering-weighted value the element's scattering
        # actually has). An unweighted mean in an element straddling a cloud edge mixed cloud g = 0.85 with Rayleigh
        # g = 0 at nodes that scatter 10^3-10^4 times less, giving g ~ 0.4 where nearly all scattering is cloud.
        g_iel = if case.g_field !== nothing
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]; l = lidx(i,j,k,ngl)
                gn[l] = case.g_field(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip)
                sn[l] = case.sigma_fn(mesh.x[gip], mesh.y[gip], mesh.z[gip], gip)
            end
            ws = sum(sn); ws > 0 ? sum(sn .* gn) / ws : sum(gn) / ngl^3
        else
            case.g_const
        end
        C_iel_all[iel] = _phase_matrix!(cur, prev, case, thetas, phis, w_ext, gridkey, Float64(g_iel), cnt)

        if nodal
            byg = Dict{Float64, Vector{Int}}(); zero_nodes = Int[]
            for l = 1:ngl^3
                sn[l] > 0 ? push!(get!(byg, round(gn[l] / PHASE_G_QUANTUM) * PHASE_G_QUANTUM, Int[]), l) : push!(zero_nodes, l)
            end
            gl = Tuple{Vector{Int},Matrix{Float64}}[]
            for (g, nodes) in byg
                push!(gl, (nodes, _phase_matrix!(cur, prev, case, thetas, phis, w_ext, gridkey, g, cnt)))
            end
            if isempty(gl)
                push!(gl, (zero_nodes, C_iel_all[iel]))
            else
                append!(gl[1][1], zero_nodes)
            end
            groups[iel] = gl; ngroups_max = max(ngroups_max, length(gl))
        end
    end
    PHASE_MATRIX_CACHE[] = cur
    PHASE_MATRIX_STATS[] = (elements = mesh.nelem, distinct = length(cur), built = cnt[1], reused_from_previous = cnt[2],
                            max_groups = ngroups_max)
    return C_iel_all, groups
end
const PHASE_MATRIX_STATS = Ref{Any}(nothing)

# lateral_cols (LWLateralColumns, lw_lateral.jl): if given, real lateral walls get the plane-parallel radiance of the
#   local column instead of the blackbody σT_back⁴/π of user_rad_bc_longwave.
# prev_override: see RadiativeCase (periodic faces treated as walls with that inflow).
function build_longwave_case(atmos_data, κ_lw::Vector{Float64}, σ_lw::Vector{Float64}, bdy;
                              ε_surface=0.97, T_space=0.0, lateral_cols=nothing, prev_override=nothing,
                              surface_reflection::Bool=false)
    lw = LWParams(ε_surface, T_space)
    kappa_fn = (x,y,z,ip) -> κ_lw[ip] + σ_lw[ip]
    sigma_fn = (x,y,z,ip) -> σ_lw[ip]
    bc_fn = lateral_cols === nothing ?
        ((x,y,z,ip,nx,ny,nz,θ,ϕ) -> user_rad_bc_longwave(x,y,z,θ,ϕ,nx,ny,nz,ip,atmos_data,bdy,lw)) :
        ((x,y,z,ip,nx,ny,nz,θ,ϕ) -> (bdy.is_top(z) || bdy.is_bottom(z)) ?
              user_rad_bc_longwave(x,y,z,θ,ϕ,nx,ny,nz,ip,atmos_data,bdy,lw) :
              lw_lateral_radiance(lateral_cols, x, y, z, cos(θ)))
    source_fn = (x,y,z,ip,θ,ϕ) -> user_rhs_longwave(x,y,z,θ,ϕ,ip,κ_lw,atmos_data)
    # LW scattering treated as isotropic per atmos_to_rad_longwave's own docs
    # (thermal-IR scattering is weak; ω₀_cloud=0.05 already folds the
    # MAGNITUDE into σ, only the ANGULAR redistribution is assumed isotropic).
    return RadiativeCase(kappa_fn=kappa_fn, sigma_fn=sigma_fn,
                          phase_fn=(θ,θ1,ϕ,ϕ1,g) -> 1/(4π), g_const=0.0, source_isotropic=true,
                          bc_fn=bc_fn, source_fn=source_fn, prev_override=prev_override,
                          # surface_reflection: Lambertian reflection of the downward LW flux with albedo 1 - ε at the
                          # bottom (bc_fn keeps the emission ε B); off by default, as the offline ICA reference in
                          # lw_lean.py (reflect=False) and the 1D lateral/initial columns (lw_lateral.jl) assume.
                          surface_reflection = (surface_reflection && ε_surface < 1.0) ?
                              (albedo = 1.0 - Float64(ε_surface), is_bottom = bdy.is_bottom) : nothing)
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
                                       sw::SWParams, F_dir::Vector{Float64}, bdy; surface_albedo=0.15,
                                       phase_norm::Symbol=:none, lateral_bc::Symbol=:plane_parallel,
                                       mesh=nothing, periodic_walls::Bool=false, comm=nothing, beam_coeff=nothing,
                                       beam_coeff_elem=nothing, diffuse_reflection::Bool=true,
                                       phase_g_mode::Symbol=:element)
    κ_ext_sw = κ_sw .+ σ_sw
    ω0_lateral = sum(σ_sw) / sum(max.(κ_ext_sw, 1e-30))   # domain-mean SSA, per user_rad_bc_shortwave_diffuse's docstring

    kappa_fn = (x,y,z,ip) -> κ_ext_sw[ip]
    sigma_fn = (x,y,z,ip) -> σ_sw[ip]
    g_field  = (x,y,z,ip) -> g_eff[ip]
    # Boundary inflow of DIFFUSE radiance. Top: none. Bottom: Lambertian reflection of the direct beam here, plus
    # (diffuse_reflection, default) the reflected diffuse downward flux via RadiativeCase.surface_reflection.
    # Lateral walls (non-periodic meshes only):
    #   :single_scatter  plane-parallel single-scattering radiance of the reference column (SWLateralSingleScatter)
    #   :none            no inflow
    #   :legacy          the old user_rad_bc_shortwave_diffuse estimate ω̄ Φ S₀ (1 - e^{-τ/μ₀})
    #   :plane_parallel  single scattering + multiply scattered part from the 1D reference solver (SWLateralPlaneParallel)
    # Built lazily on the first lateral call, so periodic meshes (no lateral faces) never pay for it.
    # The column profile is gathered over MPI (collective), so it is built here on every rank; only the per-rank
    # radiance object (which runs the 1D solver) is built lazily.
    lat_ref = Ref{Any}(nothing)
    prof = nothing
    if lateral_bc in (:single_scatter, :plane_parallel) || periodic_walls
        mesh === nothing && error("build_shortwave_diffuse_case: lateral_bc = $(lateral_bc) needs mesh=")
        prof = build_sw_column_optics(mesh, κ_ext_sw, σ_sw, g_eff, mesh.ngl)
    end
    make_lat() = begin
        SW_TOP[] !== nothing && return SW_TOP[].lat      # extended column (air above the lid + domain column), true S₀
        lateral_bc == :single_scatter ?
            SWLateralSingleScatter(prof, sw.S₀_flux, sw.μ₀, sw.φ₀, surface_albedo) :
            build_sw_lateral_plane_parallel(prof, sw.S₀_flux, sw.μ₀, sw.φ₀, surface_albedo)
    end
    if !(lateral_bc in (:single_scatter, :plane_parallel, :none, :legacy))
        error("build_shortwave_diffuse_case: unknown lateral_bc $(lateral_bc)")
    end
    bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> begin
        if bdy.is_bottom(z)
            surface_albedo * F_dir[ip] / π
        elseif bdy.is_top(z)
            st = SW_TOP[]                        # :RT_sw_top_bc => :column: diffuse light from the air above the lid
            st === nothing ? 0.0 : st.lat(st.z_lid, θ, ϕ)
        elseif lateral_bc in (:single_scatter, :plane_parallel)
            lat_ref[] === nothing && (lat_ref[] = make_lat())
            lat_ref[](z, θ, ϕ)
        elseif lateral_bc == :none
            0.0
        else
            τ_z = interp_optical_depth(z, sw.z_prof, sw.τ_from_TOA)
            user_rad_bc_shortwave_diffuse(x,y,z,θ,ϕ,bdy,sw,F_dir[ip],τ_z,ω0_lateral,g_eff[ip]; surface_albedo=surface_albedo)
        end
    end
    # Direct beam scattered into the diffuse field: σ(x) I_dir(x) Φ(Ω, Ω_sun; g(x)), I_dir = F_dir/μ₀, with the same
    # phase function as the scattering operator (phase_function_hg) and, unless phase_norm = :none, normalised per node.
    θ_sun, ϕ_sun = rt_sun_angles(sw)
    beam = sw.μ₀ > 0.0 ?
        (θ = θ_sun, ϕ = ϕ_sun, coeff = beam_coeff === nothing ? σ_sw .* F_dir ./ sw.μ₀ : beam_coeff, g = g_eff, normalize = phase_norm != :none,
         coeff_elem = beam_coeff_elem) : nothing

    # periodic_walls: periodic faces become walls whose diffuse inflow is the lateral_bc radiance (prev_override)
    ovr = nothing
    if periodic_walls
        lateral_bc in (:single_scatter, :plane_parallel, :none) ||
            error("build_shortwave_diffuse_case: periodic_walls needs lateral_bc :single_scatter, :plane_parallel or :none")
        lat = lateral_bc == :none ? nothing : make_lat()
        ovr = lateral_prev_override(mesh, (x,y,z,θ,ϕ) -> lat === nothing ? 0.0 : lat(z, θ, ϕ), comm)
    end
    return RadiativeCase(kappa_fn=kappa_fn, sigma_fn=sigma_fn,
                          phase_fn=phase_function_hg, g_field=g_field,
                          bc_fn=bc_fn, phase_norm=phase_norm, beam_source=beam, prev_override=ovr,
                          phase_g_mode=phase_g_mode,
                          surface_reflection = (diffuse_reflection && surface_albedo > 0.0) ?
                              (albedo = Float64(surface_albedo), is_bottom = bdy.is_bottom) : nothing)
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

    # Always reduce over global node ids, also on ONE rank: node_div counts every element sharing the global node,
    # and on a periodic mesh the two faces' nodes share one global id even within a rank. Skipping this on one rank
    # halved G on periodic-face nodes and quartered it on corner nodes (local uniform-cloud test: -3647 / -5473 K/day
    # there vs -1 K/day inside).
    begin
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
    if get(ENV, "JEXPRESSO_RT_DEBUG_TOPO", "") == "1"
        # Elements that are their own face neighbour (e.g. a one-element-thick periodic slab)
        n_self_adj = count(neighbor_elem[iel][fid] == iel for iel = 1:mesh.nelem, fid = 1:6)
        n_self_partner = count(iel in p for ((iel, fid), p) in same_rank_partner)
        n_bdy = count(neighbor_elem[iel][fid] == 0 for iel = 1:mesh.nelem, fid = 1:6)
        kinds = Dict{Symbol,Int}(); for (_, k) in same_rank_kind; kinds[k] = get(kinds, k, 0) + 1; end
        println("Rank $rank [topo]: nelem=$(mesh.nelem) self-adjacent faces=$n_self_adj  boundary/periodic faces=$n_bdy  " *
                "suspected=$(length(suspected_faces)) same-rank periodic partners=$(length(same_rank_partner)) " *
                "(self=$n_self_partner, kinds=$kinds) cross-rank=$(length(neighbor_rank))")
        test_periodic_face_correspondence_same_rank(mesh, same_rank_kind, same_rank_partner, ngl)
        for iel = 1:min(mesh.nelem, 2), fid = 1:6
            println("Rank $rank [topo]: iel=$iel fid=$fid neighbor=$(neighbor_elem[iel][fid]) nfid=$(neighbor_fid[iel][fid]) " *
                    "partner=$(get(same_rank_partner, (iel,fid), nothing))")
        end
    end
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
    gid_exch = compact_gid_exchange(augment_gid_exchange(build_ncf_gid_value_exchange(mesh, neighbor_rank, n_dirs, comm), ngl))
    periodic_only_faces = Dict{Tuple{TInt,TInt},TInt}(k => v for (k,v) in neighbor_rank if get(neighbor_is_periodic, k, false))
    #periodic_gid_exch = build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm)
    periodic_gid_exch = compact_gid_exchange(augment_gid_exchange(build_ncf_gid_value_exchange(mesh, periodic_only_faces, n_dirs, comm), ngl))

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

const LIMITED_ELEMS = Ref{Vector{Int32}}(Int32[])   # per local element: directions limited in the last final sweep
# The GMRES solution BEFORE the final limited sweep (nothing if no limited sweep ran). It solves the linear system, so
# it is the right warm start for the next solve; the limited field is not (warm-starting from it cost 102 matvecs
# instead of 21 on a same-state re-solve, local 2x2O3).
const UNLIMITED_SOLUTION = Ref{Any}(nothing)

# Jacobi step statistics in one pass and without allocation: max |xn - x|, ‖xn - x‖², ‖xn‖². A separate function so the
# loop is type-stable (inline in run_unified_case_solve, whose locals are boxed by its closures, it allocated per element
# with GC disabled: +34 GB live heap per rank on bubble 40x40x60 6x6O6, job 1378743 OOM).
function _jacobi_step_stats(xn::Vector{Float64}, x::Vector{Float64})
    dm = 0.0; ss = 0.0; sx = 0.0
    @inbounds @simd for i in eachindex(xn, x)
        d = xn[i] - x[i]; dm = max(dm, abs(d)); ss += d * d; sx += xn[i] * xn[i]
    end
    return dm, ss, sx
end

function run_unified_case_solve(setup, metrics, comm, case::RadiativeCase;
                                 max_iter::Int=200, tol::Float64=1e-5,
                                 I_cg_all_init::Union{Dict{Tuple{Int,Int},Vector{Float64}},Nothing}=nothing,
                                 solver::Symbol=:jacobi,
                                 gmres_restart::Int=20, gmres_max_restarts::Int=10,
                                 gmres_basis::Type{<:AbstractFloat}=Float64,
                                 final_limited_sweep::Bool=false, limited_sweeps::Int=1, limiter_mode::Symbol=:posthoc,
                                 rtol::Float64=0.0,
                                 monitor::Union{Nothing,Function}=nothing,
                                 opts::SweepOptions=SweepOptions(),
                                 report::Union{Nothing,Dict{Symbol,Any}}=nothing)
    rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    mesh = setup.mesh
    ψ, dψ, ω = setup.sem.basis.ψ, setup.sem.basis.dψ, setup.sem.ω
    ngl = mesh.ngl; n_dirs = setup.n_dirs

    ro = apply_sweep_options!(opts, np)          # sets the switches; returns the keyword values

    t_setup0 = time_ns()
    κ, σ, bc_fn_all, q_ext_all = build_case_fields(mesh, case)
    t_s1 = time_ns()
    C_iel_all, C_groups = build_phase_coupling_matrices_elementwise(mesh, case)
    t_s2 = time_ns()
    refl = surface_reflection_setup!(bc_fn_all, mesh, case)        # nothing unless case.surface_reflection

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
            mb, lines = auto_cache_budget_mb(setup, ni_, comm; krylov_vectors=krylov, krylov_bytes=sizeof(gmres_basis))
            MPI.Comm_rank(comm) == 0 && (foreach(println, lines); flush(stdout))
            mb
        end
        cache_mb = 0.0                              # the sweep sets the real budget after discovery
    end
    t_s3 = time_ns()
    cache = make_elem_cache(mesh, metrics, setup.n_dirs, ngl; budget_mb=cache_mb)
    memprobe("after make_elem_cache", comm)
    t_setup = (time_ns() - t_setup0) * 1e-9          # case fields, phase matrices, element cache (before the solve)
    # setup breakdown (max over ranks) and phase-matrix sharing
    tb = MPI.Allreduce(Float64[t_s1 - t_setup0, t_s2 - t_s1, t_s3 - t_s2, time_ns() - t_s3] .* 1e-9, MPI.MAX, comm)
    ps = PHASE_MATRIX_STATS[]
    pd = ps === nothing ? zeros(4) : MPI.Allreduce(Float64[ps.distinct, ps.built, ps.reused_from_previous, ps.max_groups], MPI.MAX, comm)
    rank == 0 && println(@sprintf("  setup: case fields %.1f s, phase matrices %.1f s (g per %s; max per rank: %d distinct, %d built, %d reused from the previous call, %d node groups per element), reflection + initial guess %.1f s, element cache %.1f s",
                                  tb[1], tb[2], case.phase_g_mode, Int(pd[1]), Int(pd[2]), Int(pd[3]), Int(pd[4]), tb[3], tb[4]))

    t_src = Float64[]; t_swp = Float64[]; t_tot = Float64[]
    # per-element RT cost prediction: the repartitioning weights (:RT_rebalance)
    cost_pred = RT_W_ON[] ? rt_predict_costs(mesh, metrics, dψ, ω, κ, ngl, collect(Int, n_dirs), cache.faceinfo) : nothing
    gc_prev = ro.gc_mode == :sync ? GC.enable(false) : true      # sync mode: the sweep collects on a schedule

    function apply_G(flat::Vector{Float64})
        t_a = time_ns()
        I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        unflatten_I!(I_cg_all, flat, order, len_each)
        refl === nothing || surface_reflection_update!(refl, I_cg_all)   # F↓ of this iterate -> bottom inflow
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]; nd = em.extra_npoin; C_iel = C_iel_all[iel]
            gips  = [mesh.connijk[iel,i,j,k] for i=1:ngl, j=1:ngl, k=1:ngl][:]      # i fastest = lidx order
            σ_loc = σ[gips]
            if all(iszero, σ_loc)                     # no scattering here: the source is the fixed one (read-only)
                for dk = 1:nd; extra_source_all[(iel,dk)] = q_ext_all[(iel,dk)]; end
                continue
            end
            Imat = Matrix{Float64}(undef, ngl^3, nd)
            for m = 1:nd; Imat[:,m] .= I_cg_all[(iel,m)]; end
            S = if C_groups === nothing || length(C_groups[iel]) == 1
                Imat * transpose(C_groups === nothing ? C_iel : C_groups[iel][1][2])
            else                                          # :nodal: each node group with its own g's matrix
                Sg = Matrix{Float64}(undef, ngl^3, nd)
                for (nodes, Cg) in C_groups[iel]; Sg[nodes, :] = Imat[nodes, :] * transpose(Cg); end
                Sg
            end
            for dk = 1:nd
                extra_source_all[(iel,dk)] = q_ext_all[(iel,dk)] .+ σ_loc .* view(S,:,dk)
            end
        end
        maybe_young_gc!()                                     # source-build temporaries (Imat, S, broadcasts)
        # Lagged iterate read across periodic faces; with prev_override some entries are replaced by fixed inflow
        I_prev = I_cg_all
        if case.prev_override !== nothing
            I_prev = copy(I_cg_all)                                   # shallow: only overridden vectors are copied
            for (key, (lids, vals)) in case.prev_override
                v = get(I_prev, key, nothing); v === nothing && continue
                w = copy(v); @inbounds for q in eachindex(lids); w[lids[q]] = vals[q]; end
                I_prev[key] = w
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
            extra_source_all=extra_source_all, I_cg_all_prev=I_prev, elem_cache=cache,
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
            # monitor (optional): monitor(kind, residual, nmv, I_dict, κ, σ) at every GMRES iterate (gmres_core)
            gmon = monitor === nothing ? nothing : (kind, res, x) -> begin
                Id = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(Id, x, order, len_each)
                monitor(kind, res, length(t_tot), Id, κ, σ)
            end
            x_init = flatten_I(I_cg_init, order); empty!(I_cg_init)         # the Dict copy is not needed any more
            x_final = gmres_fixed_point_solve(apply_G, x_init, comm;
                                               restart=gmres_restart, max_restarts=gmres_max_restarts, tol=tol, rtol=rtol,
                                               monitor=gmon, basis=gmres_basis)
        elseif solver == :gmres_precond
            κ_mean = sum(κ)/max(length(κ),1)
            periodic_axes = detect_periodic_axes(mesh)
            isempty(periodic_axes) && error("gmres_precond requested but mesh has no periodic-tagged faces at all")
            scales_flat = build_scales_flat(compute_direction_preconditioner_scales(mesh, κ_mean, periodic_axes), order, len_each)
            x_final = gmres_fixed_point_solve_precond(apply_G, flatten_I(I_cg_init, order), scales_flat, comm;
                                                   restart=gmres_restart, max_restarts=gmres_max_restarts, tol=tol)
        elseif solver == :jacobi_extrap || solver == :jacobi_extrap_block
            mode = solver == :jacobi_extrap ? :global : :block
            x_final = jacobi_with_extrapolation(apply_G, flatten_I(I_cg_init, order), comm;
                                             len_each=len_each, max_iter=max_iter, tol=tol, mode=mode)
        elseif solver == :jacobi                                               # plain Jacobi (source iteration)
            # x ← G(x) on flat vectors: two solution-sized vectors, no Krylov basis (memory-lean alternative to GMRES).
            # x_{k+1} - x_k = G(x_k) - x_k is exactly the residual of (I - G)x = b that GMRES tests, so the stop is the
            # same: ‖x_{k+1} - x_k‖₂ < max(tol, rtol·‖b‖) with b = G(0) (one extra sweep, as in GMRES) when rtol > 0.
            # ρ = ratio of successive residuals; err_est = res·ρ/(1-ρ) is the estimated distance to the fixed point.
            x = flatten_I(I_cg_init, order); empty!(I_cg_init)
            gnorm(v) = sqrt(MPI.Allreduce(sum(abs2, v), MPI.SUM, comm))
            tol_eff = tol
            if rtol > 0
                bn = gnorm(apply_G(zeros(length(x))))
                tol_eff = max(tol, rtol * bn)
                rank == 0 && println("  Jacobi tolerance: max(abs $(tol), rtol $(rtol) × ‖b‖ $(bn)) = $(tol_eff)")
            end
            res_prev = NaN; conv = false
            for it = 1:max_iter
                xn = apply_G(x)
                dm, ss, sx = _jacobi_step_stats(xn, x)           # function barrier: no allocation (GC is off here)
                x = xn; xn = Float64[]                           # old iterate is garbage from here
                dmax = MPI.Allreduce(dm, MPI.MAX, comm)
                rs = MPI.Allreduce([ss, sx], MPI.SUM, comm); res = sqrt(rs[1]); xnrm = sqrt(rs[2])
                ρ = isnan(res_prev) ? NaN : res / res_prev; res_prev = res
                err_est = (isnan(ρ) || ρ >= 1) ? NaN : res * ρ / (1 - ρ)
                rank == 0 && println(@sprintf("  Jacobi iter %d: residual %.4e (rel ‖x‖ %.3e), max |Δ| %.4e, ρ %.4f, est. error %.3e (rel %.3e)",
                                              it, res, res / max(xnrm, 1e-300), dmax, ρ, err_est, err_est / max(xnrm, 1e-300)))
                if monitor !== nothing
                    Id = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(Id, x, order, len_each)
                    monitor(:jacobi, res, length(t_tot), Id, κ, σ)
                end
                if res < tol_eff
                    rank == 0 && println("  Jacobi converged after $it iterations ($(length(t_tot)) sweeps)")
                    conv = true; break
                end
            end
            conv || (rank == 0 && println("  Jacobi NOT converged after $max_iter iterations"))
            x_final = x
        else
            error("unknown RT solver $(solver)")
        end
        if final_limited_sweep
            # One more sweep with the converged scattering source and the positivity limiter on in every
            # (element, direction) solve, so downstream elements receive the limited inflow. GMRES itself stays
            # linear. Reports how active the limiter was (to show it fades with resolution).
            # limited_sweeps > 1 (diagnostic, :RT_sw_limited_sweeps): repeat the limited sweep from its own result,
            # so periodic wrap inflow and scattering source also come from the limited field.
            x_lim = x_final; st = zeros(7)
            if limiter_mode == :posthoc
                # Post-hoc: the same Zhang-Shu scaling applied to every (element, direction) block of the converged
                # GMRES solution, without a sweep. A horizontally uniform solution stays uniform (the in-sweep
                # limiter feeds limited inflow, which is not a discrete transport solution, to the next element).
                cache.lim_pairs = 0; cache.lim_zero = 0; cache.lim_moved = 0.0; cache.lim_total = 0.0
                fill!(cache.lim_elem, 0)
                x_lim = copy(x_final); blk = zeros(len_each)
                for (q, key) in enumerate(order)
                    r = (q-1)*len_each+1 : q*len_each
                    copyto!(blk, 1, x_lim, first(r), len_each)
                    positivity_limit!(blk, cache, key[1], ngl, ω, metrics.Je)
                    copyto!(x_lim, first(r), blk, 1, len_each)
                end
                st = MPI.Allreduce(Float64[cache.lim_pairs, cache.lim_zero, cache.lim_moved, cache.lim_total,
                                           length(order), count(>(0), cache.lim_elem), mesh.nelem], MPI.SUM, comm)
                dn = MPI.Allreduce(Float64[sum(abs2, x_lim .- x_final), sum(abs2, x_final)], MPI.SUM, comm)
                rank == 0 && println(@sprintf("  positivity limiter (post-hoc): %d of %d pairs limited (%.2e), %d with mean <= 0; elements touched %d of %d (%.2e); radiance moved %.2e of total; limited vs solver solution L2 rel %.2e",
                    Int(st[1]), Int(st[5]), st[1] / max(st[5], 1), Int(st[2]), Int(st[6]), Int(st[7]), st[6] / max(st[7], 1),
                    st[3] / max(st[4], 1e-300), sqrt(dn[1] / max(dn[2], 1e-300))))
            end
            limiter_mode == :sweep || limiter_mode == :posthoc || error("unknown limiter_mode $(limiter_mode)")
            for ls = 1:(limiter_mode == :sweep ? max(limited_sweeps, 1) : 0)
                cache.limit = true; cache.lim_pairs = 0; cache.lim_zero = 0; cache.lim_moved = 0.0; cache.lim_total = 0.0
                fill!(cache.lim_elem, 0)
                x_prev = x_lim
                x_lim = apply_G(x_prev)
                cache.limit = false
                st = MPI.Allreduce(Float64[cache.lim_pairs, cache.lim_zero, cache.lim_moved, cache.lim_total,
                                           length(order), count(>(0), cache.lim_elem), mesh.nelem], MPI.SUM, comm)
                dn = MPI.Allreduce(Float64[sum(abs2, x_lim .- x_final), sum(abs2, x_final), sum(abs2, x_lim .- x_prev)], MPI.SUM, comm)
                rank == 0 && println(@sprintf("  positivity limiter (final sweep %d): %d of %d pairs limited (%.2e), %d with mean <= 0; elements touched %d of %d (%.2e); radiance moved %.2e of total; final sweep vs GMRES iterate L2 rel %.2e; vs previous sweep %.2e",
                    ls, Int(st[1]), Int(st[5]), st[1] / max(st[5], 1), Int(st[2]), Int(st[6]), Int(st[7]), st[6] / max(st[7], 1),
                    st[3] / max(st[4], 1e-300), sqrt(dn[1] / max(dn[2], 1e-300)), sqrt(dn[3] / max(dn[2], 1e-300))))
            end
            LIMITED_ELEMS[] = copy(cache.lim_elem)
            Iu = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(Iu, x_final, order, len_each)
            UNLIMITED_SOLUTION[] = Iu
            x_final = x_lim
        end
        I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}(); unflatten_I!(I_cg_all, x_final, order, len_each)
        result = (I_cg_all, κ, σ)
    finally
        GC.enable(gc_prev)
    end
    t_solve = (time_ns() - t_solve0) * 1e-9
    if USE_FASTDIAG[]
        fs = MPI.Allreduce(Float64[FD_STATS.n_fd, FD_STATS.n_fallback, FD_STATS.n_skip, FD_STATS.its], MPI.SUM, comm)
        fc = MPI.Allreduce(Float64[FD_STATS.max_check, FD_STATS.n_checked], MPI.MAX, comm)
        rank == 0 && println(@sprintf("  fast-diag element solves: %d fast (mean %.1f GMRES its), %d not converged -> dense, %d not eligible -> dense; self-check max rel diff %.1e over %d pairs/rank%s",
            Int(fs[1]), fs[4] / max(fs[1], 1), Int(fs[2]), Int(fs[3]), fc[1], Int(fc[2]), FD_DISABLED[] ? " (DISABLED after a failed check)" : ""))
        fw = MPI.Allreduce(Float64.(FD_SKIP_WHY), MPI.SUM, comm)
        rank == 0 && println("    not eligible, by reason: " * join(["$(FDW_NAMES[i]) $(Int(fw[i]))" for i in eachindex(fw) if fw[i] > 0], ", "))
        ftn = MPI.Allreduce(Float64[FD_STATS.n_tan, FD_STATS.its_tan, FD_STATS.n_fallback_tan], MPI.SUM, comm)
        ftc = MPI.Allreduce(Float64[FD_STATS.n_checked_tan], MPI.MIN, comm)
        rank == 0 && println(@sprintf("    of which tangential (direction parallel to element faces): %d fast (mean %.1f GMRES its), %d not converged; self-checked >= %d pairs/rank",
            Int(ftn[1]), ftn[2] / max(ftn[1], 1), Int(ftn[3]), Int(ftc[1])))
        fd_stats_reset!()
    end
    nyg = MPI.Allreduce(Float64(GC_YOUNG_COUNT[]), MPI.MAX, comm); GC_YOUNG_COUNT[] = 0
    rank == 0 && GC_ALLOC_LIMIT[] > 0 && println("  allocation-volume young GCs this solve (max over ranks): $(Int(nyg)), limit $(GC_ALLOC_LIMIT[]/1e9) GB")
    gsc = MPI.Allreduce(Float64[GC_SYNC_COUNT[], GC_SYNC_SKIP[]], MPI.MAX, comm); GC_SYNC_COUNT[] = 0; GC_SYNC_SKIP[] = 0
    rank == 0 && println("  scheduled full GCs this solve (max over ranks): $(Int(gsc[1])) done, $(Int(gsc[2])) skipped (< $(GC_SYNC_MIN[]/1e9) GB allocated since the previous)")

    # ---- report (collective) ----
    nmv = length(t_tot)
    # repartitioning weights: predicted element work of this solve (per sweep) times its sweeps
    # dense pairs held in the LU cache cost a full solve once, then RTC_DENSE_HIT per later sweep
    cached_cnt = zeros(Int, Int(mesh.nelem))
    if USE_FASTDIAG[]
        for k in keys(cache.factors); cached_cnt[k[1]] += 1; end
    end
    if RT_W_ON[] && cost_pred !== nothing && length(RT_W_ELEM) == length(cost_pred.w)
        RT_W_ELEM .+= nmv .* cost_pred.w .- max(nmv - 1, 0) .* cached_cnt .* (RTC_DENSE[] - RTC_DENSE_HIT[])
    end
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
    if report !== nothing                                                      # for the RT timing record (rt_timing.jl)
        report[:nmv] = nmv; report[:solve_time] = t_solve; report[:setup_time] = t_setup
        report[:matvec_times] = tot_dummy(T, nmv)                              # max over ranks, per matvec
        report[:src_times] = T[1:nmv]; report[:sweep_times] = T[nmv+1:2nmv]
        # load balance of the sweeps: each rank's total sweep time over the steady matvecs, min/max over ranks
        st_ = nmv >= 3 ? (3:nmv) : (1:nmv)
        my_sw = isempty(t_swp) ? 0.0 : sum(t_swp[st_])
        report[:sweep_rank_min] = MPI.Allreduce(my_sw, MPI.MIN, comm)
        report[:sweep_rank_max] = MPI.Allreduce(my_sw, MPI.MAX, comm)
        report[:cache_mb] = MPI.Allreduce(cache.bytes/1e6, MPI.MAX, comm)
        report[:rss_max_gb] = rss
    end
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


# ---- repartitioning by cost (:RT_rebalance => true) ------------------------------------------------------------
# Each RT call accumulates the predicted element work of its solves (rt_predict_costs × sweeps) per element; with a
# dynamics term (wall time between RT calls, per element) this gives the cost weight of every element over one RT
# interval. When the predicted max/mean over ranks exceeds :RT_rebalance_threshold (default 1.15), and at least
# :RT_rebalance_min_interval RT calls (default 1) have passed since the last repartition, RT_REBALANCE_DUE is set:
# the time loop stops the integrator and calls rt_rebalance_strategy! (p4est weighted partition, state transfer).
const RT_W_ON = Ref(false)
const RT_W_ELEM = Float64[]
const RT_REBALANCE_DUE = Ref(false)
const RT_CALLS_SINCE_REBAL = Ref(0)
const RT_T_END_PREV = Ref{UInt64}(0)
const RT_T_DYN = Ref(0.0)

function rt_weights_begin!(mesh, inputs)
    RT_W_ON[] = get(inputs, :RT_rebalance, false) == true
    RT_W_ON[] || return nothing
    RT_T_DYN[] = RT_T_END_PREV[] == 0 ? 0.0 : (time_ns() - RT_T_END_PREV[]) * 1e-9
    resize!(RT_W_ELEM, Int(mesh.nelem)); fill!(RT_W_ELEM, 0.0)
    return nothing
end

function rt_rebalance_check!(mesh, inputs, comm)          # collective
    RT_W_ON[] || return nothing
    rank = MPI.Comm_rank(comm)
    RT_CALLS_SINCE_REBAL[] += 1
    nel = Int(mesh.nelem)
    # dynamics: the wall time between RT calls is set by the rank with the most elements
    td = MPI.Allreduce(RT_T_DYN[], MPI.MAX, comm); nmax = MPI.Allreduce(nel, MPI.MAX, comm)
    cdyn = td > 0 ? td * 1e6 / max(nmax, 1) : 0.0                            # µs per element per RT interval
    w = RT_W_ELEM .+ cdyn
    S = MPI.Allgather(sum(w) * 1e-6, comm); R = MPI.Allgather(sum(RT_W_ELEM) * 1e-6, comm); np = length(S)
    imb = maximum(S) / max(sum(S) / np, 1e-300)
    thr = Float64(get(inputs, :RT_rebalance_threshold, 1.15)); PART_THRESHOLD[] = thr
    minint = Int(get(inputs, :RT_rebalance_min_interval, 1))
    PART_W_OLD[] = w; PART_W_MESH[] = mesh
    RT_REBALANCE_DUE[] = imb > thr && RT_CALLS_SINCE_REBAL[] >= minint
    rank == 0 && println(@sprintf(" # RT load balance (predicted): cost per RT interval max/mean %.3f (threshold %.2f); RT element work max %.1f s, mean %.1f s; dynamics %.1f s per interval (%.0f µs/element)%s",
                                  imb, thr, maximum(R), sum(R) / np, td, cdyn, RT_REBALANCE_DUE[] ? " -> repartition requested" : ""))
    RT_T_END_PREV[] = time_ns()
    return nothing
end

# Before a mesh rebuild (AMR step or cost repartition): drop everything RT keeps for the old mesh, so it is not held
# while the new mesh and its RT infrastructure are built (the plan caches are keyed by identity of old objects).
function rt_release_state!(params)
    mp = params.mp
    hasproperty(mp, :rt_infra) && (mp.rt_infra = nothing)
    hasproperty(mp, :rt_sol_sw) && (mp.rt_sol_sw = [])
    hasproperty(mp, :rt_sol_lw) && (mp.rt_sol_lw = [])
    invalidate_angular_caches!()
    FD_AFFINE[] = nothing
    LIMITED_ELEMS[] = Int32[]; UNLIMITED_SOLUTION[] = nothing
    return nothing
end

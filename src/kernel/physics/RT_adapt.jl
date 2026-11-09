# Regrid-aware RT solve loop.  NOT compiled or run.  Needs angular_tail_indicator.jl and adaptive_angular_regrid_patched.jl.
#
# Flow inside get_rt_heat_fluxes_sweeps!, once per RT call (when inputs[:RT_adaptive_angular] == true):
#   1. (not the first call) indicator on the PREVIOUS shortwave solution, which lives on the current angular mesh
#   2. adaptive_angular_regrid!  (a-priori criterion from g, σ, F_dir  +  indicator levels/radius, halo, per-class radius)
#   3. remap the shortwave AND longwave warm starts onto the new angular mesh (nearest old direction)
#   4. drop every global cache keyed by the old angular mesh
#   5. build the cases and solve on the new infrastructure as before
# The longwave solve shares the shortwave angular mesh; the indicator only looks at the shortwave solution (the sun lobe).

mutable struct RTAdaptState
    levels::Dict{Int,Int}          # gel => level after the last regrid (global on every rank)
    radius::Dict{Int,Float64}      # gel => radius used (local gels)
    ncalls::Int
end
const RT_ADAPT = RTAdaptState(Dict{Int,Int}(), Dict{Int,Float64}(), 0)
reset_rt_adapt!() = (empty!(RT_ADAPT.levels); empty!(RT_ADAPT.radius); RT_ADAPT.ncalls = 0; nothing)

# Every global table that is keyed by, or derived from, the angular mesh.  If you have other module-level caches, add them here:
#   grep -n "^const .*\(Dict\|Ref\)" on the sweep/ncf files.
function invalidate_angular_state!()
    clear_face_plans!()
    for nm in (:ASYNC_LAYOUTS, :ASYNC_DEPS, :GS_CYCLE_CACHE, :ANGULAR_SAME_PLANS, :ANGULAR_REMOTE_PLANS)
        if isdefined(@__MODULE__, nm)
            x = getfield(@__MODULE__, nm)
            x isa AbstractDict && empty!(x)
        end
    end
    return nothing
end

# Warm start onto the new angular mesh: every new direction takes the vector of the nearest OLD direction (unit-vector dot product),
# so unchanged directions are copied exactly and new directions inside a refined element are first-order accurate.
function remap_warm_start(old::Dict{Tuple{Int,Int},Vector{Float64}}, snapshots, mesh)
    new = Dict{Tuple{Int,Int},Vector{Float64}}()
    n3 = mesh.ngl^3
    for iel = 1:mesh.nelem
        snap = snapshots[iel]; em = mesh.extra_mesh[iel]
        no = snap.npoin; nn = em.extra_npoin
        ox = Matrix{Float64}(undef, 3, no)
        for m = 1:no
            θ = snap.coords[1,m]; ϕ = snap.coords[2,m]
            ox[1,m] = sin(θ)*cos(ϕ); ox[2,m] = sin(θ)*sin(ϕ); ox[3,m] = cos(θ)
        end
        for d = 1:nn
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            ux = sin(θ)*cos(ϕ); uy = sin(θ)*sin(ϕ); uz = cos(θ)
            best = 1; bd = -Inf
            for m = 1:no
                c = ox[1,m]*ux + ox[2,m]*uy + ox[3,m]*uz
                c > bd && (bd = c; best = m)
            end
            new[(iel,d)] = haskey(old, (iel,best)) ? copy(old[(iel,best)]) : zeros(n3)
        end
    end
    return new
end

# Convention hook: (θ,ϕ) of the direction the solar beam TRAVELS, i.e. the centre of the forward-scattering lobe, in the same angular
# coordinates as mesh.extra_mesh (θ polar from +z, ϕ azimuth).  In the manufactured tests θs=2.4 (Ωz<0, downward).
# Replace the body with whatever compute_direct_beam / your sun_direction helper uses, so the refined cap sits on the beam.
rt_sun_angles(sw) = error("define rt_sun_angles(sw): return (θ_beam, ϕ_beam) consistent with compute_direct_beam")

function rt_adaptive_regrid_step!(mesh, params, infra, comm, inputs, phase_fn, g_node, σ_node, F_node, θs, ϕs,
                                   sols::Vector; have_solution::Bool)
    rank = MPI.Comm_rank(comm)
    st = RT_ADAPT; st.ncalls += 1
    (st.ncalls - 1) % get(inputs, :RT_regrid_every, 1) == 0 || return infra, sols

    ind_tol   = get(inputs, :RT_ind_tol, 1e-3)       # ≈ e_target / (2 I_ref); do not go below the spatial error floor
    max_level = get(inputs, :RT_max_level, 2)
    tol_quad  = get(inputs, :RT_tol_quad, 1e-3)
    halo      = get(inputs, :RT_halo, 1)
    lobe_k    = get(inputs, :RT_lobe_k, 2.0)
    coarsen   = get(inputs, :RT_adapt_coarsen, true)

    t0 = time_ns()
    level_floor = nothing; radius_floor = nothing
    if have_solution
        level_floor, dgn = angular_tail_levels(mesh, sols[1], st.levels, comm; tol_rel=ind_tol, max_level=max_level,
                                               allow_coarsen=coarsen, sun=(θs, ϕs), radius_cur=st.radius)
        radius_floor = dgn.radius
    end
    t_ind = (time_ns() - t0)*1e-9

    t1 = time_ns()
    builder = make_angular_mesh_builder(inputs, inputs[:backend])
    infra_new, snapshots, info = adaptive_angular_regrid!(mesh, params, infra, comm, phase_fn, g_node, σ_node, F_node, θs, ϕs, builder;
        tol_quad=tol_quad, max_level=max_level, lobe_k=lobe_k, level_floor=level_floor, radius_floor=radius_floor,
        radius_mode=:per_class, halo=halo, check_periodic=get(inputs, :RT_adapt_check, false))
    t_regrid = (time_ns() - t1)*1e-9

    t2 = time_ns()
    new_sols = Any[remap_warm_start(s, snapshots, mesh) for s in sols]
    t_remap = (time_ns() - t2)*1e-9

    st.levels = info.levels; st.radius = info.radius
    invalidate_angular_state!()

    tmax(t) = MPI.Allreduce(t, MPI.MAX, comm)
    ti, tr, tm = tmax(t_ind), tmax(t_regrid), tmax(t_remap)
    if rank == 0
        hist = zeros(Int, max_level+1)
        for l in values(info.levels); hist[min(l, max_level)+1] += 1; end
        @info "[RT adaptive] regrid $(st.ncalls): dirs $(info.n_dirs_before) -> $(info.n_dirs_after)  levels(0..$max_level)=$hist  " *
              "indicator $(round(ti; digits=2)) s, regrid+infra $(round(tr; digits=2)) s, warm-start remap $(round(tm; digits=2)) s"
    end
    return infra_new, new_sols
end

# ------------------------------------------------------------------------------------------------
# get_rt_heat_fluxes_sweeps!: your function with the adaptive block added (marked ### NEW).  Everything else unchanged.
function get_rt_heat_fluxes_sweeps!(q, qe, mesh, micro, metrics, atmos_data, params, dψ, ψ, ω, PhysConst, inputs;
                                     rt_max_iter::Int=200, rt_tol::Float64=1e-5)
    comm = get_mpi_comm()
    bdy = make_boundary_predicates(mesh)

    dycore_to_atmos_data!(q, qe, mesh.npoin, micro.Tabs, micro.qc, micro.qi, inputs[:lmoist],
                           atmos_data, PhysConst, inputs[:SOL_VARS_TYPE])

    if micro.rt_nelem_saved != 0
        local_changed = micro.rt_nelem_saved != mesh.nelem ? 1 : 0
        any_changed = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(local_changed, MPI.MAX, comm) : local_changed
        if any_changed > 0
            @info "[RT warm-start] Spatial mesh changed: rebuilding NCF infrastructure, discarding warm start"
            micro.rt_infra = build_unified_ncf_infrastructure(mesh, params, comm)
            micro.rt_sol_sw_available = false
            micro.rt_sol_lw_available = false
            reset_rt_adapt!()                                                     ### NEW: levels are keyed by gel, invalid now
        end
    else
        micro.rt_infra = build_unified_ncf_infrastructure(mesh, params, comm)
        reset_rt_adapt!()                                                         ### NEW
    end
    infra = micro.rt_infra

    κ_sw, σ_sw, g_eff = atmos_to_rad_shortwave(atmos_data, mesh.npoin)
    κ_ext_sw = κ_sw .+ σ_sw
    z_prof, τ_from_TOA = build_sw_lateral_bc_profile(mesh, κ_ext_sw, mesh.ngl)
    sw = SWParams(inputs[:RT_S0_flux], inputs[:RT_μ0], inputs[:RT_ϕ0], 0.35, z_prof, τ_from_TOA)
    F_dir, G_dir, Q_dir = compute_direct_beam(mesh, κ_sw, σ_sw, sw)

    ### NEW: adaptive angular regrid at the start of the RT solve, warm starts remapped, caches invalidated
    if get(inputs, :RT_adaptive_angular, false)
        θs, ϕs = rt_sun_angles(sw)
        have_sol = micro.rt_sol_sw_available && micro.rt_sol_lw_available
        sols = have_sol ? Any[micro.rt_sol_sw, micro.rt_sol_lw] : Any[]
        micro.rt_infra, new_sols = rt_adaptive_regrid_step!(mesh, params, infra, comm, inputs, phase_function_hg,
                                       g_eff, σ_sw, F_dir, θs, ϕs, sols; have_solution=have_sol)
        infra = micro.rt_infra
        if have_sol
            micro.rt_sol_sw = new_sols[1]
            micro.rt_sol_lw = new_sols[2]
        end
    end
    ### END NEW

    #=sw_case = build_shortwave_diffuse_case(κ_sw, σ_sw, g_eff, sw, F_dir, bdy; surface_albedo=inputs[:RT_surface_albedo])

    I_cg_sw, κ_sw_used, σ_sw_used = run_unified_case_solve(infra, metrics, comm, sw_case;
        max_iter=rt_max_iter, tol=rt_tol,
        solver=:gmres,
        gmres_restart=60, gmres_max_restarts=10,
        I_cg_all_init = micro.rt_sol_sw_available ? micro.rt_sol_sw : nothing)
    micro.rt_sol_sw = I_cg_sw

    Q_sw, dTdt_sw, F_net_sw, G_sw = compute_rt_radiative_heating_unified(
        mesh, I_cg_sw, atmos_data, κ_sw_used, σ_sw_used, false, comm;
        G_dir=G_dir, Q_dir=Q_dir, sw_μ₀=sw.μ₀, sw_φ₀=sw.φ₀)

    if !(inputs[:energy_equation] == "theta") || inputs[:lmoist]
        micro.flux_sw .= PhysConst.cp .* dTdt_sw
    else
        micro.flux_sw .= dTdt_sw
    end=#

    κ_lw, σ_lw = atmos_to_rad_longwave(atmos_data, mesh.npoin)
    lw = LWParams(inputs[:RT_ϵ_surface], inputs[:RT_T_space])
    lw_case = build_longwave_case(atmos_data, κ_lw, σ_lw, bdy; ε_surface=lw.ε_surface, T_space=lw.T_space)

    I_cg_lw, κ_lw_used, σ_lw_used = run_unified_case_solve(infra, metrics, comm, lw_case;
        max_iter=rt_max_iter, tol=rt_tol,
        solver=:gmres,
        gmres_restart=60, gmres_max_restarts=10,
        I_cg_all_init = micro.rt_sol_lw_available ? micro.rt_sol_lw : nothing)
    micro.rt_sol_lw = I_cg_lw

    Q_lw, dTdt_lw, F_net_lw, G_lw = compute_rt_radiative_heating_unified(
        mesh, I_cg_lw, atmos_data, κ_lw_used, σ_lw_used, true, comm)

    if !(inputs[:energy_equation] == "theta") || inputs[:lmoist]
        micro.flux_lw .= PhysConst.cp .* dTdt_lw
    else
        micro.flux_lw .= dTdt_lw
    end

    micro.rt_sol_sw_available = true
    micro.rt_sol_lw_available = true
    micro.rt_nelem_saved = mesh.nelem
end

# A posteriori angular-resolution indicator. NOT compiled or run.
# Uses the previous solution (the warm start) so refinement follows the solution, not only local g.
#
# For every angular element of every spatial element it takes the nodal values at each of the ngl^3 spatial
# nodes, converts them to a tensor Legendre expansion in the element's local (theta,phi) coordinates, and
# measures the energy in the highest-degree modes (a==n or b==n).  Persson-Peraire style, but the tail is
# divided by the GLOBAL peak |I|, so a lobe of small amplitude does not trigger refinement (amplitude-aware).
#
#   eta_el = max over angular elements and spatial nodes of  sqrt(sum_top c^2) / I_ref
#   p      = polynomial degree (n-1);  one more level reduces the top-mode tail by about 2^p asymptotically
#   level  = clamp( cur + ceil( log(eta/tol) / log(2^p) ), 0, max_level )   if eta > tol
#          = cur - 1 if eta < tol / 2^p  (only if allow_coarsen)           else cur
#
# I_cg_all   : Dict{(iel,dir) => Vector (ngl^3, lidx order)} on the CURRENT angular mesh
# levels_cur : Dict{gel => level} of the current mesh (info.levels from the previous regrid; empty = all 0)
# returns    : Dict{gel => level}, plus a named tuple of diagnostics
function legendre_vander(ξ::Vector{Float64})
    n = length(ξ); V = zeros(n, n)
    for i = 1:n
        V[i,1] = 1.0
        n > 1 && (V[i,2] = ξ[i])
        for k = 2:n-1
            V[i,k+1] = ((2k-1)*ξ[i]*V[i,k] - (k-1)*V[i,k-1]) / k
        end
    end
    # normalise so that the coefficients are L2-orthonormal on [-1,1]
    for k = 1:n; V[:,k] .*= sqrt((2(k-1)+1)/2); end
    return V
end

function angular_tail_levels(mesh, I_cg_all, levels_cur::Dict{Int,Int}, comm;
                              tol_rel::Float64=1e-3, max_level::Int=3, allow_coarsen::Bool=false,
                              sun::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                              radius_cur::Dict{Int,Float64}=Dict{Int,Float64}())
    ngl = mesh.ngl; n3 = ngl^3
    I_ref = 0.0
    for v in values(I_cg_all); I_ref = max(I_ref, maximum(abs, v; init=0.0)); end
    I_ref = MPI.Allreduce(I_ref, MPI.MAX, comm)
    I_ref == 0.0 && return Dict{Int,Int}(Int(mesh.el2gel[i]) => 0 for i = 1:mesh.nelem), (eta=zeros(mesh.nelem), I_ref=0.0, radius=Dict{Int,Float64}())

    eta = zeros(mesh.nelem)
    out = Dict{Int,Int}(); rad = Dict{Int,Float64}()
    Ωs = sun === nothing ? (0.0,0.0,0.0) : (sin(sun[1])*cos(sun[2]), sin(sun[1])*sin(sun[2]), cos(sun[1]))
    vinv_cache = Dict{Vector{Float64},Matrix{Float64}}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for e = 1:em.extra_nelem
            n = em.extra_nop[e] + 1
            θ = [em.extra_coords[1, em.extra_connijk[e,a,1]] for a = 1:n]
            ϕ = [em.extra_coords[2, em.extra_connijk[e,1,b]] for b = 1:n]
            for b = 2:n; ϕ[b] < ϕ[1] - 1e-8 && (ϕ[b] += 2π); end          # unwrap the 2π seam
            ξθ = 2 .* (θ .- θ[1]) ./ (θ[end] - θ[1]) .- 1
            ξϕ = 2 .* (ϕ .- ϕ[1]) ./ (ϕ[end] - ϕ[1]) .- 1
            Vθi = get!(() -> inv(legendre_vander(ξθ)), vinv_cache, ξθ)
            Vϕi = get!(() -> inv(legendre_vander(ξϕ)), vinv_cache, ξϕ)
            Iab = Matrix{Float64}(undef, n, n); te = 0.0
            for p = 1:n3
                for a = 1:n, b = 1:n
                    Iab[a,b] = I_cg_all[(iel, em.extra_connijk[e,a,b])][p]
                end
                c = Vθi * Iab * transpose(Vϕi)                                # c[k,l]: degree (k-1,l-1)
                t = 0.0
                for k = 1:n, l = 1:n
                    (k == n || l == n) && (t += c[k,l]^2)
                end
                eta[iel] = max(eta[iel], sqrt(t) / I_ref); te = max(te, sqrt(t) / I_ref)
            end
            if sun !== nothing && te > tol_rel                             # this angular element needs resolution: the sun-centred radius must reach its farthest node
                for a = 1:n, b = 1:n
                    ip = em.extra_connijk[e,a,b]; θn = em.extra_coords[1,ip]; ϕn = em.extra_coords[2,ip]
                    dd = acos(clamp(Ωs[1]*sin(θn)*cos(ϕn) + Ωs[2]*sin(θn)*sin(ϕn) + Ωs[3]*cos(θn), -1.0, 1.0))
                    rad[Int(mesh.el2gel[iel])] = max(get(rad, Int(mesh.el2gel[iel]), 0.0), dd)
                end
            end
        end
        gel = Int(mesh.el2gel[iel]); cur = get(levels_cur, gel, 0)
        p = mesh.extra_mesh[iel].extra_nop[1]; rate = 2.0^p                 # top-mode coefficient scales like h^p (checked numerically: asymptotic ratio 8 for p=3)
        lvl = cur
        if eta[iel] > tol_rel
            lvl = cur + ceil(Int, log(eta[iel]/tol_rel) / log(rate))
        elseif allow_coarsen && eta[iel] < tol_rel / rate
            lvl = cur - 1
        end
        out[gel] = clamp(lvl, 0, max_level)
        # keep the previous radius only while the element stays refined: on an already refined mesh the tails are small, so a radius
        # recomputed from them alone would shrink the refined cap and the regrid would oscillate.  Level 0 forgets the radius.
        rad[gel] = out[gel] > 0 ? max(get(rad, gel, 0.0), get(radius_cur, gel, 0.0)) : 0.0
    end
    return out, (eta = eta, I_ref = I_ref, radius = rad)
end

# Patched adaptive_angular_regrid!. NOT compiled or run.
# Changes vs. your version (everything else is identical):
#   1. new keyword  level_floor::Union{Nothing,Dict{Int,Int}}  (gel => absolute target level, from angular_tail_levels,
#      computed on the CURRENT mesh BEFORE this call).  nothing => behaviour is exactly the old one.
#   2. an element is "flagged" if it is a-priori active OR has a positive floor.  Flagged elements get
#      L_loc = max(a-priori level (active only), floor).
#   3. refinement radius: old behaviour (own g) when level_floor===nothing.  With a floor, a-priori-active elements keep
#      their own hw; every other refined element (flagged by the indicator, or raised by balancing) uses hw_ref,
#      the narrowest lobe among the flagged elements, i.e. the lobe that can actually be advected there.
#   6. radius_mode (:per_class default / :uniform / :per_element) and check_periodic: see the comment above the radius block.
#   5. new keyword radius_floor (gel => radians, from angular_tail_levels(...; sun=(θsun,ϕsun))): the refinement radius is
#      never smaller than the farthest node of an angular element whose tail exceeded the tolerance.
#   4. hw is looked up lazily (hw_of); the old hw[bucket(g)] raised KeyError for an element that balance_angular_levels
#      raised from 0 although its g bucket was not among the active buckets.
function adaptive_angular_regrid!(mesh, sem, infra, comm, phase_fn, g_node, σ_node, F_node,
                                   θsun::Float64, ϕsun::Float64, builder::Function;
                                   tol_quad::Float64=1e-3, max_level::Int=3,
                                   strength_tol::Float64=0.05, lobe_k::Float64=2.0,
                                   level_floor::Union{Nothing,Dict{Int,Int}}=nothing,
                                   radius_floor::Union{Nothing,Dict{Int,Float64}}=nothing,
                                   radius_mode::Symbol=:per_class, check_periodic::Bool=false, periodic_fatal::Bool=true,
                                   halo::Int=0)
    ngl = mesh.ngl
    snapshots = Dict{Int,NamedTuple}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        snapshots[iel] = (coords = copy(em.extra_coords[:, 1:em.extra_npoin]), npoin = em.extra_npoin)
    end
    n_before = MPI.Allreduce(sum((s.npoin for s in values(snapshots)); init=0), +, comm)

    s_el = zeros(mesh.nelem); g_el = zeros(mesh.nelem)
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        gip = mesh.connijk[iel,i,j,k]
        s_el[iel] += σ_node[gip]*F_node[gip]/ngl^3; g_el[iel] += g_node[gip]/ngl^3
    end
    s_ref = MPI.Allreduce(maximum(s_el; init=0.0), MPI.MAX, comm)
    active = [s_ref > 0 && s_el[iel] > strength_tol*s_ref for iel = 1:mesh.nelem]
    bucket(g) = round(g; digits=2)

    floor_of(iel) = level_floor === nothing ? 0 : get(level_floor, Int(mesh.el2gel[iel]), 0)
    from_ind = [!active[iel] && floor_of(iel) > 0 for iel = 1:mesh.nelem]     # flagged only by the indicator

    g_buckets = sort(unique(Float64[bucket(g_el[iel]) for iel = 1:mesh.nelem if active[iel]]))
    hw_cache = Dict{Float64,Float64}()
    hw_of(g) = get!(() -> phase_lobe_halfwidth(phase_fn, bucket(g)), hw_cache, bucket(g))
    # narrowest lobe among all flagged elements (a-priori or indicator), reduced over ranks
    hw_loc = minimum((hw_of(g_el[iel]) for iel = 1:mesh.nelem if active[iel] || from_ind[iel]); init=Inf)
    hw_ref = MPI.Allreduce(hw_loc, MPI.MIN, comm)

    reset_angular_mesh!(mesh, builder)
    dθ0 = mesh.nelem > 0 ? base_theta_extent(mesh.extra_mesh[1]) : 0.0
    table = (mesh.nelem > 0 && !isempty(g_buckets)) ?
            sun_row_error_table!(mesh, builder, phase_fn, g_buckets, θsun, ϕsun, lobe_k, max_level, dθ0) :
            Dict{Float64,Vector{Float64}}()
    reset_angular_mesh!(mesh, builder)

    L_loc = Dict{Int,Int}()
    for iel = 1:mesh.nelem
        la = active[iel] ? levels_from_table(table, bucket(g_el[iel]), tol_quad) : 0
        L_loc[Int(mesh.el2gel[iel])] = min(max(la, floor_of(iel)), max_level)
    end

    sp = Set{Int}()
    for (iel,_) in keys(infra.ncf_as_child);  push!(sp, iel); end
    for (iel,_) in keys(infra.ncf_as_parent); push!(sp, iel); end
    gel_to_class = build_periodic_equivalence_classes(mesh, infra.same_rank_kind, infra.same_rank_partner,
        infra.same_rank_is_periodic, infra.neighbor_kind, infra.neighbor_gel_group, infra.neighbor_is_periodic, comm)
    global_adj, global_refined = build_global_element_adjacency_and_refined_set(mesh, infra.neighbor_elem, sp,
        infra.neighbor_rank, infra.neighbor_gel_group, gel_to_class, comm)
    forced_zero = Set{Int}(Int(g) for g in global_refined)
    for g in global_refined, n in get(global_adj, g, ()); push!(forced_zero, Int(n)); end
    L = balance_angular_levels(L_loc, global_adj, gel_to_class, forced_zero, comm)
    # halo: the coarse->fine inflow error at the edge of a refined region is ~8% of the local lobe peak.  Moving the edge outward by
    # `halo` spatial layers puts it where the amplitude is lower.  Each pass: level = max(own, neighbours'), then re-balance
    # (restores periodic-class equality, the one-level jump limit and forced_zero).  L is global on every rank.
    for _ = 1:halo
        L_h = Dict{Int,Int}()
        for iel = 1:mesh.nelem
            g = Int(mesh.el2gel[iel]); lv = L[g]
            for n in get(global_adj, g, ()); lv = max(lv, get(L, Int(n), 0)); end
            L_h[g] = lv
        end
        L = balance_angular_levels(L_h, global_adj, gel_to_class, forced_zero, comm)
    end

    # radius_mode:
    #   :per_element  each element uses its own radius (own g / own indicator).  Periodic partners can then get DIFFERENT meshes.
    #   :per_class    own radius, but equalised (max) over each periodic equivalence class, so periodic partners get identical meshes
    #                 (levels are already equal across a class).  Non-periodic neighbours keep their own radii.  DEFAULT.
    #   :uniform      one radius for every refined element (conservative; mesh = function of level only).
    R_flag = radius_floor === nothing ? 0.0 : MPI.Allreduce(maximum(values(radius_floor); init=0.0), MPI.MAX, comm)
    R_el = Dict{Int,Float64}()
    for iel = 1:mesh.nelem
        hw_use = (level_floor !== nothing && !active[iel] && isfinite(hw_ref)) ? hw_ref : hw_of(g_el[iel])
        r = lobe_k*hw_use
        if radius_floor !== nothing
            rf = get(radius_floor, Int(mesh.el2gel[iel]), 0.0)
            # an element refined only because of the halo has no indicator radius of its own: give it the largest flagged radius
            rf == 0.0 && L[Int(mesh.el2gel[iel])] > 0 && (rf = R_flag)
            r = max(r, rf)
        end
        R_el[Int(mesh.el2gel[iel])] = r
    end
    if radius_mode == :per_class
        R_el = equalize_over_periodic_classes(R_el, gel_to_class, comm)
    elseif radius_mode == :uniform
        Rg = MPI.Allreduce(maximum(values(R_el); init=0.0), MPI.MAX, comm)
        for k in keys(R_el); R_el[k] = Rg; end
    elseif radius_mode != :per_element
        error("radius_mode must be :per_class, :uniform or :per_element")
    end
    for iel = 1:mesh.nelem
        lvl = L[Int(mesh.el2gel[iel])]; lvl == 0 && continue
        R_use = R_el[Int(mesh.el2gel[iel])]
        refine_lobe_to_level!(mesh, iel, θsun, ϕsun, R_use, lvl, dθ0)
    end
    check_periodic && check_periodic_mesh_consistency(mesh, gel_to_class, comm; fatal=periodic_fatal)
    infra_new = build_unified_ncf_infrastructure(mesh, sem, comm)
    n_after = MPI.Allreduce(sum((mesh.extra_mesh[i].extra_npoin for i = 1:mesh.nelem); init=0), +, comm)
    return infra_new, snapshots, (levels=L, table=table, n_dirs_before=n_before, n_dirs_after=n_after, radius=R_el)
end


# gel_to_class :: Dict{TInt,Set{TInt}} maps a gel to the Set of ALL members of its periodic class (global); only gels with a periodic
# edge are keys.  Class id = minimum member.  Elements that are not keys are singleton classes and need no exchange.
periodic_class_id(gel_to_class, gel::Int) = Int(minimum(gel_to_class[gel]))

# Periodic partners must have IDENTICAL angular meshes (periodic angular NCF is unsupported).  Signature = hash of the sorted
# (θmin,θmax,ϕmin,ϕmax) of every angular element; every member of a class must share one signature.  Only elements that are keys of
# gel_to_class (i.e. on periodic faces) are examined, and their (class, signature) pairs are gathered on rank 0.  NOT compiled or run.
function check_periodic_mesh_consistency(mesh, gel_to_class, comm; fatal::Bool=true)
    rank = MPI.Comm_rank(comm)
    loc = Int[]
    for iel = 1:mesh.nelem
        gel = Int(mesh.el2gel[iel]); haskey(gel_to_class, gel) || continue
        em = mesh.extra_mesh[iel]; bnds = Tuple{Float64,Float64,Float64,Float64}[]
        for e = 1:em.extra_nelem
            b = get_angular_element_bounds(em.extra_coords, em.extra_connijk, em.extra_nop, e)
            push!(bnds, (round(b[1]; digits=9), round(b[2]; digits=9), round(b[3]; digits=9), round(b[4]; digits=9)))
        end
        push!(loc, periodic_class_id(gel_to_class, gel)); push!(loc, Int(hash(sort(bnds)) >> 12))
    end
    nmax = MPI.Allreduce(length(loc), MPI.MAX, comm)
    nbad = 0
    if nmax > 0
        buf = zeros(Int, nmax); buf[1:length(loc)] .= loc
        all_buf = MPI.Gather(buf, comm; root=0)
        if rank == 0
            sig = Dict{Int,Int}(); nmem = Dict{Int,Int}()
            for k = 1:2:length(all_buf)-1
                c = all_buf[k]; sg = all_buf[k+1]; (c == 0 && sg == 0) && continue
                nmem[c] = get(nmem, c, 0) + 1
                haskey(sig, c) ? (sig[c] == sg || (nbad += 1)) : (sig[c] = sg)
            end
            println("  periodic mesh consistency: $(length(nmem)) classes on periodic faces, $(nbad) element(s) whose angular mesh differs from their class")
        end
        nbad = MPI.bcast(nbad, comm; root=0)
    elseif rank == 0
        println("  periodic mesh consistency: no periodic faces found")
    end
    nbad > 0 && fatal && error("periodic partners have different angular meshes ($nbad); periodic angular NCF is not supported")
    return nbad
end

# max of a per-gel value over each periodic class (all members, all ranks).  Only elements on periodic faces take part.
function equalize_over_periodic_classes(vals::Dict{Int,Float64}, gel_to_class, comm)
    pg = [g for g in keys(vals) if haskey(gel_to_class, g)]
    cls = Int[periodic_class_id(gel_to_class, g) for g in pg]
    v   = Float64[vals[g] for g in pg]
    nmax = MPI.Allreduce(length(pg), MPI.MAX, comm)
    nmax == 0 && return vals
    cb = zeros(Int, nmax); vb = fill(-1.0, nmax); cb[1:length(pg)] .= cls; vb[1:length(pg)] .= v
    call = MPI.Allgather(cb, comm); vall = MPI.Allgather(vb, comm)
    cmax = Dict{Int,Float64}()
    for k in eachindex(call)
        vall[k] < 0 && continue
        cmax[call[k]] = max(get(cmax, call[k], -Inf), vall[k])
    end
    out = copy(vals)
    for (g, c) in zip(pg, cls); out[g] = cmax[c]; end
    return out
end



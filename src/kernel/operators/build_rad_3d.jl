using SparseArrays


# ──────────────────────────────────────────────────────────────────────────────
# Logging helper: only rank 0 prints to avoid redundant multi-process output.
# ──────────────────────────────────────────────────────────────────────────────
macro rankinfo(rank, msgs...)
    quote
        if $(esc(rank)) == 0
            @info $(map(esc, msgs)...)
        end
    end
end

"""
    build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je,
        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        nx, ny, nz, elem_to_face, extra_mesh, κ, σ,
        QT::Inexact, SD::NSD_3D, AD::ContGal)

Top-level driver for the 5D (3D spatial × 2D angular) radiative transfer problem
with optional adaptive angular refinement.

# Overview

Assembles and solves the discrete radiative transfer equation:

    (Ω·∇ + κ) I(x,Ω) = σ ∫ Φ(Ω,Ω') I(x,Ω') dΩ' + Q(x,Ω)

using a continuous Galerkin spectral element discretisation on a tensor-product
spatial × angular mesh.  When `inputs[:adaptive_extra_meshes]` is `true` the
angular mesh is locally refined using an h-adaptivity criterion and the resulting
non-conforming DOFs are handled via a constraint (prolongation/restriction) matrix
approach.

# Adaptive path (inputs[:adaptive_extra_meshes] == true)

1. **Initial numbering** — `adaptive_spatial_angular_numbering_3D_2D!` builds a
   combined spatial-angular connectivity array `connijk_spa` and identifies hanging
   (non-conforming) nodes.

2. **Global numbering** — `setup_global_numbering_adaptive_angular_scalable`
   assigns compact global IDs to all DOFs and builds the owner map for MPI
   partitioning.

3. **Ghost layer** — `build_nonconforming_ghost_layer` identifies off-rank parent
   nodes needed for the constraint application.

4. **Adaptivity** — An element-wise smoothness criterion drives
   `adapt_angular_grid_3Dby2D!`, which refines angular elements where the
   radiation field is insufficiently resolved.  After refinement the numbering
   and ghost layer are rebuilt for the adapted mesh.

5. **Constraint matrices** — `build_restriction_matrices_local_and_ghost` builds:
   - `nc_mat`  (R): restriction  n_free × n_spa  (maps all DOFs to free DOFs)
   - `P`       (P = R'): prolongation  n_spa × n_free
   - Ghost-side equivalents and constraint data for MPI boundary hanging nodes.

6. **Matrix assembly** — `sparse_lhs_assembly_3Dby2D_adaptive` and
   `assemble_mass_diagonal_3Dby2D_adaptive` build the local LHS and lumped mass
   matrix.  `M⁻¹ LHS` is formed using the diagonal mass inverse.

7. **Parallel restriction/prolongation** — The operator
   `A_free = R (M⁻¹ LHS) P`
   is applied in parallel across MPI ranks.  Interface hanging nodes require
   explicit communication of row effects (before R) and column effects (before P)
   via `exchange_hanging_effects`.

8. **RHS assembly and restriction** — `RHS_red = R · b` with inflow boundary
   conditions applied in-place.

9. **Solve** — `solve_parallel_lsqr` wraps Krylov.jl's LSQR with a parallel
   matrix-vector product defined via `create_parallel_linear_operator`.

10. **Prolongation of solution** — `solution_new = P · x_free`, with ghost-parent
    contributions exchanged across ranks.

11. **Post-processing** — Angular integration of the solution and reference field,
    L² error computation, and VTK output.

# Non-adaptive path (inputs[:adaptive_extra_meshes] == false)

Assembles `M⁻¹ LHS` on a uniform angular mesh, applies boundary conditions, and
solves with `solve_parallel_lsqr`.

# Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `mesh` | `Mesh3D` | Spatial mesh (connectivity, coords, boundary info) |
| `inputs` | `Dict` | Problem configuration (see keys below) |
| `neqs` | `Int` | Number of equations (= 1 for scalar RT) |
| `ngl` | `Int` | Number of Gauss–Lobatto points per spatial direction |
| `dψ` | `Array` | Spatial basis derivative matrix `dψ[m,i]` |
| `ψ` | `Array` | Spatial basis matrix `ψ[i,m]` |
| `ω` | `Vector` | 1D Gauss–Lobatto quadrature weights |
| `Je` | `Array` | Spatial Jacobian determinant `Je[iel,i,j,k]` |
| `dξdx ... dζdz` | `Array` | Spatial metric terms `[iel,i,j,k]` |
| `nx, ny, nz` | `Array` | Outward boundary face normals `[iface,i,j]` |
| `elem_to_face` | `Array` | Map from `(iel,i,j,k)` to boundary face index and local coords |
| `extra_mesh` | `Array` or struct | Angular mesh data per spatial element |
| `κ` | `Array` | Absorption coefficient (used when `inputs[:lRT_from_data]` is true) |
| `σ` | `Array` | Scattering coefficient (used when `inputs[:lRT_from_data]` is true) |
| `QT` | `Inexact` | Quadrature type dispatch tag |
| `SD` | `NSD_3D` | Spatial dimension dispatch tag |
| `AD` | `ContGal` | Approximation type dispatch tag |

## Relevant `inputs` keys

| Key | Type | Description |
|-----|------|-------------|
| `:adaptive_extra_meshes` | `Bool` | Enable h-adaptive angular refinement |
| `:rad_HG_g` | `Float64` | Henyey–Greenstein asymmetry parameter g ∈ (-1,1) |
| `:output_dir` | `String` | Directory for VTK output |
| `:lRT_from_data` | `Bool` | Read κ,σ from arrays rather than user functions |

# MPI parallelism

The spatial mesh is partitioned across MPI ranks by the caller.  Each rank owns
a contiguous subset of spatial elements.  Angular DOFs inherit the spatial
partitioning: a spatial-angular DOF `(ip, ip_ang)` is owned by the rank that
owns spatial node `ip`.

Ghost layers are built explicitly for non-conforming angular interfaces that
cross MPI partition boundaries.  All communication uses non-blocking point-to-
point sends (via `exchange_buffers`) or flat `MPI.Alltoallv` collectives.

# Output

Writes a VTK file of the angle-integrated solution and reference field to
`inputs[:output_dir]`.  Prints the L² error to the log on rank 0.

# See also

- [`solve_parallel_lsqr`](@ref) — parallel LSQR wrapper
- [`create_parallel_linear_operator`](@ref) — MPI matvec operator
- [`build_restriction_matrices_local_and_ghost`](@ref) — constraint matrix assembly
- [`setup_global_numbering_adaptive_angular_scalable`](@ref) — compact global numbering
"""
function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        nx, ny, nz, elem_to_face,
        extra_mesh, κ, σ, atmos_data, z_prof, τ_from_TOA, QT::Inexact, SD::NSD_3D, AD::ContGal;
        rt_sol_lw = Float64[], rt_sol_sw = Float64[], rt_sol_lw_available = false, rt_sol_sw_available = false)

    bdy = make_boundary_predicates(mesh)
    sw = []
    lw = []
    connijk_spa = []
    #=sw_base = SWParams(1361.0, 0.6, 0.0, 0.0)
    for δ in [0.25, 0.3, 0.4, 0.5, 0.6]
        sw_test = SWParams(sw_base.S₀_flux, sw_base.μ₀, sw_base.φ₀, δ)
        result  = check_beam_flux(extra_mesh, sw_test, extra_mesh.extra_nop[1])
        
    end=#
    F_dir  :: Union{Vector{Float64}, Nothing} = nothing
    G_dir  :: Union{Vector{Float64}, Nothing} = nothing
    Q_dir  :: Union{Vector{Float64}, Nothing} = nothing
    if(inputs[:lRT_from_data])

        sw  = SWParams(inputs[:RT_S0_flux], inputs[:RT_μ0], inputs[:RT_ϕ0], 0.35, z_prof, τ_from_TOA)#inputs[:RT_δ_beam])
        lw  = LWParams(inputs[:RT_ϵ_surface], inputs[:RT_T_space])
    end


    if (inputs[:RT_shortwave])
        bounds = (xmin=mesh.xmin, xmax=mesh.xmax, ymin=mesh.ymin, ymax = mesh.ymax, zmin = mesh.zmin, zmax = mesh.zmax)
        #=ray_grid = build_ray_grid(sw, bounds, 2*5, 2*5, 4*15)
        # Compute τ and F_dir in parallel
        κ_ext_sw      = κ .+ σ
        F_dir, τ_nodes = compute_tau_ray_parallel(ray_grid, mesh, κ_ext_sw, sw, MPI.COMM_WORLD)=#
        κ_ext_sw      = κ .+ σ
        F_dir, τ_nodes = compute_tau_direct_3D(mesh, κ_ext_sw, sw, MPI.COMM_WORLD)
        G_dir, Q_dir = compute_direct_radiation(F_dir, τ_nodes, κ, σ, sw, mesh)
        ip_toa = argmax(mesh.z)
        ip_sfc = argmin(mesh.z)
        mfp = 1.0 / mean(κ_ext_sw)
        sw_ω₀_lateral  = mean(σ ./ max.(κ_ext_sw, 1e-30))
        @info "SW domain-mean ω₀: $(round(sw_ω₀_lateral, digits=4))"
        @info "Mean free path: $(round(mfp/1000, digits=2)) km"
        @info "Domain width:   $(round((mesh.xmax-mesh.xmin)/1000, digits=2)) km"
        @info "τ at TOA (z=$(round(mesh.z[ip_toa],digits=0))m): $(round(τ_nodes[ip_toa],digits=4))"
        @info "τ at sfc (z=$(round(mesh.z[ip_sfc],digits=0))m): $(round(τ_nodes[ip_sfc],digits=4))"
    end
    nc_mat = zeros(Float64,1)
    P = zeros(Float64,1)
    MLHS_pre = zeros(Float64,1)
    rest = zeros(Float64,1)
    nc_non_global_nodes = []
    n_non_global_nodes = 0
    n_spa = 0
    npoin = mesh.npoin
    nelem = mesh.nelem
    npoin_ang_total = 0
    nprocs = MPI.Comm_size(MPI.COMM_WORLD)
    rank   = MPI.Comm_rank(MPI.COMM_WORLD)
    comm   = MPI.COMM_WORLD

    # ── Stage 1b: Initialize Spatial AMR Cache ──────────────────────────────────
    # Stage 1: Spatial Mesh Infrastructure Integration
    # Initialize cache structure to hold spatial AMR data for constraint assembly.
    spatial_amr_cache = SpatialAMRCache(
        element_refinement_levels = Vector{Int}(mesh.ad_lvl),
        num_spatial_hanging_facets = mesh.num_ncf
    )

    # Compute GLOBAL max to ensure all ranks make same decision (prevent MPI divergence)
    local_has_spatial = mesh.num_ncf > 0 ? 1 : 0
    global_has_spatial = MPI.Allreduce(local_has_spatial, MPI.MAX, MPI.COMM_WORLD)
    has_spatial_hanging_nodes = global_has_spatial == 1
    @info rank, has_spatial_hanging_nodes
    # Stage 2-3 will be called after angular mesh setup (either adaptive or uniform)
    # Cache initialized but constraint building deferred until mesh fully available
    R_spatial = nothing  # Will hold spatial restriction matrix if needed
    P_spatial = nothing  # Will hold spatial prolongation matrix if needed
    spatial_hanging_nodes_all_angular = Set{Int}()  # set for spatial hanging nodes
    if inputs[:adaptive_extra_meshes]
        gip2owner_extra = []
        local_parent_indices = Set{Int}()
        local_non_owned_parents = Set{Int}()
        nonowned_parent_gids = Set{Int}()
        gid_to_extended_parents = Dict{Int, Int}()
        extended_parents_to_gid = Int[]
        extended_parents_x = Float64[]
        extended_parents_y = Float64[]
        extended_parents_z = Float64[]
        extended_parents_θ = Float64[]
        extended_parents_ϕ = Float64[]
        extended_parents_ip = Int[]
        all_hanging_nodes = Set{Int}()
        ghost_constraint_data = Dict{Int, Vector{Tuple{Int, Float64}}}()
        ghost_constraint_data_rhs = Dict{Int, Vector{Tuple{Int, Float64}}}()
        reverse_ghost_map = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()
        # ── Allocate per-element angular mesh arrays ──────────────────────────
        extra_meshes_coords      = [Array{Float64}(undef, size(extra_mesh[e].extra_coords,1),   size(extra_mesh[e].extra_coords,2))   for e in 1:nelem]
        extra_meshes_connijk     = [Array{Int}(undef,     extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_Je    = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dξdx  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dxdξ  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dξdy  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dydξ  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dηdx  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dxdη  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dηdy  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dydη  = [Array{Float64}(undef, extra_mesh[e].extra_nelem,             extra_mesh[e].extra_nop[1]+1,          extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_nops  = [Array{Int}(undef,     extra_mesh[e].extra_nelem)             for e in 1:nelem]
        extra_meshes_extra_npoins = zeros(Int, nelem)
        extra_meshes_extra_nelems = zeros(Int, nelem)
        extra_meshes_ref_level   = [Array{Int}(undef,     extra_mesh[e].extra_nelem)             for e in 1:nelem]

        npoin_ang_total = 0
        for e = 1:nelem
            extra_meshes_coords[e]       = extra_mesh[e].extra_coords[:,:]
            extra_meshes_connijk[e]      = extra_mesh[e].extra_connijk
            extra_meshes_extra_Je[e]     = extra_mesh[e].extra_metrics.Je[:,:,:]
            extra_meshes_extra_dξdx[e]   = extra_mesh[e].extra_metrics.dξdx[:,:,:]
            extra_meshes_extra_dxdξ[e]   = extra_mesh[e].extra_metrics.dxdξ[:,:,:]
            extra_meshes_extra_dξdy[e]   = extra_mesh[e].extra_metrics.dξdy[:,:,:]
            extra_meshes_extra_dydξ[e]   = extra_mesh[e].extra_metrics.dydξ[:,:,:]
            extra_meshes_extra_dηdx[e]   = extra_mesh[e].extra_metrics.dηdx[:,:,:]
            extra_meshes_extra_dxdη[e]   = extra_mesh[e].extra_metrics.dxdη[:,:,:]
            extra_meshes_extra_dηdy[e]   = extra_mesh[e].extra_metrics.dηdy[:,:,:]
            extra_meshes_extra_dydη[e]   = extra_mesh[e].extra_metrics.dydη[:,:,:]
            extra_meshes_extra_npoins[e] = extra_mesh[e].extra_npoin
            extra_meshes_extra_nelems[e] = extra_mesh[e].extra_nelem
            extra_meshes_extra_nops[e]   = extra_mesh[e].extra_nop
            extra_meshes_ref_level[e]    = extra_mesh[e].ref_level
            npoin_ang_total += mesh.ngl * mesh.ngl * extra_mesh[e].extra_npoin
        end

        connijk_spa = [Array{Int}(undef, ngl, ngl, ngl,
                          extra_meshes_extra_nelems[iel],
                          extra_meshes_extra_nops[iel][1]+1,
                          extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]

        neighbors = zeros(Int, nelem, 26, 2)
        adapted   = false

        # ── Spatial coincident node deduplication for global numbering ────────
        # Spatial AMR creates coincident nodes: child element corners at the same
        # (x,y,z) as a parent corner.  Without remapping them to share the parent's
        # GIP, setup_global_numbering_adaptive_angular_scalable creates two compact
        # GIDs for the same location but only one rank claims the DOF, causing the
        # "Unowned global DOFs detected" assertion failure.
        ip2gip_dedup  = Int.(mesh.ip2gip)
        _remapped_ips = Set{Int}()
        if has_spatial_hanging_nodes
            _cache_pre = build_spatial_constraint_matrices(
                mesh, spatial_amr_cache,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                ngl, rank
            )
            _cache_pre = exchange_spatial_ghosts(
                mesh, _cache_pre,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                rank, comm
            )
            for (child_ip, parent_ip) in _cache_pre.coincident_nodes
                ip2gip_dedup[child_ip] = mesh.ip2gip[parent_ip]
                push!(_remapped_ips, child_ip)
            end
            for (child_ip, parent_global_spa) in _cache_pre.cross_rank_coincident_nodes
                ip2gip_dedup[child_ip] = parent_global_spa
                push!(_remapped_ips, child_ip)
            end
        end

        # ── Pre-adaptivity connectivity and numbering ─────────────────────────
        @rankinfo rank "Building initial adaptive spatial-angular connectivity..."
        flush(stdout)
        nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa =
            adaptive_spatial_angular_numbering_3D_2D!(
                connijk_spa, nelem, ngl, mesh.connijk,
                extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                extra_meshes_coords, mesh.x, mesh.y, mesh.z,
                extra_meshes_ref_level, neighbors, adapted, extra_meshes_extra_Je)
        @rankinfo rank "✓ Connectivity built: n_spa=$n_spa, n_non_global=$n_non_global_nodes"
        flush(stdout)

        @rankinfo rank "Building pre-adaptivity global numbering..."
        flush(stdout)
        ip2gip_spa, gip2ip, gip2owner_spa, gnpoin =
            setup_global_numbering_adaptive_angular_scalable(
                ip2gip_dedup, mesh.gip2owner, mesh, connijk_spa,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                n_spa, n_non_global_nodes, nc_non_global_nodes;
                remapped_ips = _remapped_ips)
        @rankinfo rank "✓ Global numbering complete: gnpoin=$gnpoin"
        flush(stdout)

        gip2owner_extra = zeros(Int, n_spa)
        for ip = 1:n_spa
            gip2owner_extra[ip] = gip2owner_spa[ip2gip_spa[ip]]
        end

            # ── Reverse GID lookup (used in many downstream functions) ────────
        gip_to_local = Dict{Int, Int}()
        sizehint!(gip_to_local, n_spa)
        for ip = 1:n_spa
            gip_to_local[ip2gip_spa[ip]] = ip
        end

        @rankinfo rank "Global DOF count (pre-adapt): $gnpoin"
        flush(stdout)

        @rankinfo rank "Building non-conforming ghost layer..."
        flush(stdout)
        ghost_layer = build_nonconforming_ghost_layer(
            mesh, connijk_spa, mesh.ip2gip, ip2gip_spa, gip2owner_spa,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
            extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
            extra_meshes_ref_level, n_spa, neighbors)
        @rankinfo rank "✓ Ghost layer built"
        flush(stdout)

        # ── Pre-adaptivity matrices for adaptivity criterion ─────────────────
        @rankinfo rank "Assembling pre-adaptivity LHS and mass matrix..."
        @info rank
        LHS = sparse_lhs_assembly_3Dby2D_adaptive(
            ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
            mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ,
            extra_meshes_connijk, extra_meshes_extra_Je,
            extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl,
            extra_meshes_extra_nelems, dξdx, dξdy, dξdz,
            dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
            extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, inputs[:lRT_from_data], κ, σ, inputs[:RT_shortwave])

        Md = assemble_mass_diagonal_3Dby2D_adaptive(
            ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
            extra_meshes_extra_Je, extra_meshes_extra_nops, n_spa, nelem,
            ngl, extra_meshes_extra_nelems, connijk_spa)
        
        M_inv = spdiagm(0 => 1.0 ./ Md)
        @info size(Md), size(LHS)
        @info "Md range" minimum(Md), maximum(Md)
        @info "Md mean" sum(Md) / length(Md)
        @info "Large Md entries" count(Md .> 1e4)
        @info "LHS diagonal range" minimum(diag(LHS)), maximum(diag(LHS))
        @info "MLHS diagonal range" minimum(diag(M_inv*LHS)), maximum(diag(M_inv*LHS))
        @info "norm(LHS, Inf)" norm(LHS, Inf)
        @info "norm(MLHS, Inf)" norm(M_inv*LHS, Inf)
        @info "norm(LHS, 1)" norm(LHS, 1)
        @info "norm(MLHS, 1)" norm(M_inv*LHS, 1)
        # ── Adaptivity criterion ──────────────────────────────────────────────
        @rankinfo rank "Computing adaptivity criterion..."
        one_vec = ones(Float64, size(LHS, 1))
        pointwise_interaction = abs.(LHS) * one_vec

        if nprocs > 1
            pointwise_interaction_g = zeros(gnpoin)
            for ip = 1:n_spa
                pointwise_interaction_g[ip2gip_spa[ip]] = pointwise_interaction[ip]
            end
            MPI.Allreduce!(pointwise_interaction_g, +, comm)
            for ip = 1:n_spa
                pointwise_interaction[ip] = pointwise_interaction_g[ip2gip_spa[ip]]
            end
        end

        criterion = compute_adaptivity_criterion3D_2D(
            pointwise_interaction, nelem, ngl, mesh.connijk,
            extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_coords, connijk_spa, extra_mesh[1].ψ, extra_mesh[1].dψ,
            extra_meshes_extra_dξdx, extra_meshes_extra_dηdx,
            extra_meshes_extra_dξdy, extra_meshes_extra_dηdy)
        #Reduction to find maximum of criterion
        crit_max = MPI.Allreduce(maximum(maximum.(criterion)), MPI.MAX, comm)
        thresholds = [crit_max * inputs[:RT_amr_threshold][1]]
        # ── Angular grid refinement ───────────────────────────────────────────
        @rankinfo rank "Adapting angular grid..."

        # Build per-element refinement records and compute which elements are
        # eligible for angular refinement.  An element that is already spatially
        # non-conforming with any neighbor (local or cross-rank) must not also
        # become angularly non-conforming with that same neighbor.
        ang_refine_records = build_element_refinement_records(
            mesh, extra_meshes_ref_level, nelem, ngl, comm)
        verify_element_refinement_records(
            ang_refine_records, nelem, comm, rank;
            uniform_spatial = !has_spatial_hanging_nodes)
        ang_refine_mask = build_angular_refinement_mask(ang_refine_records)
        @rankinfo rank "Angular refinement eligible: $(count(ang_refine_mask))/$nelem elements"

        #Force periodic angular elements on the same spatial element to conform after refinement
        max_conformity_iters = 2
        for _ = 1:max_conformity_iters
            n_forced_before = count(c -> c > thresholds[1],
                                    [criterion[iel][e] for iel=1:nelem
                                    for e=1:extra_meshes_extra_nelems[iel]])
            enforce_periodic_phi_conformity!(criterion, thresholds, extra_meshes_ref_level,
                                            nelem, extra_meshes_extra_nelems, extra_meshes_extra_nops,
                                            extra_meshes_connijk, extra_meshes_coords)
            n_forced_after = count(c -> c > thresholds[1],
                                [criterion[iel][e] for iel=1:nelem
                                    for e=1:extra_meshes_extra_nelems[iel]])
            n_forced_after == n_forced_before && break
        end

        adapt_angular_grid_3Dby2D!(
            criterion, thresholds, extra_meshes_ref_level, nelem, ngl,
            extra_meshes_extra_nelems, extra_meshes_extra_nops, neighbors,
            extra_meshes_extra_npoins, extra_meshes_connijk, extra_meshes_coords,
            extra_meshes_extra_Je,
            extra_meshes_extra_dξdx, extra_meshes_extra_dxdξ,
            extra_meshes_extra_dξdy, extra_meshes_extra_dydξ,
            extra_meshes_extra_dηdy, extra_meshes_extra_dydη,
            extra_meshes_extra_dηdx, extra_meshes_extra_dxdη,
            mesh.connijk, mesh.x, mesh.y, mesh.z,
            mesh.xmin, mesh.ymin, mesh.zmin,
            mesh.xmax, mesh.ymax, mesh.zmax,
            extra_mesh[1].ψ, extra_mesh[1].dψ, ang_refine_mask)

        nonowned_parent_indices = Set{Int}()
        nc_mat = sparse(I,n_spa,n_spa)
        nc_mat_rhs = sparse(I,n_spa,n_spa)
        P = nc_mat'
        P_vec = nc_mat_rhs'
        
        adapted = false   # set to true inside the if block when adaptation occurs
        # CRITICAL: Must use global maximum across ALL ranks to avoid MPI divergence
        # If one rank has adapted and another hasn't, they end up in different code blocks
        # with different MPI calls that can't be satisfied collectively
        local_max_ref = maximum(maximum.(extra_meshes_ref_level))
        global_max_ref = MPI.Allreduce(local_max_ref, MPI.MAX, MPI.COMM_WORLD)
        if !(global_max_ref == 0) || has_spatial_hanging_nodes
            if !(global_max_ref == 0)
                # ── Post-adaptivity connectivity ──────────────────────────────────
                connijk_spa = [Array{Int}(undef, ngl, ngl, ngl,
                                  extra_meshes_extra_nelems[iel],
                                  extra_meshes_extra_nops[iel][1]+1,
                                  extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]
    
                @rankinfo rank "Rebuilding connectivity on adapted mesh..."
                nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa =
                    adaptive_spatial_angular_numbering_3D_2D!(
                        connijk_spa, nelem, ngl, mesh.connijk,
                        extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                        extra_meshes_coords, mesh.x, mesh.y, mesh.z,
                        extra_meshes_ref_level, neighbors, adapted, extra_meshes_extra_Je)
    
                adapted = true
                @rankinfo rank "Hanging nodes: $n_non_global_nodes"
    
                # ── Post-adaptivity global numbering ──────────────────────────────
                @rankinfo rank "Building post-adaptivity global numbering..."
                ip2gip_spa, gip2ip, gip2owner_spa, gnpoin =
                    setup_global_numbering_adaptive_angular_scalable(
                        ip2gip_dedup, mesh.gip2owner, mesh, connijk_spa,
                        extra_meshes_coords, extra_meshes_connijk,
                        extra_meshes_extra_nops, extra_meshes_extra_nelems,
                        n_spa, n_non_global_nodes, nc_non_global_nodes;
                        remapped_ips = _remapped_ips)
    
                @rankinfo rank "Global DOF count (post-adapt): $gnpoin  (free: $(n_spa - n_non_global_nodes), hanging: $n_non_global_nodes)"
    
                # ── Owner map in local indexing ───────────────────────────────────
                gip2owner_extra = zeros(Int, n_spa)
                for ip = 1:n_spa
                    gip2owner_extra[ip] = gip2owner_spa[ip2gip_spa[ip]]
                end
    
                # ── Reverse GID lookup (used in many downstream functions) ────────
                gip_to_local = Dict{Int, Int}()
                sizehint!(gip_to_local, n_spa)
                for ip = 1:n_spa
                    gip_to_local[ip2gip_spa[ip]] = ip
                end
    
                # ── Ghost layer ───────────────────────────────────────────────────
                @rankinfo rank "Building ghost layer..."
                ghost_layer = build_nonconforming_ghost_layer(
                    mesh, connijk_spa, mesh.ip2gip, ip2gip_spa, gip2owner_spa,
                    extra_meshes_coords, extra_meshes_connijk,
                    extra_meshes_extra_nops, extra_meshes_extra_nelems,
                    extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
                    extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
                    extra_meshes_ref_level, n_spa, neighbors)
    
                @rankinfo rank "Building extended local numbering for ghost parents..."
                gid_to_extended_local, extended_local_to_gid, n_total =
                    build_extended_local_numbering(n_spa, ghost_layer, ip2gip_spa, rank)
    
                # ── Constraint (restriction/prolongation) matrices ────────────────
                @rankinfo rank "Building restriction and prolongation matrices..."
                nc_mat, P, nc_mat_rhs, P_vec,
                    ghost_constraint_data, ghost_constraint_data_rhs, all_hanging_nodes,
                    gid_to_extended_parents, extended_parents_to_gid,
                    extended_parents_x, extended_parents_y, extended_parents_z,
                    extended_parents_θ, extended_parents_ϕ, extended_parents_ip,
                    local_parent_indices, nonowned_parent_indices, nonowned_parent_gids =
                    build_restriction_matrices_local_and_ghost(
                        connijk_spa, nc_non_global_nodes, n_spa,
                        ghost_layer, extra_meshes_coords, extra_meshes_connijk,
                        extra_meshes_extra_nops, extra_meshes_extra_nelems,
                        extra_meshes_extra_Je, mesh, ngl, nelem, neighbors,
                        ip2gip_spa, gip2owner_extra, gid_to_extended_local,
                        extended_local_to_gid, rank)
    
                # Build reverse ghost map via AllGather so every rank learns about
                # hanging nodes on OTHER ranks that depend on parents it owns locally.
                # The local-only approach misses: rank B has hanging node H depending on
                # parent P owned by rank A; A's ghost_constraint_data never mentions H,
                # so A never sends its contribution → H stays at P * solution locally only.
                @rankinfo rank "Building reverse ghost constraint map..."
                _ext_gid_set_ang = Set(extended_parents_to_gid)
                _local_rev_ang = Float64[]
                for (ip_hanging, parent_constraints) in ghost_constraint_data
                    hanging_gid   = ip2gip_spa[ip_hanging]
                    owner_hanging = get(gip2owner_spa, hanging_gid, rank)
                    for (parent_gid, weight) in parent_constraints
                        parent_gid in _ext_gid_set_ang || continue
                        push!(_local_rev_ang, Float64(parent_gid), Float64(hanging_gid),
                              Float64(owner_hanging), weight)
                    end
                end
                _n_rev_loc_ang = Int32(length(_local_rev_ang) ÷ 4)
                _n_rev_all_ang = MPI.Allgather([_n_rev_loc_ang], comm)
                _all_rev_ang   = MPI.Allgatherv(_local_rev_ang, Int32.(_n_rev_all_ang .* 4), comm)
    
                reverse_ghost_map = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()
                for k = 1:length(_all_rev_ang) ÷ 4
                    s = 4*(k-1)
                    parent_gid    = Int(round(_all_rev_ang[s+1]))
                    hanging_gid   = Int(round(_all_rev_ang[s+2]))
                    owner_hanging = Int(round(_all_rev_ang[s+3]))
                    weight        = _all_rev_ang[s+4]
                    local_parent  = get(gip_to_local, parent_gid, 0)
                    local_parent == 0 && continue
                    get(gip2owner_spa, parent_gid, -1) != rank && continue
                    push!(get!(reverse_ghost_map, local_parent, Tuple{Int,Int,Float64}[]),
                          (hanging_gid, owner_hanging, weight))
                end
    
                # ── LHS and mass matrix on adapted mesh ───────────────────────────
                @rankinfo rank "Assembling LHS on adapted mesh..."
                LHS = sparse_lhs_assembly_3Dby2D_adaptive(
                    ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                    mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ,
                    extra_meshes_connijk, extra_meshes_extra_Je,
                    extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl,
                    extra_meshes_extra_nelems, dξdx, dξdy, dξdz,
                    dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                    extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, inputs[:lRT_from_data], κ, σ, inputs[:RT_shortwave])
    
                P    = nc_mat'
                rest = nc_mat
                
                Md = assemble_mass_diagonal_3Dby2D_adaptive(
                    ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                    extra_meshes_extra_Je, extra_meshes_extra_nops, n_spa, nelem,
                    ngl, extra_meshes_extra_nelems, connijk_spa)
                M_inv = spdiagm(0 => 1.0 ./ Md)
                @info size(Md), size(LHS)
                @info "Md range" minimum(Md), maximum(Md)
                @info "Md mean" sum(Md) / length(Md)
                @info "Large Md entries" count(Md .> 1e4)
                @info "LHS diagonal range" minimum(diag(LHS)), maximum(diag(LHS))
                @info "MLHS diagonal range" minimum(diag(M_inv*LHS)), maximum(diag(M_inv*LHS))
                @info "norm(LHS, Inf)" norm(LHS, Inf)
                @info "norm(MLHS, Inf)" norm(M_inv*LHS, Inf)
                @info "norm(LHS, 1)" norm(LHS, 1)
                @info "norm(MLHS, 1)" norm(M_inv*LHS, 1)
            else  # global_max_ref == 0, has_spatial_hanging_nodes: spatial-only combined path
                # connijk_spa, n_spa, ip2gip_spa, gip2owner_spa, gip2owner_extra, gip_to_local,
                # nc_mat (identity), nc_non_global_nodes, n_non_global_nodes, LHS, Md
                # are all already correct from pre-adaptation (lines ~270-429).
                # Only the ghost/constraint structures are missing — initialize empty.
                ghost_constraint_data     = Dict{Int, Vector{Tuple{Int,Float64}}}()
                ghost_constraint_data_rhs = Dict{Int, Vector{Tuple{Int,Float64}}}()
                all_hanging_nodes         = Set{Int}()
                gid_to_extended_parents   = Dict{Int,Int}()
                extended_parents_to_gid   = Int[]
                extended_parents_x        = Float64[]
                extended_parents_y        = Float64[]
                extended_parents_z        = Float64[]
                extended_parents_θ        = Float64[]
                extended_parents_ϕ        = Float64[]
                extended_parents_ip       = Int[]
                local_parent_indices      = Set{Int}()
                nonowned_parent_gids      = Set{Int}()
            end  # if !(global_max_ref == 0) / else (spatial-only)

            # ── Stage 2-3: Spatial Constraints with Adaptive Angular Mesh ─────────────
            if has_spatial_hanging_nodes
                @info "[$rank] Building spatial constraints with adaptive angular mesh (Stage 2+)..."
                spatial_amr_cache = build_spatial_constraint_matrices(
                    mesh, spatial_amr_cache,
                    extra_meshes_coords, extra_meshes_connijk,
                    extra_meshes_extra_nops, extra_meshes_extra_nelems,
                    ngl, rank
                )
                @info "[$rank] ✓ Spatial constraint matrices built with adaptive angular"
                @info "[$rank]   - $(length(spatial_amr_cache.parent_weights)) hanging nodes"
    
                verify_spatial_constraints(spatial_amr_cache, rank, npoin, extra_meshes_extra_nelems, extra_meshes_extra_nops)
                @info "[$rank] ✓ Verification passed"
    
                spatial_amr_cache = exchange_spatial_ghosts(
                    mesh, spatial_amr_cache,
                    extra_meshes_extra_nops, extra_meshes_extra_nelems,
                    rank, comm
                )
                @info "[$rank] ✓ Ghost exchange complete (Stage 3)"
                verify_spatial_ghost_exchange(spatial_amr_cache, rank)
            end
    
            # ── Set up combined spatial+angular restriction/prolongation ─────────────
            # Both R operators are set up HERE before application:
            #   - Angular R: already in nc_mat (built above by build_restriction_matrices_local_and_ghost)
            #   - Spatial R: augmented into nc_mat by combine_spatial_angular_restrictions
            #
            # Mathematical basis (disjoint hanging sets by can_refine_angular mask):
            #   nc_mat_combined = nc_mat + nc_mat_spatial_local - I_{spatial_hanging}
            #
            # The off-diagonals sum without interference; the subtracted identity
            # removes the identity rows nc_mat has for nodes free-in-angular but
            # spatial-hanging. One combined R is then applied with the same parallel
            # restriction/prolongation machinery (compute_hanging_row_effects_before_restriction,
            # exchange_hanging_effects, etc.).
            #
            # Guard: requires adapted=true (connijk_spa, n_spa available).
            # If no spatial hanging nodes, nc_mat is unchanged and all_hanging_nodes
            # stays as the pure angular set.
            spatial_hanging_combined = Set{Int}()
    
            if has_spatial_hanging_nodes
                # Build (x,y,z,θ,ϕ) → combined DOF index for the adapted mesh.
                point_dict_combined_adapted = Dict{NTuple{5,Float64}, Int}()
                for iel_d = 1:nelem
                    for kk_d = 1:ngl, jj_d = 1:ngl, ii_d = 1:ngl
                        ip_s_d = mesh.connijk[iel_d, ii_d, jj_d, kk_d]
                        x_d = mesh.x[ip_s_d]; y_d = mesh.y[ip_s_d]; z_d = mesh.z[ip_s_d]
                        for e_d = 1:extra_meshes_extra_nelems[iel_d]
                            nop_d = extra_meshes_extra_nops[iel_d][e_d]
                            for jθ_d = 1:nop_d+1, iθ_d = 1:nop_d+1
                                i_ang_d = extra_meshes_connijk[iel_d][e_d, iθ_d, jθ_d]
                                θ_d = extra_meshes_coords[iel_d][1, i_ang_d]
                                ϕ_d = extra_meshes_coords[iel_d][2, i_ang_d]
                                key_d = (round(x_d, digits=12), round(y_d, digits=12),
                                         round(z_d, digits=12), round(θ_d, digits=12),
                                         round(ϕ_d, digits=12))
                                get!(point_dict_combined_adapted, key_d,
                                     connijk_spa[iel_d][ii_d, jj_d, kk_d, e_d, iθ_d, jθ_d])
                            end
                        end
                    end
                end
                if global_max_ref == 0
                    # Spatial-only case: no angular adaptation, nc_mat is identity.
                    # Use the proven uniform-path approach directly (mirrors lines 952-1429)
                    # instead of combine_spatial_angular_restrictions, which has parallel bugs.
                    (nc_mat, P, nc_mat_rhs,
                     ghost_constraint_data, ghost_constraint_data_rhs,
                     gid_to_extended_parents, extended_parents_to_gid,
                     gip_to_local, all_hanging_nodes, _) =
                        build_spatial_constraints_for_combined_path(
                            spatial_amr_cache, mesh,
                            ip2gip_spa, gip2owner_extra, gip2owner_spa,
                            point_dict_combined_adapted, n_spa, rank, comm
                        )
                    spatial_hanging_combined = all_hanging_nodes
                    P_vec = sparse(nc_mat_rhs')
                else
                    # Combined angular+spatial: build spatial constraints with the proven
                    # uniform-path approach, then merge them with the angular nc_mat.
                    (R_spatial_ext_comb, _, R_spatial_rhs_ext_comb,
                     ghost_cdata_spa_c, ghost_cdata_spa_rhs_c,
                     _, ext_to_gid_spa_c,
                     _, spatial_hanging_combined_raw, _) =
                        build_spatial_constraints_for_combined_path(
                            spatial_amr_cache, mesh,
                            ip2gip_spa, gip2owner_extra, gip2owner_spa,
                            point_dict_combined_adapted, n_spa, rank, comm
                        )
                    nc_mat, P, nc_mat_rhs, P_vec, all_hanging_nodes, spatial_hanging_combined =
                        combine_spatial_angular_restrictions(
                            nc_mat, nc_mat_rhs,
                            ghost_constraint_data, ghost_constraint_data_rhs,
                            gid_to_extended_parents, extended_parents_to_gid,
                            extended_parents_x, extended_parents_y, extended_parents_z,
                            extended_parents_θ, extended_parents_ϕ, extended_parents_ip,
                            all_hanging_nodes,
                            R_spatial_ext_comb, R_spatial_rhs_ext_comb,
                            ghost_cdata_spa_c, ghost_cdata_spa_rhs_c,
                            ext_to_gid_spa_c,
                            spatial_hanging_combined_raw,
                            n_spa, rank
                        )
                end
            end

            # ── Rebuild reverse ghost map after spatial constraints are finalized ──
            # For the spatial-only combined path and the combined angular+spatial path,
            # ghost_constraint_data and gid_to_extended_parents were set/augmented inside
            # the has_spatial_hanging_nodes block. The reverse_ghost_map built earlier
            # (lines 547-575) used only the angular extended parents and is now stale.
            # Rebuild it here from the finalized ghost_constraint_data so that solution
            # prolongation at line ~2334 correctly exchanges cross-rank parent contributions.
            if has_spatial_hanging_nodes
                _ext_gid_set_all = Set(extended_parents_to_gid)
                _local_rev_all = Float64[]
                for (ip_hanging, parent_constraints) in ghost_constraint_data
                    hanging_gid   = ip2gip_spa[ip_hanging]
                    owner_hanging = get(gip2owner_spa, hanging_gid, rank)
                    for (parent_gid, weight) in parent_constraints
                        parent_gid in _ext_gid_set_all || continue
                        push!(_local_rev_all, Float64(parent_gid), Float64(hanging_gid),
                              Float64(owner_hanging), weight)
                    end
                end
                _n_rev_loc_all = Int32(length(_local_rev_all) ÷ 4)
                _n_rev_all_all = MPI.Allgather([_n_rev_loc_all], comm)
                _all_rev_all   = MPI.Allgatherv(_local_rev_all, Int32.(_n_rev_all_all .* 4), comm)

                reverse_ghost_map = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()
                for k = 1:length(_all_rev_all) ÷ 4
                    s = 4*(k-1)
                    parent_gid    = Int(round(_all_rev_all[s+1]))
                    hanging_gid   = Int(round(_all_rev_all[s+2]))
                    owner_hanging = Int(round(_all_rev_all[s+3]))
                    weight        = _all_rev_all[s+4]
                    local_parent  = get(gip_to_local, parent_gid, 0)
                    local_parent == 0 && continue
                    get(gip2owner_spa, parent_gid, -1) != rank && continue
                    push!(get!(reverse_ghost_map, local_parent, Tuple{Int,Int,Float64}[]),
                          (hanging_gid, owner_hanging, weight))
                end
                @info "[$rank] reverse_ghost_map rebuilt with $(length(reverse_ghost_map)) parent entries (spatial+angular extended)"
            end

            # ── Mass matrix assembly and inversion ────────────────────────────────
            pM = setup_assembler(SD, Md, ip2gip_spa, gip2owner_extra)
            if pM !== nothing
                assemble_mpi!(Md, pM)
            end
            M_inv = spdiagm(0 => 1.0 ./ Md)
            MLHS  = sparse(M_inv * LHS)
            # ── All-reduce parent–parent entries across ranks ─────────────────────
            # Hanging-node constraint handling requires that (parent, parent) blocks
            # of M⁻¹LHS are globally consistent before restriction/prolongation.
            if nprocs > 1
                @rankinfo rank "All-reducing parent-parent matrix entries..."
                #=MLHS = allreduce_parent_parent_entries(
                    MLHS, nc_mat, gip2owner_extra, n_spa, rank, comm,
                    ip2gip_spa, local_parent_indices, nonowned_parent_indices,
                    nonowned_parent_gids, gip_to_local)=#
            end
    
            # ── Parallel restriction: apply R from the left ───────────────────────
            # Interface hanging nodes on other ranks contribute row effects that must
            # be communicated before the local nc_mat application.
            @rankinfo rank "Applying parallel restriction (left multiply by R)..."
            row_effects_to_send, MLHS_effects =
                compute_hanging_row_effects_before_restriction(
                    ghost_constraint_data, MLHS, ip2gip_spa, gip2owner_spa, rank,
                    gid_to_extended_parents, extended_parents_to_gid, gip_to_local)
    
            A_left_restricted = nc_mat * MLHS_effects
    
            received_row_effects =
                exchange_hanging_effects(row_effects_to_send, rank, comm)
    
            A_with_row_effects = add_hanging_row_effects(
                A_left_restricted, received_row_effects, ip2gip_spa,
                size(MLHS_effects, 1), rank, gip_to_local)
    
            # ── Parallel prolongation: apply P from the right ─────────────────────
            @rankinfo rank "Applying parallel prolongation (right multiply by P)..."
            col_effects_to_send, A_ghost_effects =
                compute_hanging_col_effects_before_prolongation(
                    ghost_constraint_data, A_with_row_effects, ip2gip_spa, gip2owner_spa,
                    n_spa, size(MLHS_effects, 1), all_hanging_nodes, rank,
                    gid_to_extended_parents, extended_parents_to_gid, gip_to_local)
    
            A_both_restricted = A_ghost_effects * P
    
            received_col_effects =
                exchange_hanging_effects(col_effects_to_send, rank, comm)
    
            n_spa_g = size(MLHS_effects, 1)
            A_with_col_effects = add_hanging_col_effects(
                A_both_restricted, received_col_effects, ip2gip_spa,
                n_spa_g, rank, gip_to_local)
    
            # ── Extract free-node submatrix ───────────────────────────────────────
            n_free = n_spa - length(all_hanging_nodes)
            @rankinfo rank "Extracting free-node submatrix (removing hanging rows/cols)..."
            A_free = extract_free_submatrix_remove_all_hanging(
                A_with_col_effects, all_hanging_nodes, n_free, n_spa_g, rank)
    
            A = sparse(A_free)
    
            # Remove non-owned parent–parent entries to avoid double-counting.
            # These were all-reduced above and must be zeroed on non-owning ranks.
            for i in nonowned_parent_indices
                for j in nonowned_parent_indices
                    if abs(A[i,j]) > 0.0
                        #A[i,j] = 0.0
                    end
                end
            end
            A = sparse(A)
    
            RHS = zeros(TFloat, n_spa_g)
            ref = zeros(TFloat, n_spa)
            BDY = zeros(TFloat, n_spa_g)
    
        else  # if !(global_max_ref == 0) || has_spatial_hanging_nodes
            # Conforming adaptive case: no refinement occurred, so Md still holds the
            # locally-assembled mass diagonal from line ~374.  In parallel the diagonal
            # is incomplete (missing neighbour-rank contributions at shared DOFs), so
            # assemble it globally before inverting.
            pM_conform = setup_assembler(SD, Md, ip2gip_spa, gip2owner_extra)
            if pM_conform !== nothing
                assemble_mpi!(Md, pM_conform)
            end
            M_inv = spdiagm(0 => 1.0 ./ Md)
            A = sparse(M_inv*LHS)
            n_spa_g = n_spa
            RHS = zeros(TFloat, n_spa)
            ref = zeros(TFloat, n_spa)
            BDY = zeros(TFloat, n_spa)
        end
    else  # ── Non-adaptive path ────────────────────────────────────────────────

        # ── Stage 2-3: Spatial Constraints with Uniform Angular Mesh ─────────────
        # Defaults used when there are no spatial hanging nodes.
        n_spa_new = 0
        n_non_global_nodes_spa = 0
        connijk_spa_uniform = Vector{Array{Int,6}}()
        point_dict_spa = Dict{NTuple{5,Float64}, Int}()
        global_sp_ang_to_gid = Dict{Tuple{Int,Int}, Int}()
        _all_needed_set = Set{Int}()

        if has_spatial_hanging_nodes
            @info "[$rank] Building spatial constraints with uniform angular mesh (Stage 2+)..."

            # Create arrays mimicking adaptive mesh structure for uniform case
            extra_meshes_coords_uniform = [extra_mesh.extra_coords[:,:] for _ in 1:nelem]
            extra_meshes_connijk_uniform = [extra_mesh.extra_connijk for _ in 1:nelem]
            extra_meshes_extra_nops_uniform = [[extra_mesh.extra_nop[1] for _ in 1:extra_mesh.extra_nelem] for _ in 1:nelem]
            extra_meshes_extra_nelems_uniform = [extra_mesh.extra_nelem for _ in 1:nelem]

            # Initialize here (outside try) so they are in scope for the assembly and
            # Stage 5 blocks below. The try block updates ip2gip_dedup for coincident
            # nodes and populates nc_non_global_nodes_spa / gip_spa_to_local_ip / n_ang.
            ip2gip_dedup           = Int.(mesh.ip2gip)
            nc_non_global_nodes_spa = Int[]
            gip_spa_to_local_ip    = Dict{Int,Int}()
            n_ang                  = extra_mesh.extra_npoin
            _remapped_ips          = Set{Int}()

            try
                spatial_amr_cache = build_spatial_constraint_matrices(
                    mesh, spatial_amr_cache,
                    extra_meshes_coords_uniform, extra_meshes_connijk_uniform,
                    extra_meshes_extra_nops_uniform, extra_meshes_extra_nelems_uniform,
                    ngl, rank
                )
                @info "[$rank] ✓ Spatial constraint matrices built with uniform angular"
                @info "[$rank]   - $(length(spatial_amr_cache.parent_weights)) hanging nodes"

                verify_spatial_constraints(spatial_amr_cache, rank, npoin, extra_meshes_extra_nelems_uniform, extra_meshes_extra_nops_uniform)
                @info "[$rank] ✓ Verification passed"

                spatial_amr_cache = exchange_spatial_ghosts(
                    mesh, spatial_amr_cache,
                    extra_meshes_extra_nops_uniform, extra_meshes_extra_nelems_uniform,
                    rank, comm
                )
                @info "[$rank] ✓ Ghost exchange complete (Stage 3)"
                verify_spatial_ghost_exchange(spatial_amr_cache, rank)

                # ── Update ip2gip_dedup BEFORE mutual-exclusivity filtering ──────────
                # Coincident nodes (child at same (x,y,z) as parent corner) must share
                # the parent's spatial GIP so setup_global_numbering assigns them the same
                # compact GID on all ranks.  This must happen for ALL coincident nodes —
                # even those that will be removed below for also acting as parents.
                # Collect remapped_ips: passed to setup_global_numbering so Phase 2 treats
                # remapped IPs as nonowned (preventing both ranks from claiming ownership of
                # the same GIP, which would break sharing detection and produce duplicate GIDs).
                _remapped_ips = Set{Int}()
                for (child_ip, parent_ip) in spatial_amr_cache.coincident_nodes
                    ip2gip_dedup[child_ip] = mesh.ip2gip[parent_ip]
                    push!(_remapped_ips, child_ip)
                end
                for (child_ip, parent_global_spa) in spatial_amr_cache.cross_rank_coincident_nodes
                    ip2gip_dedup[child_ip] = parent_global_spa
                    push!(_remapped_ips, child_ip)
                end

                # ── Mutual-exclusivity post-processing ────────────────────────
                # Corner nodes at junctions of multiple NCFs can appear simultaneously
                # as parents of interpolated hanging nodes AND as coincident children of
                # another NCF. Such nodes must NOT be classified as coincident for
                # constraint purposes (coincident treatment zeroes their rows), but their
                # ip2gip_dedup remapping above is kept so global numbering remains unique.
                local_parent_ip_set = Set{Int}()
                for (_, plist) in spatial_amr_cache.parent_weights
                    for (pip, _) in plist
                        push!(local_parent_ip_set, pip)
                    end
                end
                cross_rank_parent_gip_set = Set{Int}()
                for (_, plist) in spatial_amr_cache.cross_rank_parent_weights
                    for (pgip, _) in plist
                        push!(cross_rank_parent_gip_set, pgip)
                    end
                end
                n_removed_local      = 0
                n_removed_cross_rank = 0
                for child_ip in collect(keys(spatial_amr_cache.coincident_nodes))
                    child_gip = Int(mesh.ip2gip[child_ip])
                    # Remove if acting as a parent, OR if also classified as interpolated child
                    # (parent_weights keys take precedence over coincident_nodes)
                    if child_ip in local_parent_ip_set ||
                       child_gip in cross_rank_parent_gip_set ||
                       haskey(spatial_amr_cache.parent_weights, child_ip) ||
                       haskey(spatial_amr_cache.cross_rank_parent_weights, child_ip)
                        delete!(spatial_amr_cache.coincident_nodes, child_ip)
                        n_removed_local += 1
                    end
                end
                for child_ip in collect(keys(spatial_amr_cache.cross_rank_coincident_nodes))
                    child_gip = Int(mesh.ip2gip[child_ip])
                    if child_ip in local_parent_ip_set ||
                       child_gip in cross_rank_parent_gip_set ||
                       haskey(spatial_amr_cache.parent_weights, child_ip) ||
                       haskey(spatial_amr_cache.cross_rank_parent_weights, child_ip)
                        delete!(spatial_amr_cache.cross_rank_coincident_nodes, child_ip)
                        n_removed_cross_rank += 1
                    end
                end
                @info "[$rank] Mutual-exclusivity: removed $n_removed_local coincident + $n_removed_cross_rank cross-rank coincident nodes that were also acting as parents or interpolated children"

                # ── Populate gip_spa_to_local_ip: reverse map global spatial ID → local spatial ip ──
                # Needed to convert cross-rank parent global spatial IDs to local DOF indices.
                # (initialized to empty Dict before this try block)
                for ip = 1:mesh.npoin
                    gip_spa_to_local_ip[Int(mesh.ip2gip[ip])] = ip
                end

                # ── Build connijk_spa_uniform with coordinate deduplication ──
                # Key: (x,y,z,θ,ϕ) → unique local DOF index ip_spa.
                # Coincident spatial nodes (same (x,y,z) but different mesh.connijk ip)
                # naturally map to the same ip_spa — no explicit coincident handling needed.
                # ip2gip_dedup ensures coincident nodes share the same gip signature in
                # setup_global_numbering_adaptive_angular_scalable → same compact GID.
                connijk_spa_uniform = [
                    Array{Int}(undef, ngl, ngl, ngl, extra_mesh.extra_nelem,
                               extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
                    for _ in 1:nelem]
                point_dict_spa = Dict{NTuple{5,Float64}, Int}()
                n_spa_counter = 0
                for iel_c = 1:nelem
                    for kc = 1:ngl, jc = 1:ngl, ic = 1:ngl
                        ip_c = mesh.connijk[iel_c, ic, jc, kc]
                        x_c = mesh.x[ip_c]; y_c = mesh.y[ip_c]; z_c = mesh.z[ip_c]
                        for e_ext_c = 1:extra_mesh.extra_nelem
                            nop_c = extra_mesh.extra_nop[1]
                            for jθ_c = 1:nop_c+1, iθ_c = 1:nop_c+1
                                i_ang_c = extra_mesh.extra_connijk[e_ext_c, iθ_c, jθ_c]
                                θ_c = extra_mesh.extra_coords[1, i_ang_c]
                                ϕ_c = extra_mesh.extra_coords[2, i_ang_c]
                                key = (round(x_c, digits=12), round(y_c, digits=12),
                                       round(z_c, digits=12), round(θ_c, digits=12),
                                       round(ϕ_c, digits=12))
                                idx = get(point_dict_spa, key, 0)
                                if idx == 0
                                    n_spa_counter += 1
                                    point_dict_spa[key] = n_spa_counter
                                    connijk_spa_uniform[iel_c][ic, jc, kc, e_ext_c, iθ_c, jθ_c] = n_spa_counter
                                else
                                    connijk_spa_uniform[iel_c][ic, jc, kc, e_ext_c, iθ_c, jθ_c] = idx
                                end
                            end
                        end
                    end
                end
                n_spa_new = n_spa_counter
                @info "[$rank] connijk_spa_uniform built: $n_spa_new unique DOFs (vs $(npoin * n_ang) without dedup)"

                # ── Build nc_non_global_nodes_spa: ip_spa indices for interpolated hanging nodes ──
                # Look up each hanging spatial node's (x,y,z,θ,ϕ) in point_dict_spa to get ip_spa.
                # Coincident nodes are already transparent (same ip_spa as parent via dedup).
                nc_non_global_nodes_spa = Int[]
                nc_set_spa = Set{Int}()
                for hanging_ip in keys(spatial_amr_cache.parent_weights)
                    x_h = mesh.x[hanging_ip]; y_h = mesh.y[hanging_ip]; z_h = mesh.z[hanging_ip]
                    for e_ext_h = 1:extra_mesh.extra_nelem
                        nop_h = extra_mesh.extra_nop[1]
                        for jθ_h = 1:nop_h+1, iθ_h = 1:nop_h+1
                            i_ang_h = extra_mesh.extra_connijk[e_ext_h, iθ_h, jθ_h]
                            θ_h = extra_mesh.extra_coords[1, i_ang_h]
                            ϕ_h = extra_mesh.extra_coords[2, i_ang_h]
                            key_h = (round(x_h, digits=12), round(y_h, digits=12),
                                     round(z_h, digits=12), round(θ_h, digits=12),
                                     round(ϕ_h, digits=12))
                            ip_spa_h = get(point_dict_spa, key_h, 0)
                            ip_spa_h == 0 && continue
                            if !(ip_spa_h in nc_set_spa)
                                push!(nc_non_global_nodes_spa, ip_spa_h)
                                push!(nc_set_spa, ip_spa_h)
                            end
                        end
                    end
                end
                for hanging_ip in keys(spatial_amr_cache.cross_rank_parent_weights)
                    x_h = mesh.x[hanging_ip]; y_h = mesh.y[hanging_ip]; z_h = mesh.z[hanging_ip]
                    for e_ext_h = 1:extra_mesh.extra_nelem
                        nop_h = extra_mesh.extra_nop[1]
                        for jθ_h = 1:nop_h+1, iθ_h = 1:nop_h+1
                            i_ang_h = extra_mesh.extra_connijk[e_ext_h, iθ_h, jθ_h]
                            θ_h = extra_mesh.extra_coords[1, i_ang_h]
                            ϕ_h = extra_mesh.extra_coords[2, i_ang_h]
                            key_h = (round(x_h, digits=12), round(y_h, digits=12),
                                     round(z_h, digits=12), round(θ_h, digits=12),
                                     round(ϕ_h, digits=12))
                            ip_spa_h = get(point_dict_spa, key_h, 0)
                            ip_spa_h == 0 && continue
                            if !(ip_spa_h in nc_set_spa)
                                push!(nc_non_global_nodes_spa, ip_spa_h)
                                push!(nc_set_spa, ip_spa_h)
                            end
                        end
                    end
                end
                n_non_global_nodes_spa = length(nc_non_global_nodes_spa)

                # Hanging DOF set for BC exclusion — coincident nodes excluded:
                # they are now transparent (same ip_spa as parent via coordinate dedup).
                spatial_hanging_nodes_all_angular = Set(nc_non_global_nodes_spa)
                @info "[$rank] Tracked $(length(spatial_hanging_nodes_all_angular)) spatial-angular hanging DOFs ($(length(spatial_amr_cache.parent_weights)) local + $(length(spatial_amr_cache.cross_rank_parent_weights)) cross-rank); coincident nodes handled via coordinate dedup"

                # ── Stage 4 prep: Build spatial constraint matrices for RHS/LHS assembly ────
                # Build xyz_ang_map from point_dict_spa: (x,y,z) → [(θ,ϕ,dof_idx),...]
                xyz_ang_map_spa = Dict{NTuple{3,Float64}, Vector{Tuple{Float64,Float64,Int}}}()
                for ((x,y,z,θ,ϕ), idx) in point_dict_spa
                    push!(get!(xyz_ang_map_spa, (x,y,z), Tuple{Float64,Float64,Int}[]), (θ, ϕ, idx))
                end
                @info "[$rank] Building spatial restriction and prolongation matrices (Stage 4 prep)..."
                R_spatial, P_spatial = build_spatial_restriction_and_prolongation(
                    spatial_amr_cache, n_spa_new, spatial_hanging_nodes_all_angular,
                    point_dict_spa, mesh, xyz_ang_map_spa
                )
                @info "[$rank] ✓ Spatial matrices built: R_spatial = $(size(R_spatial)), P_spatial = $(size(P_spatial))"

            catch err
                rethrow(err)
            end
        end

        npoin_ang_total = npoin * extra_mesh.extra_npoin
        @info "total number of points DOFs", npoin_ang_total, npoin, extra_mesh.extra_npoin

        if has_spatial_hanging_nodes
            # Use coordinate-dedup assembly so coincident spatial nodes accumulate
            # into the same DOF slot (n_spa_new ≤ npoin_ang_total).
            n_dofs = n_spa_new
            @rankinfo rank "Assembling LHS ($n_spa_new deduplicated DOF, spatial AMR)..."
            LHS = sparse_lhs_assembly_3Dby2D_spatial_amr(
                ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
                mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh.ψ,
                extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
                extra_mesh.extra_coords, extra_mesh.extra_nop,
                n_spa_new, nelem, ngl, extra_mesh.extra_nelem,
                dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                extra_mesh.extra_npoin, inputs[:rad_HG_g],
                connijk_spa_uniform, inputs[:lRT_from_data], κ, σ, inputs[:RT_shortwave])

            Md = assemble_mass_diagonal_3Dby2D_spatial_amr(
                ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
                extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
                extra_mesh.extra_nop, n_spa_new, nelem, ngl,
                extra_mesh.extra_nelem, connijk_spa_uniform)

            # Compact global numbering: coincident nodes share the parent's GID (via
            # ip2gip_dedup), interpolated hanging nodes appear after free nodes in GID space.
            _ip2gip_spa_full, _gip2ip_spa, _gip2owner_spa_full, gnpoin =
                setup_global_numbering_adaptive_angular_scalable(
                    ip2gip_dedup, mesh.gip2owner, mesh, connijk_spa_uniform,
                    extra_meshes_coords_uniform, extra_meshes_connijk_uniform,
                    extra_meshes_extra_nops_uniform, extra_meshes_extra_nelems_uniform,
                    n_spa_new, n_non_global_nodes_spa, nc_non_global_nodes_spa;
                    remapped_ips = _remapped_ips)
            ip2gip_spa = _ip2gip_spa_full
            gip2owner_extra = zeros(Int, n_spa_new)
            for ip = 1:n_spa_new
                gip2owner_extra[ip] = _gip2owner_spa_full[ip2gip_spa[ip]]
            end

            # ── Build global (spatial_GIP, angular_node_index) → compact_GID map ───────
            # Cross-rank NCF parents are ghost nodes — not in this rank's local elements,
            # so not in point_dict_spa.  Targeted two-round AllGather:
            #   Round 1: gather all needed parent spatial GIPs (from cross_rank_parent_weights)
            #   Round 2: each owning rank provides (sp_gip, ang_ip, compact_gid) triples
            _needed_gips = Set{Int}()
            for (_, _cw) in spatial_amr_cache.cross_rank_parent_weights
                for (_pgip, _) in _cw
                    push!(_needed_gips, _pgip)
                end
            end
            _local_needed = collect(_needed_gips)
            _nn = Int32(length(_local_needed))
            _nn_all = MPI.Allgather([_nn], comm)
            _all_needed = MPI.Allgatherv(_local_needed, _nn_all, comm)
            _all_needed_set = Set(_all_needed)

            # Each rank provides triples for needed GIPs it has in its local elements
            _local_triples = Int[]
            _seen_spa_ang = Set{Tuple{Int,Int}}()
            for _iel = 1:nelem
                for _kk = 1:ngl, _jj = 1:ngl, _ii = 1:ngl
                    _ip_s = mesh.connijk[_iel, _ii, _jj, _kk]
                    _sp_gip = ip2gip_dedup[_ip_s]
                    _sp_gip in _all_needed_set || continue
                    for _e = 1:extra_mesh.extra_nelem
                        _nop = extra_mesh.extra_nop[_e]
                        for _jt = 1:_nop+1, _it = 1:_nop+1
                            _i_ang = extra_mesh.extra_connijk[_e, _it, _jt]
                            _key_t = (_sp_gip, _i_ang)
                            _key_t in _seen_spa_ang && continue
                            push!(_seen_spa_ang, _key_t)
                            _ip_spa = connijk_spa_uniform[_iel][_ii, _jj, _kk, _e, _it, _jt]
                            push!(_local_triples, _sp_gip, _i_ang, ip2gip_spa[_ip_spa])
                        end
                    end
                end
            end
            _nt = Int32(length(_local_triples) ÷ 3)
            _nt_all = MPI.Allgather([_nt], comm)
            _all_triples = MPI.Allgatherv(_local_triples, Int32.(_nt_all .* 3), comm)
            global_sp_ang_to_gid = Dict{Tuple{Int,Int}, Int}()
            sizehint!(global_sp_ang_to_gid, length(_all_triples) ÷ 3)
            for _k = 1:length(_all_triples) ÷ 3
                _sg = _all_triples[3_k-2]; _ai = _all_triples[3_k-1]; _cg = _all_triples[3_k]
                get!(global_sp_ang_to_gid, (_sg, _ai), _cg)
            end
            @info "[$rank] cross-rank parent GID map: $(length(global_sp_ang_to_gid)) entries for $(length(_all_needed_set)) needed spatial GIPs"

            # AllGather (x,y,z,θ,ϕ) for every cross-rank parent compact GID.
            # Each rank contributes coords for GIPs it provided in _local_triples.
            # Used to resolve extended-parent row coordinates in the debug comparison.
            _local_cdata = Float64[]
            _seen_cg_c   = Set{Int}()
            for _k = 1:length(_local_triples) ÷ 3
                _sg = _local_triples[3_k-2]; _ai = _local_triples[3_k-1]; _cg = _local_triples[3_k]
                _cg in _seen_cg_c && continue
                _lip = get(gip_spa_to_local_ip, _sg, 0)
                _lip == 0 && continue
                push!(_seen_cg_c, _cg)
                push!(_local_cdata, Float64(_cg),
                      mesh.x[_lip], mesh.y[_lip], mesh.z[_lip],
                      extra_mesh.extra_coords[1, _ai], extra_mesh.extra_coords[2, _ai])
            end
            _nc_local  = Int32(length(_local_cdata) ÷ 6)
            _nc_all_c  = MPI.Allgather([_nc_local], comm)
            _all_cdata = MPI.Allgatherv(_local_cdata, Int32.(_nc_all_c .* 6), comm)
            global_gid_to_coords = Dict{Int, NTuple{5,Float64}}()
            sizehint!(global_gid_to_coords, length(_all_cdata) ÷ 6)
            for _k = 1:length(_all_cdata) ÷ 6
                _s  = 6*(_k-1)
                _cg = Int(round(_all_cdata[_s+1]))
                haskey(global_gid_to_coords, _cg) && continue
                global_gid_to_coords[_cg] = (_all_cdata[_s+2], _all_cdata[_s+3],
                                              _all_cdata[_s+4], _all_cdata[_s+5],
                                              _all_cdata[_s+6])
            end
            _local_cdata = nothing; _all_cdata = nothing; _seen_cg_c = nothing
        else
            n_dofs = npoin_ang_total
            @rankinfo rank "Assembling LHS ($npoin_ang_total DOF)..."
            LHS = sparse_lhs_assembly_3Dby2D(
                ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
                mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh.ψ,
                extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
                extra_mesh.extra_coords, extra_mesh.extra_nop,
                npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                extra_mesh.extra_npoin, inputs[:rad_HG_g],
                inputs[:lRT_from_data], κ, σ, inputs[:RT_shortwave])

            Md = assemble_mass_diagonal_3Dby2D(
                ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
                extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
                extra_mesh.extra_nop, npoin_ang_total, nelem, ngl,
                extra_mesh.extra_nelem, extra_mesh.extra_npoin)

            ip2gip_spa, gip2owner_extra, gnpoin =
                setup_global_numbering_extra_dim(
                    mesh.ip2gip, mesh.gip2owner, npoin,
                    extra_mesh.extra_npoin, npoin_ang_total)
        end
        n_spa = n_dofs

        pM = setup_assembler(SD, Md, ip2gip_spa, gip2owner_extra)
        if pM !== nothing
            assemble_mpi!(Md, pM)
        end

        M_inv = spdiagm(0 => 1.0 ./ Md)
        MLHS     = sparse(M_inv * LHS)
        A = copy(MLHS)
        @info "[$rank] A after assembly: size=$(size(A)), nnz=$(nnz(A))"
        M_inv = nothing; LHS = nothing
        GC.gc()

        BDY = zeros(TFloat, n_dofs)
        RHS = zeros(TFloat, n_dofs)
        ref = zeros(TFloat, n_dofs)

        # Extended parent structures for spatial AMR (populated inside the constraint block)
        extended_parents_to_gid_spa = Int[]
        gid_to_extended_parents_spa = Dict{Int,Int}()
        n_ext_spa = n_dofs
        A_with_rows_debug = nothing   # captured post-restriction, pre-prolongation for debug comparison

        # ── Stage 5: Apply spatial constraints to LHS matrix (uniform angular mesh) ────
        if has_spatial_hanging_nodes && R_spatial !== nothing
            @info "[$rank] Applying spatial constraints to LHS matrix (Stage 5, MPI-aware)..."

            # Build helper data structures for MPI pattern
            # gip_to_local: global DOF ID → local DOF index
            gip_to_local_spa = Dict{Int,Int}()
            sizehint!(gip_to_local_spa, n_spa_new)
            for ip = 1:n_spa_new
                gip_to_local_spa[ip2gip_spa[ip]] = ip
            end

            # gip2owner_spa_gid: global DOF ID → owner rank (for use in row/col effects functions)
            gip2owner_spa_gid = Dict{Int,Int}()
            sizehint!(gip2owner_spa_gid, n_spa_new)
            for ip = 1:n_spa_new
                gip2owner_spa_gid[ip2gip_spa[ip]] = gip2owner_extra[ip]
            end
            # Add ghost parent DOF compact GIDs for cross-rank NCFs.
            # Use global_sp_ang_to_gid (built above) — works even when the parent is not
            # in this rank's local elements (and thus absent from point_dict_spa).
            for k = 1:length(mesh.pgip_ghost)
                parent_gip_spa = Int(mesh.pgip_ghost[k])
                parent_owner   = Int(mesh.pgip_owner[k])
                parent_gip_spa <= 0 && continue
                parent_gip_spa in _all_needed_set || continue
                for ip_a = 1:n_ang
                    compact_gid = get(global_sp_ang_to_gid, (parent_gip_spa, ip_a), 0)
                    compact_gid == 0 && continue
                    gip2owner_spa_gid[compact_gid] = parent_owner
                end
            end

            # Build old→new DOF mapping: (spatial_ip-1)*n_ang + ang_ip → dedup ip_spa.
            # Used to convert old-formula DOF indices to the new coordinate-dedup indices.
            old_to_new_dof = zeros(Int, npoin * n_ang)
            for ip_s = 1:npoin
                x_s = mesh.x[ip_s]; y_s = mesh.y[ip_s]; z_s = mesh.z[ip_s]
                for ip_a = 1:n_ang
                    θ_a = extra_mesh.extra_coords[1, ip_a]
                    ϕ_a = extra_mesh.extra_coords[2, ip_a]
                    key_o = (round(x_s, digits=12), round(y_s, digits=12), round(z_s, digits=12),
                             round(θ_a, digits=12), round(ϕ_a, digits=12))
                    old_to_new_dof[(ip_s - 1) * n_ang + ip_a] = get(point_dict_spa, key_o, 0)
                end
            end

            # Find all interpolated parent spatial IPs (coincident nodes are transparent)
            parent_spatial_ips = Set{Int}()
            for (_, weights) in spatial_amr_cache.parent_weights
                for (parent_ip, _) in weights
                    push!(parent_spatial_ips, parent_ip)
                end
            end

            nonowned_parent_gids_spa = Set{Int}()
            nonowned_parent_dofs_spa  = Int[]
            for spatial_ip in parent_spatial_ips
                for ip_a = 1:n_ang
                    dof = old_to_new_dof[(spatial_ip - 1) * n_ang + ip_a]
                    (dof == 0 || dof > n_spa_new) && continue
                    if gip2owner_extra[dof] != rank
                        push!(nonowned_parent_dofs_spa, dof)
                        push!(nonowned_parent_gids_spa, ip2gip_spa[dof])
                    end
                end
            end

            # cross_rank_parent_weights may also contain parents with local representation
            # (shared boundary nodes whose NCF face is owned by another rank).
            # Directly map (raw_GIP, ang_ip) → compact_GID via global_sp_ang_to_gid,
            # then invert ip2gip_spa to find the local dedup DOF — no raw spatial GIPs needed.
            let _gid_to_local_dof = Dict{Int,Int}(Int(ip2gip_spa[ip]) => ip for ip = 1:n_spa_new),
                _seen = Set(nonowned_parent_dofs_spa)
                for (_, cross_weights) in spatial_amr_cache.cross_rank_parent_weights
                    for (parent_global_spa, _) in cross_weights
                        for ip_a = 1:n_ang
                            compact_gid = get(global_sp_ang_to_gid, (parent_global_spa, ip_a), 0)
                            compact_gid == 0 && continue
                            local_dof = get(_gid_to_local_dof, compact_gid, 0)
                            local_dof == 0 && continue
                            gip2owner_extra[local_dof] == rank && continue
                            local_dof in _seen && continue
                            push!(_seen, local_dof)
                            push!(nonowned_parent_dofs_spa, local_dof)
                            push!(nonowned_parent_gids_spa, ip2gip_spa[local_dof])
                        end
                    end
                end
            end

            # Build two ghost constraint dicts, mirroring the angular adaptive pattern:
            #
            # ghost_constraint_data_spa     → used for MATRIX row/col effects (Stage 5)
            #   Contains ALL local parent GIDs (owned + non-owned) + cross-rank GIDs.
            #   The matrix effects functions skip owned parents via `owner == rank && continue`,
            #   so including them is harmless but keeps the pattern consistent with nc_mat.
            #
            # ghost_constraint_data_spa_rhs → used for RHS effects (Stage 4)
            #   Contains ONLY non-owned local parent GIDs + cross-rank GIDs.
            #   compute_hanging_rhs_effects_before_restriction has NO ownership skip,
            #   so owned parents must be excluded (they're already in R_spatial_rhs * RHS).
            #
            # Also builds R_spatial_rhs: the RHS restriction matrix with ONLY owned-parent
            # constraints (mirrors nc_mat_rhs in the angular adaptive case).

            ghost_constraint_data_spa     = Dict{Int, Vector{Tuple{Int, Float64}}}()
            ghost_constraint_data_spa_rhs = Dict{Int, Vector{Tuple{Int, Float64}}}()

            # Pre-compute which dedup DOFs are already handled via parent_weights →
            # R_spatial_rhs (local restriction).  With coordinate dedup, two different
            # spatial IPs at the same (x,y,z) can map to the same hanging_dof — one in
            # parent_weights (local, → R_spatial_rhs) and one in cross_rank_parent_weights
            # (cross-rank, → ghost exchange).  Including such a DOF in both paths would
            # apply the restriction twice.  Exclude from ghost_constraint_data_spa_rhs any
            # hanging_dof already covered by R_spatial_rhs.
            _rhs_handled_dofs = Set{Int}()
            for (hsp, _) in spatial_amr_cache.parent_weights
                for _ipa = 1:n_ang
                    _d = old_to_new_dof[(hsp - 1) * n_ang + _ipa]
                    _d == 0 && continue
                    push!(_rhs_handled_dofs, _d)
                end
            end

            # ── Source 1: fully cross-rank parents ─────────────────────────────
            # hanging_dof uses old_to_new_dof for dedup index.
            # Parent compact GID from global_sp_ang_to_gid: works for parents that are ghost
            # nodes (not in local elements) and thus absent from point_dict_spa.
            # ghost_constraint_data_spa_rhs gets the entry ONLY if the hanging_dof has no
            # local representation in parent_weights (otherwise R_spatial_rhs handles it).
            gip_to_local_r = Dict{Int,Int}(Int(ip2gip_spa[ip]) => ip for ip = 1:n_spa_new)
            for (hanging_local_spa, cross_weights) in spatial_amr_cache.cross_rank_parent_weights
                for ip_a = 1:n_ang
                    hanging_dof = old_to_new_dof[(hanging_local_spa - 1) * n_ang + ip_a]
                    (hanging_dof == 0 || hanging_dof > n_spa_new) && continue
                    parent_constraints = Tuple{Int, Float64}[]
                    for (parent_global_spa, weight) in cross_weights
                        parent_gid = get(global_sp_ang_to_gid, (parent_global_spa, ip_a), 0)
                        parent_gid == 0 && continue
                        push!(parent_constraints, (parent_gid, weight))
                        local_ip = get(gip_to_local_r, parent_gid, 0)
                        local_ip == 0 && continue
                        R_spatial[local_ip, hanging_dof] = weight
                    end
                    isempty(parent_constraints) && continue
                    ghost_constraint_data_spa[hanging_dof] = copy(parent_constraints)
                    if !(hanging_dof in _rhs_handled_dofs)
                        ghost_constraint_data_spa_rhs[hanging_dof] = copy(parent_constraints)
                    end
                end
            end
            P_spatial = sparse(R_spatial')
            # ── Source 2: local parents split by ownership ────────────────────
            # Build R_spatial_rhs (owned parents only) and populate both ghost dicts.
            # All DOF indices use old_to_new_dof for dedup conversion.
            I_rhs_sp = Int[]; J_rhs_sp = Int[]; V_rhs_sp = Float64[]
            for dof = 1:n_spa_new
                if !(dof in spatial_hanging_nodes_all_angular)
                    push!(I_rhs_sp, dof); push!(J_rhs_sp, dof); push!(V_rhs_sp, 1.0)
                end
            end
            n_nonowned_local_added = 0
            for (hanging_local_spa, local_weights) in spatial_amr_cache.parent_weights
                for ip_a = 1:n_ang
                    hanging_dof = old_to_new_dof[(hanging_local_spa - 1) * n_ang + ip_a]
                    (hanging_dof == 0 || hanging_dof > n_spa_new) && continue
                    for (parent_local_spa, weight) in local_weights
                        parent_dof = old_to_new_dof[(parent_local_spa - 1) * n_ang + ip_a]
                        (parent_dof == 0 || parent_dof > n_spa_new) && continue
                        parent_gid = ip2gip_spa[parent_dof]
                        # Matrix dict gets ALL local parents (owned ones skipped at use time)
                        mat_entry = get!(ghost_constraint_data_spa, hanging_dof, Tuple{Int,Float64}[])
                        if !any(g == parent_gid for (g, _) in mat_entry)
                            push!(mat_entry, (parent_gid, weight))
                        end
                        if gip2owner_extra[parent_dof] == rank
                            # Owned parent → goes into R_spatial_rhs for local RHS multiply
                            push!(I_rhs_sp, parent_dof)
                            push!(J_rhs_sp, hanging_dof)
                            push!(V_rhs_sp, weight)
                        else
                            # Non-owned local parent → RHS dict (will be exchanged via MPI)
                            rhs_entry = get!(ghost_constraint_data_spa_rhs, hanging_dof, Tuple{Int,Float64}[])
                            if !any(g == parent_gid for (g, _) in rhs_entry)
                                push!(rhs_entry, (parent_gid, weight))
                                n_nonowned_local_added += 1
                            end
                        end
                    end
                end
            end
            # Sources 3 & 4 (coincident nodes) are omitted: with coordinate dedup,
            # coincident child and parent share the same ip_spa → same compact GID →
            # no explicit constraint row needed; they're transparent in assembly.

            R_spatial_rhs = sparse(I_rhs_sp, J_rhs_sp, V_rhs_sp, n_spa_new, n_spa_new)

            @info "[$rank] ghost_constraint_data_spa (matrix): $(length(ghost_constraint_data_spa)) hanging DOFs"
            @info "[$rank] ghost_constraint_data_spa_rhs: $(length(ghost_constraint_data_spa_rhs)) hanging DOFs"
            @info "[$rank] non-owned local parent RHS entries: $n_nonowned_local_added"
            @info "[$rank] R_spatial_rhs: $(size(R_spatial_rhs)), nnz=$(nnz(R_spatial_rhs))"
            @info "maxima of R_spatial_rhs", maximum(R_spatial_rhs), minimum(R_spatial_rhs)
            @info "maxima of R_spatial", maximum(R_spatial), minimum(R_spatial)
            # ── Build extended parent numbering for cross-rank parents ─────────
            # For each cross-rank parent GID that appears in ghost_constraint_data_spa,
            # assign an extended local index > npoin_ang_total on this rank.
            # This mirrors the angular adaptive case's gid_to_extended_parents mechanism:
            # rank 1 (hanging rank) computes cross-rank parent P0's row/col effects and
            # stores them in extended rows/cols of the local matrix.  The parallel matvec
            # then AllReduces these extended-row contributions to P0's global residual.
            n_ghost_ext_spa = 0
            for (_, parent_list) in ghost_constraint_data_spa
                for (parent_gid, _) in parent_list
                    owner = get(gip2owner_spa_gid, parent_gid, rank)
                    if owner != rank &&
                       !haskey(gid_to_extended_parents_spa, parent_gid) &&
                       !haskey(gip_to_local_spa, parent_gid)
                        n_ghost_ext_spa += 1
                        gid_to_extended_parents_spa[parent_gid] = n_spa_new + n_ghost_ext_spa
                        push!(extended_parents_to_gid_spa, parent_gid)
                    end
                end
            end
            n_ext_spa = n_spa_new + n_ghost_ext_spa
            @info "[$rank] Extended parent structures: $n_ghost_ext_spa cross-rank parent GIDs → n_ext=$n_ext_spa"

            # Extended gip→local map that also resolves extended parent GIDs
            gip_to_local_spa_ext = copy(gip_to_local_spa)
            for (pgid, ext_idx) in gid_to_extended_parents_spa
                gip_to_local_spa_ext[pgid] = ext_idx
            end

            # Extended R and P: same as R_spatial/P_spatial but with identity rows/cols
            # at every extended parent index so that extended rows survive R * A_eff.
            if n_ghost_ext_spa > 0
                I_re, J_re, V_re = findnz(R_spatial)
                for k = 1:n_ghost_ext_spa
                    ext_idx = n_spa_new + k
                    push!(I_re, ext_idx); push!(J_re, ext_idx); push!(V_re, 1.0)
                end
                R_spatial_ext = sparse(I_re, J_re, V_re, n_ext_spa, n_ext_spa)
                P_spatial_ext = sparse(R_spatial_ext')
            else
                R_spatial_ext = R_spatial
                P_spatial_ext = P_spatial
            end

            # ── Step 1: AllReduce parent–parent entries ───────────────────────
            if nprocs > 1
                @info "[$rank] AllReducing spatial parent-parent matrix entries..."
                #=A = allreduce_parent_parent_entries(
                    A, R_spatial_ext, gip2owner_extra, n_spa_new, rank, comm,
                    ip2gip_spa, Int[], nonowned_parent_dofs_spa,
                    nonowned_parent_gids_spa, gip_to_local_spa)=#
                @info "[$rank] A after AllReduce: nnz=$(nnz(A))"
            end

            # ── Step 2: Parallel restriction (R from left) ───────────────────
            # Cross-rank parent effects go into extended rows of A_eff (via extended
            # parent structures); R_spatial_ext has identity at those rows so they survive
            # the left-multiply.  No MPI exchange needed for the extended-parent effects.
            @info "[$rank] ghost_constraint_data_spa: $(length(ghost_constraint_data_spa)) hanging DOFs, parent_weights: $(length(spatial_amr_cache.parent_weights)), cross_rank_parent_weights: $(length(spatial_amr_cache.cross_rank_parent_weights))"
            @info "[$rank] Applying parallel spatial restriction (left multiply by R_ext)..."
            row_effects_to_send, A_eff =
                compute_hanging_row_effects_before_restriction(
                    ghost_constraint_data_spa, A, ip2gip_spa, gip2owner_spa_gid, rank,
                    gid_to_extended_parents_spa, extended_parents_to_gid_spa,
                    gip_to_local_spa_ext)
            @info "[$rank] A_eff after compute_row_effects: size=$(size(A_eff)), nnz=$(nnz(A_eff))"
            for (dest, eff) in row_effects_to_send
                @info "[$rank] → sending $(length(eff)) row effects to rank $dest"
            end
            if isempty(row_effects_to_send)
                @info "[$rank] → no row effects to send (all cross-rank parents in extended set)"
            end

            # Left-multiply by extended R: local hanging rows get constraints applied;
            # extended parent rows (identity in R_spatial_ext) pass through unchanged.
            A_left = R_spatial_ext * A_eff
            @info "[$rank] A after R_ext * A_eff: size=$(size(A_left)), nnz=$(nnz(A_left))"

            received_row_effects = exchange_hanging_effects(row_effects_to_send, rank, comm)
            for (src, eff) in received_row_effects
                @info "[$rank] ← received $(length(eff)) row effects from rank $src"
            end

            A_with_rows = add_hanging_row_effects(
                A_left, received_row_effects, ip2gip_spa,
                n_ext_spa, rank, gip_to_local_spa_ext)
            @info "[$rank] A after add_row_effects: nnz=$(nnz(A_with_rows))"
            A_with_rows_debug = copy(A_with_rows)   # snapshot for debug comparison

            # ── Step 3: Parallel prolongation (P from right) ─────────────────
            @info "[$rank] Applying parallel spatial prolongation (right multiply by P_ext)..."
            col_effects_to_send, A_col_eff =
                compute_hanging_col_effects_before_prolongation(
                    ghost_constraint_data_spa, A_with_rows, ip2gip_spa, gip2owner_spa_gid,
                    n_spa_new, n_ext_spa,
                    spatial_hanging_nodes_all_angular, rank,
                    gid_to_extended_parents_spa, extended_parents_to_gid_spa,
                    gip_to_local_spa_ext)
            @info "[$rank] A after compute_col_effects: nnz=$(nnz(A_col_eff))"
            for (dest, eff) in col_effects_to_send
                @info "[$rank] → sending $(length(eff)) col effects to rank $dest"
            end
            if isempty(col_effects_to_send)
                @info "[$rank] → no col effects to send (all cross-rank parents in extended set)"
            end

            # Right-multiply by extended P: hanging columns get prolongation applied;
            # extended parent columns (identity in P_spatial_ext) pass through unchanged.
            A_both = A_col_eff * P_spatial_ext
            @info "[$rank] A after *P_spatial_ext: size=$(size(A_both)), nnz=$(nnz(A_both)), nnz(P_spatial_ext)=$(nnz(P_spatial_ext))"

            received_col_effects = exchange_hanging_effects(col_effects_to_send, rank, comm)
            for (src, eff) in received_col_effects
                @info "[$rank] ← received $(length(eff)) col effects from rank $src"
            end

            A_with_cols = add_hanging_col_effects(
                A_both, received_col_effects, ip2gip_spa,
                n_ext_spa, rank, gip_to_local_spa_ext)
            @info "[$rank] A after add_col_effects: nnz=$(nnz(A_with_cols))"

            A = sparse(A_with_cols)  # n_ext_spa × n_ext_spa
            
            # ── Step 4: Zero non-owned parent entries (prevent GMRES double-count) ──
            if nprocs > 1 && !isempty(nonowned_parent_dofs_spa)
                @info "[$rank] Zeroing $(length(nonowned_parent_dofs_spa)) non-owned parent DOF entries..."
                n_zeroed = 0
                for i in nonowned_parent_dofs_spa
                    i <= size(A, 1) || continue
                    for j in nonowned_parent_dofs_spa
                        j <= size(A, 2) || continue
                        if abs(A[i,j]) > 0.0
                            #A[i,j] = 0.0
                            #A_with_rows_debug[i,j] = 0.0
                            n_zeroed += 1
                        end
                    end
                end
                A_with_rows_debug = sparse(A_with_rows_debug)
                A = sparse(A)
                @info "[$rank] A after Step4 zero (n_zeroed=$n_zeroed): nnz=$(nnz(A))"
            end

            @info "[$rank] ✓ Spatial LHS constraints applied (MPI-aware): A = $(size(A)), final nnz=$(nnz(A))"
        end
    end  # adaptive / non-adaptive

    # ── Common setup: nc_mat row structure, free DOF count ────────────────────
    nc_rows  = size(nc_mat, 1) > 2 ? rowvals(nc_mat) : zeros(Int, 1, 1)
    # Use all_hanging_nodes (the full combined set for all adaptive cases) rather than
    # n_non_global_nodes (angular-only) so that the combined angular+spatial path
    # correctly excludes both types of hanging DOFs from the free-DOF count.
    n_free   = inputs[:adaptive_extra_meshes] ? n_spa - length(all_hanging_nodes) :
               (has_spatial_hanging_nodes ? n_spa_new - n_non_global_nodes_spa : npoin_ang_total)

    # ── Precompute spatial factors (shared by RHS and BC loops) ──────────────
    spatial_factor = Vector{Float64}(undef, npoin)
    for ip = 1:npoin
        gip = exp(-((1/3) * (mesh.x[ip] - 1.0))^2)
        hip = exp(-4.0  * (2 - mesh.y[ip]) / 2)
        spatial_factor[ip] = gip * hip
    end

    # Precompute boundary face lookup: spatial node → list of (iface, fi, fj)
    node_to_bdy_faces = Dict{Int, Vector{Tuple{Int,Int,Int}}}()
    for iface = 1:mesh.nfaces_bdy
        for iter_i = 1:ngl, iter_j = 1:ngl
            ip1  = mesh.poin_in_bdy_face[iface, iter_i, iter_j]
            list = get!(node_to_bdy_faces, ip1, Tuple{Int,Int,Int}[])
            push!(list, (iface, iter_i, iter_j))
        end
    end

    boundary_dict = Dict{Int, Float64}()
    bdy_normals   = Dict{Int, Vector{NTuple{3,Float64}}}()
    bdy_normals_for_ghosts = Dict{Int, NTuple{3,Float64}}()
    # ── RHS and boundary condition assembly ──────────────────────────────────
    @rankinfo rank "Assembling RHS and boundary conditions..."
   
    for iel = 1:nelem
        for i = 1:ngl, j = 1:ngl, k = 1:ngl
            ip  = mesh.connijk[iel, i, j, k]
            x   = mesh.x[ip]; y = mesh.y[ip]; z = mesh.z[ip]
            sf  = spatial_factor[ip]
            is_boundary = ip in mesh.poin_in_bdy_face

            if is_boundary && !haskey(bdy_normals, ip)
                iface  = elem_to_face[iel, i, j, k, 1]
                face_i = elem_to_face[iel, i, j, k, 2]
                face_j = elem_to_face[iel, i, j, k, 3]

                # Collect all distinct face normals directly from node_to_bdy_faces.
                # This avoids the coordinate-equality fragility of the nmatches approach
                # entirely — we just gather every distinct normal from every face that
                # owns this node.
                collected = NTuple{3,Float64}[]
                for (iface2, fi2, fj2) in get(node_to_bdy_faces, ip, [])
                    nxyz = (nx[iface2, fi2, fj2],
                            ny[iface2, fi2, fj2],
                            nz[iface2, fi2, fj2])
                    is_dup = any(abs(n[1]-nxyz[1]) < 1e-12 &&
                                 abs(n[2]-nxyz[2]) < 1e-12 &&
                                 abs(n[3]-nxyz[3]) < 1e-12 for n in collected)
                    is_dup || push!(collected, nxyz)
                end

                # Fallback: if node_to_bdy_faces didn't cover this node's first face
                if isempty(collected)
                    push!(collected, (nx[iface, face_i, face_j],
                                    ny[iface, face_i, face_j],
                                    nz[iface, face_i, face_j]))
                end

                bdy_normals[ip] = collected
            end

            face_normals = is_boundary ? bdy_normals[ip] : NTuple{3,Float64}[]

            if inputs[:adaptive_extra_meshes]
                for e_ext = 1:extra_meshes_extra_nelems[iel]
                    nop = extra_meshes_extra_nops[iel][e_ext]
                    for jθ = 1:nop+1, iθ = 1:nop+1
                        ip_ext  = extra_meshes_connijk[iel][e_ext, iθ, jθ]
                        θ       = extra_meshes_coords[iel][1, ip_ext]
                        ϕ       = extra_meshes_coords[iel][2, ip_ext]
                        ip_g    = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        is_owned = (gip2owner_extra[ip_g] == rank)

                        sip = exp(-((6/(2π)) * (θ - 3π/5))^2)
                        bip = exp(-((6/(2π)) * (ϕ - 2π/3))^2)
                        ref[ip_g] = sf * sip * bip

                        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                        if is_boundary && !(ip_g in all_hanging_nodes)
                            # Find the most-inflow face normal for this direction:
                            # the one giving the most negative Ω·n.
                            # If the minimum dot product is < -1e-13 the direction
                            # is inflow through at least one face → apply BC.
                            best_dot = 0.0
                            best_nx  = face_normals[1][1]
                            best_ny  = face_normals[1][2]
                            best_nz  = face_normals[1][3]
                            
                            for (fnx, fny, fnz) in face_normals
                                
                                d = fnx*Ωx + fny*Ωy + fnz*Ωz
                                if d < best_dot
                                    best_dot = d
                                    best_nx  = fnx; best_ny = fny; best_nz = fnz
                                end
                            end
                            
                            
                            if best_dot < -1e-10
                                    val = 0.0
                                    if inputs[:RT_shortwave]
                                        val = user_rad_bc_shortwave_diffuse(x, y, z, θ, ϕ, bdy, sw, F_dir[ip], τ_nodes[ip], sw_ω₀_lateral, inputs[:rad_HG_g][ip])
                                    elseif inputs[:RT_longwave]
                                        val = user_rad_bc_longwave(x, y, z, θ, ϕ, best_nx, best_ny, best_nz, ip, atmos_data, bdy, lw)
                                    else
                                        val = user_rad_bc(x, y, z, θ, ϕ)
                                    end
                                    BDY[ip_g] = val
                                    if !haskey(boundary_dict, ip_g)
                                        boundary_dict[ip_g] = is_owned ? val : 0.0
                                    end
                            else
                                if is_owned
                                    RHS[ip_g] = if inputs[:RT_shortwave]
                                        user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip, F_dir, σ, sw, inputs[:rad_HG_g][ip])
                                    elseif inputs[:RT_longwave]
                                        user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data)
                                    else
                                        user_rhs(x, y, z, θ, ϕ)
                                    end
                                end
                            end
                        else
                            if is_owned
                                RHS[ip_g] = if inputs[:RT_shortwave]
                                    user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip, F_dir, σ, sw, inputs[:rad_HG_g][ip])
                                elseif inputs[:RT_longwave]
                                    user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data)
                                else
                                    user_rhs(x, y, z, θ, ϕ)
                                end
                            end
                        end
                    end
                end
            else
                for e_ext = 1:extra_mesh.extra_nelem
                    nop    = extra_mesh.extra_nop[e_ext]
                    for iθ = 1:nop+1, iϕ = 1:nop+1
                        ip_ext  = extra_mesh.extra_connijk[e_ext, iϕ, iθ]
                        θ       = extra_mesh.extra_coords[1, ip_ext]
                        ϕ       = extra_mesh.extra_coords[2, ip_ext]
                        # Use coordinate-dedup DOF index when spatial AMR is active,
                        # old formula otherwise (gip2owner_extra size matches accordingly).
                        ip_g = has_spatial_hanging_nodes ?
                               connijk_spa_uniform[iel][i, j, k, e_ext, iϕ, iθ] :
                               (ip-1) * extra_mesh.extra_npoin + ip_ext
                        is_owned = (gip2owner_extra[ip_g] == rank)

                        sip = exp(-((6/(2π)) * (θ - 3π/5))^2)
                        bip = exp(-((6/(2π)) * (ϕ - 2π/3))^2)
                        ref[ip_g] = sf * sip * bip

                        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                        if is_boundary && !(ip_g in spatial_hanging_nodes_all_angular)
                            best_dot = 0.0
                            best_nx  = face_normals[1][1]
                            best_ny  = face_normals[1][2]
                            best_nz  = face_normals[1][3]
                            for (fnx, fny, fnz) in face_normals
                                d = fnx*Ωx + fny*Ωy + fnz*Ωz
                                
                                if d < best_dot
                                    best_dot = d
                                    best_nx  = fnx; best_ny = fny; best_nz = fnz
                                end
                            end
                            
                            if best_dot < -1e-10
                               
                            
                                val = 0.0
                                if inputs[:RT_shortwave]
                                    val = user_rad_bc_shortwave_diffuse(x, y, z, θ, ϕ, bdy, sw, F_dir[ip], τ_nodes[ip], sw_ω₀_lateral, inputs[:rad_HG_g][ip])
                                elseif inputs[:RT_longwave]
                                    val = user_rad_bc_longwave(x, y, z, θ, ϕ, best_nx, best_ny, best_nz, ip, atmos_data, bdy, lw)
                                else
                                    val = user_rad_bc(x, y, z, θ, ϕ)
                                end
                                BDY[ip_g] = val
                                if !haskey(boundary_dict, ip_g)
                                    boundary_dict[ip_g] = is_owned ? val : 0.0
                                end
                                
                            else
                                if is_owned
                                    RHS[ip_g] = if inputs[:RT_shortwave]
                                        user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip, F_dir, σ, sw, inputs[:rad_HG_g][ip])
                                    elseif inputs[:RT_longwave]
                                        user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data)
                                    else
                                        user_rhs(x, y, z, θ, ϕ)
                                    end
                                end
                            end
                        else
                            if is_owned
                                RHS[ip_g] = if inputs[:RT_shortwave]
                                    user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip, F_dir, σ, sw, inputs[:rad_HG_g][ip])
                                elseif inputs[:RT_longwave]
                                    user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data)
                                else
                                    user_rhs(x, y, z, θ, ϕ)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    #=for i=1:npoin_ang_total
        gip = ip2gip_spa[i]
        if (gip == 56710)
            @info "before 1", RHS[i]
        elseif (gip == 48347)
            @info "before 4", RHS[i]
        end
    end=#
    # Verify neighbour / NCF consistency for the purely-spatial-AMR path
    # (adaptive_extra_meshes=false, uniform angular mesh).
    # extra_meshes_ref_level is not available here; pass nothing → all angular levels = 0.
    if !inputs[:adaptive_extra_meshes]
        begin
            _ang_ref_records_spa_only = build_element_refinement_records(
                mesh, nothing, nelem, ngl, comm)
            verify_element_refinement_records(
                _ang_ref_records_spa_only, nelem, comm, rank;
                uniform_spatial = !has_spatial_hanging_nodes)
        end
    end

    # ── Stage 4: Apply spatial constraints to RHS (uniform angular mesh, MPI-aware) ──
    if !inputs[:adaptive_extra_meshes] && has_spatial_hanging_nodes && R_spatial !== nothing
        @info "[$rank] Applying spatial constraints to RHS (Stage 4, MPI-aware)..."
        n_ang = extra_mesh.extra_npoin
        RHS_original = copy(RHS)
        # Step 1: compute weighted RHS contributions from cross-rank hanging nodes
        # and route them to the owning rank of each ghost parent.
        # Uses ghost_constraint_data_spa_rhs (non-owned local + cross-rank only),
        # since owned local parents are already handled by R_spatial_rhs * RHS.
        rhs_effects_to_send =
            compute_hanging_rhs_effects_before_restriction(
                ghost_constraint_data_spa_rhs, RHS, ip2gip_spa, gip2owner_spa_gid, rank)
        
        # Step 2: local restriction using ownership-filtered R_spatial_rhs.
        # Only owned-parent constraints are applied locally; non-owned local and
        # cross-rank parents were sent in Step 1 via ghost_constraint_data_spa_rhs.
        RHS = R_spatial_rhs*RHS#R_spatial_rhs * RHS

        for i=1:n_dofs
            gip = ip2gip_spa[i]
            if (gip == 4181)
                @info "restricted 1", RHS[i]
            elseif (gip == 4417)
                @info "restricted 4", RHS[i]
            end
        end
        #=for i=1:npoin_ang_total
            gip = ip2gip_spa[i]
            if (gip == 56710)
                @info "restricted 1", RHS[i]
            elseif (gip == 48347)
                @info "restricted 4", RHS[i]
            end
        end=#
        # Step 3: exchange cross-rank effects (collective — all ranks participate)
        received_rhs_effects =
            exchange_hanging_effects_vector(rhs_effects_to_send, rank, comm)

        # Step 4: accumulate received effects into parent rows
        RHS = add_hanging_rhs_effects(
            RHS, received_rhs_effects, ip2gip_spa, n_spa_new, rank, gip_to_local_spa)

        for i=1:n_spa_new
            gip = ip2gip_spa[i]
            if (gip == 4181)
                @info "added 1", RHS[i]
            elseif (gip == 4417)
                @info "added 4", RHS[i]
            end
        end
        #=for i=1:npoin_ang_total
            gip = ip2gip_spa[i]
            if (gip == 56710)
                @info "added 1", RHS[i]
            elseif (gip == 48347)
                @info "added 4", RHS[i]
            end
        end=#
        @info "[$rank] ✓ Spatial RHS constraints applied (MPI-aware)"
    end

    
    # ── Ghost boundary nodes (extended ghost-parent DOFs) ─────────────────────
    # extended_parents_ip/x/y/z/θ/ϕ are populated only by the angular adaptive path.
    # Spatial extended parents carry no angular coordinate data here, so skip this
    # block when there is no angular adaptation (extended_parents_ip is empty).
    if inputs[:adaptive_extra_meshes] && n_spa_g > n_spa && !isempty(extended_parents_ip)
        n_extended = n_spa_g - n_spa
        for i = 1:n_extended
            ip   = extended_parents_ip[i]
            x    = extended_parents_x[i]; y = extended_parents_y[i]; z = extended_parents_z[i]
            θ    = extended_parents_θ[i]; ϕ = extended_parents_ϕ[i]
            ip_g = n_spa + i
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            if ip in mesh.poin_in_bdy_face
                face_normals = bdy_normals[ip]
                best_dot = 0.0
                best_nx  = face_normals[1][1]
                best_ny  = face_normals[1][2]
                best_nz  = face_normals[1][3]
                for (fnx, fny, fnz) in face_normals
                    d = fnx*Ωx + fny*Ωy + fnz*Ωz
                    
                    if d < best_dot
                        best_dot = d
                        best_nx  = fnx; best_ny = fny; best_nz = fnz
                    end
                end
                if best_dot < -1e-13
                    val = user_rad_bc(x, y, z, θ, ϕ)
                    BDY[ip_g] = val
                    boundary_dict[ip_g] = val
                end
            end
        end
    end

    # ── Boundary condition application (modify matrix rows) ───────────────────
        
        boundary_set = Set(keys(boundary_dict))
        rows_A = rowvals(A)
        vals_A = nonzeros(A)
        
        

        for col = 1:size(A, 2)
            for ptr = nzrange(A, col)
                row = rows_A[ptr]
                row in boundary_set || continue
                if row == col
                    vals_A[ptr] = (row <= n_spa && gip2owner_extra[row] == rank) ? 1.0 : 0.0
                else
                    vals_A[ptr] = 0.0
                end
            end
        end
        A = dropzeros!(A)
        
        
        @info "[$rank] A after BC application + dropzeros: nnz=$(nnz(A))"


    # ── RHS restriction ───────────────────────────────────────────────────────
    if inputs[:adaptive_extra_meshes]
        #RHS_mass_weighted = Md .* RHS
        @rankinfo rank "Restricting RHS..."
        @info maximum(RHS), minimum(RHS)
        rhs_effects_to_send =
            compute_hanging_rhs_effects_before_restriction(
                ghost_constraint_data_rhs, RHS, ip2gip_spa, gip2owner_spa, rank)
                #ghost_constraint_data_rhs, RHS_mass_weighted, ip2gip_spa, gip2owner_spa, rank)

        RHS_restricted = nc_mat_rhs * RHS#_mass_weighted

        received_rhs_effects =
            exchange_hanging_effects_vector(rhs_effects_to_send, rank, comm)

        RHS_with_effects = add_hanging_rhs_effects(
            RHS_restricted, received_rhs_effects, ip2gip_spa, n_spa, rank, gip_to_local)

        RHS_red = extract_free_rhs_subvector(
            RHS_with_effects, all_hanging_nodes, n_free, n_spa, n_spa_g, rank)

        for (node, val) in boundary_dict
            if val != 0.0
                RHS_red[node] = val
                
            end
        end
        
        B = RHS_red
    
        @info "maxima of restricted RHS"
        @info maximum(RHS), minimum(RHS)
    else
        for (node, val) in boundary_dict
            RHS[node] = val
        end
        B = RHS
    end

    # ── Coincident node enforcement (uniform angular path) ────────────────────
    # Coincident nodes sit at exactly the same (x,y,z) as a parent node but have
    # a different local spatial IP. Their assembled rows duplicate the parent's
    # contribution. Zero their rows in A (identity diagonal) and zero their RHS;
    # after the solve, copy the parent solution into them.
    @info "maxima of A right before coicidence handling", maximum(A), minimum(A)
    # Coincident node enforcement is no longer needed: coordinate deduplication in
    # connijk_spa_uniform ensures child and parent nodes share the same ip_spa, so
    # their assembled contributions naturally merge into a single DOF slot. No
    # explicit row-zeroing is required.

    #-Symmetrize boundary conditions:
    As  = sparse(A)
    rows_A = rowvals(As)
    vals_A = nonzeros(As)
    @info "maxima of A right before comparison", maximum(As), minimum(As)
    #=if !inputs[:adaptive_extra_meshes] && has_spatial_hanging_nodes && R_spatial !== nothing
        _outdir = "."
        _spa_amr_ref = joinpath(_outdir, "spa_amr_serial_reference.jld2")

        # Build coordinate array for dedup DOFs: invert point_dict_spa.
        # dedup_coords[ip_spa, :] = (x, y, z, θ, ϕ) for DOF ip_spa.
        _dedup_coords = Matrix{Float64}(undef, n_spa_new, 5)
        for (key, ip_spa) in point_dict_spa
            _dedup_coords[ip_spa, 1] = key[1]; _dedup_coords[ip_spa, 2] = key[2]
            _dedup_coords[ip_spa, 3] = key[3]; _dedup_coords[ip_spa, 4] = key[4]
            _dedup_coords[ip_spa, 5] = key[5]
        end

        # ── Targeted RHS diagnostic ──────────────────────────────────────────────
        # Print B[ip_spa] for a specific (x,y,z,θ,ϕ) to diagnose serial/parallel mismatch.
        _tgt = (1.000000000001, 1.0, 1.333333333334, 1.230305638774, 2.094395102393)
        _tol = 1e-7
        for (key, ip_spa) in point_dict_spa
            if abs(key[1]-_tgt[1]) < _tol && abs(key[2]-_tgt[2]) < _tol &&
               abs(key[3]-_tgt[3]) < _tol && abs(key[4]-_tgt[4]) < _tol &&
               abs(key[5]-_tgt[5]) < _tol
                is_hanging = ip_spa in spatial_hanging_nodes_all_angular
                owner      = gip2owner_extra[ip_spa]
                gid        = ip2gip_spa[ip_spa]
                @info "[$rank] TARGET coord: ip_spa=$ip_spa gid=$gid owner=$owner " *
                      "B=$(B[ip_spa]) original RHS=$(RHS_original[ip_spa]) hanging=$is_hanging"
            end
        end

        # Build coordinate lookup for extended parent rows in A_with_rows_debug.
        # global_gid_to_coords was AllGathered alongside global_sp_ang_to_gid:
        # every rank that owns a cross-rank parent contributed its (x,y,z,θ,ϕ).
        _ext_parent_coords = NTuple{5,Float64}[]
        for gid in extended_parents_to_gid_spa
            push!(_ext_parent_coords, get(global_gid_to_coords, gid, (NaN, NaN, NaN, NaN, NaN)))
        end
        _n_nan_ext = count(c -> any(isnan, c), _ext_parent_coords)
        _n_nan_ext > 0 && @warn "[$rank] $_n_nan_ext / $(length(_ext_parent_coords)) extended parent coords are NaN — those entries will be skipped in comparison"
        
        if nprocs == 1
            @info "[$rank] Saving serial spatial-AMR reference → $_spa_amr_ref"
            save_serial_spatial_amr(
                B, As, As, mesh, extra_mesh, n_ang, n_spa_new,
                spatial_hanging_nodes_all_angular, _spa_amr_ref;
                dedup_coords = _dedup_coords)
        elseif isfile(_spa_amr_ref)
            @info "[$rank] Comparing parallel spatial-AMR assembly against serial reference..."
            reduce_and_compare_parallel_spatial_amr(
                B, As, As, mesh, extra_mesh, n_ang, n_spa_new,
                spatial_hanging_nodes_all_angular,
                ip2gip_spa, gnpoin, extended_parents_to_gid_spa,
                comm, rank, _spa_amr_ref;
                dedup_coords      = _dedup_coords,
                gip2owner_extra   = gip2owner_extra,
                ext_parent_coords = _ext_parent_coords)
        else
            @info "[$rank] No serial reference found at $_spa_amr_ref — run with 1 rank first to generate it"
        end
        MPI.Barrier(comm)

        # ── Targeted row diagnostic ───────────────────────────────────────────
        # Target: a Δ=-4 row from the As nnz comparison (serial=71, parallel=67).
        # Run serial (1 rank) → note column list. Run parallel → compare.
        _diag_target = (2.172673164644, 0.0, 0.333333333333, 2.611012560984, 0.0)
        diagnose_row(A_with_rows_debug, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "A_with_rows")
        MPI.Barrier(comm)
        diagnose_row(MLHS, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "MLHS")
        MPI.Barrier(comm)
        diagnose_row(As, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "As")
        MPI.Barrier(comm)
         _diag_target = (1.173, 0.666666666, 1.3333333, 2.611, 2.455)
        diagnose_row(A_with_rows_debug, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "A_with_rows")
        MPI.Barrier(comm)
        diagnose_row(MLHS, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "MLHS")
        MPI.Barrier(comm)
        diagnose_row(As, _diag_target, _dedup_coords, n_spa_new,
                     _ext_parent_coords, ip2gip_spa, rank, comm, "As")
        MPI.Barrier(comm)
        if rank ==1
            @info rank, nnz(R_spatial[299423,:])
            @info "finding local column target for prints on rank", rank
            for i=1:n_spa_new
                if (ip2gip_spa[i] == 20406)
                    @info rank, "local target col ip =", i
                end
            end
            @info rank, MLHS[299423,46930]
            @info rank, A_eff[299423,46930]
            @info rank, A_left[299423,46930]
            @info rank, A_with_rows[299423,46930]
            @info rank, A_with_rows_debug[299423,46930]
            @info rank, A_col_eff[299423,46930]
            @info rank, A_both[299423,46930]
            @info rank, A_with_cols[299423,46930]
            @info rank, As[299423,46930]
        else
            @info rank, nnz(R_spatial[53326,:])
            @info "finding local column target for prints on rank", rank
            for i=1:n_spa_new
                if (ip2gip_spa[i] == 20406)
                    @info rank, "local target col ip =", i
                end
            end
            @info rank, MLHS[53326,53482]
            @info rank, A_eff[53326,53482]
            @info rank, A_left[53326,53482]
            @info rank, A_with_rows[53326,53482]
            @info rank, A_with_rows_debug[53326,53482]
            @info rank, A_col_eff[53326,53482]
            @info rank, A_both[53326,53482]
            @info rank, A_with_cols[53326,53482]
            @info rank, As[53326,53482]
        end

    end=#

    #=for col in boundary_set
        val = BDY[col]
        
        for ptr in nzrange(As, col)
            row = rows_A[ptr]
            row in boundary_set && continue  # skip BC-BC interactions
            row > n_free && continue
            #if(gip2owner_extra[row] == rank)
                B[row] -= vals_A[ptr] * val
            #end
            vals_A[ptr] = 0.0
        end
    end=#
    As = dropzeros!(As)
    if inputs[:adaptive_extra_meshes]
        for ip in all_hanging_nodes
            As[ip,ip] = 1.0
            B[ip] = 0.0
        end
    elseif has_spatial_hanging_nodes && R_spatial !== nothing
        @info length(spatial_hanging_nodes_all_angular)
        for ip in spatial_hanging_nodes_all_angular
            
            As[ip,ip] = 1.0
            B[ip] = 0.0
        end
    end

    n_fixed = 0
    n_small = 0
    n_big = 0
    if (inputs[:adaptive_extra_meshes] && (inputs[:RT_shortwave] || inputs[:RT_longwave]))
    for ip = 1:n_spa
        ip in boundary_set && continue
        if d[ip] < 0
            # Get original physics diagonal from MLHS
            original_diag = MLHS[ip,ip]
            if (original_diag > 1e-6)
                n_big +=1
            end
            if (original_diag < 1e-6)
                n_small +=1
            end
            if original_diag > 0
                As[ip, ip] = original_diag
                #@info "Fixed negative diagonal at node $ip: $(diag(As)[ip]) → $original_diag"
                n_fixed +=1
            end
        end
    end
    end
    @info "n_fixed=$n_fixed, number of large diag values used=$n_big, number of small diag values used=$n_small"

    # ── Parallel correctness check: compare assembled (As, B) against serial ──
    # To use: run once with 1 rank to generate the reference file, then run with
    # multiple ranks.  Set the path below to a writable location.
    

    # ── Linear solve ──────────────────────────────────────────────────────────
    @rankinfo rank "Solving system ($(size(A,1)) DOF)..."

    A   = nothing; GC.gc()
    diagnose_rt_matrix(As, ip2gip_spa, gnpoin, "Longwave")
    d = diag(As)

    #solution=As\B    
    

    npoin_ang_total = size(B, 1)
    @info maximum(As), minimum(As), maximum(B), minimum(B)
    
    solution = if (inputs[:adaptive_extra_meshes] && inputs[:RT_shortwave])
        x_warm = zeros(Float64,n_spa)
        if (rt_sol_sw_available)
            x_warm = rt_sol_sw
        #=else
            diag_vals = diag(As)
            for ip = 1:n_free
                dv = abs(diag_vals[ip])
                x_warm[ip] = dv > 1e-14 ? B[ip] / dv : 0.0
            end=#

        end
        solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_spa, x_warm;
            npoin_g       = n_spa_g,
            g_ip2gip      = extended_parents_to_gid,
            g_gip2ip      = gid_to_extended_parents,
            precond       = :jacobi,
            restart       = 100,
            tol           = 1e-7,
            connijk_spa   = connijk_spa,
            extra_nelem   = extra_meshes_extra_nelems,
            extra_nops    = extra_meshes_extra_nops,
            extra_connijk = extra_meshes_connijk,
            nelem         = nelem,
            ngl           = ngl,
            ladaptive     = true,
            npoin_ang     = extra_meshes_extra_npoins)

    elseif (inputs[:adaptive_extra_meshes] && inputs[:RT_longwave])
        x_warm = zeros(Float64,n_spa)
        if (rt_sol_lw_available)
            x_warm = rt_sol_lw
        #=else
            diag_vals = diag(As)
            for ip = 1:n_free
                dv = abs(diag_vals[ip])
                x_warm[ip] = dv > 1e-14 ? B[ip] / dv : 0.0
            end=#
        end
        solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_spa, x_warm;
            npoin_g  = n_spa_g,
            g_ip2gip = extended_parents_to_gid,
            g_gip2ip = gid_to_extended_parents,
            precond  = :none,
            restart  = 200,
            tol      = 1e-4)

    elseif inputs[:adaptive_extra_meshes]
        # Row scaling (left) — normalize by row norm
        x_warm = Float64[]
        if (inputs[:lmanufactured_solution])
            x_warm = ref
            
            for i in all_hanging_nodes
                x_warm[i] = 0.0
            end
           
        end
        
        #=solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_spa, x_warm;
            npoin_g  = n_spa_g,
            g_ip2gip = extended_parents_to_gid,
            g_gip2ip = gid_to_extended_parents,
            precond  = :ilu,
            restart  = 500,
            tol      = 1e-6)=#
        #diagnose_ghost_rows(As, ip2gip_spa, gip2owner_extra, gnpoin, n_spa_g,
        #            extended_parents_to_gid, MPI.COMM_WORLD)
        
        solve_parallel_gmres_asm(ip2gip_spa, gip2owner_extra, As, B, gnpoin,
            n_spa, x_warm;
            precond     = :global_ilu,
            asm_solver  = :rcmilu,
            asm_ilu_tau = 0.1,
            npoin_g     = n_spa_g,
            g_ip2gip    = extended_parents_to_gid,
            g_gip2ip    = gid_to_extended_parents,
            restart     = 100,
            tol         = 1e-7)

    elseif (inputs[:RT_shortwave])
        x_warm = Float64[]
        if (rt_sol_sw_available)
            x_warm = rt_sol_sw
        end
        solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_dofs, x_warm;
            npoin_g       = n_dofs,
            precond       = :none,
            restart       = 60,
            tol           = 1e-7,
            extra_nelem   = extra_mesh.extra_nelem,
            extra_nops    = extra_mesh.extra_nop,
            extra_connijk = extra_mesh.extra_connijk,
            nelem         = nelem,
            ngl           = ngl,
            ladaptive     = false,
            npoin_ang     = extra_mesh.extra_npoin,
            npoin_space   = npoin)

    elseif (inputs[:RT_longwave])
        x_warm = Float64[]
        if (rt_sol_lw_available)
            x_warm = rt_sol_lw
        end
        solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_dofs, x_warm;
            npoin_g = n_dofs,
            precond = :none,
            restart = 60,
            tol     = 1e-7)
    else
        @info maximum(ip2gip_spa), minimum(ip2gip_spa)
        x_warm = Float64[]
        if (inputs[:lmanufactured_solution])
            x_warm = ref
            if has_spatial_hanging_nodes && R_spatial !== nothing
                for i in spatial_hanging_nodes_all_angular
                    x_warm[i] = 0.0
                end
            end
        end
        #=nonowned_indices = Set{Int}()
        for i=1:npoin_ang_total
            if gip2owner_extra[i] != rank
                push!(nonowned_indices, ip2gip_spa[i])
            end
        end
        gip_to_local = Dict{Int, Int}()
        sizehint!(gip_to_local, npoin_ang_total)
        for ip = 1:npoin_ang_total
            gip_to_local[ip2gip_spa[ip]] = ip
        end
        As = allreduce_shared_entries(
            As, gip2owner_extra, npoin_ang_total, comm, ip2gip_spa,
            nonowned_indices,
            gip_to_local
        )
        As = sparse(As)=#
        #=solve_parallel_gmres(ip2gip_spa, gip2owner_extra, As, B, gnpoin, npoin_ang_total, x_warm;
        npoin_g  = n_ext_spa,
        g_ip2gip = extended_parents_to_gid_spa,
        g_gip2ip = gid_to_extended_parents_spa,
        precond  = :none,
        restart  = 50,
        tol      = 1e-7)=#
        #=target_gid = 5761
        for ip = 1:npoin_ang_total
            @info ip, ip2gip_spa[ip]
            if ip2gip_spa[ip] == target_gid
                diag_val = As[ip, ip]
                row_nnz  = nnz(As[ip:ip, :])
                @info "Rank $rank: GID $target_gid found at local ip=$ip diagonal=$diag_val row_nnz=$row_nnz"
            end
        end=#
        solve_parallel_gmres_asm(ip2gip_spa, gip2owner_extra, As, B, gnpoin, npoin_ang_total, x_warm;
        npoin_g  = n_ext_spa,
        g_ip2gip = extended_parents_to_gid_spa,
        g_gip2ip = gid_to_extended_parents_spa,
        #precond  = :inner_gmres,
        precond  = :global_ilu,
        asm_solver = :rcmsplu,
        asm_ilu_tau = 0.1,
        restart  = 50,
        tol      = 1e-7)

    end
    
    @rankinfo rank "Solve complete."
    #A = nothing; RHS = nothing; GC.gc()
    @info maximum(solution), minimum(solution)
    # ── Solution prolongation ─────────────────────────────────────────────────
    solution_new = zeros(Float64, n_spa)
    if inputs[:adaptive_extra_meshes]
        @rankinfo rank "Prolonging solution to full mesh..."
        solution_contributions_to_send =
            compute_solution_prolongation_contributions(
                reverse_ghost_map, solution, ip2gip_spa, n_free, rank)

        solution_local = P * solution

        received_solution_contributions =
            exchange_hanging_effects_vector(solution_contributions_to_send, rank, comm)

        solution_new = add_solution_prolongation_contributions(
            @view(solution_local[1:n_spa]), received_solution_contributions,
            ip2gip_spa, n_spa, rank, gip_to_local)
        @info "prolonged solution maxima"
        @info maximum(solution_new), minimum(solution_new)
        @info maximum(P_vec), minimum(P_vec), maximum(nc_mat_rhs), minimum(nc_mat_rhs)
    end
    
    if !inputs[:adaptive_extra_meshes] && has_spatial_hanging_nodes && R_spatial !== nothing
        @info "[$rank] Applying spatial prolongation to solution..."

        # Build reverse ghost map via AllGather: each rank broadcasts its
        # (parent_gid, hanging_gid, owner_hanging, weight) tuples for extended parents;
        # owning ranks collect entries where they hold the parent locally.
        _ext_gid_set_prol = Set(extended_parents_to_gid_spa)
        _local_rev = Float64[]
        for (hanging_dof, constraints) in ghost_constraint_data_spa
            hanging_dof > n_spa_new && continue
            hanging_gid   = ip2gip_spa[hanging_dof]
            owner_hanging = get(gip2owner_spa_gid, hanging_gid, rank)
            for (parent_gid, weight) in constraints
                parent_gid in _ext_gid_set_prol || continue
                push!(_local_rev, Float64(parent_gid), Float64(hanging_gid),
                      Float64(owner_hanging), weight)
            end
        end
        _n_rev_loc = Int32(length(_local_rev) ÷ 4)
        _n_rev_all = MPI.Allgather([_n_rev_loc], comm)
        _all_rev   = MPI.Allgatherv(_local_rev, Int32.(_n_rev_all .* 4), comm)

        spa_reverse_ghost_map = Dict{Int, Vector{Tuple{Int,Int,Float64}}}()
        for k = 1:length(_all_rev) ÷ 4
            s = 4*(k-1)
            parent_gid    = Int(round(_all_rev[s+1]))
            hanging_gid   = Int(round(_all_rev[s+2]))
            owner_hanging = Int(round(_all_rev[s+3]))
            weight        = _all_rev[s+4]
            local_parent  = get(gip_to_local_spa, parent_gid, 0)
            local_parent == 0 && continue
            get(gip2owner_spa_gid, parent_gid, -1) != rank && continue
            push!(get!(spa_reverse_ghost_map, local_parent, Tuple{Int,Int,Float64}[]),
                  (hanging_gid, owner_hanging, weight))
        end

        solution_contribs_spa = compute_solution_prolongation_contributions(
            spa_reverse_ghost_map, solution, ip2gip_spa, n_spa_new, rank)

        solution_local_spa = P_spatial_ext * solution

        received_spa_contribs = exchange_hanging_effects_vector(
            solution_contribs_spa, rank, comm)

        solution_new_spa = add_solution_prolongation_contributions(
            @view(solution_local_spa[1:n_spa_new]), received_spa_contribs,
            ip2gip_spa, n_spa_new, rank, gip_to_local_spa)
        solution[1:n_spa_new] .= solution_new_spa

        @info "[$rank] ✓ Spatial prolongation applied"
    end
    # ── Angular integration and error computation ─────────────────────────────
    @rankinfo rank "Integrating solution in angle..."
    int_sol      = zeros(TFloat, npoin, 1)
    int_ref      = zeros(TFloat, npoin, 1)
    int_sol_accum = zeros(npoin)
    int_ref_accum = zeros(npoin)
    L2_err = 0.0
    L2_ref = 0.0

    # Compute node sharing divisors for correct accumulation at shared nodes.
    # Local count: how many local elements reference each spatial node.
    node_div_local = zeros(Int, npoin)
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j, k]
            node_div_local[ip] += 1
        end
    end

    # For spatial AMR, parent nodes at the coarse-fine interface appear in local
    # elements on MULTIPLE ranks (fine-element owner and coarse-element owner).
    # We need the GLOBAL total of contributions so that all ranks' shares sum to
    # exactly one angular integral per spatial node. Use spatial GIDs (mesh.ip2gip)
    # which are consistent with the g_int_sol AllReduce below.
    #==# 
    #if has_spatial_hanging_nodes
        g_count = zeros(Int, mesh.gnpoin)
        for ip = 1:npoin
            g_count[mesh.ip2gip[ip]] += node_div_local[ip]
        end
        MPI.Allreduce!(g_count, MPI.SUM, comm)
        node_div_global = zeros(Int, npoin)
        for ip = 1:npoin
            node_div_global[ip] = g_count[mesh.ip2gip[ip]]
        end
        node_div_global
    #else
    #    node_div_local
    #end
    node_div = node_div_global

        if inputs[:adaptive_extra_meshes]
            for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[iel, i, j, k]
                x   = mesh.x[ip]; y = mesh.y[ip]; z = mesh.z[ip]
                for e_ext = 1:extra_meshes_extra_nelems[iel]
                    for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                        for iϕ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                            ip_ext = extra_meshes_connijk[iel][e_ext, iθ, iϕ]
                            θ      = extra_meshes_coords[iel][1, ip_ext]
                            ϕ      = extra_meshes_coords[iel][2, ip_ext]
                            ip_g   = connijk_spa[iel][i, j, k, e_ext, iθ, iϕ]
                            weight = extra_meshes_extra_Je[iel][e_ext, iθ, iϕ] *
                                     extra_mesh[iel].ωθ[iθ] * extra_mesh[iel].ωθ[iϕ]
                            int_sol_accum[ip] += solution_new[ip_g] * weight / node_div[ip]
                            if (abs(x -1.5) < 1e-3 && abs(y - 4/3) < 1e-3 && abs(z - 0.0) < 1e-2)
                                @info "1", rank, solution_new[ip_g], int_sol_accum[ip], e_ext, solution[ip_g], solution_local[ip_g]
                            end
                            if (inputs[:lmanufactured_solution])
                                int_ref_accum[ip] += ref[ip_g] * weight / node_div[ip]
                                L2_ref += ref[ip_g]^2 * weight * ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                                L2_err += (ref[ip_g] - solution_new[ip_g])^2 *
                                      weight * ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                            end
                        end
                    end
                end
            end
        else
            for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[iel, i, j, k]
                for e_ext = 1:extra_mesh.extra_nelem
                    for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                        for iϕ = 1:extra_mesh.extra_nop[e_ext]+1
                            ip_ext = extra_mesh.extra_connijk[e_ext, iθ, iϕ]
                            θ      = extra_mesh.extra_coords[1, ip_ext]
                            ϕ      = extra_mesh.extra_coords[2, ip_ext]
                            ip_g   = has_spatial_hanging_nodes ?
                                     connijk_spa_uniform[iel][i, j, k, e_ext, iθ, iϕ] :
                                     (ip-1) * extra_mesh.extra_npoin + ip_ext
                            int_sol_accum[ip] += solution[ip_g] *
                                extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ] *
                                extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ] / node_div[ip]
                            
                            if (inputs[:lmanufactured_solution])
                                int_ref_accum[ip] += ref[ip_g] *
                                extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ] *
                                extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ] / node_div[ip]
                                L2_ref += ref[ip_g]^2 *
                                extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ] *
                                extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ] *
                                ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                                L2_err += (ref[ip_g] - solution[ip_g])^2 *
                                extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ] *
                                extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ] *
                                ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                            end
                        end
                    end
                end
            end
        end

        for iel = 1:nelem
            for i = 1:ngl, j = 1:ngl, k = 1:ngl
                ip  = mesh.connijk[iel, i, j, k]
                x   = mesh.x[ip]; y = mesh.y[ip]; z = mesh.z[ip]
                if (abs(x -1.5) < 1e-3 && abs(y - 4/3) < 1e-3 && abs(z - 0.0) < 1e-2)
                    @info "1", mesh.ip2gip[ip], ip, int_sol_accum[ip], node_div[ip], int_sol_accum[ip] * node_div[ip]
                end
                if (abs(x -0.0) < 1e-3 && abs(y - 4/3) < 1e-2 && abs(z - 0.0) < 1e-2)
                    @info "2", mesh.ip2gip[ip], ip, int_sol_accum[ip], node_div[ip], int_sol_accum[ip] * node_div[ip]
                end
            end
        end
        gnpoin_spa  = mesh.gnpoin
        g_int_sol   = zeros(Float64, gnpoin_spa)
        g_int_ref   = zeros(Float64, gnpoin_spa)

        for ip = 1:npoin
            if int_sol_accum[ip] != 0.0 || int_ref_accum[ip] != 0.0
                gip = mesh.ip2gip[ip]
                g_int_sol[gip] += int_sol_accum[ip]
                g_int_ref[gip] += int_ref_accum[ip]
            end
        end
        MPI.Allreduce!(g_int_sol, MPI.SUM, comm)
        MPI.Allreduce!(g_int_ref, MPI.SUM, comm)
        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            int_sol[ip] = g_int_sol[gip]
            int_ref[ip] = g_int_ref[gip]
            x   = mesh.x[ip]; y = mesh.y[ip]; z = mesh.z[ip]
            if (abs(x -1.5) < 1e-3 && abs(y - 4/3) < 1e-3 && abs(z - 0.0) < 1e-2)
                @info "1", gip, ip, int_sol[ip], g_int_sol[gip]
            end
            if (abs(x -0.0) < 1e-3 && abs(y - 4/3) < 1e-2 && abs(z - 0.0) < 1e-2)
                @info "2", gip, ip, int_sol[ip], g_int_sol[gip]
            end
        end
        if (inputs[:lmanufactured_solution])
            L2_ref_g = MPI.Allreduce(L2_ref, MPI.SUM, comm)
            L2_err_g = MPI.Allreduce(L2_err, MPI.SUM, comm)
        end

    if (inputs[:lmanufactured_solution])
        @info "local errors", rank, sqrt(L2_err), sqrt(L2_ref), sqrt(L2_err / L2_ref) 
        @rankinfo rank @sprintf("L2 error: ‖e‖ = %.6e  ‖u‖ = %.6e  relative = %.6e",
        sqrt(L2_err_g), sqrt(L2_ref_g), sqrt(L2_err_g / L2_ref_g))
    end

    if (inputs[:RT_radiative_heating]) && (inputs[:RT_longwave] || inputs[:RT_shortwave])
        if (inputs[:adaptive_extra_meshes])
            if (inputs[:RT_shortwave])
                Q, dTdt, F_net, G = compute_rt_radiative_heating(
                solution_new, atmos_data, mesh,
                connijk_spa, extra_meshes_connijk, extra_meshes_coords,
                extra_meshes_extra_Je, extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_extra_npoins,
                nelem, ngl, κ, σ, extra_mesh[1].ωθ, extra_mesh[1].ωθ, node_div, true, inputs[:RT_longwave];
                G_dir = G_dir,
                Q_dir = Q_dir,
                sw_μ₀ = sw.μ₀,
                sw_φ₀ = sw.φ₀)
            else
                Q, dTdt, F_net, G = compute_rt_radiative_heating(
                solution_new, atmos_data, mesh,
                connijk_spa, extra_meshes_connijk, extra_meshes_coords,
                extra_meshes_extra_Je, extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_extra_npoins,
                nelem, ngl, κ, σ, extra_mesh[1].ωθ, extra_mesh[1].ωθ, node_div, true, inputs[:RT_longwave])
            end

        else

            if (inputs[:RT_shortwave])
                Q, dTdt, F_net, G = compute_rt_radiative_heating(
                solution, atmos_data, mesh,
                connijk_spa, extra_mesh.extra_connijk, extra_mesh.extra_coords,
                extra_mesh.extra_metrics.Je, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_npoin,
                nelem, ngl, κ, σ, extra_mesh.ωθ, extra_mesh.ωθ, node_div, false, inputs[:RT_longwave];
                G_dir = G_dir,
                Q_dir = Q_dir,
                sw_μ₀ = sw.μ₀,
                sw_φ₀ = sw.φ₀)
            else
                Q, dTdt, F_net, G = compute_rt_radiative_heating(
                solution, atmos_data, mesh,
                connijk_spa, extra_mesh.extra_connijk, extra_mesh.extra_coords,
                extra_mesh.extra_metrics.Je, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_npoin,
                nelem, ngl, κ, σ, extra_mesh.ωθ, extra_mesh.ωθ, node_div, false, inputs[:RT_longwave])
        
            end
        end
    end


    if (inputs[:RT_atmos_coupling])
        return Q, dTdt, solution
    else
        # ── VTK output ────────────────────────────────────────────────────────────
        if (inputs[:RT_radiative_heating]) && (inputs[:RT_longwave] || inputs[:RT_shortwave])
            out_vectors = zeros(mesh.npoin,12)
            if (inputs[:adaptive_extra_meshes])
                out_vectors = zeros(mesh.npoin,13)
            end
            for i=1:mesh.npoin
                out_vectors[i,1] = int_sol[i]
                out_vectors[i,2] = Q[i]
                out_vectors[i,3] = dTdt[i]
                out_vectors[i,4] = F_net[i]
                out_vectors[i,5] = G[i]
                out_vectors[i,6] = atmos_data.t_lev[i]
                out_vectors[i,7] = atmos_data.q_liq[i]
                out_vectors[i,8] = atmos_data.q_ice[i]
                out_vectors[i,9] = κ[i]
                out_vectors[i,10] = σ[i]
                if (inputs[:RT_shortwave])
                    out_vectors[i,11] = F_dir[i]
                    out_vectors[i,12] = τ_nodes[i]
                end
            end
            if (inputs[:adaptive_extra_meshes])
                for iel = 1:nelem
                    for i = 1:ngl, j = 1:ngl, k = 1:ngl
                        ip = mesh.connijk[iel,i,j,k]
                        out_vectors[ip,13] = extra_meshes_extra_nelems[iel]
                    end
                end
            end
            @rankinfo rank "Writing output"
            title = @sprintf "Solution-Radiation"
            if (inputs[:outformat] == VTK())
                if (inputs[:adaptive_extra_meshes])
                    write_vtk(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes", "nelem_ang"], ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes", "nelem_ang"]; iout=1, nvar=13)
                    return
                else
                    write_vtk(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes"], ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes"]; iout=1, nvar=12)
                    return
                end
            elseif (inputs[:outformat] == NETCDF())
                if (inputs[:adaptive_extra_meshes])
                    write_NetCDF(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes", "nelem_ang"], ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes", "nelem_ang"]; iout=1, nvar=13)
                    return
                else
                    write_NetCDF(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes"], ["Ang_int","Q","dTdt","F_net","G", "T", "q_liq", "q_ice", "kappa", "sigma", "F_dir", "τ_nodes"]; iout=1, nvar=12)
                    return
                end
            end
        else
            
            out_vectors = zeros(mesh.npoin,1)
            if (inputs[:adaptive_extra_meshes])
                out_vectors = zeros(mesh.npoin,2)
            end
            for i=1:mesh.npoin
                out_vectors[i,1] = int_sol[i]
            end
            if (inputs[:adaptive_extra_meshes])
                for iel = 1:nelem
                    for i = 1:ngl, j = 1:ngl, k = 1:ngl
                        ip = mesh.connijk[iel,i,j,k]
                        out_vectors[ip,2] = extra_meshes_extra_nelems[iel]
                    end
                end
            end
            @rankinfo rank "Writing output"
            title = @sprintf "Solution-Radiation"
            if (inputs[:outformat] == VTK())
                if (inputs[:adaptive_extra_meshes])
                    write_vtk(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int","nelem_ang"], ["Ang_int","nelem_ang"]; iout=1, nvar=2)
                    return
                else
                    write_vtk(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int"], ["Ang_int"]; iout=1, nvar=1)
                    return
                end
            elseif (inputs[:outformat] == NETCDF())
                if (inputs[:adaptive_extra_meshes])
                    write_NetCDF(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int", "nelem_ang"], ["Ang_int","nelem_ang"]; iout=1, nvar=2)
                    return
                else
                    write_NetCDF(SD, mesh, out_vectors, out_vectors, nothing, nothing, nothing,
                        0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                        ["Ang_int"], ["Ang_int"]; iout=1, nvar=1)
                    return
                end
            end
        end
    end

end

# ──────────────────────────────────────────────────────────────────────────────
"""
    sparse_lhs_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, z,
        ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
        npoin_ang_total, nelem, ngl, nelem_ang,
        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        npoin_ang, rad_HG_g, rad_data, κ_data, σ_data)

Assemble the sparse LHS matrix for the **non-adaptive** 3D×2D radiative transfer
problem using a continuous Galerkin spectral element formulation.

The assembled matrix represents the discretisation of:

    (Ω·∇ + κ - σ ∫Φ dΩ') I

on a tensor-product spatial × angular mesh with uniform angular elements.

# Arguments (key ones)

| Argument | Description |
|----------|-------------|
| `ω` | 1D spatial GL quadrature weights |
| `Je` | Spatial Jacobian `[iel,i,j,k]` |
| `connijk` | Spatial connectivity `[iel,i,j,k]` → local node index |
| `ωθ, ωϕ` | Angular quadrature weights in θ and ϕ |
| `ψ, dψ` | Spatial basis and derivative matrices |
| `connijk_ang` | Angular connectivity `[e,iθ,jθ]` → angular node index |
| `Je_ang` | Angular Jacobian `[e,iθ,jθ]` |
| `coords_ang` | Angular coordinates `[1,:] = θ`, `[2,:] = ϕ` |
| `nop_ang` | Polynomial order per angular element |
| `rad_HG_g` | Henyey–Greenstein g parameter |
| `rad_data` | If `true`, read κ and σ from arrays; else call user functions |

# Returns

Sparse matrix of size `npoin_ang_total × npoin_ang_total`.

# Notes

Scattering integrals `∫ Φ(Ω,Ω') dΩ'` are precomputed once per angular node
before the element loop to avoid O(n²) redundant evaluations.

The matrix is built using a pre-allocated COO triplet with dynamic resizing
via `_grow_if_needed!` to avoid `push!` overhead.
"""
function sparse_lhs_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
                                    npoin_ang_total, nelem, ngl, nelem_ang,
                                    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                                    npoin_ang, g_eff, rad_data, κ_data, σ_data, lshortwave)

    intΦ = zeros(Float64, npoin_ang)
    g_table    = [0.0, 0.85]
    intΦ_table = zeros(Float64, npoin_ang, length(g_table))
    
    if (lshortwave)
        for (ig, g_val) in enumerate(g_table)
            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ip_ext = connijk_ang[e_ext, iθ, jθ]
                    θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                    val = 0.0
                    for e_ext_s = 1:nelem_ang
                        for nθ = 1:nop_ang[e_ext_s]+1, mθ = 1:nop_ang[e_ext_s]+1
                            ipθ = connijk_ang[e_ext_s, mθ, nθ]
                            Φ   = user_scattering_functions(θ, coords_ang[1,ipθ],
                                                            ϕ, coords_ang[2,ipθ],
                                                            g_val)
                            val += ωθ[mθ] * ωϕ[nθ] * Je_ang[e_ext_s, mθ, nθ] * Φ
                        end
                    end
                    intΦ_table[ip_ext, ig] = val
                end
            end
        end
    else
        # Precompute scattering integrals once per angular node
        for e_ext = 1:nelem_ang
            for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                ip_ext = connijk_ang[e_ext, iθ, jθ]
                θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                val = 0.0
                for e_ext_s = 1:nelem_ang
                    for nθ = 1:nop_ang[e_ext_s]+1, mθ = 1:nop_ang[e_ext_s]+1
                        ipθ = connijk_ang[e_ext_s, mθ, nθ]
                        Φ   = user_scattering_functions(θ, coords_ang[1,ipθ], ϕ, coords_ang[2,ipθ], g_eff)
                        val += ωθ[mθ] * ωϕ[nθ] * Je_ang[e_ext_s, mθ, nθ] * Φ
                    end
                end
                intΦ[ip_ext] = val
            end
        end
    end

    ngl_stencil = 3 * (ngl - 1) + 1
    estimated   = round(Int, npoin_ang_total * ngl_stencil * 1.2)
    I_vec = Vector{Int}(undef, estimated)
    J_vec = Vector{Int}(undef, estimated)
    V_vec = Vector{Float64}(undef, estimated)
    pos   = 1

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip       = connijk[iel,i,j,k]
            ωJac     = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            dξdx_ij  = dξdx[iel,i,j,k]; dξdy_ij = dξdy[iel,i,j,k]; dξdz_ij = dξdz[iel,i,j,k]
            dηdx_ij  = dηdx[iel,i,j,k]; dηdy_ij = dηdy[iel,i,j,k]; dηdz_ij = dηdz[iel,i,j,k]
            dζdx_ij  = dζdx[iel,i,j,k]; dζdy_ij = dζdy[iel,i,j,k]; dζdz_ij = dζdz[iel,i,j,k]
            κ = rad_data ? κ_data[ip] : user_extinction(x[ip], y[ip], z[ip])
            σ = rad_data ? σ_data[ip] : user_scattering_coef(x[ip], y[ip], z[ip])

            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ip_ext     = connijk_ang[e_ext, iθ, jθ]
                    ωJac_rad   = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext, iθ, jθ]
                    ωJac_full  = ωJac * ωJac_rad
                    θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                    idx_ip = (ip-1)*npoin_ang + ip_ext
                    intϕ_eff = 0.0
                    if (lshortwave)
                        α = g_eff[ip] / 0.85
                        intΦ_eff = (1 - α) * intΦ_table[ip_ext, 1] + α * intΦ_table[ip_ext, 2]
                    else
                        intΦ_eff = intΦ[ip_ext]
                    end
                    _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                    I_vec[pos] = idx_ip; J_vec[pos] = idx_ip
                    if (rad_data)
                        V_vec[pos] = (κ + σ*(1 - intΦ_eff)) * ωJac_full
                    else
                        V_vec[pos] = (κ - σ * intΦ_eff) * ωJac_full
                    end
                    pos += 1

                    prop_i = (dξdx_ij*Ωx + dξdy_ij*Ωy + dξdz_ij*Ωz) * ωJac_full
                    for m = 1:ngl
                        abs(dψ[m,i]) < eps(Float64) && continue
                        jp    = connijk[iel,m,j,k]
                        idx_m = (jp-1)*npoin_ang + ip_ext
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_m; V_vec[pos] = dψ[m,i]*prop_i
                        pos += 1
                    end

                    prop_j = (dηdx_ij*Ωx + dηdy_ij*Ωy + dηdz_ij*Ωz) * ωJac_full
                    for n = 1:ngl
                        abs(dψ[n,j]) < eps(Float64) && continue
                        jp    = connijk[iel,i,n,k]
                        idx_n = (jp-1)*npoin_ang + ip_ext
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_n; V_vec[pos] = dψ[n,j]*prop_j
                        pos += 1
                    end

                    prop_k = (dζdx_ij*Ωx + dζdy_ij*Ωy + dζdz_ij*Ωz) * ωJac_full
                    for o = 1:ngl
                        abs(dψ[o,k]) < eps(Float64) && continue
                        jp    = connijk[iel,i,j,o]
                        idx_o = (jp-1)*npoin_ang + ip_ext
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_o; V_vec[pos] = dψ[o,k]*prop_k
                        pos += 1
                    end
                end
            end
        end
    end
    return sparse(view(I_vec,1:pos-1), view(J_vec,1:pos-1), view(V_vec,1:pos-1))
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    sparse_lhs_assembly_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ, x, y, z,
        ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
        npoin_ang_total, nelem, ngl, nelem_ang,
        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        npoin_ang, rad_HG_g, connijk_spa)

Assemble the sparse LHS matrix for the **adaptive** 3D×2D radiative transfer
problem.  Identical physics to `sparse_lhs_assembly_3Dby2D` but uses the
spatially-varying angular mesh described by `connijk_spa` and
element-local `connijk_ang[iel]`, `coords_ang[iel]`, etc.

Scattering integrals are precomputed per spatial element (since the angular
mesh may differ between elements after adaptive refinement).

# Returns

Sparse matrix of size `n_spa × n_spa` where `n_spa` is the total number of
unique spatial-angular DOFs.
"""
function sparse_lhs_assembly_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
                                    npoin_ang_total, nelem, ngl, nelem_ang,
                                    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                                    npoin_ang, g_eff, connijk_spa, rad_data, κ_data, σ_data, lshortwave)

    ngl_stencil = 3 * (ngl - 1) + 1
    estimated   = round(Int, npoin_ang_total * ngl_stencil * 1.2)
    I_vec = Vector{Int}(undef, estimated)
    J_vec = Vector{Int}(undef, estimated)
    V_vec = Vector{Float64}(undef, estimated)
    pos   = 1

    # Precompute per-element scattering integrals
    intΦ_per_iel = Vector{Vector{Float64}}(undef, nelem)
    g_table    = [0.0, 0.85]
    intΦ_table_per_iel = Vector{Matrix{Float64}}(undef, nelem)
    
    if (lshortwave)
        for iel = 1:nelem
            n_ang  = npoin_ang[iel]
            intΦ   = zeros(Float64, n_ang, length(g_table))
            for (ig, g_val) in enumerate(g_table)
                for e_ext = 1:nelem_ang[iel]
                    for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                        ip_ext = connijk_ang[iel][e_ext, iθ, jθ]
                        θ = coords_ang[iel][1, ip_ext]; ϕ = coords_ang[iel][2, ip_ext]
                        val = 0.0
                        for e_ext_s = 1:nelem_ang[iel]
                            for nθ = 1:nop_ang[iel][e_ext_s]+1, mθ = 1:nop_ang[iel][e_ext_s]+1
                                ipθ = connijk_ang[iel][e_ext_s, mθ, nθ]
                                Φ   = user_scattering_functions(θ, coords_ang[iel][1,ipθ],
                                                                ϕ, coords_ang[iel][2,ipθ],
                                                                g_val)

                                val += ωθ[mθ] * ωϕ[nθ] * Je_ang[iel][e_ext_s, mθ, nθ] * Φ
                            end
                        end
                        intΦ[ip_ext, ig] = val
                    end
                end
            end
            intΦ_table_per_iel[iel] = intΦ
        end
    else
        for iel = 1:nelem
            n_ang  = npoin_ang[iel]
            intΦ   = zeros(Float64, n_ang)
            for e_ext = 1:nelem_ang[iel]
                for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                    ip_ext = connijk_ang[iel][e_ext, iθ, jθ]
                    θ = coords_ang[iel][1, ip_ext]; ϕ = coords_ang[iel][2, ip_ext]
                    val = 0.0
                    for e_ext_s = 1:nelem_ang[iel]
                        for nθ = 1:nop_ang[iel][e_ext_s]+1, mθ = 1:nop_ang[iel][e_ext_s]+1
                            ipθ = connijk_ang[iel][e_ext_s, mθ, nθ]
                            Φ   = user_scattering_functions(θ, coords_ang[iel][1,ipθ],
                                                            ϕ, coords_ang[iel][2,ipθ], g_eff)
                            val += ωθ[mθ] * ωϕ[nθ] * Je_ang[iel][e_ext_s, mθ, nθ] * Φ
                        end
                    end
                    intΦ[ip_ext] = val
                end
            end
            intΦ_per_iel[iel] = intΦ
        end
    end

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip       = connijk[iel,i,j,k]
            ωJac     = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            dξdx_ij  = dξdx[iel,i,j,k]; dξdy_ij = dξdy[iel,i,j,k]; dξdz_ij = dξdz[iel,i,j,k]
            dηdx_ij  = dηdx[iel,i,j,k]; dηdy_ij = dηdy[iel,i,j,k]; dηdz_ij = dηdz[iel,i,j,k]
            dζdx_ij  = dζdx[iel,i,j,k]; dζdy_ij = dζdy[iel,i,j,k]; dζdz_ij = dζdz[iel,i,j,k]
            κ = rad_data ? κ_data[ip] : user_extinction(x[ip], y[ip], z[ip])
            σ = rad_data ? σ_data[ip] : user_scattering_coef(x[ip], y[ip], z[ip])

            for e_ext = 1:nelem_ang[iel]
                for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                    ip_ext    = connijk_ang[iel][e_ext, iθ, jθ]
                    ωJac_rad  = ωθ[iθ]*ωϕ[jθ]*Je_ang[iel][e_ext, iθ, jθ]
                    ωJac_full = ωJac * ωJac_rad
                    θ = coords_ang[iel][1, ip_ext]; ϕ = coords_ang[iel][2, ip_ext]
                    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                    idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                    if (lshortwave)
                        α = g_eff[ip] / 0.85
                        intΦ_eff = (1 - α) * intΦ_table_per_iel[iel][ip_ext, 1] + α * intΦ_table_per_iel[iel][ip_ext, 2]
                    else
                        intΦ_eff = intΦ_per_iel[iel][ip_ext]
                    end
                    _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                    I_vec[pos] = idx_ip; J_vec[pos] = idx_ip
                    if (rad_data)
                        V_vec[pos] = (κ + σ * (1 - intΦ_eff)) * ωJac_full
                    else
                        V_vec[pos] = (κ - σ*intΦ_eff) * ωJac_full
                    end
                    pos += 1

                    prop_i = (dξdx_ij*Ωx + dξdy_ij*Ωy + dξdz_ij*Ωz) * ωJac_full
                    for m = 1:ngl
                        abs(dψ[m,i]) < eps(Float64) && continue
                        idx_m = connijk_spa[iel][m,j,k,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_m; V_vec[pos] = dψ[m,i]*prop_i
                        pos += 1
                    end

                    prop_j = (dηdx_ij*Ωx + dηdy_ij*Ωy + dηdz_ij*Ωz) * ωJac_full
                    for n = 1:ngl
                        abs(dψ[n,j]) < eps(Float64) && continue
                        idx_n = connijk_spa[iel][i,n,k,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_n; V_vec[pos] = dψ[n,j]*prop_j
                        pos += 1
                    end

                    prop_k = (dζdx_ij*Ωx + dζdy_ij*Ωy + dζdz_ij*Ωz) * ωJac_full
                    for o = 1:ngl
                        abs(dψ[o,k]) < eps(Float64) && continue
                        idx_o = connijk_spa[iel][i,j,o,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_o; V_vec[pos] = dψ[o,k]*prop_k
                        pos += 1
                    end
                end
            end
        end
    end
    return sparse(view(I_vec,1:pos-1), view(J_vec,1:pos-1), view(V_vec,1:pos-1))
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    _grow_if_needed!(I, J, V, pos)

Double the capacity of COO triplet arrays `I`, `J`, `V` when the write cursor
`pos` is about to overflow.  Called inline inside assembly loops to avoid
pre-allocating a worst-case buffer.
"""
@inline function _grow_if_needed!(I, J, V, pos)
    if pos > length(I)
        n = length(I)
        resize!(I, 2n); resize!(J, 2n); resize!(V, 2n)
    end
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    assemble_mass_diagonal_3Dby2D(ω, Je, connijk, ωθ, ωϕ, connijk_ang,
        Je_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang)

Assemble the lumped (diagonal) mass matrix for the non-adaptive problem using
inexact (Gauss–Lobatto) integration.  Returns a vector `Md` of length
`npoin_ang_total` with the diagonal entries.
"""
function assemble_mass_diagonal_3Dby2D(ω, Je, connijk, ωθ, ωϕ, connijk_ang,
        Je_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang)

    Md = zeros(Float64, npoin_ang_total)
    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            ip   = connijk[iel,i,j,k]
            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                    ip_ext   = connijk_ang[e_ext,iθ,jθ]
                    idx      = (ip-1)*npoin_ang + ip_ext
                    Md[idx] += ωJac * ωJac_rad
                end
            end
        end
    end
    return Md
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    assemble_mass_diagonal_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ,
        Je_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, connijk_spa)

Adaptive variant of `assemble_mass_diagonal_3Dby2D`.  Uses `connijk_spa` for
the combined spatial-angular connectivity and element-local angular metrics.

Returns a vector `Md` of length `npoin_ang_total = n_spa`.
"""
function assemble_mass_diagonal_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ,
        Je_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, connijk_spa)

    Md = zeros(Float64, npoin_ang_total)
    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            for e_ext = 1:nelem_ang[iel]
                for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                    ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[iel][e_ext,iθ,jθ]
                    idx      = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                    Md[idx] += ωJac * ωJac_rad
                end
            end
        end
    end
    return Md
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    sparse_lhs_assembly_3Dby2D_spatial_amr(ω, Je, connijk, ωθ, ωϕ, x, y, z,
        ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
        n_spa_new, nelem, ngl, nelem_ang,
        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        npoin_ang, rad_HG_g, connijk_spa, rad_data, κ_data, σ_data, lshortwave)

Assemble the sparse LHS matrix for spatial AMR with **uniform** angular mesh.
Uses coordinate-deduplicated `connijk_spa[iel][i,j,k,e_ext,iθ,jθ]` for DOF indexing.
The angular mesh is uniform across spatial elements (same `connijk_ang`, `Je_ang`, etc.).
Output size: `n_spa_new × n_spa_new`.
"""
function sparse_lhs_assembly_3Dby2D_spatial_amr(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang,
                                    n_spa_new, nelem, ngl, nelem_ang,
                                    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                                    npoin_ang, g_eff, connijk_spa, rad_data, κ_data, σ_data, lshortwave)

    intΦ = zeros(Float64, npoin_ang)
    g_table    = [0.0, 0.85]
    intΦ_table = zeros(Float64, npoin_ang, length(g_table))

    if lshortwave
        for (ig, g_val) in enumerate(g_table)
            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ip_ext = connijk_ang[e_ext, iθ, jθ]
                    θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                    val = 0.0
                    for e_ext_s = 1:nelem_ang
                        for nθ = 1:nop_ang[e_ext_s]+1, mθ = 1:nop_ang[e_ext_s]+1
                            ipθ = connijk_ang[e_ext_s, mθ, nθ]
                            Φ   = user_scattering_functions(θ, coords_ang[1,ipθ],
                                                            ϕ, coords_ang[2,ipθ], g_val)
                            val += ωθ[mθ] * ωϕ[nθ] * Je_ang[e_ext_s, mθ, nθ] * Φ
                        end
                    end
                    intΦ_table[ip_ext, ig] = val
                end
            end
        end
    else
        for e_ext = 1:nelem_ang
            for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                ip_ext = connijk_ang[e_ext, iθ, jθ]
                θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                val = 0.0
                for e_ext_s = 1:nelem_ang
                    for nθ = 1:nop_ang[e_ext_s]+1, mθ = 1:nop_ang[e_ext_s]+1
                        ipθ = connijk_ang[e_ext_s, mθ, nθ]
                        Φ   = user_scattering_functions(θ, coords_ang[1,ipθ], ϕ, coords_ang[2,ipθ], g_eff)
                        val += ωθ[mθ] * ωϕ[nθ] * Je_ang[e_ext_s, mθ, nθ] * Φ
                    end
                end
                intΦ[ip_ext] = val
            end
        end
    end

    ngl_stencil = 3 * (ngl - 1) + 1
    estimated   = round(Int, n_spa_new * ngl_stencil * 1.2)
    I_vec = Vector{Int}(undef, estimated)
    J_vec = Vector{Int}(undef, estimated)
    V_vec = Vector{Float64}(undef, estimated)
    pos   = 1

    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip       = connijk[iel,i,j,k]
            ωJac     = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            dξdx_ij  = dξdx[iel,i,j,k]; dξdy_ij = dξdy[iel,i,j,k]; dξdz_ij = dξdz[iel,i,j,k]
            dηdx_ij  = dηdx[iel,i,j,k]; dηdy_ij = dηdy[iel,i,j,k]; dηdz_ij = dηdz[iel,i,j,k]
            dζdx_ij  = dζdx[iel,i,j,k]; dζdy_ij = dζdy[iel,i,j,k]; dζdz_ij = dζdz[iel,i,j,k]
            κ = rad_data ? κ_data[ip] : user_extinction(x[ip], y[ip], z[ip])
            σ = rad_data ? σ_data[ip] : user_scattering_coef(x[ip], y[ip], z[ip])

            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ip_ext     = connijk_ang[e_ext, iθ, jθ]
                    ωJac_rad   = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext, iθ, jθ]
                    ωJac_full  = ωJac * ωJac_rad
                    θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
                    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                    idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                    if lshortwave
                        α = g_eff[ip] / 0.85
                        intΦ_eff = (1 - α) * intΦ_table[ip_ext, 1] + α * intΦ_table[ip_ext, 2]
                    else
                        intΦ_eff = intΦ[ip_ext]
                    end
                    _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                    I_vec[pos] = idx_ip; J_vec[pos] = idx_ip
                    if rad_data
                        V_vec[pos] = (κ + σ*(1 - intΦ_eff)) * ωJac_full
                    else
                        V_vec[pos] = (κ - σ * intΦ_eff) * ωJac_full
                    end
                    pos += 1

                    prop_i = (dξdx_ij*Ωx + dξdy_ij*Ωy + dξdz_ij*Ωz) * ωJac_full
                    for m = 1:ngl
                        abs(dψ[m,i]) < eps(Float64) && continue
                        idx_m = connijk_spa[iel][m,j,k,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_m; V_vec[pos] = dψ[m,i]*prop_i
                        pos += 1
                    end

                    prop_j = (dηdx_ij*Ωx + dηdy_ij*Ωy + dηdz_ij*Ωz) * ωJac_full
                    for n = 1:ngl
                        abs(dψ[n,j]) < eps(Float64) && continue
                        idx_n = connijk_spa[iel][i,n,k,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_n; V_vec[pos] = dψ[n,j]*prop_j
                        pos += 1
                    end

                    prop_k = (dζdx_ij*Ωx + dζdy_ij*Ωy + dζdz_ij*Ωz) * ωJac_full
                    for o = 1:ngl
                        abs(dψ[o,k]) < eps(Float64) && continue
                        idx_o = connijk_spa[iel][i,j,o,e_ext,iθ,jθ]
                        _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                        I_vec[pos] = idx_ip; J_vec[pos] = idx_o; V_vec[pos] = dψ[o,k]*prop_k
                        pos += 1
                    end
                end
            end
        end
    end
    return sparse(view(I_vec,1:pos-1), view(J_vec,1:pos-1), view(V_vec,1:pos-1), n_spa_new, n_spa_new)
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    assemble_mass_diagonal_3Dby2D_spatial_amr(ω, Je, connijk, ωθ, ωϕ, connijk_ang,
        Je_ang, nop_ang, n_spa_new, nelem, ngl, nelem_ang, connijk_spa)

Assemble the lumped mass matrix for spatial AMR with **uniform** angular mesh.
Uses coordinate-deduplicated `connijk_spa[iel][i,j,k,e_ext,iθ,jθ]` for DOF indexing.
Returns a vector of length `n_spa_new`.
"""
function assemble_mass_diagonal_3Dby2D_spatial_amr(ω, Je, connijk, ωθ, ωϕ, connijk_ang,
        Je_ang, nop_ang, n_spa_new, nelem, ngl, nelem_ang, connijk_spa)

    Md = zeros(Float64, n_spa_new)
    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
            for e_ext = 1:nelem_ang
                for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
                    ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                    idx      = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                    Md[idx] += ωJac * ωJac_rad
                end
            end
        end
    end
    return Md
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    compute_adaptivity_criterion3D_2D(pointwise_interaction, nelem, ngl, connijk,
        connijk_ang, nop_ang, nelem_ang, coords_ang, connijk_spa,
        ψ_ang, dψ_ang, dξdθ, dηdθ, dξdϕ, dηdϕ)

Compute an element-wise h-adaptivity refinement criterion for the angular mesh.

The criterion for each angular element `e_ext` of spatial element `iel` is the
maximum over all quadrature points of `|∇_Ω (pointwise_interaction)|`, where
the gradient is computed in the angular (θ,ϕ) directions using the provided
metric terms.

# Returns

`criterion[iel][e_ext]` — maximum smoothness indicator for each element.
Elements with large values (above the threshold in `adapt_angular_grid_3Dby2D!`)
are marked for refinement.
"""
function compute_adaptivity_criterion3D_2D(pointwise_interaction, nelem, ngl, connijk,
        connijk_ang, nop_ang, nelem_ang, coords_ang, connijk_spa,
        ψ_ang, dψ_ang, dξdθ, dηdθ, dξdϕ, dηdϕ)

    criterion = [Vector{Float64}(undef, nelem_ang[e]) for e = 1:nelem]
    for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            for e_ext = 1:nelem_ang[iel]
                criterion[iel][e_ext] = 0.0
                for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                    ip_ext = connijk_ang[iel][e_ext, iθ, jθ]
                    idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                    gradξ  = 0.0; gradη = 0.0
                    for kθ = 1:nop_ang[iel][e_ext]+1
                        gradξ += pointwise_interaction[idx_ip] * dψ_ang[iθ,kθ] * ψ_ang[jθ,kθ]
                        gradη += pointwise_interaction[idx_ip] * ψ_ang[iθ,kθ]  * dψ_ang[jθ,kθ]
                    end
                    gradθ = gradξ*dξdθ[iel][e_ext,iθ,jθ] + gradη*dηdθ[iel][e_ext,iθ,jθ]
                    gradϕ = gradξ*dξdϕ[iel][e_ext,iθ,jθ] + gradη*dηdϕ[iel][e_ext,iθ,jθ]
                    criterion[iel][e_ext] = max(criterion[iel][e_ext], abs(gradθ + gradϕ))
                end
            end
        end
    end
    return criterion
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    adapt_angular_grid_3Dby2D!(criterion, thresholds, ref_level, nelem, ngl,
        nelem_ang, nop_ang, neighbors, npoin_ang, connijk_ang, coords_ang,
        Je_ang, dξdx_ang, dxdξ_ang, dξdy_ang, dydξ_ang, dηdy_ang, dydη_ang,
        dηdx_ang, dxdη_ang, connijk, x, y, z,
        xmin_grid, ymin_grid, zmin_grid, xmax_grid, ymax_grid, zmax_grid, ψ, dψ)

Perform one pass of h-adaptive refinement on the angular mesh.

Each angular element whose criterion exceeds `thresholds[1]` is split into four
child elements (bisection in both θ and ϕ).  Hanging nodes at the interface
between refined and unrefined elements are identified by the subsequent call to
`adaptive_spatial_angular_numbering_3D_2D!`.

After refinement, ϕ-periodicity is enforced by merging points at ϕ=2π with
their counterparts at ϕ=0.

The `neighbors` array is updated to track which element pairs have non-
conforming angular interfaces, which drives the constraint matrix construction
in `build_restriction_matrices_local_and_ghost`.

# Modifies in-place

`connijk_ang`, `coords_ang`, `Je_ang`, all metric arrays, `nop_ang`,
`nelem_ang`, `npoin_ang`, `ref_level`, `criterion`, `neighbors`.
"""
function adapt_angular_grid_3Dby2D!(criterion, thresholds, ref_level, nelem, ngl, nelem_ang, nop_ang, neighbors, npoin_ang,
        connijk_ang, coords_ang, Je_ang, dξdx_ang, dxdξ_ang, dξdy_ang, dydξ_ang, dηdy_ang, dydη_ang, dηdx_ang, dxdη_ang, connijk,
        x, y, z, xmin_grid, ymin_grid, zmin_grid, xmax_grid, ymax_grid, zmax_grid, ψ, dψ,
        can_refine::Union{BitVector,Nothing}=nothing)

    lgl = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    adapted_ang = zeros(Int, nelem)
    ang_adapted  = zeros(Int, nelem)

    for iel = 1:nelem
        # Skip angular refinement for elements that are spatially non-conforming
        # with at least one neighbor.  Allowing both types of non-conformity at
        # the same interface is not supported.
        can_refine !== nothing && !can_refine[iel] && continue

        e_ext = 1

        while e_ext <= nelem_ang[iel]

            if (abs(criterion[iel][e_ext]) > thresholds[1])
               #@info iel, e_ext
               #@info rank, criterion[iel][e_ext], thresholds[1], e_ext, iel
                adapted_ang[iel] = 1
                ref_level[iel][e_ext] += 1
                ang_adapted[iel] = 1
                nelem_ang[iel] += 3
                npoin_ang[iel] += (4*(nop_ang[iel][e_ext]+1)^2 - 4*(nop_ang[iel][e_ext]+1))

                θmax = coords_ang[iel][1, connijk_ang[iel][e_ext, nop_ang[iel][e_ext]+1, nop_ang[iel][e_ext]+1]]
                θmin = coords_ang[iel][1, connijk_ang[iel][e_ext, 1, 1]]
                ϕmax = coords_ang[iel][2, connijk_ang[iel][e_ext, nop_ang[iel][e_ext]+1, nop_ang[iel][e_ext]+1]]
                ϕmin = coords_ang[iel][2, connijk_ang[iel][e_ext, 1, 1]]
                ϕmax == 0 && (ϕmax = 2π)
                θ12 = (θmax + θmin) / 2
                ϕ12 = (ϕmax + ϕmin) / 2

                connijk_ang_new = zeros(Int, nelem_ang[iel], nop_ang[iel][1]+1, nop_ang[iel][1]+1)
                coords_new      = zeros(Float64, 2, npoin_ang[iel])
                metrics         = allocate_metrics(NSD_2D(), nelem_ang[iel], 0, nop_ang[iel][e_ext], TFloat, CPU())
                nop_ang_new     = zeros(Int, nelem_ang[iel])
                nop_ang_new[1:nelem_ang[iel]-1] .= nop_ang[iel][e_ext]
                nop_ang_new[nelem_ang[iel]]      = nop_ang[iel][1]
                criterion_new   = zeros(Float64, nelem_ang[iel])
                ref_level_new   = zeros(Int,     nelem_ang[iel])

                point_dict = Dict{NTuple{2,Float64}, Int}()
                iter = 1
                function insert_point!(θ, ϕ, target_arr, e, i, j)
                    key = (θ, ϕ)
                    idx = get(point_dict, key, 0)
                    if idx != 0
                        target_arr[e, i, j] = idx
                    else
                        point_dict[key]    = iter
                        target_arr[e,i,j]  = iter
                        coords_new[1,iter] = θ
                        coords_new[2,iter] = ϕ
                        iter += 1
                    end
                end

                # Copy elements before the refined one
                if e_ext > 1
                    for e_ext1 = 1:e_ext-1
                        for i = 1:nop_ang[iel][e_ext1]+1, j = 1:nop_ang[iel][e_ext1]+1
                            θ = coords_ang[iel][1, connijk_ang[iel][e_ext1,i,j]]
                            ϕ = coords_ang[iel][2, connijk_ang[iel][e_ext1,i,j]]
                            metrics.dxdξ[e_ext1,i,j] = dxdξ_ang[iel][e_ext1,i,j]
                            metrics.Je[e_ext1,i,j]   = Je_ang[iel][e_ext1,i,j]
                            metrics.dξdx[e_ext1,i,j] = dξdx_ang[iel][e_ext1,i,j]
                            metrics.dxdη[e_ext1,i,j] = dxdη_ang[iel][e_ext1,i,j]
                            metrics.dηdx[e_ext1,i,j] = dηdx_ang[iel][e_ext1,i,j]
                            metrics.dηdy[e_ext1,i,j] = dηdy_ang[iel][e_ext1,i,j]
                            metrics.dξdy[e_ext1,i,j] = dξdy_ang[iel][e_ext1,i,j]
                            metrics.dydη[e_ext1,i,j] = dydη_ang[iel][e_ext1,i,j]
                            metrics.dydξ[e_ext1,i,j] = dydξ_ang[iel][e_ext1,i,j]
                            criterion_new[e_ext1] = criterion[iel][e_ext1]
                            ref_level_new[e_ext1] = ref_level[iel][e_ext1]
                            insert_point!(θ, ϕ, connijk_ang_new, e_ext1, i, j)
                        end
                    end
                end

                lgl = basis_structs_ξ_ω!(LGL(), nop_ang[iel][1], CPU())

                # Build four child elements
                child_ranges = [
                    (θmin, θ12, ϕmin, ϕ12, e_ext),
                    (θ12, θmax, ϕmin, ϕ12, e_ext+1),
                    (θmin, θ12, ϕ12, ϕmax, e_ext+2),
                    (θ12, θmax, ϕ12, ϕmax, e_ext+3),
                ]
                for (θlo, θhi, ϕlo, ϕhi, e_child) in child_ranges
                    criterion_new[e_child] = 0.0
                    ref_level_new[e_child] = ref_level[iel][e_ext]
                    for i = 1:nop_ang_new[e_child]+1
                        ξi = lgl.ξ[i]
                        θ  = θlo*(1-ξi)*0.5 + θhi*(1+ξi)*0.5
                        for j = 1:nop_ang_new[e_child]+1
                            ξj = lgl.ξ[j]
                            ϕ  = ϕlo*(1-ξj)*0.5 + ϕhi*(1+ξj)*0.5
                            insert_point!(θ, ϕ, connijk_ang_new, e_child, i, j)
                        end
                    end
                end

                # Copy elements after the refined one
                if e_ext < nelem_ang[iel] - 3
                    for e_ext1 = e_ext+4:nelem_ang[iel]
                        src = e_ext1 - 3
                        for i = 1:nop_ang_new[e_ext1]+1, j = 1:nop_ang_new[e_ext1]+1
                            θ = coords_ang[iel][1, connijk_ang[iel][src,i,j]]
                            ϕ = coords_ang[iel][2, connijk_ang[iel][src,i,j]]
                            metrics.dxdξ[e_ext1,i,j] = dxdξ_ang[iel][src,i,j]
                            metrics.Je[e_ext1,i,j]   = Je_ang[iel][src,i,j]
                            metrics.dξdx[e_ext1,i,j] = dξdx_ang[iel][src,i,j]
                            metrics.dxdη[e_ext1,i,j] = dxdη_ang[iel][src,i,j]
                            metrics.dηdx[e_ext1,i,j] = dηdx_ang[iel][src,i,j]
                            metrics.dηdy[e_ext1,i,j] = dηdy_ang[iel][src,i,j]
                            metrics.dξdy[e_ext1,i,j] = dξdy_ang[iel][src,i,j]
                            metrics.dydη[e_ext1,i,j] = dydη_ang[iel][src,i,j]
                            metrics.dydξ[e_ext1,i,j] = dydξ_ang[iel][src,i,j]
                            criterion_new[e_ext1] = criterion[iel][src]
                            ref_level_new[e_ext1] = ref_level[iel][src]
                            insert_point!(θ, ϕ, connijk_ang_new, e_ext1, i, j)
                        end
                    end
                end

                npoin_ang[iel] = iter - 1
                coords_new     = coords_new[:, 1:npoin_ang[iel]]

                # Build metrics for the four new child elements
                for e_ext1 = e_ext:e_ext+3
                    for i = 1:nop_ang_new[e_ext1]+1, j = 1:nop_ang_new[e_ext1]+1
                        ip  = connijk_ang_new[e_ext1,i,j]
                        θij = coords_new[1,ip]; ϕij = coords_new[2,ip]
                        @turbo for l = 1:nop_ang_new[e_ext1]+1
                            for k = 1:nop_ang_new[e_ext1]+1
                                a = dψ[i,k]*ψ[j,l]; b = ψ[i,k]*dψ[j,l]
                                metrics.dxdξ[e_ext1,k,l] += a*θij
                                metrics.dxdη[e_ext1,k,l] += b*θij
                                metrics.dydξ[e_ext1,k,l] += a*ϕij
                                metrics.dydη[e_ext1,k,l] += b*ϕij
                            end
                        end
                    end
                    @inbounds for l = 1:nop_ang_new[e_ext1]+1, k = 1:nop_ang_new[e_ext1]+1
                        dxdξ_v = metrics.dxdξ[e_ext1,k,l]; dydη_v = metrics.dydη[e_ext1,k,l]
                        dydξ_v = metrics.dydξ[e_ext1,k,l]; dxdη_v = metrics.dxdη[e_ext1,k,l]
                        ip  = connijk_ang_new[e_ext1,k,l]
                        θ = coords_new[1,ip]
                        Je_v   = dxdξ_v*dydη_v - dydξ_v*dxdη_v
                        Jinv   = 1.0 / Je_v
                        metrics.Je[e_ext1,k,l]   = sin(θ) * Je_v
                        metrics.dξdx[e_ext1,k,l] =  dydη_v*Jinv
                        metrics.dξdy[e_ext1,k,l] = -dxdη_v*Jinv
                        metrics.dηdx[e_ext1,k,l] = -dydξ_v*Jinv
                        metrics.dηdy[e_ext1,k,l] =  dxdξ_v*Jinv
                    end
                end

                # Enforce ϕ-periodicity: merge ϕ=2π nodes with ϕ=0 counterparts
                zero_ϕ_dict = Dict{Float64, Int}()
                for iper = 1:iter-1   # use iter-1, not npoin_ang[iel] which isn't updated yet
                    if coords_new[2, iper] < 1e-10   # ϕ ≈ 0
                        zero_ϕ_dict[coords_new[1, iper]] = iper
                    end
                end
                merge_pairs = Dict{Int, Int}()
                for iper = 1:iter-1
                    ϕ = coords_new[2, iper]
                    if abs(ϕ / π - 2.0) < 1e-10   # ϕ ≈ 2π
                        θ = coords_new[1, iper]
                        canonical = get(zero_ϕ_dict, θ, 0)
                        if canonical != 0 && canonical != iper
                            merge_pairs[iper] = canonical
                        end
                    end
                end
    
                if !isempty(merge_pairs)
                    for ip_old in sort(collect(keys(merge_pairs)), rev=true)
                        ip_new = merge_pairs[ip_old]
                        for e = 1:nelem_ang[iel]
                            for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                                connijk_ang_new[e,i,j] == ip_old && (connijk_ang_new[e,i,j] = ip_new)
                            end
                        end
                        for e = 1:nelem_ang[iel]
                            for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                                connijk_ang_new[e,i,j] > ip_old && (connijk_ang_new[e,i,j] -= 1)
                            end
                        end

                        # Rebuild merge_pairs with updated keys AND values
                        merge_pairs = Dict(
                            (k > ip_old ? k - 1 : k) => (v > ip_old ? v - 1 : v)
                                for (k, v) in merge_pairs if k != ip_old
                            )

                        for i = ip_old+1:npoin_ang[iel]
                            coords_new[1,i-1] = coords_new[1,i]
                            coords_new[2,i-1] = coords_new[2,i]
                        end
                        npoin_ang[iel] -= 1
                    end
                end

                connijk_ang[iel] = connijk_ang_new
                coords_ang[iel]  = coords_new
                dxdξ_ang[iel]    = metrics.dxdξ[:,:,:]
                dxdη_ang[iel]    = metrics.dxdη[:,:,:]
                dydξ_ang[iel]    = metrics.dydξ[:,:,:]
                dydη_ang[iel]    = metrics.dydη[:,:,:]
                Je_ang[iel]      = metrics.Je[:,:,:]
                dξdx_ang[iel]    = metrics.dξdx[:,:,:]
                dξdy_ang[iel]    = metrics.dξdy[:,:,:]
                dηdx_ang[iel]    = metrics.dηdx[:,:,:]
                dηdy_ang[iel]    = metrics.dηdy[:,:,:]
                nop_ang[iel]     = nop_ang_new
                criterion[iel]   = criterion_new
                ref_level[iel]   = ref_level_new
                
            
            end
            e_ext += 1
            
        end
    end

    # ── Neighbor detection ────────────────────────────────────────────────────
    for iel = 1:nelem
        
        xmin = xmax = x[connijk[iel,1,1,1]]; ymin = ymax = y[connijk[iel,1,1,1]]
        zmin = zmax = z[connijk[iel,1,1,1]]
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[iel,i,j,k]
            xmin = min(xmin,x[ip]); xmax = max(xmax,x[ip])
            ymin = min(ymin,y[ip]); ymax = max(ymax,y[ip])
            zmin = min(zmin,z[ip]); zmax = max(zmax,z[ip])
        end
        match_bdy = (xmin==xmin_grid)+(xmax==xmax_grid)+(ymin==ymin_grid)+
                    (ymax==ymax_grid)+(zmin==zmin_grid)+(zmax==zmax_grid)
        iter = 1; found_neighbors = 0
        while iter <= nelem && found_neighbors < 26
            xmin_i = xmax_i = x[connijk[iter,1,1,1]]; ymin_i = ymax_i = y[connijk[iter,1,1,1]]
            zmin_i = zmax_i = z[connijk[iter,1,1,1]]
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = connijk[iter,i,j,k]
                xmin_i = min(xmin_i,x[ip]); xmax_i = max(xmax_i,x[ip])
                ymin_i = min(ymin_i,y[ip]); ymax_i = max(ymax_i,y[ip])
                zmin_i = min(zmin_i,z[ip]); zmax_i = max(zmax_i,z[ip])
            end
            m1 = (abs(xmin-xmin_i)<1e-5)+(abs(xmax-xmax_i)<1e-5)+
                 (abs(ymin-ymin_i)<1e-5)+(abs(ymax-ymax_i)<1e-5)+
                 (abs(zmin-zmin_i)<1e-5)+(abs(zmax-zmax_i)<1e-5)
            m2 = (abs(xmin-xmax_i)<1e-5)+(abs(xmax-xmin_i)<1e-5)+
                 (abs(ymin-ymax_i)<1e-5)+(abs(ymax-ymin_i)<1e-5)+
                 (abs(zmin-zmax_i)<1e-5)+(abs(zmax-zmin_i)<1e-5)
            if (m1 > 2 && m2 > 0) || (m1 > 0 && m2 > 1) || m2 > 2
                found_neighbors += 1
                neighbors[iel, found_neighbors, 1] = iter
                if !(adapted_ang[iel]==0 && adapted_ang[iter]==0)
                    if nelem_ang[iel] != nelem_ang[iter] || coords_ang[iel] != coords_ang[iter]
                        neighbors[iel, found_neighbors, 2] = 1
                    end
                end
            end
            if (match_bdy==1 && found_neighbors==17) ||
               (match_bdy==2 && found_neighbors==11) ||
               (match_bdy==3 && found_neighbors==7)
                found_neighbors = 26
            end
            iter += 1
        end
        
    end
end


# ──────────────────────────────────────────────────────────────────────────────
"""
    adaptive_spatial_angular_numbering_3D_2D!(connijk_spa, nelem, ngl, connijk,
        connijk_ang, nop_ang, nelem_ang, coords_ang, x, y, z,
        ref_level, neighbors, adapted, extra_meshes_extra_Je)

Build the combined spatial-angular connectivity array `connijk_spa` and
identify hanging (non-conforming) nodes at adaptive refinement interfaces.

# Algorithm

1. **Point deduplication** — each unique `(x,y,z,θ,ϕ)` tuple is assigned a
   single integer index using a `Dict`.  Shared DOFs (e.g. nodes at spatial
   element boundaries with the same angular coordinates) receive the same index.

2. **Hanging node identification** — for each pair of neighbouring spatial
   elements with different angular refinement levels, interpolation matrices
   (Lagrange) are built between the coarse and fine angular grids.  Nodes on the
   fine side that are not images of coarse nodes (i.e. interpolation weights ≠
   {0,1}) are marked as hanging.

3. **Renumbering** — free nodes are numbered `1:n_free`, hanging nodes
   `n_free+1:n_spa` so that the free block appears first in all subsequent
   matrices.

4. **Constraint matrix** — the sparse matrix `nc_mat` of size `n_free × n_spa`
   encodes the interpolation weights.  Column `j` (a hanging node) contains the
   weights of the free (parent) nodes that interpolate to it.  Free nodes have
   unit diagonal entries.

# Returns

`(nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa)`

| Return value | Description |
|---|---|
| `nc_mat` | Sparse constraint matrix R, size `n_free × n_spa` |
| `nc_non_global_nodes` | Sorted list of hanging node local indices |
| `n_non_global_nodes` | Number of hanging nodes |
| `n_spa` | Total number of unique spatial-angular DOFs |
"""
function adaptive_spatial_angular_numbering_3D_2D!(connijk_spa, nelem, ngl, connijk,
        connijk_ang, nop_ang, nelem_ang, coords_ang, x, y, z,
        ref_level, neighbors, adapted, extra_meshes_extra_Je)

    iter = 1
    interp_sourcesθ = zeros(Float64, nop_ang[1][1]+1)
    interp_targetsθ = zeros(Float64, nop_ang[1][1]+1)
    interp_sourcesϕ = zeros(Float64, nop_ang[1][1]+1)
    interp_targetsϕ = zeros(Float64, nop_ang[1][1]+1)
    ωθ = zeros(Float64, nop_ang[1][1]+1)
    ωϕ = zeros(Float64, nop_ang[1][1]+1)
    Lθ = zeros(Float64, nop_ang[1][1]+1, nop_ang[1][1]+1)
    Lϕ = zeros(Float64, nop_ang[1][1]+1, nop_ang[1][1]+1)
    nc_non_global_nodes = []
    n_non_global_nodes  = 0

    # Phase 1: assign unique indices to all (x,y,z,θ,ϕ) points
    point_dict = Dict{NTuple{5,Float64}, Int}()
    @inbounds for iel = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip  = connijk[iel,i,j,k]
            x_p = x[ip]; y_p = y[ip]; z_p = z[ip]
            for e_ext = 1:nelem_ang[iel]
                for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                    ip_ang = connijk_ang[iel][e_ext,iθ,jθ]
                    θ_p    = coords_ang[iel][1, ip_ang]
                    ϕ_p    = coords_ang[iel][2, ip_ang]
                    key    = (round(x_p,digits=12), round(y_p,digits=12),
                              round(z_p,digits=12), round(θ_p,digits=12),
                              round(ϕ_p,digits=12))
                    idx = get(point_dict, key, 0)
                    if idx == 0
                        point_dict[key] = iter
                        connijk_spa[iel][i,j,k,e_ext,iθ,jθ] = iter
                        iter += 1
                    else
                        connijk_spa[iel][i,j,k,e_ext,iθ,jθ] = idx
                    end
                end
            end
        end
    end

    n_spa = iter - 1

    # Phase 2: identify hanging nodes at non-conforming interfaces
    interpolation_cache = Dict{NTuple{4,Int}, Tuple{Matrix{Float64}, Matrix{Float64}, Int, Int}}()
    nc_non_global_set   = Set{Int}()

    for iel = 1:nelem
        1 ∈ view(neighbors, iel, :, 2) || continue
        for ineighbor = 1:26
            neighbors[iel, ineighbor, 2] != 1 && continue
            adapted = true
            iel1    = neighbors[iel, ineighbor, 1]
            matching_nodes = find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)

            for (i, j, k, i1, j1, k1) in matching_nodes
                for e_ext = 1:nelem_ang[iel]
                    has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang) && continue
                    θmin, θmax, ϕmin, ϕmax = get_element_bounds_fast(iel, e_ext, coords_ang, connijk_ang, nop_ang)
                    for e_ext1 = 1:nelem_ang[iel1]
                        θmin1, θmax1, ϕmin1, ϕmax1 = get_element_bounds_fast(iel1, e_ext1, coords_ang, connijk_ang, nop_ang)
                        is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1) || continue
                        cache_key = (iel, e_ext, iel1, e_ext1)
                        if !haskey(interpolation_cache, cache_key)
                            Lθ_c, Lϕ_c = build_interpolation_matrices!(
                                iel, e_ext, iel1, e_ext1, coords_ang, connijk_ang, nop_ang,
                                θmin, θmax, ϕmin, ϕmax, ϕmin1, ϕmax1,
                                interp_sourcesθ, interp_targetsθ, interp_sourcesϕ, interp_targetsϕ,
                                ωθ, ωϕ, Lθ, Lϕ)
                            nop_p = nop_ang[iel][e_ext]; nop_c = nop_ang[iel1][e_ext1]
                            interpolation_cache[cache_key] = (copy(Lθ_c), copy(Lϕ_c), nop_p, nop_c)
                        end
                        Lθ_use, Lϕ_use, _, nop_child = interpolation_cache[cache_key]
                        for jθ = 1:nop_child+1, iθ = 1:nop_child+1
                            is_vertex_θ = 1 in Lθ_use[iθ,:]
                            is_vertex_ϕ = 1 in Lϕ_use[jθ,:]
                            if !is_vertex_θ || !is_vertex_ϕ
                                jp_spa = connijk_spa[iel1][i1,j1,k1,e_ext1,iθ,jθ]
                                push!(nc_non_global_set, jp_spa)
                            end
                        end
                    end
                end
            end
        end
    end
    nc_non_global_nodes = sort(collect(nc_non_global_set))
    n_non_global_nodes  = length(nc_non_global_nodes)

    # Phase 3: renumber — free nodes first, hanging nodes last
    if n_non_global_nodes > 0
        node_map      = Dict{Int,Int}()
        free_counter  = 1
        for old_idx = 1:n_spa
            old_idx in nc_non_global_set && continue
            node_map[old_idx] = free_counter; free_counter += 1
        end
        hanging_counter = n_spa - n_non_global_nodes + 1
        for old_idx in nc_non_global_nodes
            node_map[old_idx] = hanging_counter; hanging_counter += 1
        end
        @inbounds for iel = 1:nelem
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                for e_ext = 1:nelem_ang[iel]
                    for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                        old_idx = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                        connijk_spa[iel][i,j,k,e_ext,iθ,jθ] = node_map[old_idx]
                    end
                end
            end
        end
        nc_non_global_nodes = [node_map[old_idx] for old_idx in nc_non_global_nodes]
        sort!(nc_non_global_nodes)
    end

    # Phase 4: build constraint matrix
    n_free           = n_spa - n_non_global_nodes
    estimated_entries = n_free + n_non_global_nodes * 20
    I_vec = sizehint!(Int[],     estimated_entries)
    J_vec = sizehint!(Int[],     estimated_entries)
    V_vec = sizehint!(Float64[], estimated_entries)

    for ip_g = 1:n_free  # identity block for free nodes
        push!(I_vec, ip_g); push!(J_vec, ip_g); push!(V_vec, 1.0)
    end

    @inbounds for iel = 1:nelem
        1 ∈ view(neighbors, iel, :, 2) || continue
        for ineighbor = 1:26
            neighbors[iel, ineighbor, 2] != 1 && continue
            iel1 = neighbors[iel, ineighbor, 1]
            matching_nodes = find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)
            for (i, j, k, i1, j1, k1) in matching_nodes
                for e_ext = 1:nelem_ang[iel]
                    has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang) && continue
                    θmin, θmax, ϕmin, ϕmax = get_element_bounds_fast(iel, e_ext, coords_ang, connijk_ang, nop_ang)
                    for e_ext1 = 1:nelem_ang[iel1]
                        θmin1, θmax1, ϕmin1, ϕmax1 = get_element_bounds_fast(iel1, e_ext1, coords_ang, connijk_ang, nop_ang)
                        is_child_element(θmin,θmax,ϕmin,ϕmax,θmin1,θmax1,ϕmin1,ϕmax1) || continue
                        cache_key = (iel, e_ext, iel1, e_ext1)
                        Lθ_use, Lϕ_use, nop_parent, nop_child = interpolation_cache[cache_key]
                        for jθ = 1:nop_parent+1, iθ = 1:nop_parent+1
                            ip_spa = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                            ip_spa > n_free && continue
                            Je_parent = extra_meshes_extra_Je[iel][e_ext,iθ,jθ]
                            for lθ = 1:nop_child+1, kθ = 1:nop_child+1
                                jp_spa = connijk_spa[iel1][i1,j1,k1,e_ext1,kθ,lθ]
                                jp_spa <= n_free && continue
                                Lθ_val = Lθ_use[kθ,iθ]; Lϕ_val = Lϕ_use[lθ,jθ]
                                (abs(Lθ_val) < 1e-14 || abs(Lϕ_val) < 1e-14) && continue
                                (abs(Lθ_val-1.0) < 1e-14 && abs(Lϕ_val-1.0) < 1e-14) && continue
                                Je_child = extra_meshes_extra_Je[iel1][e_ext1,kθ,lθ]
                                push!(I_vec, ip_spa)
                                push!(J_vec, jp_spa)
                                push!(V_vec, Lθ_val*Lϕ_val*(Je_parent/Je_child))
                            end
                        end
                    end
                end
            end
        end
    end

    nc_mat = sparse(I_vec, J_vec, V_vec, n_free, n_spa)

    # Normalise columns so each hanging node's weights sum to 1
    @inbounds for j = 1:n_spa
        col_sum = 0.0
        for idx in nzrange(nc_mat, j)
            col_sum += nonzeros(nc_mat)[idx]
        end
        if abs(col_sum - 1.0) > 1e-13 && abs(col_sum) > 1e-14
            for idx in nzrange(nc_mat, j)
                nonzeros(nc_mat)[idx] /= col_sum
            end
        end
    end

    return nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa
end


# ──────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ──────────────────────────────────────────────────────────────────────────────

"""
    find_face_node_match(ngl, iel, ip, connijk)

Search the face nodes of element `iel` for node `ip`.
Returns `(i, j, k, ip)` where `(i,j,k)` are the local indices, or `(*, *, *, 0)`
if not found.

Iterates over the 3D surface (faces) of a hexahedral element in a fixed order.
"""
function find_face_node_match(ngl, iel, ip, connijk)
    found = false; ip1 = 0; iter = 0
    i = 0; j = 0; k = 0
    while !found && iter <= ngl^3-(ngl-2)^3
        if iter < ngl^2
            i = 1; j = (iter ÷ ngl)+1; k = iter % ngl+1
        elseif iter < (2*ngl^2)-ngl
            i = (iter ÷ ngl)%(ngl-1)+1; if i==1 i=ngl end
            j = ngl; k = iter % ngl+1
        elseif iter < (3*ngl^2)-(2*ngl)
            i = ngl; j = (iter ÷ ngl)%(ngl-1)+1; k = iter % ngl+1
        elseif iter < (4*ngl^2)-(4*ngl)
            i = (iter ÷ ngl)%(ngl-1)+1; j = 1; k = iter % ngl+1
        elseif iter < (4*ngl^2)-(4*ngl)+(ngl-2)^2
            d = ((4*ngl^2)-(4*ngl)) ÷ (ngl-2); r = ((4*ngl^2)-(4*ngl)) - d*(ngl-2)
            i = ((iter-r) ÷ (ngl-2))-d + 2; j = iter % (ngl-2)+2; k = 1
        else
            d = ((4*ngl^2)-(4*ngl)+(ngl-2)^2) ÷ (ngl-2)
            r = ((4*ngl^2)-(4*ngl)+(ngl-2)^2) - d*(ngl-2)
            i = ((iter-r) ÷ (ngl-2))-d + 2; j = iter % (ngl-2)+2; k = ngl
        end
        if connijk[iel,i,j,k] == ip
            ip1 = connijk[iel,i,j,k]; found = true
        end
        iter += 1
    end
    return i, j, k, ip1
end


"""
    find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)

Return all face node pairs `(i,j,k, i1,j1,k1)` shared between elements `iel`
and `iel1`.  Uses a `Set` for O(1) lookup of `iel1` nodes.
"""
function find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)
    matching  = Tuple{Int,Int,Int,Int,Int,Int}[]
    iel1_nodes = Set(connijk[iel1,:,:,:])
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        is_face = (i==1 || i==ngl || j==1 || j==ngl || k==1 || k==ngl)
        is_face || continue
        ip = connijk[iel,i,j,k]
        ip in iel1_nodes || continue
        i1, j1, k1, ip1 = find_face_node_match(ngl, iel1, ip, connijk)
        ip1 != 0 && push!(matching, (i,j,k,i1,j1,k1))
    end
    return matching
end


"""
    has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang)

Return `true` if angular element `e_ext` of `iel` has an identical counterpart
in `iel1` (same coordinates at all quadrature points).
"""
function has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang)
    for e_check = 1:nelem_ang[iel1]
        if (coords_ang[iel][1, connijk_ang[iel][e_ext,:,:]] ==
            coords_ang[iel1][1, connijk_ang[iel1][e_check,:,:]] &&
            coords_ang[iel][2, connijk_ang[iel][e_ext,:,:]] ==
            coords_ang[iel1][2, connijk_ang[iel1][e_check,:,:]])
            return true
        end
    end
    return false
end


"""
    get_element_bounds_fast(iel, e_ext, coords_ang, connijk_ang, nop_ang)

Return `(θmin, θmax, ϕmin, ϕmax)` for angular element `e_ext` of spatial
element `iel`.  ϕmax=0 is remapped to 2π for the periodicity boundary.
"""
function get_element_bounds_fast(iel, e_ext, coords_ang, connijk_ang, nop_ang)
    nop  = nop_ang[iel][e_ext]
    θmin = coords_ang[iel][1, connijk_ang[iel][e_ext,1,1]]
    θmax = coords_ang[iel][1, connijk_ang[iel][e_ext,nop+1,nop+1]]
    ϕmin = coords_ang[iel][2, connijk_ang[iel][e_ext,1,1]]
    ϕmax = coords_ang[iel][2, connijk_ang[iel][e_ext,nop+1,nop+1]]
    ϕmax == 0.0 && (ϕmax = 2π)
    return minmax(θmin,θmax)..., minmax(ϕmin,ϕmax)...
end


"""
    is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1; tol=1e-14)

Return `true` if the angular element with bounds `(θmin1,θmax1,ϕmin1,ϕmax1)` is
contained within `(θmin,θmax,ϕmin,ϕmax)` (i.e. is a child after refinement).
"""
function is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1; tol=1e-14)
    return θmin1 >= θmin-tol && θmax1 <= θmax+tol &&
           ϕmin1 >= ϕmin-tol && ϕmax1 <= ϕmax+tol
end


"""
    build_interpolation_matrices!(iel, e_ext, iel1, e_ext1, coords_ang, connijk_ang, nop_ang,
        θmin, θmax, ϕmin, ϕmax, ϕmin1, ϕmax1,
        interp_sourcesθ, interp_targetsθ, interp_sourcesϕ, interp_targetsϕ,
        ωθ, ωϕ, Lθ, Lϕ)

Build and row-normalise the Lagrange interpolation matrices `Lθ` and `Lϕ` that
map from the parent angular element `(iel, e_ext)` to the child element
`(iel1, e_ext1)`.

Interpolation is performed via barycentric weights in θ and ϕ independently.
Rows are normalised to sum to 1 (partition of unity) to ensure exact
reproduction of constants.

Returns views into `Lθ` and `Lϕ` sized to `(nop_child+1, nop_parent+1)`.
"""
function build_interpolation_matrices!(iel, e_ext, iel1, e_ext1, coords_ang, connijk_ang, nop_ang,
                                      θmin, θmax, ϕmin, ϕmax, ϕmin1, ϕmax1,
                                      interp_sourcesθ, interp_targetsθ, interp_sourcesϕ, interp_targetsϕ,
                                      ωθ, ωϕ, Lθ, Lϕ)
    nop_parent = nop_ang[iel][e_ext]
    nop_child  = nop_ang[iel1][e_ext1]

    for iθ = 1:nop_parent+1
        interp_sourcesθ[iθ] = coords_ang[iel][1, connijk_ang[iel][e_ext,iθ,iθ]]
        interp_sourcesϕ[iθ] = coords_ang[iel][2, connijk_ang[iel][e_ext,iθ,iθ]]
        ϕmax == 2π && iθ == nop_parent+1 && (interp_sourcesϕ[iθ] = ϕmax)
    end
    for iθ = 1:nop_child+1
        interp_targetsθ[iθ] = coords_ang[iel1][1, connijk_ang[iel1][e_ext1,iθ,iθ]]
        interp_targetsϕ[iθ] = coords_ang[iel1][2, connijk_ang[iel1][e_ext1,iθ,iθ]]
        ϕmax1 == 2π && iθ == nop_child+1 && (interp_targetsϕ[iθ] = ϕmax1)
    end

    fill!(Lθ, 0.0); fill!(Lϕ, 0.0); fill!(ωθ, 0.0); fill!(ωϕ, 0.0)
    BarycentricWeights!(view(interp_sourcesθ, 1:nop_parent+1), view(ωθ, 1:nop_parent+1))
    BarycentricWeights!(view(interp_sourcesϕ, 1:nop_parent+1), view(ωϕ, 1:nop_parent+1))
    PolynomialInterpolationMatrix!(
        view(interp_sourcesθ, 1:nop_parent+1), view(ωθ, 1:nop_parent+1),
        view(interp_targetsθ, 1:nop_child+1),  view(Lθ, 1:nop_child+1, 1:nop_parent+1))
    PolynomialInterpolationMatrix!(
        view(interp_sourcesϕ, 1:nop_parent+1), view(ωϕ, 1:nop_parent+1),
        view(interp_targetsϕ, 1:nop_child+1),  view(Lϕ, 1:nop_child+1, 1:nop_parent+1))

    for i = 1:nop_child+1
        rs_θ = sum(Lθ[i, 1:nop_parent+1]); abs(rs_θ) > 1e-14 && (Lθ[i,1:nop_parent+1] ./= rs_θ)
        rs_ϕ = sum(Lϕ[i, 1:nop_parent+1]); abs(rs_ϕ) > 1e-14 && (Lϕ[i,1:nop_parent+1] ./= rs_ϕ)
    end

    return view(Lθ, 1:nop_child+1, 1:nop_parent+1),
           view(Lϕ, 1:nop_child+1, 1:nop_parent+1)
end

function enforce_periodic_phi_conformity!(criterion, thresholds, ref_level,
                                           nelem, nelem_ang, nop_ang,
                                           connijk_ang, coords_ang)
    @info "periodic_phi_conformity entered"
    for iel = 1:nelem
        for e_ext = 1:nelem_ang[iel]
            nop = nop_ang[iel][e_ext]

            # ϕ at first and last node columns of this element
            ϕ_first = coords_ang[iel][2, connijk_ang[iel][e_ext, 1,     1]]
            ϕ_last  = coords_ang[iel][2, connijk_ang[iel][e_ext, 1, nop+1]]

            # An element wraps the periodic boundary if its ϕ sequence
            # increases but the last column is less than the first —
            # i.e. the sequence went e.g. 5.7 → 5.9 → 6.1 → 0.0
            touches_2pi = ϕ_last < ϕ_first

            # ϕ=0 side: first column is near 0 and last column is larger
            # (the periodic partner of the above)
            touches_0   = ϕ_first < 1e-10 && ϕ_last > ϕ_first
            
            (touches_2pi || touches_0) || continue

            # θ range of this element — used to find the paired element
            θ_lo = coords_ang[iel][1, connijk_ang[iel][e_ext, 1,     1]]
            θ_hi = coords_ang[iel][1, connijk_ang[iel][e_ext, nop+1, 1]]

            for e_pair = 1:nelem_ang[iel]
                e_pair == e_ext && continue
                nop_p = nop_ang[iel][e_pair]

                ϕ_pair_first = coords_ang[iel][2, connijk_ang[iel][e_pair, 1,       1]]
                ϕ_pair_last  = coords_ang[iel][2, connijk_ang[iel][e_pair, 1, nop_p+1]]

                pair_2pi = ϕ_pair_last < ϕ_pair_first
                pair_0   = ϕ_pair_first < 1e-10 && ϕ_pair_last > ϕ_pair_first

                θ_lo_p = coords_ang[iel][1, connijk_ang[iel][e_pair, 1,       1]]
                θ_hi_p = coords_ang[iel][1, connijk_ang[iel][e_pair, nop_p+1, 1]]

                # Counterpart: same θ range, opposite ϕ boundary
                if abs(θ_lo_p - θ_lo) < 1e-10 && abs(θ_hi_p - θ_hi) < 1e-10 &&
                   ((touches_2pi && pair_0) || (touches_0 && pair_2pi))
                    
                    if abs(criterion[iel][e_ext]) > abs(criterion[iel][e_pair])
                        criterion[iel][e_pair] = criterion[iel][e_ext]
                    end
                    break
                end
            end
        end
    end
end

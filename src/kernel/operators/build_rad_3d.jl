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
        extra_mesh, κ, σ, atmos_data, QT::Inexact, SD::NSD_3D, AD::ContGal)

    bdy = make_boundary_predicates(mesh)
    sw = []
    lw = []
    
    #=sw_base = SWParams(1361.0, 0.6, 0.0, 0.0)
    for δ in [0.1, 0.2, 0.3, 0.4, 0.5]
        sw_test = SWParams(sw_base.S₀_flux, sw_base.μ₀, sw_base.φ₀, δ)
        result  = check_beam_flux(extra_mesh, sw_test, extra_mesh.extra_nop[1])
        
    end=#
    
    if(inputs[:lRT_from_data])

        sw  = SWParams(inputs[:RT_S0_flux], inputs[:RT_μ0], inputs[:RT_ϕ0], 0.25)#inputs[:RT_δ_beam])
        lw  = LWParams(inputs[:RT_ϵ_surface], inputs[:RT_T_space])
    end
    
    nc_mat = zeros(Float64,1)
    P = zeros(Float64,1)
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

        # ── Pre-adaptivity connectivity and numbering ─────────────────────────
        @rankinfo rank "Building initial adaptive spatial-angular connectivity..."
        nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa =
            adaptive_spatial_angular_numbering_3D_2D!(
                connijk_spa, nelem, ngl, mesh.connijk,
                extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                extra_meshes_coords, mesh.x, mesh.y, mesh.z,
                extra_meshes_ref_level, neighbors, adapted, extra_meshes_extra_Je)

        @rankinfo rank "Building pre-adaptivity global numbering..."
        ip2gip_spa, gip2ip, gip2owner_spa, gnpoin =
            setup_global_numbering_adaptive_angular_scalable(
                mesh.ip2gip, mesh.gip2owner, mesh, connijk_spa,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                n_spa, n_non_global_nodes, nc_non_global_nodes)

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

        ghost_layer = build_nonconforming_ghost_layer(
            mesh, connijk_spa, mesh.ip2gip, ip2gip_spa, gip2owner_spa,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
            extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
            extra_meshes_ref_level, n_spa, neighbors)

        # ── Pre-adaptivity matrices for adaptivity criterion ─────────────────
        @rankinfo rank "Assembling pre-adaptivity LHS and mass matrix..."
        LHS = sparse_lhs_assembly_3Dby2D_adaptive(
            ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
            mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ,
            extra_meshes_connijk, extra_meshes_extra_Je,
            extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl,
            extra_meshes_extra_nelems, dξdx, dξdy, dξdz,
            dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
            extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, inputs[:lRT_from_data], κ, σ)

        Md = assemble_mass_diagonal_3Dby2D_adaptive(
            ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
            extra_meshes_extra_Je, extra_meshes_extra_nops, n_spa, nelem,
            ngl, extra_meshes_extra_nelems, connijk_spa)

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

        # ── Angular grid refinement ───────────────────────────────────────────
        @rankinfo rank "Adapting angular grid..."
        
        adapt_angular_grid_3Dby2D!(
            criterion, inputs[:RT_amr_threshold], extra_meshes_ref_level, nelem, ngl,
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
            extra_mesh[1].ψ, extra_mesh[1].dψ)

            nonowned_parent_indices = Set{Int}()
            nc_mat = sparse(I,n_spa,n_spa)
            nc_mat_rhs = sparse(I,n_spa,n_spa)
            P = nc_mat'
            P_vec = nc_mat_rhs'

        if !(maximum(maximum.(extra_meshes_ref_level)) == 0)

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
                    mesh.ip2gip, mesh.gip2owner, mesh, connijk_spa,
                    extra_meshes_coords, extra_meshes_connijk,
                    extra_meshes_extra_nops, extra_meshes_extra_nelems,
                    n_spa, n_non_global_nodes, nc_non_global_nodes)

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

            # Precompute reverse ghost map (reused for solution prolongation)
            @rankinfo rank "Building reverse ghost constraint map..."
            reverse_ghost_map = build_reverse_ghost_constraint_map(
                ghost_constraint_data, ip2gip_spa, gip2owner_spa, rank, gip_to_local)

            # ── LHS and mass matrix on adapted mesh ───────────────────────────
            @rankinfo rank "Assembling LHS on adapted mesh..."
            LHS = sparse_lhs_assembly_3Dby2D_adaptive(
                ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ,
                extra_meshes_connijk, extra_meshes_extra_Je,
                extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl,
                extra_meshes_extra_nelems, dξdx, dξdy, dξdz,
                dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
                extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, inputs[:lRT_from_data], κ, σ)

            P    = nc_mat'
            rest = nc_mat

            Md = assemble_mass_diagonal_3Dby2D_adaptive(
                ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                extra_meshes_extra_Je, extra_meshes_extra_nops, n_spa, nelem,
                ngl, extra_meshes_extra_nelems, connijk_spa)
        end  # if adapted

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
            MLHS = allreduce_parent_parent_entries(
                MLHS, nc_mat, gip2owner_extra, n_spa, rank, comm,
                ip2gip_spa, local_parent_indices, nonowned_parent_indices,
                nonowned_parent_gids, gip_to_local)
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
        n_free = n_spa - length(nc_non_global_nodes)
        @rankinfo rank "Extracting free-node submatrix (removing hanging rows/cols)..."
        A_free = extract_free_submatrix_remove_all_hanging(
            A_with_col_effects, all_hanging_nodes, n_free, n_spa_g, rank)

        A = sparse(A_free)

        # Remove non-owned parent–parent entries to avoid double-counting.
        # These were all-reduced above and must be zeroed on non-owning ranks.
        for i in nonowned_parent_indices
            for j in nonowned_parent_indices
                if abs(A[i,j]) > 0.0
                    A[i,j] = 0.0
                end
            end
        end
        A = sparse(A)

        RHS = zeros(TFloat, n_spa_g)
        ref = zeros(TFloat, n_spa)
        BDY = zeros(TFloat, n_spa_g)

    else  # ── Non-adaptive path ────────────────────────────────────────────────

        npoin_ang_total = npoin * extra_mesh.extra_npoin
        n_spa = npoin_ang_total
        @rankinfo rank "Assembling LHS ($npoin_ang_total DOF)..."
        LHS = sparse_lhs_assembly_3Dby2D(
            ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
            mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh.ψ,
            extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
            extra_mesh.extra_coords, extra_mesh.extra_nop,
            npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
            dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
            extra_mesh.extra_npoin, inputs[:rad_HG_g],
            inputs[:lRT_from_data], κ, σ)

        Md = assemble_mass_diagonal_3Dby2D(
            ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ,
            extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je,
            extra_mesh.extra_nop, npoin_ang_total, nelem, ngl,
            extra_mesh.extra_nelem, extra_mesh.extra_npoin)

        ip2gip_spa, gip2owner_extra, gnpoin =
            setup_global_numbering_extra_dim(
                mesh.ip2gip, mesh.gip2owner, npoin,
                extra_mesh.extra_npoin, npoin_ang_total)

        pM = setup_assembler(SD, Md, ip2gip_spa, gip2owner_extra)
        if pM !== nothing
            assemble_mpi!(Md, pM)
        end

        M_inv = spdiagm(0 => 1.0 ./ Md)
        A     = sparse(M_inv * LHS)
        M_inv = nothing; LHS = nothing
        GC.gc()

        BDY = zeros(TFloat, npoin_ang_total)
        RHS = zeros(TFloat, npoin_ang_total)
        ref = zeros(TFloat, npoin_ang_total)
    end  # adaptive / non-adaptive

    # ── Common setup: nc_mat row structure, free DOF count ────────────────────
    nc_rows  = size(nc_mat, 1) > 2 ? rowvals(nc_mat) : zeros(Int, 1, 1)
    n_free   = inputs[:adaptive_extra_meshes] ? n_spa - n_non_global_nodes : npoin_ang_total

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
                matchx = (x == mesh.xmax || x == mesh.xmin)
                matchy = (y == mesh.ymax || y == mesh.ymin)
                matchz = (z == mesh.zmax || z == mesh.zmin)
                nmatches = matchx + matchy + matchz

                # Collect all distinct face normals for this node
                # (one per face the node belongs to)
                collected = NTuple{3,Float64}[]
                push!(collected, (nx[iface, face_i, face_j],
                                  ny[iface, face_i, face_j],
                                  nz[iface, face_i, face_j]))

                if nmatches >= 2
                    added  = 0
                    needed = nmatches - 1
                    for (iface2, fi2, fj2) in get(node_to_bdy_faces, ip, [])
                        iface2 == iface && continue
                        if nmatches == 2
                            if abs(collected[1][1] - nx[iface2,fi2,fj2]) > 1e-12 ||
                               abs(collected[1][2] - ny[iface2,fi2,fj2]) > 1e-12 ||
                               abs(collected[1][3] - nz[iface2,fi2,fj2]) > 1e-12
                                push!(collected, (nx[iface2,fi2,fj2],
                                                  ny[iface2,fi2,fj2],
                                                  nz[iface2,fi2,fj2]))
                                added += 1
                            end
                        else
                            push!(collected, (nx[iface2,fi2,fj2],
                                              ny[iface2,fi2,fj2],
                                              nz[iface2,fi2,fj2]))
                            added += 1
                        end
                        added >= needed && break
                    end
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

                        if is_boundary
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
                            bdy_normals_for_ghosts[ip] = (best_nx, best_ny, best_nz)

                            if best_dot < -1e-13
                                if ip_g <= n_free && !(ip_g in all_hanging_nodes)
                                    val = 0.0
                                    if inputs[:RT_shortwave]
                                        val = user_rad_bc_shortwave(x, y, z, θ, ϕ, best_nx, best_ny, best_nz, bdy, sw)
                                    elseif inputs[:RT_longwave]
                                        val = user_rad_bc_longwave(x, y, z, θ, ϕ, best_nx, best_ny, best_nz, ip, atmos_data, bdy, lw)
                                    else
                                        val = user_rad_bc(x, y, z, θ, ϕ)
                                    end
                                    BDY[ip_g] = val
                                    if !haskey(boundary_dict, ip_g)
                                        boundary_dict[ip_g] = is_owned ? val : 0.0
                                    end
                                end
                            else
                                if is_owned
                                    RHS[ip_g] = if inputs[:RT_shortwave]
                                        user_rhs_shortwave(x, y, z, θ, ϕ)
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
                                    user_rhs_shortwave(x, y, z, θ, ϕ)
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
                        ip_g    = (ip-1) * extra_mesh.extra_npoin + ip_ext
                        is_owned = (gip2owner_extra[ip_g] == rank)

                        sip = exp(-((6/(2π)) * (θ - 3π/5))^2)
                        bip = exp(-((6/(2π)) * (ϕ - 2π/3))^2)
                        ref[ip_g] = sf * sip * bip

                        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                        if is_boundary
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
                            bdy_normals_for_ghosts[ip] = (best_nx, best_ny, best_nz)
                            if best_dot < -1e-13
                                val = 0.0
                                if inputs[:RT_shortwave]
                                    val = user_rad_bc_shortwave(x, y, z, θ, ϕ, best_nx, best_ny, best_nz, bdy, sw)
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
                                        user_rhs_shortwave(x, y, z, θ, ϕ)
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
                                    user_rhs_shortwave(x, y, z, θ, ϕ)
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

    # ── Ghost boundary nodes (extended ghost-parent DOFs) ─────────────────────
    if inputs[:adaptive_extra_meshes] && n_spa_g > n_spa
        n_extended = n_spa_g - n_spa
        for i = 1:n_extended
            ip   = extended_parents_ip[i]
            x    = extended_parents_x[i]; y = extended_parents_y[i]; z = extended_parents_z[i]
            θ    = extended_parents_θ[i]; ϕ = extended_parents_ϕ[i]
            ip_g = n_spa + i

            if ip in mesh.poin_in_bdy_face
                nx_new = bdy_normals_for_ghosts[ip][1]
                ny_new = bdy_normals_for_ghosts[ip][2]
                nz_new = bdy_normals_for_ghosts[ip][3]
                if nx_new*sin(θ)*cos(ϕ) + ny_new*sin(θ)*sin(ϕ) + nz_new*cos(θ) < -1e-13
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
    

    # ── RHS restriction ───────────────────────────────────────────────────────
    if inputs[:adaptive_extra_meshes]
        @rankinfo rank "Restricting RHS..."
        rhs_effects_to_send =
            compute_hanging_rhs_effects_before_restriction(
                ghost_constraint_data_rhs, RHS, ip2gip_spa, gip2owner_spa, rank)

        RHS_restricted = nc_mat_rhs * RHS

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
    else
        for (node, val) in boundary_dict
            RHS[node] = val
        end
        B = RHS
    end

    # ── Linear solve ──────────────────────────────────────────────────────────
    @rankinfo rank "Solving system ($(size(A,1)) DOF)..."
    As  = sparse(A)
    A   = nothing; GC.gc()

    npoin_ang_total = size(B, 1)
    @info maximum(As), minimum(As), maximum(B), minimum(B)
    solution = if inputs[:adaptive_extra_meshes]
        solve_parallel_lsqr(ip2gip_spa, gip2owner_extra, As, B, gnpoin, n_spa, pM;
            npoin_g = n_spa_g,
            g_ip2gip = extended_parents_to_gid,
            g_gip2ip = gid_to_extended_parents)
    else
        solve_parallel_lsqr(ip2gip_spa, gip2owner_extra, As, B, gnpoin, npoin_ang_total, pM;
        npoin_g = npoin_ang_total)
    end

    @rankinfo rank "Solve complete."
    A = nothing; RHS = nothing; GC.gc()
    @info maximum(solution), minimum(solution)
    # ── Solution prolongation ─────────────────────────────────────────────────
    solution_new = zeros(Float64, n_spa)
    if inputs[:adaptive_extra_meshes]
        @rankinfo rank "Prolonging solution to full mesh..."
        solution_contributions_to_send =
            compute_solution_prolongation_contributions(
                reverse_ghost_map, solution, ip2gip_spa, n_free, rank)

        solution_local = P_vec * solution

        received_solution_contributions =
            exchange_hanging_effects_vector(solution_contributions_to_send, rank, comm)

        solution_new = add_solution_prolongation_contributions(
            @view(solution_local[1:n_spa]), received_solution_contributions,
            ip2gip_spa, n_spa, rank, gip_to_local)
    end

    # ── Angular integration and error computation ─────────────────────────────
    @rankinfo rank "Integrating solution in angle..."
    int_sol      = zeros(TFloat, npoin, 1)
    int_ref      = zeros(TFloat, npoin, 1)
    int_sol_accum = zeros(npoin)
    int_ref_accum = zeros(npoin)
    L2_err = 0.0
    L2_ref = 0.0

    # Compute node sharing divisors for correct accumulation at shared nodes
    node_div = ones(Int, npoin)
    
        for iel = 1:nelem
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip    = mesh.connijk[iel, i, j, k]
                x_ip  = mesh.x[ip]; y_ip = mesh.y[ip]; z_ip = mesh.z[ip]
                matchx  = abs(x_ip - mesh.xmin) < 1e-10 || abs(x_ip - mesh.xmax) < 1e-10
                matchy  = abs(y_ip - mesh.ymin) < 1e-10 || abs(y_ip - mesh.ymax) < 1e-10
                matchz  = abs(z_ip - mesh.zmin) < 1e-10 || abs(z_ip - mesh.zmax) < 1e-10
                nmatches = matchx + matchy + matchz
                on_bdy   = ip in mesh.poin_in_bdy_face

                is_corner = (i ∈ (1,ngl)) && (j ∈ (1,ngl)) && (k ∈ (1,ngl))
                is_edge   = !is_corner && (
                    ((i ∈ (1,ngl)) && (j ∈ (1,ngl))) ||
                    ((i ∈ (1,ngl)) && (k ∈ (1,ngl))) ||
                    ((j ∈ (1,ngl)) && (k ∈ (1,ngl))))
                is_face   = !is_corner && !is_edge &&
                            (i ∈ (1,ngl) || j ∈ (1,ngl) || k ∈ (1,ngl))

                div = 1
                if is_corner
                    if !on_bdy; div = 8
                    elseif nmatches == 1; div = 4
                    elseif nmatches == 2; div = 2
                    end
                elseif is_edge
                    if !on_bdy; div = 4
                    elseif nmatches == 1; div = 2
                    end
                elseif is_face
                    if !on_bdy; div = 2 end
                end
                node_div[ip] = max(node_div[ip], div)
            end
        end

        if inputs[:adaptive_extra_meshes]
            for iel = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[iel, i, j, k]
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
                            ip_g   = (ip-1) * extra_mesh.extra_npoin + ip_ext
                            int_sol_accum[ip] += solution[ip_g] *
                                extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ] *
                                extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ] / node_div[ip]
                            if (inputs[:lmanufactured_solution])    
                                int_ref_accum[ip] += (ref[ip_g] - solution[ip_g]) *
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
    

    # ── Global reduction of integrals and norms ───────────────────────────────
    
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
        end
        if (inputs[:lmanufactured_solution])
            L2_ref_g = MPI.Allreduce(L2_ref, MPI.SUM, comm)
            L2_err_g = MPI.Allreduce(L2_err, MPI.SUM, comm)
        end
    

    if (inputs[:lmanufactured_solution])
        @rankinfo rank @sprintf("L2 error: ‖e‖ = %.6e  ‖u‖ = %.6e  relative = %.6e",
        sqrt(L2_err_g), sqrt(L2_ref_g), sqrt(L2_err_g / L2_ref_g))
    end

    # ── VTK output ────────────────────────────────────────────────────────────
    @rankinfo rank "Writing output"
    title = @sprintf "Solution-Radiation"
    write_vtk(SD, mesh, int_sol, int_sol, nothing, nothing, nothing,
              0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
              ["Ang_int"], ["Ang_int"]; iout=1, nvar=1)

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
                                    npoin_ang, rad_HG_g, rad_data, κ_data, σ_data)

    intΦ = zeros(Float64, npoin_ang)
    HG, _ = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)),
                   0, 2π, rtol=1e-13, atol=1e-13)

    # Precompute scattering integrals once per angular node
    for e_ext = 1:nelem_ang
        for jθ = 1:nop_ang[e_ext]+1, iθ = 1:nop_ang[e_ext]+1
            ip_ext = connijk_ang[e_ext, iθ, jθ]
            θ = coords_ang[1, ip_ext]; ϕ = coords_ang[2, ip_ext]
            val = 0.0
            for e_ext_s = 1:nelem_ang
                for nθ = 1:nop_ang[e_ext_s]+1, mθ = 1:nop_ang[e_ext_s]+1
                    ipθ = connijk_ang[e_ext_s, mθ, nθ]
                    Φ   = user_scattering_functions(θ, coords_ang[1,ipθ], ϕ, coords_ang[2,ipθ], HG)
                    val += ωθ[mθ] * ωϕ[nθ] * Je_ang[e_ext_s, mθ, nθ] * Φ
                end
            end
            intΦ[ip_ext] = val
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

                    _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                    I_vec[pos] = idx_ip; J_vec[pos] = idx_ip
                    V_vec[pos] = (κ - σ*intΦ[ip_ext]) * ωJac_full
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
                                    npoin_ang, rad_HG_g, connijk_spa, rad_data, κ_data, σ_data)

    ngl_stencil = 3 * (ngl - 1) + 1
    estimated   = round(Int, npoin_ang_total * ngl_stencil * 1.2)
    I_vec = Vector{Int}(undef, estimated)
    J_vec = Vector{Int}(undef, estimated)
    V_vec = Vector{Float64}(undef, estimated)
    pos   = 1

    HG, _ = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)),
                   0, 2π, rtol=1e-13, atol=1e-13)

    # Precompute per-element scattering integrals
    intΦ_per_iel = Vector{Vector{Float64}}(undef, nelem)
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
                                                        ϕ, coords_ang[iel][2,ipθ], HG)
                        val += ωθ[mθ] * ωϕ[nθ] * Je_ang[iel][e_ext_s, mθ, nθ] * Φ
                    end
                end
                intΦ[ip_ext] = val
            end
        end
        intΦ_per_iel[iel] = intΦ
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

                    _grow_if_needed!(I_vec, J_vec, V_vec, pos)
                    I_vec[pos] = idx_ip; J_vec[pos] = idx_ip
                    V_vec[pos] = (κ - σ*intΦ_per_iel[iel][ip_ext]) * ωJac_full
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
        x, y, z, xmin_grid, ymin_grid, zmin_grid, xmax_grid, ymax_grid, zmax_grid, ψ, dψ)

    lgl = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    adapted_ang = zeros(Int, nelem)
    ang_adapted  = zeros(Int, nelem)

    for iel = 1:nelem
        e_ext = 1
        while e_ext <= nelem_ang[iel]
            if abs(criterion[iel][e_ext]) > thresholds[1]

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
                for iper = 1:npoin_ang[iel]
                    if coords_new[2, iper] <= eps(Float64)
                        zero_ϕ_dict[coords_new[1,iper]] = iper
                    end
                end
                merge_pairs = Dict{Int, Int}()
                for iper = 1:npoin_ang[iel]
                    ϕ = coords_new[2, iper]
                    if abs(ϕ/π - 2.0) <= eps(Float64)
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
                        for k in keys(merge_pairs)
                            merge_pairs[k] > ip_old && (merge_pairs[k] -= 1)
                        end
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
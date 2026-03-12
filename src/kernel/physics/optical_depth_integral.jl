
using NearestNeighbors

# ── 1. Ray grid definition ────────────────────────────────────────────────────

"""
    RayGrid

A structured grid of parallel rays in direction Ω_sun, used to compute
slant-path optical depth τ in a general unstructured parallel mesh.

Rays are parameterised by their entry point on the TOA plane, sampled
on a uniform 2D grid in the (e₁, e₂) basis perpendicular to Ω_sun.
"""
struct RayGrid
    n_rays_u  :: Int              # number of rays in first transverse direction
    n_rays_v  :: Int              # number of rays in second transverse direction
    n_steps   :: Int              # number of integration steps along each ray
    ds        :: Float64          # step size along ray (m)
    origins   :: Matrix{Float64}  # (3, n_rays) entry points at TOA (z=zmax)
    Ω_sun     :: NTuple{3,Float64}# unit beam direction (downward)
    e1        :: NTuple{3,Float64}# first basis vector perpendicular to Ω_sun
    e2        :: NTuple{3,Float64}# second basis vector perpendicular to Ω_sun
end

"""
    build_ray_grid(sw, mesh_global_bounds, n_rays_u, n_rays_v, n_steps_per_domain)

Construct the ray grid for solar direction defined by sw.μ₀, sw.φ₀.

# Arguments
- `mesh_global_bounds`: NamedTuple with xmin,xmax,ymin,ymax,zmin,zmax
- `n_rays_u`, `n_rays_v`: ray density in the two transverse directions.
  Should be fine enough to resolve κ_ext variations. A good default is
  2× the number of spatial elements per direction.
- `n_steps_per_domain`: number of integration steps along the full domain
  height. Should resolve the vertical κ_ext profile — 4× the number of
  vertical elements is usually sufficient.
"""
function build_ray_grid(sw, bounds, n_rays_u, n_rays_v, n_steps_per_domain)

    μ₀    = sw.μ₀
    φ₀    = sw.φ₀
    θ_sun = π - acos(clamp(μ₀, 0.0, 1.0))

    # Beam direction (downward into atmosphere)
    Ωx = sin(θ_sun)*cos(φ₀)
    Ωy = sin(θ_sun)*sin(φ₀)
    Ωz = cos(θ_sun)   # negative (downward)
    Ω  = (Ωx, Ωy, Ωz)

    # Two basis vectors perpendicular to Ω_sun
    # e1: horizontal, perpendicular to Ω projected onto xy plane
    # e2: completes right-hand system
    norm_xy = sqrt(Ωx^2 + Ωy^2)
    if norm_xy > 1e-10
        e1 = (-Ωy/norm_xy, Ωx/norm_xy, 0.0)
    else
        e1 = (1.0, 0.0, 0.0)   # beam is vertical
    end
    e2x = Ωy*e1[3] - Ωz*e1[2]
    e2y = Ωz*e1[1] - Ωx*e1[3]
    e2z = Ωx*e1[2] - Ωy*e1[1]
    e2  = (e2x, e2y, e2z)

    # Domain geometry
    Lx = bounds.xmax - bounds.xmin
    Ly = bounds.ymax - bounds.ymin
    Lz = bounds.zmax - bounds.zmin

    # Step size: resolve vertical profile with enough steps
    # Slant path length = Lz / |μ₀|
    slant_length = Lz / abs(μ₀)
    ds           = slant_length / n_steps_per_domain

    # Ray origins at TOA (z = zmax), distributed over the TOA plane.
    # The TOA plane in (e1,e2) coordinates spans the full domain footprint
    # plus a margin to account for the beam displacement at low sun angles.
    # Margin = horizontal displacement of beam across full domain height
    margin_x = abs(Ωx/abs(Ωz)) * Lz
    margin_y = abs(Ωy/abs(Ωz)) * Lz

    u_min = bounds.xmin - margin_x; u_max = bounds.xmax + margin_x
    v_min = bounds.ymin - margin_y; v_max = bounds.ymax + margin_y

    du = (u_max - u_min) / (n_rays_u - 1)
    dv = (v_max - v_min) / (n_rays_v - 1)

    n_rays   = n_rays_u * n_rays_v
    origins  = zeros(Float64, 3, n_rays)
    ir       = 1
    for iu = 1:n_rays_u, iv = 1:n_rays_v
        u = u_min + (iu-1)*du
        v = v_min + (iv-1)*dv
        # Origin in 3D: start at TOA (z=zmax), horizontal position (u,v)
        # mapped back from (e1,e2) to (x,y,z)
        origins[1, ir] = u
        origins[2, ir] = v
        origins[3, ir] = bounds.zmax
        ir += 1
    end

    return RayGrid(n_rays_u, n_rays_v, n_steps_per_domain,
                   ds, origins, Ω, e1, e2)
end


# ── 2. κ_ext interpolation from unstructured mesh onto ray points ─────────────

"""
    build_mesh_kd_tree(mesh_local)

Build a KD-tree for fast nearest-neighbour lookup of mesh nodes.
Used to interpolate κ_ext from the unstructured mesh onto ray points.

In practice use NearestNeighbors.jl KDTree. For a structured mesh this
degenerates to direct index lookup but the interface remains general.
"""
function build_mesh_kd_tree(mesh_local)
    pts = Matrix{Float64}(undef, 3, mesh_local.npoin)
    for ip = 1:mesh_local.npoin
        pts[1,ip] = mesh_local.x[ip]
        pts[2,ip] = mesh_local.y[ip]
        pts[3,ip] = mesh_local.z[ip]
    end
    return KDTree(pts)
end

"""
    interp_kappa_to_point(tree, κ_ext, x, y, z, k_neighbours=8)

Interpolate κ_ext from the unstructured mesh to an arbitrary point (x,y,z)
using inverse-distance weighting from the k nearest mesh nodes.

Returns 0.0 if the point is outside the local subdomain (all neighbours
are far away), allowing the ray integration to skip points not owned by
this rank.
"""
function interp_kappa_to_point(tree, κ_ext, mesh_local, x, y, z;
                                 k_neighbours=8, max_dist=Inf)
    idxs, dists = knn(tree, [x, y, z], k_neighbours, true)

    # Check if point is within this rank's subdomain
    if dists[1] > max_dist
        return 0.0, false   # point not owned by this rank
    end

    # Inverse distance weighting
    if dists[1] < 1e-10
        return κ_ext[idxs[1]], true   # exactly on a node
    end

    w_sum = 0.0
    κ_sum = 0.0
    for (idx, d) in zip(idxs, dists)
        w      = 1.0 / d^2
        w_sum += w
        κ_sum += w * κ_ext[idx]
    end
    return κ_sum / w_sum, true
end


# ── 3. Parallel ray integration ───────────────────────────────────────────────

"""
    compute_tau_ray_parallel(ray_grid, mesh_local, κ_ext, sw, comm)

Compute slant-path optical depth τ at every mesh node via parallel
ray marching.

Algorithm:
1. Each rank owns a subdomain. For each ray, determine which steps fall
   within this rank's subdomain (bounding box check).
2. Integrate κ_ext along those steps using the local mesh interpolation.
3. Communicate cumulative τ along beam ordering via MPI point-to-point.
4. Interpolate τ from ray grid back onto mesh nodes.

This is general for any mesh decomposition — structured or unstructured.
"""
function compute_tau_ray_parallel(ray_grid, mesh_local, κ_ext, sw, comm)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    n_rays = ray_grid.n_rays_u * ray_grid.n_rays_v

    # ── Build KD-tree for κ interpolation ────────────────────────────────────
    tree = build_mesh_kd_tree(mesh_local)

    # Max distance for a point to be considered "owned" by this rank:
    # use half the minimum element size as a conservative threshold
    min_elem_size = estimate_min_element_size(mesh_local)
    max_dist      = 2.0 * min_elem_size

    # ── Determine beam rank ordering ─────────────────────────────────────────
    _, my_pos, upbeam_rank, downbeam_rank =
        compute_beam_rank_ordering(mesh_local, sw, comm)

    # ── For each ray: integrate κ_ext over steps owned by this rank ──────────
    # τ_local[i_ray, i_step] = cumulative Δτ from start of this rank's
    #                           segment to step i_step on ray i_ray
    # We store per-step values so we can interpolate τ back to mesh nodes.

    # First pass: determine which steps each rank owns for each ray
    # A step is "owned" by this rank if the step point is within max_dist
    # of any local mesh node.

    # Allocate: for each ray store (step_index, x, y, z, κ_at_step)
    # Only for steps owned by this rank.
    ray_segments = Vector{Vector{NTuple{5,Float64}}}(undef, n_rays)
    for ir = 1:n_rays
        ray_segments[ir] = NTuple{5,Float64}[]
    end

    Ωx, Ωy, Ωz = ray_grid.Ω_sun

    for ir = 1:n_rays
        ox = ray_grid.origins[1, ir]
        oy = ray_grid.origins[2, ir]
        oz = ray_grid.origins[3, ir]

        for is = 0:ray_grid.n_steps-1
            # Position along ray
            px = ox + is*ray_grid.ds*Ωx
            py = oy + is*ray_grid.ds*Ωy
            pz = oz + is*ray_grid.ds*Ωz

            # Skip if outside global domain
            (pz < mesh_local.zmin || pz > mesh_local.zmax) && continue
            (px < mesh_local.xmin || px > mesh_local.xmax) && continue
            (py < mesh_local.ymin || py > mesh_local.ymax) && continue

            κ_val, owned = interp_kappa_to_point(tree, κ_ext, mesh_local,
                                                   px, py, pz;
                                                   max_dist=max_dist)
            if owned
                push!(ray_segments[ir], (Float64(is), px, py, pz, κ_val))
            end
        end
    end

    # ── Compute local Δτ per ray (trapezoidal integration) ───────────────────
    # τ_local_exit[ir] = total Δτ along this rank's segment of ray ir
    # τ_local_profile[ir] = vector of (step_index, cumulative_Δτ) pairs
    τ_local_exit    = zeros(Float64, n_rays)
    τ_local_profile = Vector{Vector{Tuple{Float64,Float64}}}(undef, n_rays)

    for ir = 1:n_rays
        seg = ray_segments[ir]
        τ_local_profile[ir] = Tuple{Float64,Float64}[]
        isempty(seg) && continue

        # Sort by step index to ensure correct integration order
        sort!(seg, by=s->s[1])

        τ_acc = 0.0
        push!(τ_local_profile[ir], (seg[1][1], 0.0))
        for i = 2:length(seg)
            κ_prev = seg[i-1][5]
            κ_curr = seg[i][5]
            τ_acc += 0.5*(κ_prev + κ_curr) * ray_grid.ds
            push!(τ_local_profile[ir], (seg[i][1], τ_acc))
        end
        τ_local_exit[ir] = τ_acc
    end

    # ── Parallel prefix sum along beam ordering ───────────────────────────────
    # Each rank receives τ_incoming[ir] from its upbeam neighbour,
    # representing the cumulative τ from TOA to the top of this rank's
    # segment for each ray.
    τ_incoming = zeros(Float64, n_rays)

    if upbeam_rank != MPI.PROC_NULL
        buf = zeros(Float64, n_rays)
        MPI.Recv!(buf, upbeam_rank, 42, comm)
        τ_incoming .= buf
    end

    # Send τ_incoming + τ_local_exit to downbeam rank
    if downbeam_rank != MPI.PROC_NULL
        MPI.Send(τ_incoming .+ τ_local_exit, downbeam_rank, 42, comm)
    end

    # ── Total τ at each ray step owned by this rank ───────────────────────────
    # τ_total at step is = τ_incoming[ir] + τ_local at that step
    # Build a lookup: (ir, step_index) → τ_total
    τ_step_lookup = Dict{Tuple{Int,Int}, Float64}()
    for ir = 1:n_rays
        for (step_idx, τ_loc) in τ_local_profile[ir]
            τ_step_lookup[(ir, round(Int, step_idx))] = τ_incoming[ir] + τ_loc
        end
    end

    # ── Interpolate τ from ray steps back onto mesh nodes ────────────────────
    τ_nodes = zeros(Float64, mesh_local.npoin)

    # For each mesh node find the nearest ray and interpolate τ
    # using the two bounding steps along that ray
    for ip = 1:mesh_local.npoin
        px = mesh_local.x[ip]
        py = mesh_local.y[ip]
        pz = mesh_local.z[ip]

        # Find which ray this node is closest to
        # Project node onto the plane perpendicular to Ω_sun
        # The ray passing through a point (px,py,pz) has its TOA origin at:
        # t = (zmax - pz) / Ωz   (parameter to reach TOA from this point)
        t_to_toa = (mesh_local.zmax - pz) / Ωz   # Ωz < 0 so t > 0
        ox_ray   = px - t_to_toa * Ωx
        oy_ray   = py - t_to_toa * Ωy
        # oz_ray = zmax by construction

        # Find nearest ray origin
        best_ir   = 1
        best_dist = Inf
        for ir = 1:n_rays
            dx = ray_grid.origins[1,ir] - ox_ray
            dy = ray_grid.origins[2,ir] - oy_ray
            d  = dx^2 + dy^2
            if d < best_dist
                best_dist = d
                best_ir   = ir
            end
        end

        # Find the step along best_ir that brackets this node's z
        # Step parameter: s = distance from TOA along beam
        s_node    = -t_to_toa   # negative because t_to_toa goes upward
        step_node = s_node / ray_grid.ds

        # Interpolate between bounding steps
        step_lo = floor(Int, step_node)
        step_hi = step_lo + 1
        frac    = step_node - step_lo

        τ_lo = get(τ_step_lookup, (best_ir, step_lo), -1.0)
        τ_hi = get(τ_step_lookup, (best_ir, step_hi), -1.0)

        if τ_lo >= 0.0 && τ_hi >= 0.0
            τ_nodes[ip] = τ_lo + frac*(τ_hi - τ_lo)
        elseif τ_lo >= 0.0
            τ_nodes[ip] = τ_lo
        elseif τ_hi >= 0.0
            τ_nodes[ip] = τ_hi
        else
            # Neither step found — this node is at a rank boundary.
            # Use τ_incoming for this ray as a lower bound.
            τ_nodes[ip] = τ_incoming[best_ir]
        end
    end

    # ── Compute F_dir from τ ─────────────────────────────────────────────────
    F_dir = sw.S₀_flux .* sw.μ₀ .* exp.(-τ_nodes)

    @info "Rank $rank: τ ∈ [$(round(minimum(τ_nodes),digits=4)), " *
          "$(round(maximum(τ_nodes),digits=4))], " *
          "F_dir ∈ [$(round(minimum(F_dir),digits=2)), " *
          "$(round(maximum(F_dir),digits=2))] W/m²"

    return F_dir, τ_nodes
end


# ── 4. Helper functions ───────────────────────────────────────────────────────

function estimate_min_element_size(mesh_local)
    # Estimate minimum element size from the range of coordinates
    # For a structured mesh this is exact; for unstructured it is approximate
    Lx = maximum(mesh_local.x)- minimum(mesh_local.x)
    Ly = maximum(mesh_local.y) - minimum(mesh_local.y)
    Lz = maximum(mesh_local.z) - minimum(mesh_local.z)
    # Assume roughly cubic elements
    n_est = mesh_local.npoin^(1/3)
    return min(Lx, Ly, Lz) / max(n_est, 1)
end

"""
    compute_direct_radiation(F_dir, τ_nodes, κ_ext, σ, sw, mesh_local)

Compute all direct beam quantities needed for the diffuse RTE source term
and heating rate computation.

# Returns
- `G_dir[ip]`      : mean intensity from direct beam = F_dir/μ₀ (W/m²)
                     This is ∫ I_dir dΩ. Since I_dir is a delta function
                     at Ω_sun, G_dir = I_dir_peak × solid_angle = F_dir/μ₀.
- `S_dir[ip, ip_ang]`: scattering source at each (spatial, angular) node
                     = σ(ip) × Φ(Ω, Ω_sun) × F_dir(ip)/μ₀
                     This is the RHS term for the diffuse RTE.
- `Q_dir[ip]`      : direct beam heating rate contribution (W/m³)
                     = κ_abs(ip) × G_dir(ip)
"""
function compute_direct_radiation(F_dir, τ_nodes, κ_abs, σ, sw, mesh_local)

    npoin     = mesh_local.npoin
    

    # ── Direct beam mean intensity ────────────────────────────────────────────
    # I_dir is a delta function at Ω_sun in angle space.
    # Its contribution to G = ∫ I dΩ is:
    #   G_dir = ∫ I_dir(Ω) dΩ = F_dir / μ₀
    # because F_dir = ∫ I_dir |cosθ| dΩ = I_dir_peak × μ₀ × solid_angle_of_beam
    # and for a delta function beam, solid_angle × I_dir_peak = F_dir/μ₀.
    G_dir = F_dir ./ sw.μ₀

    # ── Direct beam heating rate ──────────────────────────────────────────────
    # Only absorption removes energy from the beam — scattering redirects it
    # into the diffuse field (handled via S_dir below) but does not heat.
    Q_dir = κ_abs .* G_dir   # W/m³

    # ── Scattering source for diffuse RTE ─────────────────────────────────────
    # The diffuse RTE has an extra source term from photons scattered out of
    # the direct beam into direction Ω:
    #
    #   S_dir(x, Ω) = σ(x) × Φ(Ω, Ω_sun) × I_dir(x)
    #               = σ(x) × Φ(Ω, Ω_sun) × F_dir(x) / μ₀
    #
    # This is evaluated at every (spatial node ip, angular node ip_ang) pair.
    # For the assembly it is added to RHS[ip_g] where ip_g = (ip-1)*npoin_ang + ip_ang.

    θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))
    φ_sun = sw.φ₀

    

    # S_dir is assembled into the RHS during the BC/RHS loop.
    # We return Φ_sun and F_dir/μ₀ separately so the assembly can compute
    # S_dir[ip_g] = σ[ip] × Φ_sun[ip_ang] × G_dir[ip]
    # without allocating the full (npoin × npoin_ang) array.

    @info "Direct beam radiation:"
    @info "  G_dir range : [$(round(minimum(G_dir),digits=2)), " *
          "$(round(maximum(G_dir),digits=2))] W/m²"
    @info "  Q_dir range : [$(round(minimum(Q_dir),sigdigits=4)), " *
          "$(round(maximum(Q_dir),sigdigits=4))] W/m³"

    return G_dir, Q_dir
end

function compute_beam_rank_ordering(mesh_local, sw, comm)
    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    x_c = sum(mesh_local.x) / mesh_local.npoin
    y_c = sum(mesh_local.y) / mesh_local.npoin
    z_c = sum(mesh_local.z) / mesh_local.npoin

    θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))
    Ωx    =  sin(θ_sun)*cos(sw.φ₀)
    Ωy    =  sin(θ_sun)*sin(sw.φ₀)
    Ωz    =  cos(θ_sun)

    proj = x_c*Ωx + y_c*Ωy + z_c*Ωz

    # Gather projections from all ranks
    # all_projs[i] = projection of rank (i-1)  (1-based indexing, 0-based ranks)
    all_projs = zeros(Float64, nranks)
    MPI.Allgather!(Ref(proj), all_projs, comm)

    # beam_order[i] = rank index (0-based) of the i-th rank in beam order
    # sortperm gives 1-based positions into all_projs, subtract 1 for 0-based ranks
    beam_order = sortperm(all_projs) .- 1   # now contains 0-based rank indices

    # my_pos: 1-based position of this rank in beam_order
    my_pos = findfirst(==(rank), beam_order)

    # Guard against findfirst returning nothing (should not happen but be safe)
    if isnothing(my_pos)
        @error "Rank $rank not found in beam_order=$(beam_order). " *
               "all_projs=$(all_projs)"
        my_pos = 1
    end

    upbeam_rank   = my_pos > 1      ? beam_order[my_pos-1] : MPI.PROC_NULL
    downbeam_rank = my_pos < nranks ? beam_order[my_pos+1] : MPI.PROC_NULL

    @info "Rank $rank: proj=$(round(proj,digits=3)), " *
          "beam_pos=$my_pos/$nranks, " *
          "upbeam=$upbeam_rank, downbeam=$downbeam_rank"

    return beam_order, my_pos, upbeam_rank, downbeam_rank
end

function compute_tau_direct_1D(mesh_local, κ_ext, sw, comm)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    ngl    = mesh_local.ngl

    zmin_global = MPI.Allreduce(minimum(mesh_local.z), MPI.MIN, comm)
    zmax_global = MPI.Allreduce(maximum(mesh_local.z), MPI.MAX, comm)
    tol_z       = (zmax_global - zmin_global) * 1e-8

    τ_nodes = zeros(Float64, mesh_local.npoin)

    x_size = abs(mesh_local.x[mesh_local.connijk[1,1,1,1]] -
                 mesh_local.x[mesh_local.connijk[1,mesh_local.ngl,1,1]])
    y_size = abs(mesh_local.y[mesh_local.connijk[1,1,1,1]] -
                 mesh_local.y[mesh_local.connijk[1,1,mesh_local.ngl,1]])
    tol_xy = sqrt(x_size^2 + y_size^2) * 1e-4

    # ── Step 1: build column map ──────────────────────────────────────────────
    col_map = Dict{Tuple{Float64,Float64}, Vector{Int}}()

    for iel = 1:mesh_local.nelem
        for i = 1:ngl, j = 1:ngl
            ip_ref = mesh_local.connijk[iel, i, j, 1]
            x_ref  = mesh_local.x[ip_ref]
            y_ref  = mesh_local.y[ip_ref]
            key    = (round(x_ref/tol_xy)*tol_xy,
                      round(y_ref/tol_xy)*tol_xy)
            for k = 1:ngl
                ip    = mesh_local.connijk[iel, i, j, k]
                nodes = get!(col_map, key, Int[])
                ip ∉ nodes && push!(nodes, ip)
            end
        end
    end

    @info "Rank $rank: found $(length(col_map)) distinct (x,y) columns"

    # ── Step 2: integrate τ locally within each column ───────────────────────
    # col_sorted[key] = (z_u, ip_u, τ_col) — deduplicated, z descending.
    # τ_nodes[ip] = local cumulative τ from the top of THIS rank's subdomain.

    col_sorted = Dict{Tuple{Float64,Float64},
                      Tuple{Vector{Float64}, Vector{Int}, Vector{Float64}}}()

    for (key, nodes) in col_map
        z_vals   = mesh_local.z[nodes]
        sort_idx = sortperm(z_vals, rev=true)
        nodes_s  = nodes[sort_idx]
        z_s      = z_vals[sort_idx]

        z_u  = Float64[]
        κ_u  = Float64[]
        ip_u = Int[]
        for (ip, zv) in zip(nodes_s, z_s)
            if !isempty(z_u) && abs(zv - z_u[end]) < tol_z
                κ_u[end] = 0.5*(κ_u[end] + κ_ext[ip])
            else
                push!(z_u, zv); push!(κ_u, κ_ext[ip]); push!(ip_u, ip)
            end
        end

        n     = length(z_u)
        τ_col = zeros(Float64, n)
        for i = 2:n
            dz       = z_u[i-1] - z_u[i]
            τ_col[i] = τ_col[i-1] + 0.5*(κ_u[i-1] + κ_u[i]) * dz / sw.μ₀
        end

        # Assign τ to deduplicated nodes
        for (i, ip) in enumerate(ip_u)
            τ_nodes[ip] = τ_col[i]
        end
        # Assign τ to z-duplicate nodes merged during deduplication
        for ip in nodes
            zv  = mesh_local.z[ip]
            idx = findfirst(z -> abs(z - zv) < tol_z, z_u)
            !isnothing(idx) && (τ_nodes[ip] = τ_col[idx])
        end

        col_sorted[key] = (z_u, ip_u, τ_col)
    end

    # ── Step 3: node-level iterative reduction ────────────────────────────────
    # Broadcast τ at every node via gip using Allreduce MAX. For L-shaped or
    # staircase rank boundaries, the shared interface node between two ranks
    # may be an interior node of one rank's column — not its top or bottom.
    # Broadcasting at every node ensures the correct upbeam τ is found
    # regardless of where the rank boundary cuts through the column.
    #
    # Each iteration:
    #   1. Each rank writes its current τ_nodes[ip] into τ_global[gip]
    #   2. Allreduce MAX — upbeam ranks have larger τ at shared nodes
    #   3. Each rank scans its columns top-to-bottom. At the first node
    #      where τ_global[gip] > τ_nodes[ip], compute δτ and shift that
    #      node and all nodes below it in the column.
    #   4. break after first shift per column — next iter handles deeper nodes
    #
    # Converges in at most nranks iterations (one rank boundary per iter).

    if nranks > 1

        gnpoin    = mesh_local.gnpoin
        max_iters = nranks

        for iter = 1:max_iters

            τ_global = zeros(Float64, gnpoin)
            for ip = 1:mesh_local.npoin
                gip = mesh_local.ip2gip[ip]
                τ_global[gip] = max(τ_global[gip], τ_nodes[ip])
            end
            MPI.Allreduce!(τ_global, MPI.MAX, comm)

            # Compute δτ per node directly — each node gets corrected exactly once.
            # δτ[ip] = how much to add to this node based on the largest δτ
            # from any node above it in its column.
            δτ_node = zeros(Float64, mesh_local.npoin)

            for (key, nodes) in col_map
                z_u, ip_u, _ = col_sorted[key]
                n = length(ip_u)

                # Find the topmost node in this column with a discrepancy
                # and propagate its δτ downward
                δτ_col = 0.0
                for i = 1:n
                    ip  = ip_u[i]
                    gip = mesh_local.ip2gip[ip]
                    δτ_from_global = τ_global[gip] - τ_nodes[ip]

                    if δτ_from_global > δτ_col + tol_z
                        # New larger discrepancy found — this is the correction
                        # that all nodes from here downward need
                        δτ_col = δτ_from_global
                    end

                    # The correction for this node is the running maximum δτ
                    # from all nodes at or above it in the column
                    δτ_node[ip] = max(δτ_node[ip], δτ_col)
                end

                # Also apply to z-duplicate nodes
                for ip_dup in nodes
                    zv  = mesh_local.z[ip_dup]
                    idx = findfirst(z -> abs(z - zv) < tol_z, z_u)
                    if !isnothing(idx)
                        δτ_node[ip_dup] = max(δτ_node[ip_dup], δτ_node[ip_u[idx]])
                    end
                end
            end

            # Apply corrections — each node updated exactly once
            local_updated = false
            for ip = 1:mesh_local.npoin
                if δτ_node[ip] > tol_z
                    τ_nodes[ip] += δτ_node[ip]
                    local_updated = true
                end
            end

            global_updated = MPI.Allreduce(local_updated, MPI.LOR, comm)
            @info "Rank $rank: τ propagation iter $iter, local_updated=$local_updated"
            !global_updated && break
        end

    end  # nranks > 1

    # ── Diagnostic ────────────────────────────────────────────────────────────
    n_violations = 0
    for (key, nodes) in col_map
        z_u, ip_u, _ = col_sorted[key]
        for i = 2:length(ip_u)
            if τ_nodes[ip_u[i]] < τ_nodes[ip_u[i-1]] - tol_z
                n_violations += 1
                break
            end
        end
    end
    @info "Rank $rank: τ monotonicity violations: " *
          "$n_violations / $(length(col_map)) columns"

    ip_top_g = argmax(mesh_local.z)
    ip_bot_g = argmin(mesh_local.z)
    @info "Rank $rank: " *
          "z_top=$(round(mesh_local.z[ip_top_g], digits=1)) m  " *
          "τ_top=$(round(τ_nodes[ip_top_g], digits=6))  |  " *
          "z_bot=$(round(mesh_local.z[ip_bot_g], digits=1)) m  " *
          "τ_bot=$(round(τ_nodes[ip_bot_g], digits=6))"

    # ── Step 4: compute F_dir ─────────────────────────────────────────────────
    F_dir = sw.S₀_flux .* sw.μ₀ .* exp.(-τ_nodes)

    @info "Rank $rank: τ ∈ [$(round(minimum(τ_nodes), digits=4)), " *
          "$(round(maximum(τ_nodes), digits=4))], " *
          "F_dir ∈ [$(round(minimum(F_dir), digits=2)), " *
          "$(round(maximum(F_dir), digits=2))] W/m²"

    return F_dir, τ_nodes
end

function compute_tau_direct_3D(mesh_local, κ_ext, sw, comm)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    ngl    = mesh_local.ngl

    zmin_global = MPI.Allreduce(minimum(mesh_local.z), MPI.MIN, comm)
    zmax_global = MPI.Allreduce(maximum(mesh_local.z), MPI.MAX, comm)
    tol_z       = (zmax_global - zmin_global) * 1e-8

    # Domain extents for periodic wrapping of slant path
    xmin_global = MPI.Allreduce(minimum(mesh_local.x), MPI.MIN, comm)
    xmax_global = MPI.Allreduce(maximum(mesh_local.x), MPI.MAX, comm)
    ymin_global = MPI.Allreduce(minimum(mesh_local.y), MPI.MIN, comm)
    ymax_global = MPI.Allreduce(maximum(mesh_local.y), MPI.MAX, comm)
    Lx = xmax_global - xmin_global
    Ly = ymax_global - ymin_global

    # Sun direction — horizontal displacement per unit z ascent
    θ_sun = acos(clamp(sw.μ₀, 0.0, 1.0))
    dx_dz = sin(θ_sun) * cos(sw.φ₀) / sw.μ₀   # dx per dz along slant path
    dy_dz = sin(θ_sun) * sin(sw.φ₀) / sw.μ₀   # dy per dz along slant path

    τ_nodes = zeros(Float64, mesh_local.npoin)

    # ── Build global κ_ext field via Allgatherv ───────────────────────────────
    # Each node's slant path passes through (x',y') positions that may be
    # owned by other ranks. We need κ_ext at arbitrary (x,y,z) positions.
    # Solution: gather all node positions and κ_ext values globally, then
    # interpolate along the slant path using nearest-node lookup.

    local_n = mesh_local.npoin
    all_n   = MPI.Allgather(local_n, comm)
    total_n = sum(all_n)

    global_x   = zeros(Float64, total_n)
    global_y   = zeros(Float64, total_n)
    global_z   = zeros(Float64, total_n)
    global_κ   = zeros(Float64, total_n)

    MPI.Allgatherv!(mesh_local.x,  VBuffer(global_x, all_n), comm)
    MPI.Allgatherv!(mesh_local.y,  VBuffer(global_y, all_n), comm)
    MPI.Allgatherv!(mesh_local.z,  VBuffer(global_z, all_n), comm)
    MPI.Allgatherv!(κ_ext,         VBuffer(global_κ, all_n), comm)

    # ── Build vertical level structure ────────────────────────────────────────
    # For interpolation along the slant path we need to know all unique z levels.
    # At each z level we have a 2D field of κ_ext that we can interpolate in x,y.
    z_levels = sort(unique(round.(global_z, digits=6)), rev=true)  # TOA first
    nz       = length(z_levels)

    # Group global nodes by z level for fast 2D lookup
    # z_level_nodes[iz] = indices into global arrays at that z level
    z_level_nodes = [Int[] for _ = 1:nz]
    z_tol_level   = tol_z * 100   # slightly looser for level matching
    for ig = 1:total_n
        iz = findfirst(z -> abs(z - global_z[ig]) < z_tol_level, z_levels)
        !isnothing(iz) && push!(z_level_nodes[iz], ig)
    end

    # ── Integrate τ along slant path for each local node ─────────────────────
    for ip = 1:mesh_local.npoin
        x0 = mesh_local.x[ip]
        y0 = mesh_local.y[ip]
        z0 = mesh_local.z[ip]

        τ = 0.0
        iz_start = findfirst(z -> abs(z - z0) < z_tol_level, z_levels)
        isnothing(iz_start) && continue

        # Walk from this node's z level up to TOA
        # At each level pair (iz, iz-1), integrate κ along the slant path
        # using trapezoidal rule. The x,y position at each level is:
        #   x(z') = x0 + (z' - z0) * dx_dz
        #   y(z') = y0 + (z' - z0) * dy_dz
        # with periodic wrapping in x and y.

        κ_prev = interpolate_κ_at_level(x0, y0, iz_start,
                                          z_level_nodes, global_x, global_y,
                                          global_κ, xmin_global, Lx,
                                          ymin_global, Ly)

        for iz = iz_start-1 : -1 : 1   # walk upward toward TOA
            z_cur  = z_levels[iz]
            dz     = z_cur - z_levels[iz+1]   # positive upward step

            # Slant path position at this z level
            x_cur  = x0 + (z_cur - z0) * dx_dz
            y_cur  = y0 + (z_cur - z0) * dy_dz

            # Periodic wrapping
            x_cur  = xmin_global + mod(x_cur - xmin_global, Lx)
            y_cur  = ymin_global + mod(y_cur - ymin_global, Ly)

            κ_cur  = interpolate_κ_at_level(x_cur, y_cur, iz,
                                              z_level_nodes, global_x, global_y,
                                              global_κ, xmin_global, Lx,
                                              ymin_global, Ly)

            # Trapezoidal integration — dτ = κ * ds = κ * dz / μ₀
            τ += 0.5 * (κ_prev + κ_cur) * dz / sw.μ₀

            κ_prev = κ_cur
        end

        τ_nodes[ip] = τ
    end

    F_dir = sw.S₀_flux .* sw.μ₀ .* exp.(-τ_nodes)

    @info "Rank $rank: τ ∈ [$(round(minimum(τ_nodes),digits=4)), " *
          "$(round(maximum(τ_nodes),digits=4))], " *
          "F_dir ∈ [$(round(minimum(F_dir),digits=2)), " *
          "$(round(maximum(F_dir),digits=2))] W/m²"

    return F_dir, τ_nodes
end

# ── 2D interpolation at a given z level ──────────────────────────────────────
# Finds the nearest node at this z level to (x_query, y_query) and returns
# its κ value. For a structured mesh this is exact (query point coincides
# with a node). For an unstructured mesh this is nearest-node interpolation.
# Periodic wrapping is already applied to x_query, y_query before calling.

function interpolate_κ_at_level(x_q, y_q, iz,
                                  z_level_nodes, global_x, global_y,
                                  global_κ, xmin, Lx, ymin, Ly)
    best_dist = Inf
    best_κ    = 0.0
    for ig in z_level_nodes[iz]
        # Distance with periodic wrapping
        dx   = abs(global_x[ig] - x_q)
        dy   = abs(global_y[ig] - y_q)
        dx   = min(dx, Lx - dx)
        dy   = min(dy, Ly - dy)
        dist = dx^2 + dy^2
        if dist < best_dist
            best_dist = dist
            best_κ    = global_κ[ig]
        end
    end
    return best_κ
end
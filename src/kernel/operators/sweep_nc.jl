using Graphs

function test_spatial_ncf_interpolation(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    # ── Apply the single-element spatial refinement, exactly as given ────────
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl  = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, $(mesh.num_hanging_facets) hanging facets, " *
            "num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets produced — nothing to test.")
        return nothing
    end

    debug_ncf_line_coords(sem.mesh, 3, sem.mesh.ngl)

    # ── Take the first NCF entry ──────────────────────────────────────────────
    ncf_idx = 3
    child_elem_id  = Int(mesh.cip[ncf_idx])
    parent_elem_id = Int(mesh.pip[ncf_idx])
    local_facet_id = Int(mesh.lfid[ncf_idx])
    IPc_nodes = mesh.IPc_list[:, ncf_idx]
    IPp_nodes = mesh.IPp_list[:, ncf_idx]
    println("Testing NCF #$ncf_idx: child=$child_elem_id parent=$parent_elem_id fid=$local_facet_id")

    L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
        child_elem_id, parent_elem_id, local_facet_id, IPc_nodes, IPp_nodes, mesh, ngl)

    # ── Independently re-derive the SAME line extraction the builder used,
    # so we know EXACTLY which physical nodes L1/L2 correspond to — rather
    # than trust any downstream consumer's indexing convention. ─────────────
    parent_line_1 = Int[]; parent_line_2 = Int[]
    child_line_1  = Int[]; child_line_2  = Int[]
    for i = 1:ngl
        idx = (i - 1) * ngl + 1
        idx <= length(IPp_nodes) || continue
        push!(parent_line_1, Int(IPp_nodes[idx])); push!(child_line_1, Int(IPc_nodes[idx]))
    end
    for j = 1:ngl
        idx = j
        idx <= length(IPp_nodes) || continue
        push!(parent_line_2, Int(IPp_nodes[idx])); push!(child_line_2, Int(IPc_nodes[idx]))
    end

    # ── For each line, empirically determine which physical coordinate
    # actually VARIES (don't assume x/y/z from local_facet_id at all) ────────
    function varying_coord(nodes)
        xs = [mesh.x[n] for n in nodes]; ys = [mesh.y[n] for n in nodes]; zs = [mesh.z[n] for n in nodes]
        spread(v) = maximum(v) - minimum(v)
        sx, sy, sz = spread(xs), spread(ys), spread(zs)
        m = max(sx, sy, sz)
        m == sx && return xs
        m == sy && return ys
        return zs
    end

    # ── Closed-form 1D test: polynomial of degree ngl-1, EXACTLY representable
    # by ngl-node Lagrange interpolation — expect machine-precision agreement. ──
    poly(t) = sum(0.3^k * ((t - 9000.0)/1000.0)^k for k = 0:ngl-1)

    function test_line(L, parent_nodes, child_nodes, label)
        pc = varying_coord(parent_nodes)
        cc = varying_coord(child_nodes)
        parent_vals = poly.(pc)
        child_computed = L * parent_vals
        child_exact = poly.(cc)
        max_err = maximum(abs.(child_computed .- child_exact))
        println("  $label: max abs error = $max_err  " * (max_err < 1e-9 ? "✓ PASS" : "✗ FAIL"))
        return max_err
    end

    e1 = test_line(L1, parent_line_1, child_line_1, "L1 (line 1)")
    e2 = test_line(L2, parent_line_2, child_line_2, "L2 (line 2)")

    return max(e1, e2)
end

function debug_ncf_line_coords(mesh, ncf_idx::Int, ngl::Int)
    child_elem_id  = Int(mesh.cip[ncf_idx])
    parent_elem_id = Int(mesh.pip[ncf_idx])
    local_facet_id = Int(mesh.lfid[ncf_idx])
    IPc_nodes = mesh.IPc_list[:, ncf_idx]
    IPp_nodes = mesh.IPp_list[:, ncf_idx]

    println("NCF #$ncf_idx: child=$child_elem_id parent=$parent_elem_id fid=$local_facet_id")
    println("─"^70)

    parent_line_1 = Int[]; parent_line_2 = Int[]
    child_line_1  = Int[]; child_line_2  = Int[]
    for i = 1:ngl
        idx = (i - 1) * ngl + 1
        idx <= length(IPp_nodes) || continue
        push!(parent_line_1, Int(IPp_nodes[idx])); push!(child_line_1, Int(IPc_nodes[idx]))
    end
    for j = 1:ngl
        idx = j
        idx <= length(IPp_nodes) || continue
        push!(parent_line_2, Int(IPp_nodes[idx])); push!(child_line_2, Int(IPc_nodes[idx]))
    end

    println("parent_line_1 (idx=(i-1)*ngl+1, i=1:$ngl):")
    for n in parent_line_1
        println("  node=$n  (x,y,z)=($(mesh.x[n]), $(mesh.y[n]), $(mesh.z[n]))")
    end
    println("child_line_1:")
    for n in child_line_1
        println("  node=$n  (x,y,z)=($(mesh.x[n]), $(mesh.y[n]), $(mesh.z[n]))")
    end
    println("parent_line_2 (idx=j, j=1:$ngl):")
    for n in parent_line_2
        println("  node=$n  (x,y,z)=($(mesh.x[n]), $(mesh.y[n]), $(mesh.z[n]))")
    end
    println("child_line_2:")
    for n in child_line_2
        println("  node=$n  (x,y,z)=($(mesh.x[n]), $(mesh.y[n]), $(mesh.z[n]))")
    end

    # Also print the FULL face node sets (all ngl² each), for context —
    # are parent/child face extents even the relationship we'd expect
    # (child ⊂ parent, roughly half the extent in each direction)?
    println("─"^70)
    xs_p = [mesh.x[Int(n)] for n in IPp_nodes if n > 0]
    ys_p = [mesh.y[Int(n)] for n in IPp_nodes if n > 0]
    zs_p = [mesh.z[Int(n)] for n in IPp_nodes if n > 0]
    xs_c = [mesh.x[Int(n)] for n in IPc_nodes if n > 0]
    ys_c = [mesh.y[Int(n)] for n in IPc_nodes if n > 0]
    zs_c = [mesh.z[Int(n)] for n in IPc_nodes if n > 0]
    println("Full parent face extent: x∈[$(minimum(xs_p)),$(maximum(xs_p))] " *
            "y∈[$(minimum(ys_p)),$(maximum(ys_p))] z∈[$(minimum(zs_p)),$(maximum(zs_p))]")
    println("Full child  face extent: x∈[$(minimum(xs_c)),$(maximum(xs_c))] " *
            "y∈[$(minimum(ys_c)),$(maximum(ys_c))] z∈[$(minimum(zs_c)),$(maximum(zs_c))]")
end

"""
    build_dg_ncf_parent_to_child(mesh, connijk_dg, x_dg,
                                  child_iel, parent_iel, local_facet_id,
                                  IPc_nodes, IPp_nodes, L1, L2, ngl)

Evaluates the PARENT (coarse) element's own already-solved DG trace at each
of the CHILD (fine) element's face nodes, via the L1⊗L2 tensor-product
interpolation weights. Reuses the EXACT index conventions already
established in the CG hanging-node code (w_1=L1[j,:], w_2=L2[i,:],
parent_face_idx=(p2-1)*ngl+p1) — not a new convention, so it stays
consistent with the already-verified L1/L2 construction.

Returns a Vector{Float64} of length(IPc_nodes), one interpolated value per
child face node, in the SAME node_idx ordering as IPc_nodes/IPp_nodes.
"""
function build_dg_ncf_parent_to_child(mesh, connijk_dg, x_dg,
                                       child_iel::Int, parent_iel::Int, local_facet_id::Int,
                                       IPc_nodes, IPp_nodes, L1::Matrix{Float64}, L2::Matrix{Float64}, ngl::Int)
    # Reverse lookup: parent's face node (raw mesh.connijk index) -> (i,j,k)
    # within parent_iel, needed to index connijk_dg/x_dg correctly.
    parent_ijk = Dict{Int,NTuple{3,Int}}()
    for k=1:ngl, j=1:ngl, i=1:ngl
        parent_ijk[mesh.connijk[parent_iel,i,j,k]] = (i,j,k)
    end

    # Gather parent's face DG values into an ngl×ngl grid, using the SAME
    # (p1,p2) convention as build_hanging_node_constraints_with_angular:
    # parent_face_idx = (p2-1)*ngl + p1.
    parent_grid = zeros(ngl, ngl)
    for idx = 1:length(IPp_nodes)
        node = Int(IPp_nodes[idx])
        node <= 0 && continue
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        (i,j,k) = parent_ijk[node]
        parent_grid[p1,p2] = x_dg[connijk_dg[parent_iel,i,j,k]]
    end

    # Evaluate at each child face node, using the SAME node_idx decode and
    # w_1/w_2 assignment as the CG code (w_1=L1[j,:], w_2=L2[i,:]).
    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        node = Int(IPc_nodes[node_idx])
        node <= 0 && continue
        ci = div(node_idx-1, ngl) + 1
        cj = mod(node_idx-1, ngl) + 1
        (ci < 1 || ci > size(L1,1) || cj < 1 || cj > size(L2,1)) && continue
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += w_1[p1]*w_2[p2]*parent_grid[p1,p2]
        end
        child_vals[node_idx] = val
    end
    return child_vals
end

function test_dg_ncf_parent_to_child(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    # ── Apply the single-element spatial refinement ──────────────────────────
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl  = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, $(mesh.num_hanging_facets) hanging facets, " *
            "num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets produced — nothing to test.")
        return nothing
    end

    # ── Build DG infrastructure on the refined mesh ───────────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    # ── Take the first NCF entry ───────────────────────────────────────────────
    ncf_idx = 1
    child_iel  = Int(mesh.cip[ncf_idx])
    parent_iel = Int(mesh.pip[ncf_idx])
    local_facet_id = Int(mesh.lfid[ncf_idx])
    IPc_nodes = mesh.IPc_list[:, ncf_idx]
    IPp_nodes = mesh.IPp_list[:, ncf_idx]
    println("Testing NCF #$ncf_idx: child=$child_iel parent=$parent_iel fid=$local_facet_id")

    L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
        child_iel, parent_iel, local_facet_id, IPc_nodes, IPp_nodes, mesh, ngl)

    # ── Well-scaled polynomial field, additive in normalized x,y,z ────────────
    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    x_dg = zeros(npoin_dg)
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        ip = mesh.connijk[iel,i,j,k]
        x_dg[connijk_dg[iel,i,j,k]] = field(mesh.x[ip], mesh.y[ip], mesh.z[ip])
    end

    child_vals = build_dg_ncf_parent_to_child(mesh, connijk_dg, x_dg,
                                               child_iel, parent_iel, local_facet_id,
                                               IPc_nodes, IPp_nodes, L1, L2, ngl)

    max_err = 0.0
    for idx = 1:length(IPc_nodes)
        node = Int(IPc_nodes[idx])
        node <= 0 && continue
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        max_err = max(max_err, abs(child_vals[idx] - exact))
    end
    println("Parent→child NCF interpolation vs exact (well-scaled) field: max abs error = $max_err")
    return max_err
end


function test_dg_ncf_child_to_parent_v3(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl  = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)


    # ── Group NCF entries by (pip, lfid) ───────────────────────────────────────
    groups = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        key = (Int(mesh.pip[ncf_idx]), Int(mesh.lfid[ncf_idx]))
        push!(get!(groups, key, Int[]), ncf_idx)
    end

    valid_group = nothing
    for (key, members) in groups
        if length(members) == 4
            valid_group = (key, members)
            break
        end
    end
    if valid_group === nothing
        println("⚠ No (parent,face) group with exactly 4 children found.")
        return nothing
    end
    (parent_iel, local_facet_id), member_ncf_idxs = valid_group
    println("Testing group: parent=$parent_iel fid=$local_facet_id, NCF indices=$member_ncf_idxs")

    ncf_entries = Tuple{Int,Int,Vector,Vector}[]
    for ncf_idx in member_ncf_idxs
        child_iel = Int(mesh.cip[ncf_idx])
        IPc_nodes = mesh.IPc_list[:, ncf_idx]
        IPp_nodes = mesh.IPp_list[:, ncf_idx]
        push!(ncf_entries, (ncf_idx, child_iel, IPc_nodes, IPp_nodes))
    end

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    x_dg = zeros(npoin_dg)
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        ip = mesh.connijk[iel,i,j,k]
        x_dg[connijk_dg[iel,i,j,k]] = field(mesh.x[ip], mesh.y[ip], mesh.z[ip])
    end

    parent_projected = build_dg_ncf_child_to_parent_v5(mesh, connijk_dg, x_dg,
                                                     parent_iel, local_facet_id,
                                                     ncf_entries, ngl)

    # Parent's own face node list — reuse IPp_nodes from any member entry
    IPp_ref = mesh.IPp_list[:, member_ncf_idxs[1]]
    max_err = 0.0
    for idx = 1:length(IPp_ref)
        node = Int(IPp_ref[idx])
        node <= 0 && continue
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        computed = parent_projected[p1,p2]
        max_err = max(max_err, abs(computed - exact))
    end
    println("Child→parent NCF L2 projection vs exact field AT PARENT NODES: max abs error = $max_err")
    return max_err
end

@inline function axis_of_max_spread(coords_by_axis)
    spreads = [maximum(c)-minimum(c) for c in coords_by_axis]
    return argmax(spreads)
end

function debug_v4_axis_detection(mesh, ncf_entries, ngl)
    for (ncf_idx, child_iel, IPc_nodes, IPp_nodes) in ncf_entries
        child_line_1 = Int[]; child_line_2 = Int[]
        for i = 1:ngl
            idx = (i-1)*ngl + 1; idx <= length(IPc_nodes) || continue
            push!(child_line_1, Int(IPc_nodes[idx]))
        end
        for j = 1:ngl
            idx = j; idx <= length(IPc_nodes) || continue
            push!(child_line_2, Int(IPc_nodes[idx]))
        end
        xs1=[mesh.x[n] for n in child_line_1]; ys1=[mesh.y[n] for n in child_line_1]; zs1=[mesh.z[n] for n in child_line_1]
        xs2=[mesh.x[n] for n in child_line_2]; ys2=[mesh.y[n] for n in child_line_2]; zs2=[mesh.z[n] for n in child_line_2]

        println("NCF $ncf_idx, child=$child_iel:")
        println("  line1 spreads: x=$(maximum(xs1)-minimum(xs1)) y=$(maximum(ys1)-minimum(ys1)) z=$(maximum(zs1)-minimum(zs1))")
        println("  line2 spreads: x=$(maximum(xs2)-minimum(xs2)) y=$(maximum(ys2)-minimum(ys2)) z=$(maximum(zs2)-minimum(zs2))")
        println("  line1 coords: x=$xs1 y=$ys1 z=$zs1")
        println("  line2 coords: x=$xs2 y=$ys2 z=$zs2")
    end
end

function test_v4_axis_detection(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl  = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    # ── Group NCF entries by (pip, lfid), same as before ─────────────────────
    groups = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        key = (Int(mesh.pip[ncf_idx]), Int(mesh.lfid[ncf_idx]))
        push!(get!(groups, key, Int[]), ncf_idx)
    end

    valid_group = nothing
    for (key, members) in groups
        if length(members) == 4
            valid_group = (key, members)
            break
        end
    end
    if valid_group === nothing
        println("⚠ No (parent,face) group with exactly 4 children found.")
        return nothing
    end
    (parent_iel, local_facet_id), member_ncf_idxs = valid_group
    println("Testing group: parent=$parent_iel fid=$local_facet_id, NCF indices=$member_ncf_idxs")

    ncf_entries = Tuple{Int,Int,Vector,Vector}[]
    for ncf_idx in member_ncf_idxs
        child_iel = Int(mesh.cip[ncf_idx])
        IPc_nodes = mesh.IPc_list[:, ncf_idx]
        IPp_nodes = mesh.IPp_list[:, ncf_idx]
        push!(ncf_entries, (ncf_idx, child_iel, IPc_nodes, IPp_nodes))
    end

    # ── Actual diagnostic ──────────────────────────────────────────────────────
    for (ncf_idx, child_iel, IPc_nodes, IPp_nodes) in ncf_entries
        child_line_1 = Int[]; child_line_2 = Int[]
        for i = 1:ngl
            idx = (i-1)*ngl + 1; idx <= length(IPc_nodes) || continue
            push!(child_line_1, Int(IPc_nodes[idx]))
        end
        for j = 1:ngl
            idx = j; idx <= length(IPc_nodes) || continue
            push!(child_line_2, Int(IPc_nodes[idx]))
        end
        xs1=[mesh.x[n] for n in child_line_1]; ys1=[mesh.y[n] for n in child_line_1]; zs1=[mesh.z[n] for n in child_line_1]
        xs2=[mesh.x[n] for n in child_line_2]; ys2=[mesh.y[n] for n in child_line_2]; zs2=[mesh.z[n] for n in child_line_2]

        println("NCF $ncf_idx, child=$child_iel:")
        println("  line1 spreads: x=$(maximum(xs1)-minimum(xs1)) y=$(maximum(ys1)-minimum(ys1)) z=$(maximum(zs1)-minimum(zs1))")
        println("  line2 spreads: x=$(maximum(xs2)-minimum(xs2)) y=$(maximum(ys2)-minimum(ys2)) z=$(maximum(zs2)-minimum(zs2))")
        println("  line1 coords: x=$xs1 y=$ys1 z=$zs1")
        println("  line2 coords: x=$xs2 y=$ys2 z=$zs2")
    end

    return ncf_entries
end

function build_dg_ncf_child_to_parent_v5(mesh, connijk_dg, x_dg,
                                          parent_iel::Int, local_facet_id::Int,
                                          ncf_entries, ngl::Int)
    parent_result = fill(NaN, ngl, ngl)

    for (ncf_idx, child_iel, IPc_nodes, IPp_nodes) in ncf_entries
        child_ijk = Dict{Int,NTuple{3,Int}}()
        for k=1:ngl, j=1:ngl, i=1:ngl
            child_ijk[mesh.connijk[child_iel,i,j,k]] = (i,j,k)
        end

        child_line_1 = Int[]; child_line_2 = Int[]
        for i = 1:ngl
            idx = (i-1)*ngl + 1; idx <= length(IPc_nodes) || continue
            push!(child_line_1, Int(IPc_nodes[idx]))
        end
        for j = 1:ngl
            idx = j; idx <= length(IPc_nodes) || continue
            push!(child_line_2, Int(IPc_nodes[idx]))
        end
        # SWAP: child_line_1 actually varies c2 (grid's SLOW index),
        # child_line_2 actually varies c1 (grid's FAST index) — derived
        # explicitly above from idx=(i-1)*ngl+j (i slow,j fast) vs
        # child_grid's idx=(c2-1)*ngl+c1 (c1 fast,c2 slow): c1=j, c2=i.

        xs2=[mesh.x[n] for n in child_line_2]; ys2=[mesh.y[n] for n in child_line_2]; zs2=[mesh.z[n] for n in child_line_2]
        axis_c1 = axis_of_max_spread((xs2,ys2,zs2))   # child_line_2 -> TRUE c1 axis
        xs1=[mesh.x[n] for n in child_line_1]; ys1=[mesh.y[n] for n in child_line_1]; zs1=[mesh.z[n] for n in child_line_1]
        axis_c2 = axis_of_max_spread((xs1,ys1,zs1))   # child_line_1 -> TRUE c2 axis
        coord(node,axis) = axis==1 ? mesh.x[node] : (axis==2 ? mesh.y[node] : mesh.z[node])

        c1_coords = [coord(n,axis_c1) for n in child_line_2]   # SWAPPED source
        c2_coords = [coord(n,axis_c2) for n in child_line_1]   # SWAPPED source
        c1_min,c1_max = extrema(c1_coords); c2_min,c2_max = extrema(c2_coords)
        tol1 = 1e-6*max(c1_max-c1_min,1.0); tol2 = 1e-6*max(c2_max-c2_min,1.0)

        child_grid = zeros(ngl,ngl)
        for idx = 1:length(IPc_nodes)
            node = Int(IPc_nodes[idx]); node<=0 && continue
            c1 = ((idx-1)%ngl)+1; c2 = ((idx-1)÷ngl)+1
            child_grid[c1,c2] = x_dg[connijk_dg[child_iel, child_ijk[node]...]]
        end

        for idx = 1:length(IPp_nodes)
            node = Int(IPp_nodes[idx]); node<=0 && continue
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            pc1 = coord(node,axis_c1); pc2 = coord(node,axis_c2)
            (pc1 < c1_min-tol1 || pc1 > c1_max+tol1 || pc2 < c2_min-tol2 || pc2 > c2_max+tol2) && continue

            L_c1 = build_1d_lagrange_matrix(c1_coords, [pc1])
            L_c2 = build_1d_lagrange_matrix(c2_coords, [pc2])
            val = 0.0
            for a=1:ngl, b=1:ngl
                val += L_c1[1,a]*L_c2[1,b]*child_grid[a,b]
            end
            parent_result[p1,p2] = val
        end
    end
    return parent_result
end

function test_ncf_local_facet_id_convention(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh
    ngl  = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    n_match_child  = 0
    n_match_parent_via_opp = 0
    n_neither      = 0

    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        stored_lfid = Int(mesh.lfid[ncf_idx])
        IPc_nodes = mesh.IPc_list[:, ncf_idx]
        IPp_nodes = mesh.IPp_list[:, ncf_idx]

        ic_full = Set(Int(n) for n in IPc_nodes if n > 0)
        child_fid_found = 0
        for face_el = 1:6
            my_fid = CONN_FACE_EL_TO_MY_FID[face_el]
            face_nodes = Set{Int}()
            irange, jrange, krange = face_ip_range(my_fid, ngl)
            for k in krange, j in jrange, i in irange
                push!(face_nodes, mesh.connijk[child_iel,i,j,k])
            end
            Set(ic_full) == face_nodes && (child_fid_found = my_fid)
        end

        ip_full = Set(Int(n) for n in IPp_nodes if n > 0)
        parent_fid_found = 0
        for face_el = 1:6
            my_fid = CONN_FACE_EL_TO_MY_FID[face_el]
            face_nodes = Set{Int}()
            irange, jrange, krange = face_ip_range(my_fid, ngl)
            for k in krange, j in jrange, i in irange
                push!(face_nodes, mesh.connijk[parent_iel,i,j,k])
            end
            Set(ip_full) == face_nodes && (parent_fid_found = my_fid)
        end

        matches_child  = (stored_lfid == child_fid_found)
        matches_parent_via_opp = (parent_fid_found != 0 && OPP_FACE[stored_lfid] == parent_fid_found)

        matches_child && (n_match_child += 1)
        matches_parent_via_opp && (n_match_parent_via_opp += 1)
        (!matches_child && !matches_parent_via_opp) && (n_neither += 1)

        println("NCF #$ncf_idx (child=$child_iel, parent=$parent_iel): stored lfid=$stored_lfid, " *
                "child_fid_found=$child_fid_found, parent_fid_found=$parent_fid_found, " *
                "matches_child=$matches_child, matches_parent_via_OPP_FACE=$matches_parent_via_opp")
    end

    println("─"^70)
    println("Summary over $(mesh.num_ncf) NCF entries:")
    println("  stored lfid == CHILD's own fid directly:        $n_match_child")
    println("  OPP_FACE[stored lfid] == PARENT's own fid:      $n_match_parent_via_opp")
    println("  matched neither (unexpected, needs investigation): $n_neither")

    return (n_match_child=n_match_child, n_match_parent_via_opp=n_match_parent_via_opp, n_neither=n_neither)
end

function find_direction_with_child_inflow(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl

    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    child_faces = unique((Int(mesh.cip[i]), Int(mesh.lfid[i])) for i = 1:mesh.num_ncf)
    println("Distinct (child_iel, child_fid) pairs: $child_faces")

    candidates = [(1.0,0.0,0.0), (-1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,-1.0,0.0),
                  (0.0,0.0,1.0), (0.0,0.0,-1.0), (0.3,0.5,0.8), (-0.3,-0.5,-0.8),
                  (0.6,-0.7,0.4), (-0.6,0.7,-0.4)]

    for (Ωx,Ωy,Ωz) in candidates
        found = Tuple{Int,Int}[]
        for (child_iel, child_fid) in child_faces
            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, child_iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)
            if is_active[child_fid] && is_inflow[child_fid]
                push!(found, (child_iel, child_fid))
            end
        end
        println("Ω=($Ωx,$Ωy,$Ωz): child-inflow NCF faces found: $found")
    end

    return mesh, metrics
end

function test_dg_sweep_spatial_ncf_child_inflow(inputs, nparts, distribute, sem, partitioned_model;
                                                 κ_val=0.0001, I0=4.2)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — aborting.")
        return nothing
    end

    Ωx, Ωy, Ωz = 0.0, 0.0, 1.0

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    # ── ONE shared dict, built once, passed to BOTH the ordering and the sweep ──
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
        L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
            child_iel, parent_iel, child_fid, IPc_nodes, IPp_nodes, mesh, ngl)
        ncf_L1L2[ncf_idx] = (L1, L2)
    end

    order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                             metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                             ncf_child_to_parent = ncf_as_child)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    bc_fn = (x,y,z) -> I0*exp(-κ_val*(z - mesh.zmin))

    x_dg = fill(NaN, npoin_dg)
    run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                  connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                  ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2)

    n_nan = count(isnan, x_dg)
    n_nan > 0 && println("⚠ $n_nan DG entries never written")

    max_err = 0.0
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        I_exact = I0*exp(-κ_val*(mesh.z[g]-mesh.zmin))
        v = x_dg[connijk_dg[iel,i,j,k]]
        max_err = max(max_err, abs(v-I_exact))
    end
    println("Spatial NCF child-inflow sweep vs exact: max abs error = $max_err")
    return max_err
end

function test_dg_sweep_spatial_ncf_parent_inflow(inputs, nparts, distribute, sem, partitioned_model;
                                                  κ_val=0.0001, I0=4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics; ngl = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    mesh.num_ncf == 0 && (println("⚠ no NCF"); return nothing)

    Ωx, Ωy, Ωz = 1.0, 0.0, 0.0

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        parent_iel = Int(mesh.pip[ncf_idx])
        parent_fid = OPP_FACE[Int(mesh.lfid[ncf_idx])]
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
    end

    order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                             metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                             ncf_parent_to_children = Dict(k=>[Int(mesh.cip[i]) for i in v] for (k,v) in ncf_as_parent))

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    bc_fn = (x,y,z) -> I0*exp(-κ_val*(x - mesh.xmin))

    x_dg = fill(NaN, npoin_dg)
    run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                  connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                  ncf_as_parent=ncf_as_parent)

    n_nan = count(isnan, x_dg)
    n_nan > 0 && println("⚠ $n_nan unwritten")

    max_err = 0.0
    for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
        g = mesh.connijk[iel,i,j,k]
        I_exact = I0*exp(-κ_val*(mesh.x[g]-mesh.xmin))
        max_err = max(max_err, abs(x_dg[connijk_dg[iel,i,j,k]] - I_exact))
    end
    println("Spatial NCF parent-inflow sweep vs exact: max abs error = $max_err")
    return max_err
end

function debug_ncf_classification_mismatch(inputs, nparts, distribute, sem, partitioned_model,
                                            Ωx, Ωy, Ωz)
    rank = MPI.Comm_rank(get_mpi_comm())
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics; ngl = mesh.ngl

    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]

        ic = Vector{Bool}(undef,6); ia_c = Vector{Bool}(undef,6)
        classify_element_faces!(ic, ia_c, child_iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        ip = Vector{Bool}(undef,6); ia_p = Vector{Bool}(undef,6)
        classify_element_faces!(ip, ia_p, parent_iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)

        child_normal  = compute_face_normal_avg(child_iel, child_fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz,metrics.dηdx,metrics.dηdy,metrics.dηdz,metrics.dζdx,metrics.dζdy,metrics.dζdz)
        parent_normal = compute_face_normal_avg(parent_iel, parent_fid, ngl,
            metrics.dξdx,metrics.dξdy,metrics.dξdz,metrics.dηdx,metrics.dηdy,metrics.dηdz,metrics.dζdx,metrics.dζdy,metrics.dζdz)

        # Expected: parent_normal ≈ -child_normal (opposite faces of the same
        # interface), and hence parent_is_inflow should be the OPPOSITE of
        # child_is_inflow — NOT equal.
        expected_opposite = (ip[parent_fid] != ic[child_fid])
        active_agree = (ia_p[parent_fid] == ia_c[child_fid])

        println("NCF $ncf_idx: child=$child_iel(fid=$child_fid) parent=$parent_iel(fid=$parent_fid)")
        println("  child:  active=$(ia_c[child_fid]) inflow=$(ic[child_fid]) normal=$child_normal")
        println("  parent: active=$(ia_p[parent_fid]) inflow=$(ip[parent_fid]) normal=$parent_normal")
        println("  expected opposite inflow: $expected_opposite  " * (expected_opposite ? "✓" : "✗ SIGN MISMATCH"))
        println("  active status agrees: $active_agree  " * (active_agree ? "✓" : "✗ ACTIVE MISMATCH"))
    end
end

function test_ncf_dependency_graph(inputs, nparts, distribute, sem, partitioned_model,
                                    Ωx::Float64, Ωy::Float64, Ωz::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Build the SAME NCF lookup tables as the full test ─────────────────────
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end

    # ── Elements of interest: the 3 parents + all their children ─────────────
    parents = unique(Int(mesh.pip[i]) for i = 1:mesh.num_ncf)
    children = unique(Int(mesh.cip[i]) for i = 1:mesh.num_ncf)
    elements_of_interest = vcat(parents, children)
    println("Parents: $parents")
    println("Children: $children")
    println("─"^70)

    for iel in elements_of_interest
        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        deps = Int[]
        dep_kinds = String[]
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem[iel][fid]
            if other != zero(TInt)
                push!(deps, Int(other)); push!(dep_kinds, "same-rank(fid=$fid)")
            elseif haskey(ncf_as_child, (iel,fid))
                ncf_idx = ncf_as_child[(iel,fid)]
                p = Int(mesh.pip[ncf_idx])
                push!(deps, p); push!(dep_kinds, "NCF-child->parent(fid=$fid)")
            elseif haskey(ncf_parent_to_children, (iel,fid))
                for c in ncf_parent_to_children[(iel,fid)]
                    push!(deps, c); push!(dep_kinds, "NCF-parent->child(fid=$fid)")
                end
            end
        end
        println("iel=$iel depends on: $(collect(zip(deps,dep_kinds)))")
    end

    return (parents=parents, children=children,
            ncf_as_child=ncf_as_child, ncf_parent_to_children=ncf_parent_to_children,
            neighbor_elem=neighbor_elem, mesh=mesh, metrics=metrics)
end

function test_ncf_boundary_breakdown_and_duplicate_deps(inputs, nparts, distribute, sem, partitioned_model,
                                                         Ωx::Float64, Ωy::Float64, Ωz::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    if mesh.num_ncf == 0
        println("⚠ No non-conforming facets — nothing to test.")
        return nothing
    end

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    tagged_face_keys = Set{NTuple{4,TInt}}()
    for iface_bdy = 1:size(mesh.bdy_face_type, 1)
        ip1 = mesh.poin_in_bdy_face[iface_bdy, 1, 1]
        ip2 = mesh.poin_in_bdy_face[iface_bdy, ngl, 1]
        ip3 = mesh.poin_in_bdy_face[iface_bdy, ngl, ngl]
        ip4 = mesh.poin_in_bdy_face[iface_bdy, 1, ngl]
        push!(tagged_face_keys, Tuple(sort(TInt[ip1, ip2, ip3, ip4])))
    end

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end

    # ── Check 1: boundary/NCF/same-rank breakdown per child ───────────────────
    println("─"^70)
    println("CHECK 1: face-type breakdown for children 1-8")
    for iel = 1:8
        n_true_boundary = 0
        n_ncf = 0
        n_same_rank = 0
        for fid = 1:6
            fe = findfirst(==(fid), CONN_FACE_EL_TO_MY_FID)
            n1=mesh.conn_face_el[1,fe,iel]; n2=mesh.conn_face_el[2,fe,iel]
            n3=mesh.conn_face_el[3,fe,iel]; n4=mesh.conn_face_el[4,fe,iel]
            key = Tuple(sort(TInt[n1,n2,n3,n4]))
            is_tagged_boundary = key in tagged_face_keys
            is_ncf = haskey(ncf_as_child, (iel,fid))
            if is_tagged_boundary
                n_true_boundary += 1
            elseif is_ncf
                n_ncf += 1
            else
                n_same_rank += 1
            end
        end
        println("  child $iel: same_rank=$n_same_rank, true_boundary=$n_true_boundary, ncf=$n_ncf")
    end

    # ── Check 2: duplicate dependency edges for the given direction ───────────
    println("─"^70)
    println("CHECK 2: duplicate dependency edges for Ω=($Ωx,$Ωy,$Ωz)")
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    for iel = 1:8
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        deps = Tuple{Int,Int,String}[]   # (fid, target, kind)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem[iel][fid]
            if other != zero(TInt)
                push!(deps, (fid, Int(other), "same-rank"))
            elseif haskey(ncf_as_child, (iel,fid))
                p = Int(mesh.pip[ncf_as_child[(iel,fid)]])
                push!(deps, (fid, p, "NCF-child->parent"))
            elseif haskey(ncf_parent_to_children, (iel,fid))
                for c in ncf_parent_to_children[(iel,fid)]
                    push!(deps, (fid, c, "NCF-parent->child"))
                end
            end
        end
        targets = [d[2] for d in deps]
        n_dup = length(targets) - length(unique(targets))
        println("  iel=$iel: deps=$deps  " * (n_dup > 0 ? "✗ $n_dup DUPLICATE TARGET(S)" : "✓ no duplicates"))
    end

    return (mesh=mesh, metrics=metrics, neighbor_elem=neighbor_elem,
            ncf_as_child=ncf_as_child, ncf_parent_to_children=ncf_parent_to_children,
            tagged_face_keys=tagged_face_keys)
end

function test_ab_dependency_comparison(inputs, nparts, distribute, sem, partitioned_model,
                                        Ωx::Float64, Ωy::Float64, Ωz::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end

    # ── Path A: freshly re-derive iel=1's deps directly (matches CHECK 2's logic) ──
    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    classify_element_faces!(is_inflow, is_active, 1, Ωx, Ωy, Ωz, ngl,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz)
    deps_A = Int[]
    for fid = 1:6
        (is_inflow[fid] && is_active[fid]) || continue
        other = neighbor_elem[1][fid]
        if other != zero(TInt)
            push!(deps_A, Int(other))
        elseif haskey(ncf_as_child, (1,fid))
            push!(deps_A, Int(mesh.pip[ncf_as_child[(1,fid)]]))
        elseif haskey(ncf_parent_to_children, (1,fid))
            append!(deps_A, ncf_parent_to_children[(1,fid)])
        end
    end
    println("Path A (direct, isolated re-derivation) — iel=1 deps: $deps_A")

    # ── Path B: run the ACTUAL build_conforming_sweep_order (NCF-aware, full
    # mesh, exactly as the failing full-direction test calls it) and extract
    # iel=1's upwind_deps from its internal state. Since the function doesn't
    # currently expose upwind_deps, we inline its EXACT loop body here,
    # over the FULL nelem range, to reproduce byte-for-byte what it computes
    # for element 1 as a SIDE EFFECT of processing the whole mesh — the only
    # difference from Path A being that this iterates over ALL elements
    # first, in the SAME order the real function does. ──────────────────────
    upwind_deps_B = [Int[] for _ = 1:mesh.nelem]
    is_inflow2 = Vector{Bool}(undef,6); is_active2 = Vector{Bool}(undef,6)
    for iel = 1:mesh.nelem
        classify_element_faces!(is_inflow2, is_active2, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        for fid = 1:6
            (is_inflow2[fid] && is_active2[fid]) || continue
            other = neighbor_elem[iel][fid]
            if other == zero(TInt)
                if haskey(ncf_as_child, (Int(iel),fid))
                    p = Int(mesh.pip[ncf_as_child[(Int(iel),fid)]])
                    push!(upwind_deps_B[iel], p)
                elseif haskey(ncf_parent_to_children, (Int(iel),fid))
                    for c in ncf_parent_to_children[(Int(iel),fid)]
                        push!(upwind_deps_B[iel], c)
                    end
                end
                continue
            end
            push!(upwind_deps_B[iel], Int(other))
        end
    end
    println("Path B (full-mesh loop, same logic) — iel=1 deps: $(upwind_deps_B[1])")

    println("─"^70)
    println("Match: $(sort(deps_A) == sort(upwind_deps_B[1]))")

    return (deps_A=deps_A, deps_B=upwind_deps_B[1])
end

function build_conforming_sweep_order_debug(Ωx::Float64, Ωy::Float64, Ωz::Float64,
                                             nelem::TInt, ngl::TInt,
                                             neighbor_elem,
                                             dξdx, dξdy, dξdz,
                                             dηdx, dηdy, dηdz,
                                             dζdx, dζdy, dζdz;
                                             ncf_child_to_parent = nothing,   # NOW expects (child_iel,child_fid) -> PARENT ELEMENT ID directly
                                             ncf_parent_to_children = nothing)
    indegree    = zeros(TInt, nelem)
    upwind_deps = [TInt[] for _ = 1:nelem]
    downstream  = [TInt[] for _ = 1:nelem]

    is_inflow = Vector{Bool}(undef, 6)
    is_active = Vector{Bool}(undef, 6)

    for iel = 1:nelem
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem[iel][fid]
            if other == zero(TInt)
                if ncf_child_to_parent !== nothing && haskey(ncf_child_to_parent, (Int(iel),fid))
                    p = TInt(ncf_child_to_parent[(Int(iel),fid)])   # already an ELEMENT ID now, no translation needed
                    push!(upwind_deps[iel], p); push!(downstream[p], iel); indegree[iel] += 1
                elseif ncf_parent_to_children !== nothing && haskey(ncf_parent_to_children, (Int(iel),fid))
                    for c in ncf_parent_to_children[(Int(iel),fid)]
                        cT = TInt(c)
                        push!(upwind_deps[iel], cT); push!(downstream[cT], iel); indegree[iel] += 1
                    end
                end
                continue
            end
            push!(upwind_deps[iel], other)
            push!(downstream[other], iel)
            indegree[iel] += 1
        end
    end

    order  = TInt[]
    levels = zeros(TInt, nelem)
    queue  = TInt[iel for iel = 1:nelem if indegree[iel] == 0]
    indeg_work = copy(indegree)

    round = 0
    while !isempty(queue)
        round += 1
        next_queue = TInt[]
        for iel in queue
            push!(order, iel)
            for dsn in downstream[iel]
                indeg_work[dsn] -= 1
                levels[dsn] = max(levels[dsn], levels[iel] + 1)
                if indeg_work[dsn] < 0
                    println("⚠ NEGATIVE indegree at element $dsn after processing $iel")
                end
                indeg_work[dsn] == 0 && push!(next_queue, dsn)
            end
        end
        queue = next_queue
    end

    if length(order) != nelem
        stuck = [iel for iel = 1:nelem if indeg_work[iel] > 0]
        println("═"^70)
        println("STALL DEBUG: $(length(stuck)) stuck elements after $round rounds")
        for iel in stuck[1:min(10,end)]
            println("  iel=$iel: original_indegree=$(indegree[iel]), remaining=$(indeg_work[iel]), " *
                     "upwind_deps=$(upwind_deps[iel]), deps_processed=$([d in order for d in upwind_deps[iel]])")
        end
        println("═"^70)
        return (order=order, levels=levels, stuck=stuck, upwind_deps=upwind_deps,
                downstream=downstream, indegree=indegree, indeg_work=indeg_work, failed=true)
    end

    return (order=order, levels=levels, upwind_deps=upwind_deps, downstream=downstream,
            indegree=indegree, indeg_work=indeg_work, failed=false)
end

function test_sweep_order_debug_full(inputs, nparts, distribute, sem, partitioned_model,
                                      Ωx::Float64, Ωy::Float64, Ωz::Float64)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── ncf_as_child: (child_iel,child_fid) -> ncf_idx — kept EXACTLY as
    # run_dg_sweep! needs it (for L1/L2/IPc_nodes lookups). ─────────────────
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end

    # ── NEW: separate, ordering-ONLY dict — (child_iel,child_fid) -> PARENT
    # ELEMENT ID directly, fixing the missing mesh.pip[ncf_idx] translation
    # that build_conforming_sweep_order_debug was silently skipping before. ──
    ncf_child_to_parent_elemid = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child)

    result = build_conforming_sweep_order_debug(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                                 metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                                 ncf_child_to_parent=ncf_child_to_parent_elemid,
                                                 ncf_parent_to_children=ncf_parent_to_children)

    if result.failed
        println("STILL FAILING — upwind_deps[1] = $(result.upwind_deps[1])")
    else
        println("✓ Ordering succeeded, $(length(result.order)) elements processed")
    end

    return result
end

function test_dg_sweep_spatial_ncf_all_directions(inputs, nparts, distribute, sem, partitioned_model;
                                                   κ_val=0.0001, I0=4.2)
    rank = MPI.Comm_rank(get_mpi_comm())
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[415] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics; ngl = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω
    extra_npoin = mesh.extra_mesh.extra_npoin

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    mesh.num_ncf == 0 && (println("⚠ no NCF"); return nothing)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
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
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end
    # FIX: ordering-only dict, mapping directly to PARENT ELEMENT ID, not ncf_idx.
    ncf_child_to_parent_elemid = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    worst_err = 0.0; worst_ip_ext = 0
    all_errs = Float64[]
    n_cycle_failures = 0

    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]; ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        bc_fn = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1]) + Ωy*(y-r_bc[2]) + Ωz*(z-r_bc[3])))

        local order
        try
            order, _ = build_conforming_sweep_order(Ωx, Ωy, Ωz, mesh.nelem, ngl, neighbor_elem,
                                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                                      metrics.dζdx,metrics.dζdy,metrics.dζdz;
                                                      ncf_child_to_parent = ncf_child_to_parent_elemid,
                                                      ncf_parent_to_children = ncf_parent_to_children)
        catch e
            n_cycle_failures += 1
            println("ip_ext=$ip_ext (Ω=($Ωx,$Ωy,$Ωz)): ✗ CYCLE/ORDERING FAILURE — $(sprint(showerror, e))")
            continue
        end

        x_dg = fill(NaN, npoin_dg)
        run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, order, neighbor_elem,
                      connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn;
                      ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2, ncf_as_parent=ncf_as_parent)

        n_nan = count(isnan, x_dg)
        if n_nan > 0
            println("ip_ext=$ip_ext: ⚠ $n_nan unwritten entries — skipping error calc")
            continue
        end

        max_err = 0.0
        for iel = 1:mesh.nelem, k=1:ngl, j=1:ngl, i=1:ngl
            g = mesh.connijk[iel,i,j,k]
            I_exact = bc_fn(mesh.x[g],mesh.y[g],mesh.z[g])
            max_err = max(max_err, abs(x_dg[connijk_dg[iel,i,j,k]] - I_exact))
        end
        push!(all_errs, max_err)
        if max_err > worst_err
            worst_err = max_err; worst_ip_ext = ip_ext
        end
    end

    println("─"^60)
    println("Cycle/ordering failures: $n_cycle_failures / $extra_npoin")
    if !isempty(all_errs)
        println("Worst direction: ip_ext=$worst_ip_ext, max abs error = $worst_err")
        println("Mean error: $(sum(all_errs)/length(all_errs))")
        println("Median error: $(sort(all_errs)[length(all_errs)÷2 + 1])")
    end
    return all_errs
end

function test_manufactured_isotropic_spatial_ncf(inputs, nparts, distribute, sem, partitioned_model;
                                                  I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                  κ_val=0.001, σ_val=0.0005,
                                                  max_iter=500, tol=1e-9)
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; metrics = sem.metrics; ngl = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω
    extra_npoin = mesh.extra_mesh.extra_npoin
    npoin_cg = length(mesh.x)

    println("Refined mesh: $(mesh.nelem) elements, num_ncf=$(mesh.num_ncf)")
    mesh.num_ncf == 0 && (println("⚠ no NCF"); return nothing)

    # ── Domain, manufactured solution — SAME as the conforming isotropic test ──
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    dirs = Vector{NTuple{3,Float64}}(undef, extra_npoin)
    for ip_ext = 1:extra_npoin
        θ = mesh.extra_mesh.extra_coords[1, ip_ext]; ϕ = mesh.extra_mesh.extra_coords[2, ip_ext]
        dirs[ip_ext] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
    end
    w_ext = zeros(Float64, extra_npoin)
    for e_ext = 1:mesh.extra_mesh.extra_nelem
        nop = mesh.extra_mesh.extra_nop[e_ext]
        for iθ=1:nop+1, iϕ=1:nop+1
            ip_ext = mesh.extra_mesh.extra_connijk[e_ext,iθ,iϕ]
            Je_ang = mesh.extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]
            w_ext[ip_ext] += Je_ang * mesh.extra_mesh.ωθ[iθ] * mesh.extra_mesh.ωθ[iϕ]
        end
    end
    W = sum(w_ext)
    C = [σ_val*w_ext[m]/W for k=1:extra_npoin, m=1:extra_npoin]

    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    q_ext_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
            grad = dI_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            Ω_dot_grad = Ωx*grad[1] + Ωy*grad[2] + Ωz*grad[3]
            I_true = I_manuf(x_,y_,z_,Ωx,Ωy,Ωz)
            scatter_true = sum(C[k,m]*I_manuf(x_,y_,z_,dirs[m]...) for m = 1:extra_npoin)
            q_ext_all[k][g] = Ω_dot_grad + κ_val*I_true - scatter_true
        end
    end
    bc_fn_all = [(x,y,z) -> I_manuf(x,y,z,dirs[k]...) for k = 1:extra_npoin]

    # ── Mesh/NCF infrastructure — built ONCE, direction- and
    # iteration-independent, reused by every call below. ────────────────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    face_neighbor_dg = build_face_neighbor_dg(mesh, connijk_dg, neighbor_elem, neighbor_fid, face_correspondence)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
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
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end
    ncf_child_to_parent_elemid = Dict{Tuple{Int,Int}, Int}(
        key => Int(mesh.pip[ncf_idx]) for (key, ncf_idx) in ncf_as_child)

    orders = [build_conforming_sweep_order(dirs[k]..., mesh.nelem, ngl, neighbor_elem,
                  metrics.dξdx,metrics.dξdy,metrics.dξdz,
                  metrics.dηdx,metrics.dηdy,metrics.dηdz,
                  metrics.dζdx,metrics.dζdy,metrics.dζdz;
                  ncf_child_to_parent = ncf_child_to_parent_elemid,
                  ncf_parent_to_children = ncf_parent_to_children)[1] for k = 1:extra_npoin]

    # ── Iterative solve — SAME structure as the conforming isotropic test,
    # just with NCF kwargs threaded into run_dg_sweep!. ──────────────────────
    I_cg_all = [zeros(npoin_cg) for _ = 1:extra_npoin]
    for it = 1:max_iter
        max_delta = 0.0
        for k = 1:extra_npoin
            Ωx,Ωy,Ωz = dirs[k]
            q_k = copy(q_ext_all[k])
            for mm = 1:extra_npoin
                C[k,mm] == 0.0 && continue
                q_k .+= C[k,mm] .* I_cg_all[mm]
            end
            x_dg = fill(NaN, npoin_dg)
            run_dg_sweep!(x_dg, mesh, metrics, ψ, dψ, ω, Ωx, Ωy, Ωz, orders[k], neighbor_elem,
                          connijk_dg, face_neighbor_dg, κ, σz, intΦ_eff, bc_fn_all[k];
                          extra_source=q_k, ncf_as_child=ncf_as_child, ncf_L1L2=ncf_L1L2,
                          ncf_as_parent=ncf_as_parent)
            I_new = zeros(npoin_cg)
            for iel = 1:mesh.nelem, kk=1:ngl, jj=1:ngl, ii=1:ngl
                I_new[mesh.connijk[iel,ii,jj,kk]] = x_dg[connijk_dg[iel,ii,jj,kk]]
            end
            max_delta = max(max_delta, maximum(abs.(I_new .- I_cg_all[k])))
            I_cg_all[k] = I_new
        end
        println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (println("Converged after $it iterations"); break)
        it == max_iter && @warn "Did not converge (max_delta=$max_delta)"
    end

    max_err = 0.0
    for k = 1:extra_npoin
        Ωx,Ωy,Ωz = dirs[k]
        for g = 1:npoin_cg
            I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
            max_err = max(max_err, abs(I_cg_all[k][g] - I_true))
        end
    end
    println("Manufactured (isotropic, spatial NCF) — max abs error: $max_err")
    return max_err
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_angular_ncf_correspondence  (DRAFT — untested)
#
#  Given ONE spatial element's angular mesh (possibly refined) and a
#  NEIGHBORING spatial element's angular mesh (assumed unrefined, per the
#  "spatial-conforming, angular-only NCF" scope agreed for this first pass),
#  finds, PURELY GEOMETRICALLY (no index convention trusted, since
#  adapt_angular_grid_3Dby2D! does not preserve one after refinement):
#
#    - for each group of ≤4 fine angular elements at `iel` that exactly
#      tile one coarse (θ,ϕ) rectangle, the matching coarse angular
#      element at `jel` (if the rectangles agree)
#
#  Returns a Vector of NamedTuples, one per matched group:
#    (coarse_e_ext_at_jel, fine_e_exts_at_iel, θmin,θmax,ϕmin,ϕmax)
#
#  Does NOT yet build any interpolation weights — pure geometric matching,
#  verified independently before any of that.
# ─────────────────────────────────────────────────────────────────────────────

"""
    get_angular_element_bounds(coords_ang, connijk_ang, nop_ang, e_ext)

Returns (θmin, θmax, ϕmin, ϕmax) for one angular element, correctly
disambiguating the ϕ=0/2π periodic seam merge performed by
adapt_angular_grid_3Dby2D! (which replaces every physical ϕ=2π node with
its merged ϕ=0 counterpart). An element that GENUINELY spans up to ϕ=2π
will show ϕ=0.0 at its "high" corner nodes — this must be treated as 2π
for bounds purposes, or the bounding box collapses incorrectly.

Disambiguation: look at the element's OWN interior ϕ values (not just the
corners). If the element's non-zero interior ϕ values are large (near 2π),
a 0.0 at the boundary is almost certainly a wrapped 2π, not a genuine
ϕ≈0 element — since a genuine low-ϕ element's interior values would be
small too, not clustered near 2π elsewhere in the SAME element.
"""
function get_angular_element_bounds(coords_ang, connijk_ang, nop_ang, e_ext)
    n = nop_ang[e_ext] + 1
    θs = Float64[]; ϕs_raw = Float64[]
    for i = 1:n, j = 1:n
        ip = connijk_ang[e_ext, i, j]
        push!(θs, coords_ang[1, ip])
        push!(ϕs_raw, coords_ang[2, ip])
    end

    # Determine if this element genuinely straddles/touches the seam: true
    # iff it contains BOTH near-zero AND large (>π, say) ϕ values among its
    # OWN nodes — the signature of a merged 2π/0 boundary, not a genuine
    # low-ϕ element (which would have ALL small ϕ values, no large ones).
    has_near_zero = any(ϕ -> ϕ < 1e-8, ϕs_raw)
    has_large     = any(ϕ -> ϕ > π, ϕs_raw)

    if has_near_zero && has_large
        # Reinterpret every ϕ≈0 as 2π for THIS element's bounds only.
        ϕs = [ϕ < 1e-8 ? 2π : ϕ for ϕ in ϕs_raw]
    else
        ϕs = ϕs_raw
    end

    return (minimum(θs), maximum(θs), minimum(ϕs), maximum(ϕs))
end


function build_angular_ncf_correspondence(coords_ang_iel, connijk_ang_iel, nop_ang_iel, nelem_ang_iel,
                                           coords_ang_jel, connijk_ang_jel, nop_ang_jel, nelem_ang_jel;
                                           tol = 1e-8)
    bounds_iel = [get_angular_element_bounds(coords_ang_iel, connijk_ang_iel, nop_ang_iel, e)
                  for e = 1:nelem_ang_iel]
    bounds_jel = [get_angular_element_bounds(coords_ang_jel, connijk_ang_jel, nop_ang_jel, e)
                  for e = 1:nelem_ang_jel]

    matches = NamedTuple[]
    used_iel = falses(nelem_ang_iel)
    n_identical_rejected = 0
    n_other_rejected = 0

    for (je, (jθmin,jθmax,jϕmin,jϕmax)) in enumerate(bounds_jel)
        candidates = Int[]
        for (ie, (iθmin,iθmax,iϕmin,iϕmax)) in enumerate(bounds_iel)
            used_iel[ie] && continue
            inside = (iθmin >= jθmin-tol) && (iθmax <= jθmax+tol) &&
                     (iϕmin >= jϕmin-tol) && (iϕmax <= jϕmax+tol)
            exact = isapprox(iθmin,jθmin;atol=tol) && isapprox(iθmax,jθmax;atol=tol) &&
                    isapprox(iϕmin,jϕmin;atol=tol) && isapprox(iϕmax,jϕmax;atol=tol)
            (inside && !exact) && push!(candidates, ie)
        end

        length(candidates) == 0 && continue

        if length(candidates) != 4
            is_identical_case = false
            if length(candidates) == 1
                ie = candidates[1]
                iθmin,iθmax,iϕmin,iϕmax = bounds_iel[ie]
                is_identical_case = isapprox(iθmin,jθmin;atol=tol) && isapprox(iθmax,jθmax;atol=tol) &&
                                     isapprox(iϕmin,jϕmin;atol=tol) && isapprox(iϕmax,jϕmax;atol=tol)
            end
            if is_identical_case
                n_identical_rejected += 1
            else
                n_other_rejected += 1
                @warn "build_angular_ncf_correspondence: found $(length(candidates)) fine " *
                      "elements inside coarse element $je's range (not identical, not 4) — " *
                      "REJECTING. candidates=$(candidates), coarse θ∈[$jθmin,$jθmax] ϕ∈[$jϕmin,$jϕmax]"
            end
            continue
        end

        # ── Coordinate-exact tiling check, replacing the area check: the 4
        # candidates' combined MIDPOINT SPLIT must exactly reproduce the
        # coarse rectangle's own bisection point (θ12, ϕ12) — the same
        # construction adapt_angular_grid_3Dby2D! itself uses, checked here
        # independently rather than inferred from area alone. ──────────────
        θ12_expected = (jθmin + jθmax) / 2
        ϕ12_expected = (jϕmin + jϕmax) / 2

        # Each candidate should have EITHER θmin==jθmin OR θmax==jθmax (i.e.
        # abuts one of the coarse element's θ-edges), and its OTHER θ-bound
        # should equal θ12_expected exactly — same for ϕ. Verify directly.
        all_exact = true
        for ie in candidates
            iθmin,iθmax,iϕmin,iϕmax = bounds_iel[ie]
            θ_ok = (isapprox(iθmin,jθmin;atol=tol) && isapprox(iθmax,θ12_expected;atol=tol)) ||
                   (isapprox(iθmin,θ12_expected;atol=tol) && isapprox(iθmax,jθmax;atol=tol))
            ϕ_ok = (isapprox(iϕmin,jϕmin;atol=tol) && isapprox(iϕmax,ϕ12_expected;atol=tol)) ||
                   (isapprox(iϕmin,ϕ12_expected;atol=tol) && isapprox(iϕmax,jϕmax;atol=tol))
            all_exact &= (θ_ok && ϕ_ok)
        end

        if !all_exact
            @warn "build_angular_ncf_correspondence: 4 candidates found for coarse element $je " *
                  "but their bounds do NOT exactly quarter it at the expected midpoint " *
                  "(θ12=$θ12_expected, ϕ12=$ϕ12_expected) — REJECTING, this is not a clean " *
                  "quadtree split of this coarse element."
            n_other_rejected += 1
            continue
        end

        for ie in candidates; used_iel[ie] = true; end
        push!(matches, (coarse_e_ext_at_jel=je, fine_e_exts_at_iel=candidates,
                         θmin=jθmin, θmax=jθmax, ϕmin=jϕmin, ϕmax=jϕmax))
    end

    if n_identical_rejected > 0
        println("Note: $n_identical_rejected coarse element(s) matched exactly one identical " *
                "(conforming) fine element — correctly excluded.")
    end
    if n_other_rejected > 0
        println("⚠ $n_other_rejected coarse element(s) had an unexplained/non-exact partial " *
                "match — see warnings above.")
    end

    return matches
end

function adapt_single_angular_element!(mesh, target_iel::Int, target_e_ext::Int)
    iel = target_iel
    e_ext = target_e_ext

    em = mesh.extra_mesh[iel]
    nop_ang     = em.extra_nop
    nelem_ang   = em.extra_nelem
    connijk_ang = em.extra_connijk
    coords_ang  = em.extra_coords

    θmax = coords_ang[1, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    θmin = coords_ang[1, connijk_ang[e_ext, 1, 1]]
    ϕmax = coords_ang[2, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    ϕmin = coords_ang[2, connijk_ang[e_ext, 1, 1]]
    ϕmax == 0 && (ϕmax = 2π)
    θ12 = (θmax + θmin) / 2
    ϕ12 = (ϕmax + ϕmin) / 2

    nelem_ang_new = nelem_ang + 3
    npoin_ang_new_upper = em.extra_npoin + (4*(nop_ang[e_ext]+1)^2 - 4*(nop_ang[e_ext]+1))

    connijk_ang_new = zeros(Int, nelem_ang_new, nop_ang[1]+1, nop_ang[1]+1)
    coords_new      = zeros(Float64, 2, npoin_ang_new_upper)
    nop_ang_new     = zeros(Int, nelem_ang_new)
    nop_ang_new[1:nelem_ang_new-1] .= nop_ang[e_ext]
    nop_ang_new[nelem_ang_new]      = nop_ang[1]

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

    lgl = basis_structs_ξ_ω!(LGL(), nop_ang[1], CPU())

    if e_ext > 1
        for e_ext1 = 1:e_ext-1
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                insert_point!(θ, ϕ, connijk_ang_new, e_ext1, i, j)
            end
        end
    end

    child_ranges = [
        (θmin, θ12, ϕmin, ϕ12, e_ext),
        (θ12, θmax, ϕmin, ϕ12, e_ext+1),
        (θmin, θ12, ϕ12, ϕmax, e_ext+2),
        (θ12, θmax, ϕ12, ϕmax, e_ext+3),
    ]
    for (θlo, θhi, ϕlo, ϕhi, e_child) in child_ranges
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

    if e_ext < nelem_ang
        for e_ext1 = e_ext+1:nelem_ang
            dst = e_ext1 + 3
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                insert_point!(θ, ϕ, connijk_ang_new, dst, i, j)
            end
        end
    end

    npoin_ang_new = iter - 1
    coords_new = coords_new[:, 1:npoin_ang_new]

    # ── RESTORED: enforce ϕ-periodicity — merge ϕ=2π nodes with their ϕ=0
    # counterparts, EXACT logic from adapt_angular_grid_3Dby2D!. Only
    # produces a geometrically meaningful result when both sides of the
    # seam are refined to matching resolution (confirmed) — that
    # precondition is now enforced by adapt_single_angular_element_conforming!
    # calling this on BOTH seam partners, and independently checked by
    # verify_angular_phi_seam_conformity. ────────────────────────────────────
    zero_ϕ_dict = Dict{Float64, Int}()
    for iper = 1:npoin_ang_new
        if coords_new[2, iper] < 1e-10
            zero_ϕ_dict[coords_new[1, iper]] = iper
        end
    end
    merge_pairs = Dict{Int, Int}()
    for iper = 1:npoin_ang_new
        ϕ = coords_new[2, iper]
        if abs(ϕ / π - 2.0) < 1e-10
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
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] == ip_old && (connijk_ang_new[e,i,j] = ip_new)
                end
            end
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] > ip_old && (connijk_ang_new[e,i,j] -= 1)
                end
            end
            merge_pairs = Dict(
                (k > ip_old ? k - 1 : k) => (v > ip_old ? v - 1 : v)
                    for (k, v) in merge_pairs if k != ip_old
                )
            for i = ip_old+1:npoin_ang_new
                coords_new[1,i-1] = coords_new[1,i]
                coords_new[2,i-1] = coords_new[2,i]
            end
            npoin_ang_new -= 1
        end
    end
    coords_new = coords_new[:, 1:npoin_ang_new]

    em.extra_connijk = connijk_ang_new
    em.extra_coords   = coords_new
    em.extra_nop      = nop_ang_new
    em.extra_nelem    = nelem_ang_new
    em.extra_npoin    = npoin_ang_new

    return nothing
end

function test_angular_ncf_correspondence(mesh, metrics)
    target_iel = 1
    target_e_ext = 4

    adapt_single_angular_element!(mesh, target_iel, target_e_ext)

    em_iel = mesh.extra_mesh[target_iel]
    println("Element $target_iel angular mesh after refinement: " *
            "$(em_iel.extra_nelem) angular elements, $(em_iel.extra_npoin) angular points")

    # ── Find a spatially-adjacent, UNREFINED neighbor of element 1 ──────────
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand)
        break
    end
    if jel == 0
        println("⚠ No same-rank spatial neighbor found for element $target_iel — cannot test.")
        return nothing
    end
    em_jel = mesh.extra_mesh[jel]
    println("Testing against spatial neighbor jel=$jel " *
            "($(em_jel.extra_nelem) angular elements, unrefined)")

    matches = build_angular_ncf_correspondence(
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, em_iel.extra_nelem,
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, em_jel.extra_nelem)

    println("Matches found: $(length(matches))")
    for m in matches
        println("  coarse e_ext=$(m.coarse_e_ext_at_jel) <-> fine e_exts=$(m.fine_e_exts_at_iel) " *
                 "θ∈[$(m.θmin),$(m.θmax)] ϕ∈[$(m.ϕmin),$(m.ϕmax)]")
    end

    return matches
end

"""
    build_angular_ncf_parent_to_child_weights(coords_ang_jel, connijk_ang_jel, nop_ang_jel, coarse_e_ext,
                                               coords_ang_iel, connijk_ang_iel, nop_ang_iel, fine_e_ext)

Builds (Lθ, Lϕ) tensor-product interpolation matrices evaluating the
COARSE (parent) angular element's nodal values at the FINE (child)
element's own node locations — direct evaluation, NOT a projection,
matching the lesson learned from the spatial case (mass-matrix projection
was the WRONG operator there; the flux integral only ever needs a genuine
point value, since quadrature weighting is already supplied separately by
the consuming code).

Handles the ϕ=0/2π periodic seam the same way as
get_angular_element_bounds: if this SPECIFIC (coarse,fine) pair straddles
the seam, coarse ϕ nodes that are genuinely 0.0-after-merge are shifted to
2π before building the interpolation, so build_1d_lagrange_matrix always
sees a single, unambiguous monotonic coordinate range.
"""
function build_angular_ncf_parent_to_child_weights(coords_ang_jel, connijk_ang_jel, nop_ang_jel, coarse_e_ext,
                                                     coords_ang_iel, connijk_ang_iel, nop_ang_iel, fine_e_ext)
    ncoarse = nop_ang_jel[coarse_e_ext] + 1
    nfine   = nop_ang_iel[fine_e_ext] + 1

    # ── Extract coarse element's own θ-line and ϕ-line (tensor-product grid,
    # so row/column extraction is trivial — no line-finding needed). ────────
    θ_coarse = [coords_ang_jel[1, connijk_ang_jel[coarse_e_ext, i, 1]] for i = 1:ncoarse]
    ϕ_coarse = [coords_ang_jel[2, connijk_ang_jel[coarse_e_ext, 1, j]] for j = 1:ncoarse]

    θ_fine = [coords_ang_iel[1, connijk_ang_iel[fine_e_ext, i, 1]] for i = 1:nfine]
    ϕ_fine = [coords_ang_iel[2, connijk_ang_iel[fine_e_ext, 1, j]] for j = 1:nfine]

    # ── Periodicity: does EITHER side's ϕ range touch the seam? Check using
    # the SAME signature as get_angular_element_bounds — near-zero AND
    # large values present together, this time checked ACROSS both the
    # coarse and fine ϕ-lines jointly, since a merged fine child might show
    # 0.0 while its coarse parent still correctly shows near-2π (or vice
    # versa), and both need to be shifted CONSISTENTLY relative to each
    # other for the interpolation to be meaningful at all. ──────────────────
    all_ϕ = vcat(ϕ_coarse, ϕ_fine)
    has_near_zero = any(ϕ -> ϕ < 1e-8, all_ϕ)
    has_large     = any(ϕ -> ϕ > π, all_ϕ)
    if has_near_zero && has_large
        ϕ_coarse = [ϕ < 1e-8 ? 2π : ϕ for ϕ in ϕ_coarse]
        ϕ_fine   = [ϕ < 1e-8 ? 2π : ϕ for ϕ in ϕ_fine]
    end

    Lθ = build_1d_lagrange_matrix(θ_coarse, θ_fine)
    Lϕ = build_1d_lagrange_matrix(ϕ_coarse, ϕ_fine)

    return Lθ, Lϕ
end

"""
    evaluate_angular_ncf_parent_to_child(coarse_grid, Lθ, Lϕ)

Applies the tensor-product weights to the coarse element's own nodal
solution grid, returning the fine element's interpolated grid. `coarse_grid`
is (ncoarse × ncoarse), indexed [i,j] matching the SAME (θ,ϕ) node
convention as connijk_ang[e_ext,i,j].
"""
function evaluate_angular_ncf_parent_to_child(coarse_grid::Matrix{Float64}, Lθ::Matrix{Float64}, Lϕ::Matrix{Float64})
    nfine_θ = size(Lθ,1); nfine_ϕ = size(Lϕ,1)
    fine_grid = zeros(nfine_θ, nfine_ϕ)
    for a = 1:nfine_θ, b = 1:nfine_ϕ
        val = 0.0
        for p = 1:size(Lθ,2), q = 1:size(Lϕ,2)
            val += Lθ[a,p]*Lϕ[b,q]*coarse_grid[p,q]
        end
        fine_grid[a,b] = val
    end
    return fine_grid
end

function test_angular_ncf_parent_to_child(mesh)
    target_iel = 1
    target_e_ext = 4   # the seam-straddling element confirmed working above

    adapt_single_angular_element!(mesh, target_iel, target_e_ext)
    em_iel = mesh.extra_mesh[target_iel]

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); break
    end
    em_jel = mesh.extra_mesh[jel]

    matches = build_angular_ncf_correspondence(
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, em_iel.extra_nelem,
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, em_jel.extra_nelem)

    if length(matches) != 1
        println("⚠ expected exactly 1 match, got $(length(matches)) — aborting")
        return nothing
    end
    m = matches[1]
    println("Testing coarse e_ext=$(m.coarse_e_ext_at_jel) <-> fine e_exts=$(m.fine_e_exts_at_iel)")

    # ── Well-scaled field, normalized against THIS coarse element's own
    # θ,ϕ range — same discipline as the spatial NCF fix. ────────────────────
    θ0, θ1 = m.θmin, m.θmax
    ϕ0, ϕ1 = m.ϕmin, m.ϕmax
    Lθ_range = θ1-θ0; Lϕ_range = ϕ1-ϕ0
    ncoarse = em_jel.extra_nop[m.coarse_e_ext_at_jel] + 1
    P(t) = sum(0.3^k * t^k for k=0:ncoarse-1)
    field(θ,ϕ) = P((θ-θ0)/Lθ_range) + P((ϕ-ϕ0)/Lϕ_range)

	
    coarse_grid = zeros(ncoarse, ncoarse)
    for i = 1:ncoarse, j = 1:ncoarse
        ip = em_jel.extra_connijk[m.coarse_e_ext_at_jel, i, j]
    	θ, ϕ = em_jel.extra_coords[1,ip], em_jel.extra_coords[2,ip]
    	ϕ_eval = (ϕ < 1e-8 && ϕ1 > π) ? 2π : ϕ   # SAME seam correction, applied HERE too
    	coarse_grid[i,j] = field(θ, ϕ_eval)
    end

    max_err = 0.0
    for fine_e_ext in m.fine_e_exts_at_iel
        Lθ, Lϕ = build_angular_ncf_parent_to_child_weights(
            em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, m.coarse_e_ext_at_jel,
            em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, fine_e_ext)
        fine_computed = evaluate_angular_ncf_parent_to_child(coarse_grid, Lθ, Lϕ)

        nfine = em_iel.extra_nop[fine_e_ext] + 1
        for a = 1:nfine, b = 1:nfine
            ip = em_iel.extra_connijk[fine_e_ext, a, b]
            θ, ϕ = em_iel.extra_coords[1,ip], em_iel.extra_coords[2,ip]
            ϕ_eval = (ϕ < 1e-8 && ϕ1 > π) ? 2π : ϕ   # same seam correction for the EXACT comparison
            exact = field(θ, ϕ_eval)
            max_err = max(max_err, abs(fine_computed[a,b] - exact))
        end
    end
    println("Angular NCF parent→child interpolation vs exact field: max abs error = $max_err")
    return max_err
end

"""
    build_angular_ncf_child_to_parent(em_iel, fine_e_exts, coarse_grid_shape,
                                       em_jel, coarse_e_ext, θmin,θmax,ϕmin,ϕmax,
                                       child_values_lookup)

Angular analogue of build_dg_ncf_child_to_parent_v5: for each of the COARSE
(parent) element's own node locations, finds the ONE covering fine (child)
element and evaluates that child's own nodal values there via DIRECT 2D
Lagrange interpolation — NOT a projection (same lesson as the spatial
case). Every coarse node must be claimed by exactly one child; a NaN
sentinel flags any gap.

`child_values_lookup(fine_e_ext)` is a caller-supplied function returning
that fine element's (nfine x nfine) solved nodal grid — kept abstract
here since the actual DG storage/indexing convention will be decided at
sweep-wiring time, same as the spatial case's design.

Seam handling: ϕ=0/2π disambiguation is applied CONSISTENTLY to every raw
coordinate read here — the coarse node's own ϕ, and each candidate child's
ϕ range — using the SAME `ϕ_range_upper > π` signal already established,
learned directly from the earlier test bug (weights and data must agree
on ONE seam convention, not each apply their own).
"""
function build_angular_ncf_child_to_parent(em_iel, fine_e_exts::Vector{Int},
                                            em_jel, coarse_e_ext::Int,
                                            θmin::Float64, θmax::Float64, ϕmin::Float64, ϕmax::Float64,
                                            child_values_lookup)
    ncoarse = em_jel.extra_nop[coarse_e_ext] + 1
    parent_result = fill(NaN, ncoarse, ncoarse)

    seam_case = ϕmax > π && (ϕmin < 1e-8 || ϕmax > π)   # this coarse element's OWN range touches/exceeds π→2π side
    # More precise, consistent signal: does the OVERALL coarse element range
    # indicate we're on the high-ϕ side of a seam-merged element at all —
    # i.e., ϕmax genuinely close to 2π (not just >π, since a normal element
    # could span up to π without touching the seam at all).
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)

    fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ

    for fine_e_ext in fine_e_exts
        nfine = em_iel.extra_nop[fine_e_ext] + 1
        θ_fine = [em_iel.extra_coords[1, em_iel.extra_connijk[fine_e_ext,i,1]] for i = 1:nfine]
        ϕ_fine_raw = [em_iel.extra_coords[2, em_iel.extra_connijk[fine_e_ext,1,j]] for j = 1:nfine]
        ϕ_fine = fix_ϕ.(ϕ_fine_raw)

        θf_min,θf_max = extrema(θ_fine); ϕf_min,ϕf_max = extrema(ϕ_fine)
        tol = 1e-6*max(θf_max-θf_min, ϕf_max-ϕf_min, 1.0)

        child_grid = child_values_lookup(fine_e_ext)   # (nfine x nfine), caller-supplied

        for pi = 1:ncoarse, pj = 1:ncoarse
            isnan(parent_result[pi,pj]) || continue   # already claimed by an earlier child
            ip = em_jel.extra_connijk[coarse_e_ext, pi, pj]
            pθ = em_jel.extra_coords[1, ip]
            pϕ = fix_ϕ(em_jel.extra_coords[2, ip])

            (pθ < θf_min-tol || pθ > θf_max+tol || pϕ < ϕf_min-tol || pϕ > ϕf_max+tol) && continue

            Lθ = build_1d_lagrange_matrix(θ_fine, [pθ])
            Lϕ = build_1d_lagrange_matrix(ϕ_fine, [pϕ])
            val = 0.0
            for a = 1:nfine, b = 1:nfine
                val += Lθ[1,a]*Lϕ[1,b]*child_grid[a,b]
            end
            parent_result[pi,pj] = val
        end
    end

    return parent_result
end

function test_angular_ncf_child_to_parent(mesh)
    target_iel = 1
    target_e_ext = 4

    adapt_single_angular_element!(mesh, target_iel, target_e_ext)
    em_iel = mesh.extra_mesh[target_iel]

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); break
    end
    em_jel = mesh.extra_mesh[jel]

    matches = build_angular_ncf_correspondence(
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, em_iel.extra_nelem,
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, em_jel.extra_nelem)
    length(matches) != 1 && (println("⚠ expected 1 match, got $(length(matches))"); return nothing)
    m = matches[1]

    θ0,θ1,ϕ0,ϕ1 = m.θmin, m.θmax, m.ϕmin, m.ϕmax
    near_seam = isapprox(ϕ1, 2π; atol=1e-6)
    fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ

    ncoarse = em_jel.extra_nop[m.coarse_e_ext_at_jel] + 1
    P(t) = sum(0.3^k * t^k for k=0:ncoarse-1)
    field(θ,ϕ) = P((θ-θ0)/(θ1-θ0)) + P((ϕ-ϕ0)/(ϕ1-ϕ0))

    # Build EXACT fine-element grids (what a real solve would have produced)
    # using field() evaluated at each child's own (seam-corrected) nodes —
    # these stand in for "already solved" child values.
    child_values_cache = Dict{Int,Matrix{Float64}}()
    for fe in m.fine_e_exts_at_iel
        nfine = em_iel.extra_nop[fe] + 1
        grid = zeros(nfine,nfine)
        for a=1:nfine, b=1:nfine
            ip = em_iel.extra_connijk[fe,a,b]
            θ,ϕ = em_iel.extra_coords[1,ip], em_iel.extra_coords[2,ip]
            grid[a,b] = field(θ, fix_ϕ(ϕ))
        end
        child_values_cache[fe] = grid
    end
    lookup(fe) = child_values_cache[fe]

    parent_computed = build_angular_ncf_child_to_parent(
        em_iel, m.fine_e_exts_at_iel, em_jel, m.coarse_e_ext_at_jel, θ0,θ1,ϕ0,ϕ1, lookup)

    n_nan = count(isnan, parent_computed)
    n_nan > 0 && println("⚠ $n_nan coarse nodes never claimed by any child")

    max_err = 0.0
    for pi = 1:ncoarse, pj = 1:ncoarse
        ip = em_jel.extra_connijk[m.coarse_e_ext_at_jel, pi, pj]
        θ,ϕ = em_jel.extra_coords[1,ip], em_jel.extra_coords[2,ip]
        exact = field(θ, fix_ϕ(ϕ))
        max_err = max(max_err, abs(parent_computed[pi,pj] - exact))
    end
    println("Angular NCF child→parent interpolation vs exact field: max abs error = $max_err")
    return max_err
end

function build_angular_dependency_graph(mesh, neighbor_elem; tol::Float64 = 1e-6)
    angular_deps = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Vector{Int}}}}()
    processed_pairs = Set{Tuple{Int,Int}}()

    for iel = 1:mesh.nelem
        em_iel = mesh.extra_mesh[iel]
        for fid = 1:6
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            jel = Int(jel)

            pair_key = (min(iel,jel), max(iel,jel))
            pair_key in processed_pairs && continue
            push!(processed_pairs, pair_key)

            em_jel = mesh.extra_mesh[jel]

            if em_iel.extra_nelem == em_jel.extra_nelem &&
               isapprox(em_iel.extra_coords, em_jel.extra_coords; atol=tol)
                continue
            end

            jel_bounds = [get_angular_element_bounds(em_jel.extra_coords, em_jel.extra_connijk,
                                                       em_jel.extra_nop, e) for e = 1:em_jel.extra_nelem]
            iel_bounds = [get_angular_element_bounds(em_iel.extra_coords, em_iel.extra_connijk,
                                                       em_iel.extra_nop, e) for e = 1:em_iel.extra_nelem]

            done_iel = falses(em_iel.extra_npoin)
            for e_i = 1:em_iel.extra_nelem
                θmin_i,θmax_i,ϕmin_i,ϕmax_i = iel_bounds[e_i]
                near_seam_i = isapprox(ϕmax_i, 2π; atol=1e-6)
                ni = em_iel.extra_nop[e_i] + 1
                for a = 1:ni, b = 1:ni
                    ip_ext_i = em_iel.extra_connijk[e_i,a,b]
                    done_iel[ip_ext_i] && continue
                    done_iel[ip_ext_i] = true

                    θi = em_iel.extra_coords[1, ip_ext_i]
                    ϕi_raw = em_iel.extra_coords[2, ip_ext_i]
                    ϕi = (near_seam_i && ϕi_raw < 1e-8) ? 2π : ϕi_raw

                    matched_jel_node = 0
                    for ip_ext_j = 1:em_jel.extra_npoin
                        θj = em_jel.extra_coords[1,ip_ext_j]
                        ϕj_raw = em_jel.extra_coords[2,ip_ext_j]
                        ϕj_candidates = ϕj_raw < 1e-8 ? (ϕj_raw, 2π) : (ϕj_raw,)
                        for ϕj in ϕj_candidates
                            if isapprox(θi,θj;atol=tol) && isapprox(ϕi,ϕj;atol=tol)
                                matched_jel_node = ip_ext_j
                                break
                            end
                        end
                        matched_jel_node != 0 && break
                    end
                    if matched_jel_node != 0
                        key = (iel, ip_ext_i)
                        push!(get!(angular_deps, key, Tuple{Int,Vector{Int}}[]), (jel, [matched_jel_node]))
                        # NEW: reverse singleton entry
                        rkey = (jel, matched_jel_node)
                        push!(get!(angular_deps, rkey, Tuple{Int,Vector{Int}}[]), (iel, [ip_ext_i]))
                        continue
                    end

                    for e_j = 1:em_jel.extra_nelem
                        θmin_j,θmax_j,ϕmin_j,ϕmax_j = jel_bounds[e_j]
                        inside = (θi >= θmin_j-tol) && (θi <= θmax_j+tol) &&
                                 (ϕi >= ϕmin_j-tol) && (ϕi <= ϕmax_j+tol)
                        inside || continue
                        nj = em_jel.extra_nop[e_j] + 1
                        jel_nodes = [em_jel.extra_connijk[e_j,p,q] for p=1:nj, q=1:nj] |> vec
                        key = (iel, ip_ext_i)
                        push!(get!(angular_deps, key, Tuple{Int,Vector{Int}}[]), (jel, jel_nodes))
                        break
                    end
                end
            end

            done_jel = falses(em_jel.extra_npoin)
            for e_j = 1:em_jel.extra_nelem
                θmin_j,θmax_j,ϕmin_j,ϕmax_j = jel_bounds[e_j]
                near_seam_j = isapprox(ϕmax_j, 2π; atol=1e-6)
                nj = em_jel.extra_nop[e_j] + 1
                for a = 1:nj, b = 1:nj
                    ip_ext_j = em_jel.extra_connijk[e_j,a,b]
                    done_jel[ip_ext_j] && continue
                    done_jel[ip_ext_j] = true

                    θj = em_jel.extra_coords[1, ip_ext_j]
                    ϕj_raw = em_jel.extra_coords[2, ip_ext_j]
                    ϕj = (near_seam_j && ϕj_raw < 1e-8) ? 2π : ϕj_raw

                    already_matched = false
                    for ip_ext_i = 1:em_iel.extra_npoin
                        θi = em_iel.extra_coords[1,ip_ext_i]
                        ϕi_raw = em_iel.extra_coords[2,ip_ext_i]
                        ϕi_candidates = ϕi_raw < 1e-8 ? (ϕi_raw, 2π) : (ϕi_raw,)
                        for ϕi in ϕi_candidates
                            if isapprox(θj,θi;atol=tol) && isapprox(ϕj,ϕi;atol=tol)
                                already_matched = true; break
                            end
                        end
                        already_matched && break
                    end
                    already_matched && continue   # reverse ALREADY recorded by Case A above

                    for e_i = 1:em_iel.extra_nelem
                        θmin_i,θmax_i,ϕmin_i,ϕmax_i = iel_bounds[e_i]
                        inside = (θj >= θmin_i-tol) && (θj <= θmax_i+tol) &&
                                 (ϕj >= ϕmin_i-tol) && (ϕj <= ϕmax_i+tol)
                        inside || continue
                        ni = em_iel.extra_nop[e_i] + 1
                        iel_nodes = [em_iel.extra_connijk[e_i,p,q] for p=1:ni, q=1:ni] |> vec
                        key = (jel, ip_ext_j)
                        push!(get!(angular_deps, key, Tuple{Int,Vector{Int}}[]), (iel, iel_nodes))
                        break
                    end
                end
            end
        end
    end

    return angular_deps
end

function test_angular_dependency_graph(mesh, metrics)
    target_iel = 1
    target_e_ext = 4

    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    em_iel = mesh.extra_mesh[target_iel]

    verify_angular_phi_seam_conformity(mesh, 1)
    
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); break
    end
    em_jel = mesh.extra_mesh[jel]
    println("Testing iel=$target_iel (refined, $(em_iel.extra_npoin) directions) " *
            "vs jel=$jel (unrefined, $(em_jel.extra_npoin) directions)")

    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)


    test_angular_scheduling(mesh, metrics, neighbor_elem, angular_deps)

    test_angular_ncf_plumbing(mesh, metrics, neighbor_elem, angular_deps, 1, 2)

    n_singleton = 0; n_fine_to_coarse = 0; n_coarse_to_fine = 0; n_other = 0
    coarse_elem_size = (em_jel.extra_nop[1]+1)^2

    for ((e, d), targets) in angular_deps
        (e == target_iel || e == jel) || continue
        for (te, tnodes) in targets
            if (e == target_iel && te == jel) || (e==jel && te == target_iel)
	        n = length(tnodes)
            	if n == 1
		    n_singleton += 1
		    @info e, te, d, tnodes, n_singleton, em_iel.extra_coords[1,d], em_iel.extra_coords[2,d], em_jel.extra_coords[1,tnodes[1]], em_jel.extra_coords[2,tnodes[1]]
            	elseif e == target_iel && n == coarse_elem_size
                    n_fine_to_coarse += 1
            	elseif e == jel
                    n_coarse_to_fine += 1
            	else
                    n_other += 1
               	    println("  UNEXPECTED shape: (iel=$e, ip_ext=$d) -> (jel=$te, $(n) nodes)")
            	end
            end
    	end
    end

    max_plausible_fine = 4 * coarse_elem_size
    println("Singleton (ordinary 1-to-1): $n_singleton")
    println("Fine(iel)->whole coarse element (jel): $n_fine_to_coarse  " *
            "(plausible upper bound: $max_plausible_fine)  " *
            (n_fine_to_coarse <= max_plausible_fine ? "✓" : "✗ EXCEEDS BOUND"))
    println("Coarse(jel)->whole fine element (iel): $n_coarse_to_fine  " *
            "(plausible upper bound: $coarse_elem_size)  " *
            (n_coarse_to_fine <= coarse_elem_size ? "✓" : "✗ EXCEEDS BOUND"))
    println("Unexpected shapes: $n_other")

    return angular_deps
end

"""
    find_phi_seam_partner(coords_ang, connijk_ang, nop_ang, nelem_ang, e_ext)

Given one angular element, determines if it touches iel's OWN internal
ϕ=0/2π seam, and if so, finds its periodic PARTNER element (same θ range,
opposite ϕ side) — geometric matching, NOT index-based, reusing the exact
touches_2pi/touches_0 pattern from the CG conformity-enforcement code
(confirmed reliable here since i->θ, j->φ is the actual CONSTRUCTION
convention used in adapt_angular_grid_3Dby2D!'s own child-building loop,
not an assumption we're making after the fact).

Returns partner_e_ext (0 if this element doesn't touch the seam at all,
or no partner is currently found).
"""
function find_phi_seam_partner(coords_ang, connijk_ang, nop_ang, nelem_ang, e_ext)
    nop = nop_ang[e_ext]
    ϕ_first = coords_ang[2, connijk_ang[e_ext, 1, 1]]
    ϕ_last  = coords_ang[2, connijk_ang[e_ext, 1, nop+1]]
    touches_2pi = ϕ_last < ϕ_first
    touches_0   = ϕ_first < 1e-10 && ϕ_last > ϕ_first
    (touches_2pi || touches_0) || return 0

    θ_lo = coords_ang[1, connijk_ang[e_ext, 1, 1]]
    θ_hi = coords_ang[1, connijk_ang[e_ext, nop+1, 1]]

    for e_pair = 1:nelem_ang
        e_pair == e_ext && continue
        nop_p = nop_ang[e_pair]
        ϕ_pair_first = coords_ang[2, connijk_ang[e_pair, 1, 1]]
        ϕ_pair_last  = coords_ang[2, connijk_ang[e_pair, 1, nop_p+1]]
        pair_2pi = ϕ_pair_last < ϕ_pair_first
        pair_0   = ϕ_pair_first < 1e-10 && ϕ_pair_last > ϕ_pair_first
        θ_lo_p = coords_ang[1, connijk_ang[e_pair, 1, 1]]
        θ_hi_p = coords_ang[1, connijk_ang[e_pair, nop_p+1, 1]]
        if abs(θ_lo_p-θ_lo) < 1e-10 && abs(θ_hi_p-θ_hi) < 1e-10 &&
           ((touches_2pi && pair_0) || (touches_0 && pair_2pi))
            return e_pair
        end
    end
    return 0
end

"""
    verify_angular_phi_seam_conformity(mesh, iel; tol=1e-8)

PRECONDITION CHECK: for spatial element `iel`, verifies every angular
element touching the φ=0/2π seam has its geometric partner refined to
matching resolution (same nelem contribution / same node count) — fails
LOUDLY (error), not silently, if an asymmetric seam refinement slipped
through.
"""
function verify_angular_phi_seam_conformity(mesh, iel::Int; tol::Float64=1e-8)
    em = mesh.extra_mesh[iel]
    n_checked = 0
    for e = 1:em.extra_nelem
        partner = find_phi_seam_partner(em.extra_coords, em.extra_connijk, em.extra_nop, em.extra_nelem, e)
        partner == 0 && continue
        n_checked += 1
        nop_e = em.extra_nop[e]; nop_p = em.extra_nop[partner]
        if nop_e != nop_p
            error("verify_angular_phi_seam_conformity: iel=$iel has ASYMMETRIC φ-seam " *
                  "refinement — element $e (nop=$nop_e) and its partner $partner (nop=$nop_p) " *
                  "don't match. Refinement must be applied to both sides of the seam together.")
        end
    end
    return n_checked
end

function adapt_single_angular_element_with_metrics!(mesh, target_iel::Int, target_e_ext::Int, ψ, dψ)
    iel = target_iel
    e_ext = target_e_ext

    em = mesh.extra_mesh[iel]
    nop_ang     = em.extra_nop
    nelem_ang   = em.extra_nelem
    connijk_ang = em.extra_connijk
    coords_ang  = em.extra_coords
    Je_ang      = em.extra_metrics.Je
    dξdx_ang    = em.extra_metrics.dξdx
    dxdξ_ang    = em.extra_metrics.dxdξ
    dξdy_ang    = em.extra_metrics.dξdy
    dydξ_ang    = em.extra_metrics.dydξ
    dηdy_ang    = em.extra_metrics.dηdy
    dydη_ang    = em.extra_metrics.dydη
    dηdx_ang    = em.extra_metrics.dηdx
    dxdη_ang    = em.extra_metrics.dxdη

    θmax = coords_ang[1, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    θmin = coords_ang[1, connijk_ang[e_ext, 1, 1]]
    ϕmax = coords_ang[2, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    ϕmin = coords_ang[2, connijk_ang[e_ext, 1, 1]]
    ϕmax == 0 && (ϕmax = 2π)
    θ12 = (θmax + θmin) / 2
    ϕ12 = (ϕmax + ϕmin) / 2

    nelem_ang_new = nelem_ang + 3
    npoin_ang_new_upper = em.extra_npoin + (4*(nop_ang[e_ext]+1)^2 - 4*(nop_ang[e_ext]+1))

    connijk_ang_new = zeros(Int, nelem_ang_new, nop_ang[1]+1, nop_ang[1]+1)
    coords_new      = zeros(Float64, 2, npoin_ang_new_upper)
    metrics_new     = allocate_metrics(NSD_2D(), nelem_ang_new, 0, nop_ang[e_ext], TFloat, CPU())
    nop_ang_new     = zeros(Int, nelem_ang_new)
    nop_ang_new[1:nelem_ang_new-1] .= nop_ang[e_ext]
    nop_ang_new[nelem_ang_new]      = nop_ang[1]

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

    lgl = basis_structs_ξ_ω!(LGL(), nop_ang[1], CPU())

    # ── Copy elements before e_ext, INCLUDING their metrics ──────────────────
    if e_ext > 1
        for e_ext1 = 1:e_ext-1
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                metrics_new.dxdξ[e_ext1,i,j] = dxdξ_ang[e_ext1,i,j]
                metrics_new.Je[e_ext1,i,j]   = Je_ang[e_ext1,i,j]
                metrics_new.dξdx[e_ext1,i,j] = dξdx_ang[e_ext1,i,j]
                metrics_new.dxdη[e_ext1,i,j] = dxdη_ang[e_ext1,i,j]
                metrics_new.dηdx[e_ext1,i,j] = dηdx_ang[e_ext1,i,j]
                metrics_new.dηdy[e_ext1,i,j] = dηdy_ang[e_ext1,i,j]
                metrics_new.dξdy[e_ext1,i,j] = dξdy_ang[e_ext1,i,j]
                metrics_new.dydη[e_ext1,i,j] = dydη_ang[e_ext1,i,j]
                metrics_new.dydξ[e_ext1,i,j] = dydξ_ang[e_ext1,i,j]
                insert_point!(θ, ϕ, connijk_ang_new, e_ext1, i, j)
            end
        end
    end

    child_ranges = [
        (θmin, θ12, ϕmin, ϕ12, e_ext),
        (θ12, θmax, ϕmin, ϕ12, e_ext+1),
        (θmin, θ12, ϕ12, ϕmax, e_ext+2),
        (θ12, θmax, ϕ12, ϕmax, e_ext+3),
    ]
    for (θlo, θhi, ϕlo, ϕhi, e_child) in child_ranges
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

    # ── Copy elements after e_ext, INCLUDING their metrics ───────────────────
    if e_ext < nelem_ang
        for e_ext1 = e_ext+1:nelem_ang
            dst = e_ext1 + 3
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                metrics_new.dxdξ[dst,i,j] = dxdξ_ang[e_ext1,i,j]
                metrics_new.Je[dst,i,j]   = Je_ang[e_ext1,i,j]
                metrics_new.dξdx[dst,i,j] = dξdx_ang[e_ext1,i,j]
                metrics_new.dxdη[dst,i,j] = dxdη_ang[e_ext1,i,j]
                metrics_new.dηdx[dst,i,j] = dηdx_ang[e_ext1,i,j]
                metrics_new.dηdy[dst,i,j] = dηdy_ang[e_ext1,i,j]
                metrics_new.dξdy[dst,i,j] = dξdy_ang[e_ext1,i,j]
                metrics_new.dydη[dst,i,j] = dydη_ang[e_ext1,i,j]
                metrics_new.dydξ[dst,i,j] = dydξ_ang[e_ext1,i,j]
                insert_point!(θ, ϕ, connijk_ang_new, dst, i, j)
            end
        end
    end

    npoin_ang_new = iter - 1
    coords_new = coords_new[:, 1:npoin_ang_new]

    # ── Build metrics for the 4 NEW child elements (EXACT logic from
    # adapt_angular_grid_3Dby2D!, restored — this is the piece that was
    # missing). ────────────────────────────────────────────────────────────
    for e_ext1 = e_ext:e_ext+3
        for i = 1:nop_ang_new[e_ext1]+1, j = 1:nop_ang_new[e_ext1]+1
            ip  = connijk_ang_new[e_ext1,i,j]
            θij = coords_new[1,ip]; ϕij = coords_new[2,ip]
            for l = 1:nop_ang_new[e_ext1]+1
                for k = 1:nop_ang_new[e_ext1]+1
                    a = dψ[i,k]*ψ[j,l]; b = ψ[i,k]*dψ[j,l]
                    metrics_new.dxdξ[e_ext1,k,l] += a*θij
                    metrics_new.dxdη[e_ext1,k,l] += b*θij
                    metrics_new.dydξ[e_ext1,k,l] += a*ϕij
                    metrics_new.dydη[e_ext1,k,l] += b*ϕij
                end
            end
        end
        for l = 1:nop_ang_new[e_ext1]+1, k = 1:nop_ang_new[e_ext1]+1
            dxdξ_v = metrics_new.dxdξ[e_ext1,k,l]; dydη_v = metrics_new.dydη[e_ext1,k,l]
            dydξ_v = metrics_new.dydξ[e_ext1,k,l]; dxdη_v = metrics_new.dxdη[e_ext1,k,l]
            ip  = connijk_ang_new[e_ext1,k,l]
            θ = coords_new[1,ip]
            Je_v   = dxdξ_v*dydη_v - dydξ_v*dxdη_v
            Jinv   = 1.0 / Je_v
            metrics_new.Je[e_ext1,k,l]   = sin(θ) * Je_v
            metrics_new.dξdx[e_ext1,k,l] =  dydη_v*Jinv
            metrics_new.dξdy[e_ext1,k,l] = -dxdη_v*Jinv
            metrics_new.dηdx[e_ext1,k,l] = -dydξ_v*Jinv
            metrics_new.dηdy[e_ext1,k,l] =  dxdξ_v*Jinv
        end
    end

    # ── ϕ-periodicity merge (unchanged from before, but now ALSO needs to
    # keep metrics arrays in sync — metrics are PER-ELEMENT, indexed by
    # (e,i,j), not by global point id, so NO remapping needed here at all;
    # only coords_new/connijk_ang_new (point-indexed) were ever touched by
    # the merge). ─────────────────────────────────────────────────────────
    zero_ϕ_dict = Dict{Float64, Int}()
    for iper = 1:npoin_ang_new
        if coords_new[2, iper] < 1e-10
            zero_ϕ_dict[coords_new[1, iper]] = iper
        end
    end
    merge_pairs = Dict{Int, Int}()
    for iper = 1:npoin_ang_new
        ϕ = coords_new[2, iper]
        if abs(ϕ / π - 2.0) < 1e-10
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
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] == ip_old && (connijk_ang_new[e,i,j] = ip_new)
                end
            end
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] > ip_old && (connijk_ang_new[e,i,j] -= 1)
                end
            end
            merge_pairs = Dict(
                (k > ip_old ? k - 1 : k) => (v > ip_old ? v - 1 : v)
                    for (k, v) in merge_pairs if k != ip_old
                )
            for i = ip_old+1:npoin_ang_new
                coords_new[1,i-1] = coords_new[1,i]
                coords_new[2,i-1] = coords_new[2,i]
            end
            npoin_ang_new -= 1
        end
    end
    coords_new = coords_new[:, 1:npoin_ang_new]

    em.extra_connijk = connijk_ang_new
    em.extra_coords   = coords_new
    em.extra_nop      = nop_ang_new
    em.extra_nelem    = nelem_ang_new
    em.extra_npoin    = npoin_ang_new
    em.extra_metrics  = metrics_new

    return nothing
end

function build_angular_scheduling_state(mesh, neighbor_elem, angular_deps)
    # ── Per-element direction lists: which ip_ext actually EXIST at iel ─────
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    # ── Precompute Ω for every (iel, ip_ext) pair — direction-specific,
    # element-specific (since a refined element's ip_ext doesn't correspond
    # to the same physical Ω as an unrefined neighbor's same-index ip_ext). ──
    Ω_of = Vector{Vector{NTuple{3,Float64}}}(undef, mesh.nelem)
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        Ω_of[iel] = Vector{NTuple{3,Float64}}(undef, em.extra_npoin)
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ω_of[iel][d] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        end
    end

    return (n_dirs=n_dirs, Ω_of=Ω_of)
end

function test_angular_scheduling(mesh, metrics, neighbor_elem, angular_deps; max_rounds::Int = 500)
    state = build_angular_scheduling_state(mesh, neighbor_elem, angular_deps)
    ngl = mesh.ngl

    solved = [falses(state.n_dirs[iel]) for iel = 1:mesh.nelem]

    total_pairs = sum(state.n_dirs)
    n_solved_total = 0
    round = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("test_angular_scheduling: exceeded max_rounds=$max_rounds — " *
                                     "$(total_pairs - n_solved_total) (iel,ip_ext) pairs still " *
                                     "unresolved, likely a real cycle or a graph-construction bug.")

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:state.n_dirs[iel]
                solved[iel][d] && continue
                Ωx,Ωy,Ωz = state.Ω_of[iel][d]

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel = neighbor_elem[iel][fid]
                    jel == zero(TInt) && continue   # true boundary — always satisfied
                    jel = Int(jel)

                    # Angular-dependency-aware: check angular_deps FIRST,
                    # fall back to ordinary same-index match only if no
                    # explicit entry exists for this (iel,d).
                    key = (iel, d)
                    if haskey(angular_deps, key)
                        found_this_jel = false
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            found_this_jel = true
                            for dn in dep_nodes
                                if !solved[jel][dn]
                                    ready = false
                                end
                            end
                        end
                        # If jel isn't even in this direction's angular_deps
                        # list (e.g. THIS face has identical meshes, but a
                        # DIFFERENT face of iel is what triggered an entry),
                        # fall back to ordinary index match for THIS face.
                        if !found_this_jel
                            (d <= state.n_dirs[jel] && !solved[jel][d]) && (ready = false)
                        end
                    else
                        (d <= state.n_dirs[jel] && !solved[jel][d]) && (ready = false)
                    end
                    ready || break
                end

                if ready
                    solved[iel][d] = true
                    n_newly_solved += 1
                end
            end
        end

        n_solved_total += n_newly_solved
        println("round $round: $n_newly_solved newly solved, $n_solved_total/$total_pairs total")

        if n_newly_solved == 0 && n_solved_total < total_pairs
            error("test_angular_scheduling: STALLED — 0 progress this round with " *
                  "$(total_pairs - n_solved_total) pairs remaining. Real cycle or bug.")
        end
    end

    println("Converged in $round rounds, all $total_pairs (element,direction) pairs solved.")
    return solved
end

function test_angular_ncf_plumbing(mesh, metrics, neighbor_elem, angular_deps, target_iel::Int, target_jel::Int)
    em_iel = mesh.extra_mesh[target_iel]
    em_jel = mesh.extra_mesh[target_jel]

    # ── Find one (iel, d) pair whose dependency on jel is a "whole coarse
    # element" group (fine depends on coarse) — i.e. a list length matching
    # some jel angular element's node count, NOT a singleton. ───────────────
    sample_d = 0
    sample_jel_nodes = Int[]
    for d = 1:em_iel.extra_npoin
        key = (target_iel, d)
        haskey(angular_deps, key) || continue
        for (dep_jel, dep_nodes) in angular_deps[key]
            dep_jel == target_jel || continue
            length(dep_nodes) > 1 || continue   # skip singletons, want the group case
            sample_d = d
            sample_jel_nodes = dep_nodes
            break
        end
        sample_d != 0 && break
    end

    if sample_d == 0
        println("⚠ No fine->coarse group dependency found from iel=$target_iel to jel=$target_jel")
        return nothing
    end
    println("Testing (iel=$target_iel, d=$sample_d) -> jel=$target_jel nodes=$sample_jel_nodes")

    # ── Independently identify which jel ANGULAR ELEMENT contains these
    # exact nodes, by searching connijk_ang directly — NOT reusing any
    # bookkeeping from angular_deps itself, so this is a genuine
    # cross-check, not a restatement of the same computation. ───────────────
    coarse_e_ext = 0
    for e = 1:em_jel.extra_nelem
        n = em_jel.extra_nop[e] + 1
        elem_nodes = Set(em_jel.extra_connijk[e,a,b] for a=1:n, b=1:n)
        if elem_nodes == Set(sample_jel_nodes)
            coarse_e_ext = e
            break
        end
    end
    if coarse_e_ext == 0
        println("✗ FAIL: angular_deps' node list for (iel=$target_iel,d=$sample_d) does NOT " *
                "exactly match any single angular element's node set at jel=$target_jel")
        return nothing
    end
    println("✓ Confirmed: jel's coarse angular element is e_ext=$coarse_e_ext")

    # ── Identify iel's OWN containing angular element for direction d ───────
    fine_e_ext = 0
    for e = 1:em_iel.extra_nelem
        n = em_iel.extra_nop[e] + 1
        for a=1:n, b=1:n
            if em_iel.extra_connijk[e,a,b] == sample_d
                fine_e_ext = e
                break
            end
        end
        fine_e_ext != 0 && break
    end
    if fine_e_ext == 0
        println("✗ FAIL: could not find direction d=$sample_d in ANY of iel's angular elements " *
                "(connijk_ang inconsistency)")
        return nothing
    end
    println("✓ Confirmed: iel's fine angular element is e_ext=$fine_e_ext")

    # ── Now actually build the interpolation weights and check the SAME
    # field-reproduction test as before, but entered through the SCHEDULER'S
    # OWN node identification rather than build_angular_ncf_correspondence's
    # element-bounds matching — a genuine independent path. ──────────────────
    Lθ, Lϕ = build_angular_ncf_parent_to_child_weights(
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, coarse_e_ext,
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, fine_e_ext)

    ncoarse = em_jel.extra_nop[coarse_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(em_jel.extra_coords, em_jel.extra_connijk,
                                                       em_jel.extra_nop, coarse_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ

    P(t) = sum(0.3^k * t^k for k=0:ncoarse-1)
    field(θ,ϕ) = P((θ-θmin)/(θmax-θmin)) + P((ϕ-ϕmin)/(ϕmax-ϕmin))

    coarse_grid = zeros(ncoarse,ncoarse)
    for i=1:ncoarse, j=1:ncoarse
        ip = em_jel.extra_connijk[coarse_e_ext,i,j]
        θ,ϕ = em_jel.extra_coords[1,ip], em_jel.extra_coords[2,ip]
        coarse_grid[i,j] = field(θ, fix_ϕ(ϕ))
    end

    fine_computed = evaluate_angular_ncf_parent_to_child(coarse_grid, Lθ, Lϕ)

    # Find (a,b) of sample_d within the fine element to index fine_computed
    nfine = em_iel.extra_nop[fine_e_ext] + 1
    ai,bi = 0,0
    for a=1:nfine, b=1:nfine
        if em_iel.extra_connijk[fine_e_ext,a,b] == sample_d
            ai,bi = a,b; break
        end
    end

    θd, ϕd = em_iel.extra_coords[1,sample_d], em_iel.extra_coords[2,sample_d]
    exact = field(θd, fix_ϕ(ϕd))
    err = abs(fine_computed[ai,bi] - exact)
    println("Plumbing round-trip error at d=$sample_d: $err  " * (err < 1e-9 ? "✓ PASS" : "✗ FAIL"))
    return err
end

# ── build_angular_group_upwind (parent->child), corrected: ─────────────────
function build_angular_group_upwind(iel::Int, fid::Int, jel::Int, coarse_e_ext::Int, fine_e_ext::Int,
                                     d::Int, mesh, em_jel, em_iel, x_dg, face_correspondence, ngl::Int)
    ncoarse = em_jel.extra_nop[coarse_e_ext] + 1
    coarse_ips = [em_jel.extra_connijk[coarse_e_ext,a,b] for a=1:ncoarse, b=1:ncoarse]

    Lθ, Lϕ = build_angular_ncf_parent_to_child_weights(
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, coarse_e_ext,
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, fine_e_ext)

    nfine = em_iel.extra_nop[fine_e_ext] + 1
    ai, bi = 0, 0
    for a=1:nfine, b=1:nfine
        em_iel.extra_connijk[fine_e_ext,a,b] == d && (ai,bi = a,b)
    end

    return (i,j,k) -> begin
        (ij,jj,kj) = face_correspondence[(iel,fid)][(i,j,k)]   # jel's LOCAL spatial position
        coarse_grid = zeros(ncoarse, ncoarse)
        for a = 1:ncoarse, b = 1:ncoarse
            ip = coarse_ips[a,b]
            coarse_grid[a,b] = x_dg[(jel,ip)][lidx(ij,jj,kj,ngl)]
        end
        fine_computed = evaluate_angular_ncf_parent_to_child(coarse_grid, Lθ, Lϕ)
        fine_computed[ai,bi]
    end
end

"""
Coarse-side (child->parent) angular closure: `iel`/`fid` here are the
COARSE element's OWN identity and OWN local face id (matching the
convention every other closure in this codebase uses — `(i,j,k)` fed to
the returned closure is the COARSE element's own face-node index). For
each such spatial face node, gathers ALL relevant fine children's angular
grids at the CORRESPONDING spatial position on each child (via the
ALREADY-VALIDATED face_correspondence — never assumes aligned local axes
across elements), direct-evaluates each child's contribution at its OWN
covering angular region (no projection — same lesson as the spatial case),
producing the coarse direction's value there.
"""
function build_angular_coarse_upwind(iel::Int, fid::Int, fine_elements::Vector{Int},
                                      coarse_e_ext::Int, d::Int,
                                      mesh, em_coarse, x_dg, face_correspondence, ngl::Int)
    ncoarse = em_coarse.extra_nop[coarse_e_ext] + 1
    pi_ = 0
    pj_ = 0
    for a = 1:ncoarse
    	for b = 1:ncoarse
            if em_coarse.extra_connijk[coarse_e_ext,a,b] == d
            	pi_ = a
            	pj_ = b
            end
    	end
    end    

    θp = em_coarse.extra_coords[1, em_coarse.extra_connijk[coarse_e_ext,pi_,pj_]]
    ϕp = em_coarse.extra_coords[2, em_coarse.extra_connijk[coarse_e_ext,pi_,pj_]]

    # Precompute, ONCE, which fine element(s) at each fel correspond to
    # THIS coarse element (geometry-only, Ω-independent — safe to build
    # outside the returned closure, reusing the ALREADY-VALIDATED
    # build_angular_ncf_correspondence).
    per_child = Tuple{Int,Vector{Int}}[]
    for fel in fine_elements
        em_fel = mesh.extra_mesh[fel]
        matches = build_angular_ncf_correspondence(
            em_fel.extra_coords, em_fel.extra_connijk, em_fel.extra_nop, em_fel.extra_nelem,
            em_coarse.extra_coords, em_coarse.extra_connijk, em_coarse.extra_nop, em_coarse.extra_nelem)
        m = findfirst(mm -> mm.coarse_e_ext_at_jel == coarse_e_ext, matches)
        m === nothing && continue
        push!(per_child, (fel, matches[m].fine_e_exts_at_iel))
    end

    return (i,j,k) -> begin
        result = NaN
        for (fel, fine_e_exts) in per_child
            (ij,jj,kj) = face_correspondence[(iel,fid)][(i,j,k)]   # fel's LOCAL spatial position

            for fine_e_ext in fine_e_exts
                em_fel = mesh.extra_mesh[fel]
                θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(
                    em_fel.extra_coords, em_fel.extra_connijk, em_fel.extra_nop, fine_e_ext)
                near_seam = isapprox(ϕmax, 2π; atol=1e-6)
                ϕp_fix = (near_seam && ϕp < 1e-8) ? 2π : ϕp
                (θp < θmin-1e-8 || θp > θmax+1e-8 || ϕp_fix < ϕmin-1e-8 || ϕp_fix > ϕmax+1e-8) && continue

                nfine = em_fel.extra_nop[fine_e_ext] + 1
                θ_fine = [em_fel.extra_coords[1, em_fel.extra_connijk[fine_e_ext,a,1]] for a=1:nfine]
                ϕ_fine_raw = [em_fel.extra_coords[2, em_fel.extra_connijk[fine_e_ext,1,b]] for b=1:nfine]
                ϕ_fine = [(near_seam && ϕ<1e-8) ? 2π : ϕ for ϕ in ϕ_fine_raw]
                Lθ = build_1d_lagrange_matrix(θ_fine, [θp])
                Lϕ = build_1d_lagrange_matrix(ϕ_fine, [ϕp_fix])

                child_grid = zeros(nfine,nfine)
                for a=1:nfine, b=1:nfine
                    ip = em_fel.extra_connijk[fine_e_ext,a,b]
                    child_grid[a,b] = x_dg[(fel,ip)][lidx(ij,jj,kj,ngl)]
                end
                val = 0.0
                for a=1:nfine, b=1:nfine
                    val += Lθ[1,a]*Lϕ[1,b]*child_grid[a,b]
                end
                result = val
                break
            end
            !isnan(result) && break
        end
        isnan(result) && error("build_angular_coarse_upwind: no covering fine element found for (i,j,k)=($i,$j,$k)")
        return result
    end
end

function test_angular_coarse_upwind_isolated(mesh, target_iel::Int, target_e_ext::Int)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    em_iel = mesh.extra_mesh[target_iel]

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0; test_fid = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); test_fid = fid; break
    end
    em_jel = mesh.extra_mesh[jel]

    matches = build_angular_ncf_correspondence(
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, em_iel.extra_nelem,
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, em_jel.extra_nelem)
    println("Found $(length(matches)) coarse/fine group(s) against jel=$jel")

    ngl = mesh.ngl
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    coarse_fid = Int(neighbor_fid[target_iel][test_fid])   # jel's OWN local face id for this interface

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly=mesh.ymax-mesh.ymin; Lz=mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k*t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    overall_max_err = 0.0

    for (mi, m) in enumerate(matches)
        coarse_e_ext = m.coarse_e_ext_at_jel
        fine_e_exts  = m.fine_e_exts_at_iel
        println("─"^50)
        println("Testing match $mi: coarse e_ext=$coarse_e_ext <-> fine e_exts=$fine_e_exts")

        # ── Fabricate x_dg for the fine element (target_iel), full ngl^3
        # spatial grid, at every angular node in fine_e_exts. ────────────────
        x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for e in fine_e_exts
            n = em_iel.extra_nop[e]+1
            for a=1:n, b=1:n
                ip = em_iel.extra_connijk[e,a,b]
                haskey(x_dg, (target_iel,ip)) && continue
                vals = zeros(ngl^3)
                for k=1:ngl,j=1:ngl,i=1:ngl
                    gip = mesh.connijk[target_iel,i,j,k]
                    vals[lidx(i,j,k,ngl)] = field(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                end
                x_dg[(target_iel,ip)] = vals
            end
        end

        # ── Build the closure with jel/coarse_fid as the COARSE element's
        # own identity — matching build_angular_coarse_upwind's convention. ──
        ncoarse = em_jel.extra_nop[coarse_e_ext] + 1
        max_err = 0.0
        for a = 1:min(3,ncoarse), b = 1:min(3,ncoarse)
            d = em_jel.extra_connijk[coarse_e_ext, a, b]
            closure = build_angular_coarse_upwind(jel, coarse_fid, [target_iel], coarse_e_ext, d,
                                                   mesh, em_jel, x_dg, face_correspondence, ngl)

            # Sample nodes from jel's OWN face convention (coarse_fid) —
            # matching what the closure now expects as (i,j,k).
            irange,jrange,krange = face_ip_range(coarse_fid, ngl)
            for (i,j,k) in Iterators.take(Iterators.product(irange,jrange,krange), 3)
                computed = closure(i,j,k)
                gip = mesh.connijk[jel,i,j,k]   # jel's OWN node — no correspondence needed here,
                                                  # (i,j,k) IS already jel's own local index
                exact = field(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                max_err = max(max_err, abs(computed - exact))
            end
        end
        println("  match $mi max abs error = $max_err")
        overall_max_err = max(overall_max_err, max_err)
    end

    println("─"^50)
    println("Overall max abs error across all matches: $overall_max_err")
    return overall_max_err
end

function test_angular_coarse_upwind_isolated_v2(mesh, target_iel::Int, target_e_ext::Int)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    em_iel = mesh.extra_mesh[target_iel]

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0; test_fid = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); test_fid = fid; break
    end
    em_jel = mesh.extra_mesh[jel]

    matches = build_angular_ncf_correspondence(
        em_iel.extra_coords, em_iel.extra_connijk, em_iel.extra_nop, em_iel.extra_nelem,
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, em_jel.extra_nelem)
    println("Found $(length(matches)) coarse/fine group(s) against jel=$jel")

    ngl = mesh.ngl
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    coarse_fid = Int(neighbor_fid[target_iel][test_fid])

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly=mesh.ymax-mesh.ymin; Lz=mesh.zmax-mesh.zmin
    Pspatial(t) = sum(0.3^k*t^k for k=0:ngl-1)

    overall_max_err = 0.0

    for (mi, m) in enumerate(matches)
        coarse_e_ext = m.coarse_e_ext_at_jel
        fine_e_exts  = m.fine_e_exts_at_iel
        println("─"^50)
        println("Testing match $mi: coarse e_ext=$coarse_e_ext <-> fine e_exts=$fine_e_exts")

        # ── Well-scaled angular field, normalized against THIS coarse
        # element's own (θ,ϕ) range — same discipline as the earlier angular
        # interpolation tests. Now genuinely varies with direction, so a
        # wrong Lθ/Lϕ weight would show up as real error, not cancel. ───────
        θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(
            em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, coarse_e_ext)
        near_seam = isapprox(ϕmax, 2π; atol=1e-6)
        fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ
        ncoarse = em_jel.extra_nop[coarse_e_ext] + 1
        Pang(t) = sum(0.3^k*t^k for k=0:ncoarse-1)
        field(x,y,z,θ,ϕ) = Pspatial((x-x0)/Lx) + Pspatial((y-y0)/Ly) + Pspatial((z-z0)/Lz) +
                            Pang((θ-θmin)/(θmax-θmin)) + Pang((fix_ϕ(ϕ)-ϕmin)/(ϕmax-ϕmin))

        # ── Fabricate x_dg: now genuinely direction-dependent per angular
        # node, not spatially-identical across the fine element's directions. ──
        x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for e in fine_e_exts
            n = em_iel.extra_nop[e]+1
            for a=1:n, b=1:n
                ip = em_iel.extra_connijk[e,a,b]
                haskey(x_dg, (target_iel,ip)) && continue
                θd = em_iel.extra_coords[1,ip]; ϕd = em_iel.extra_coords[2,ip]
                vals = zeros(ngl^3)
                for k=1:ngl,j=1:ngl,i=1:ngl
                    gip = mesh.connijk[target_iel,i,j,k]
                    vals[lidx(i,j,k,ngl)] = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θd, ϕd)
                end
                x_dg[(target_iel,ip)] = vals
            end
        end

        max_err = 0.0
        for a = 1:ncoarse, b = 1:ncoarse   # ALL coarse nodes now, not just a sample
            d = em_jel.extra_connijk[coarse_e_ext, a, b]
            θd = em_jel.extra_coords[1,d]; ϕd = em_jel.extra_coords[2,d]

            closure = build_angular_coarse_upwind(jel, coarse_fid, [target_iel], coarse_e_ext, d,
                                                   mesh, em_jel, x_dg, face_correspondence, ngl)

            irange,jrange,krange = face_ip_range(coarse_fid, ngl)
            for (i,j,k) in Iterators.product(irange,jrange,krange)   # ALL face nodes, not just 3
                computed = closure(i,j,k)
                gip = mesh.connijk[jel,i,j,k]
                exact = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θd, ϕd)
                max_err = max(max_err, abs(computed - exact))
            end
        end
        println("  match $mi max abs error = $max_err")
        overall_max_err = max(overall_max_err, max_err)
    end

    println("─"^50)
    println("Overall max abs error across all matches: $overall_max_err")
    return overall_max_err
end

function find_containing_angular_element(em, node_set::Set{Int})
    found_e = 0
    for e = 1:em.extra_nelem
        n = em.extra_nop[e] + 1
        elem_nodes = Set{Int}()
        for a = 1:n
            for b = 1:n
                push!(elem_nodes, em.extra_connijk[e,a,b])
            end
        end
        elem_nodes == node_set && (found_e = e)
    end
    return found_e
end

"""
UNIFIED angular upwind closure: given ONE source angular element (at jel,
identified from angular_deps' dep_nodes — always exactly one element,
whether coarse-serving-fine or fine-serving-coarse), interpolates that
element's full grid to iel's single target direction d, evaluated PER
SPATIAL FACE NODE via face_correspondence (never assumes aligned local
axes across elements).
"""
function build_angular_ncf_upwind_unified(iel::Int, fid::Int, jel::Int, d::Int,
                                           src_e_ext::Int, mesh, em_jel,
                                           x_dg, face_correspondence, ngl::Int)
    em_iel = mesh.extra_mesh[iel]
    θd = 0.0; ϕd = 0.0
    for e = 1:em_iel.extra_nelem
        n = em_iel.extra_nop[e] + 1
        for a = 1:n
            for b = 1:n
                if em_iel.extra_connijk[e,a,b] == d
                    θd = em_iel.extra_coords[1, d]
                    ϕd = em_iel.extra_coords[2, d]
                end
            end
        end
    end

    nsrc = em_jel.extra_nop[src_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd

    θ_src = Vector{Float64}(undef, nsrc)
    for a = 1:nsrc
        θ_src[a] = em_jel.extra_coords[1, em_jel.extra_connijk[src_e_ext,a,1]]
    end
    ϕ_src = Vector{Float64}(undef, nsrc)
    for b = 1:nsrc
        raw = em_jel.extra_coords[2, em_jel.extra_connijk[src_e_ext,1,b]]
        ϕ_src[b] = (near_seam && raw < 1e-8) ? 2π : raw
    end

    Lθ = build_1d_lagrange_matrix(θ_src, [θd])
    Lϕ = build_1d_lagrange_matrix(ϕ_src, [ϕd_fix])

    src_ips = Matrix{Int}(undef, nsrc, nsrc)
    for a = 1:nsrc
        for b = 1:nsrc
            src_ips[a,b] = em_jel.extra_connijk[src_e_ext,a,b]
        end
    end

    return (i,j,k) -> begin
        corr = face_correspondence[(iel,fid)][(i,j,k)]
        ii2 = corr[1]; jj2 = corr[2]; kk2 = corr[3]
        val = 0.0
        for a = 1:nsrc
            for b = 1:nsrc
                ip = src_ips[a,b]
                val += Lθ[1,a]*Lϕ[1,b]*x_dg[(jel,ip)][lidx(ii2,jj2,kk2,ngl)]
            end
        end
        val
    end
end

#=function run_angular_ncf_sweep!(mesh, metrics, ψ, dψ, ω,
                                 neighbor_elem, neighbor_fid, face_correspondence,
                                 angular_deps,
                                 κ_val::Float64, σ_val::Float64, I0::Float64;
                                 max_rounds::Int = 500, max_cycle_iter::Int = 100, cycle_tol::Float64 = 1e-9)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    round = 0
    n_cycles_resolved = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("run_angular_ncf_sweep!: exceeded max_rounds=$max_rounds — " *
                                     "$(total_pairs-n_solved_total) unresolved")

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                upwind_fns = Dict{Int, Function}()
                r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]
                    if jel_raw == zero(TInt)
                        upwind_fns[fid] = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]
                                    ready = false
                                else
                                    dj_c = dj; fid_c = fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved
                                    ready = false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c = d; fid_c = fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                npoin_cg = length(mesh.x)
                κ = fill(κ_val, npoin_cg); σz = fill(σ_val, npoin_cg); intΦ_eff = zeros(npoin_cg)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σz, intΦ_eff, mesh.connijk)
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns, fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end

        n_solved_total += n_newly_solved
        println("round $round: $n_newly_solved solved, $n_solved_total/$total_pairs")

        # ── STALL HANDLING: detect and resolve cycles, then let the loop
        # continue rather than error. ────────────────────────────────────────
        if n_newly_solved == 0 && n_solved_total < total_pairs
            cyclic_sccs = find_cyclic_sccs(mesh, metrics, neighbor_elem, angular_deps, solved, n_dirs)
            if isempty(cyclic_sccs)
                error("run_angular_ncf_sweep!: STALLED with $(total_pairs-n_solved_total) unsolved " *
                      "and NO cycle detected — genuine deadlock, likely a real bug.")
            end
            println("  Stall at round $round: found $(length(cyclic_sccs)) cyclic SCC(s), resolving...")
            for scc in cyclic_sccs
                resolve_angular_cycle!(mesh, metrics, ψ, dψ, ω, neighbor_elem, neighbor_fid,
                                       face_correspondence, angular_deps, x_dg, solved, n_dirs,
                                       scc, κ_val, σ_val, I0; max_iter=max_cycle_iter, tol=cycle_tol)
                n_cycles_resolved += 1
                n_solved_total += length(scc)
            end
            println("  Resolved $(length(cyclic_sccs)) cycle(s), resuming round-based solve")
        end
    end

    println("Converged in $round rounds, $n_cycles_resolved cycle(s) resolved locally")
    return x_dg, n_dirs
end=#

function run_angular_ncf_sweep!(mesh, metrics, ψ, dψ, ω,
                                 neighbor_elem, neighbor_fid, face_correspondence,
                                 angular_deps, κ, σ, intΦ_eff, bc_fn_all;
                                 extra_source_all = nothing,
                                 max_rounds::Int = 500, max_cycle_iter::Int = 100, cycle_tol::Float64 = 1e-9)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    round = 0
    n_cycles_resolved = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("run_angular_ncf_sweep!: exceeded max_rounds=$max_rounds — " *
                                     "$(total_pairs-n_solved_total) unresolved")

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                upwind_fns = Dict{Int, Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]
                    if jel_raw == zero(TInt)
                        bc = bc_fn_all[(iel,d)]
                        upwind_fns[fid] = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]
                                    ready = false
                                else
                                    dj_c = dj; fid_c = fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved
                                    ready = false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c = d; fid_c = fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns, fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end

        n_solved_total += n_newly_solved

        if n_newly_solved == 0 && n_solved_total < total_pairs
            cyclic_sccs = find_cyclic_sccs(mesh, metrics, neighbor_elem, angular_deps, solved, n_dirs)
            if isempty(cyclic_sccs)
                error("run_angular_ncf_sweep!: STALLED with $(total_pairs-n_solved_total) unsolved " *
                      "and NO cycle detected — genuine deadlock.")
            end
            for scc in cyclic_sccs
                resolve_angular_cycle!(mesh, metrics, ψ, dψ, ω, neighbor_elem, neighbor_fid,
                                       face_correspondence, angular_deps, x_dg, solved, n_dirs,
                                       scc, κ, σ, intΦ_eff, bc_fn_all;
                                       extra_source_all=extra_source_all,
                                       max_iter=max_cycle_iter, tol=cycle_tol)
                n_cycles_resolved += 1
                n_solved_total += length(scc)
            end
        end
    end

    return x_dg, n_dirs
end

function test_angular_ncf_exponential(mesh, metrics, basis, ω_gll; target_iel=442, target_e_ext=4,
                                       κ_val=0.0001, I0=4.2)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    ψ, dψ = basis.ψ, basis.dψ
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            bc_fn_all[(iel,d)] = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1])+Ωy*(y-r_bc[2])+Ωz*(z-r_bc[3])))
        end
    end

    x_dg, n_dirs = run_angular_ncf_sweep!(mesh, metrics, ψ, dψ, ω_gll,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          angular_deps, κ_val, 0.0, 0.0, bc_fn_all)

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    max_err = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                computed = x_dg[(iel,d)][lidx(i,j,k,mesh.ngl)]
                max_err = max(max_err, abs(computed-exact))
            end
        end
    end
    println("Angular NCF exponential test — max abs error: $max_err")
    return max_err
end

function debug_angular_ncf_error_location(mesh, x_dg, n_dirs, κ_val, I0, target_iel, jel_list)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    ngl = mesh.ngl

    for iel in vcat([target_iel], jel_list)
        em = mesh.extra_mesh[iel]
        max_err_iel = 0.0
        worst_d = 0
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                computed = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                err = abs(computed-exact)
                if err > max_err_iel
                    max_err_iel = err; worst_d = d
                end
            end
        end
        println("iel=$iel: max error = $max_err_iel (worst direction d=$worst_d)")
    end
end

function test_angular_ncf_error_location(mesh, metrics, basis, ω_gll; target_iel=1, target_e_ext=4,
                                          κ_val=0.0001, I0=4.2)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    ψ, dψ = basis.ψ, basis.dψ
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    x_dg, n_dirs = run_angular_ncf_sweep!(mesh, metrics, ψ, dψ, ω_gll,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          angular_deps, κ_val, 0.0, I0)

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    ngl = mesh.ngl

    # ── Same-rank spatial neighbors of target_iel ────────────────────────────
    same_rank_neighbors = Int[]
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        push!(same_rank_neighbors, Int(cand))
    end

    # ── A handful of elements FAR from target_iel — just pick evenly spaced
    # indices across the whole mesh, excluding target_iel/its neighbors. ────
    far_elements = Int[]
    step = max(1, mesh.nelem ÷ 5)
    for e = 1:step:mesh.nelem
        (e == target_iel || e in same_rank_neighbors) && continue
        push!(far_elements, e)
        length(far_elements) >= 5 && break
    end

    println("─"^70)
    println("Refined element: $target_iel")
    println("Same-rank neighbors: $same_rank_neighbors")
    println("Far elements (control group): $far_elements")
    println("─"^70)

    function report(iel)
        em = mesh.extra_mesh[iel]
        max_err_iel = 0.0
        worst_d = 0
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                computed = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                err = abs(computed-exact)
                if err > max_err_iel
                    max_err_iel = err; worst_d = d
                end
            end
        end
        println("iel=$iel ($(n_dirs[iel]) directions): max error = $max_err_iel (worst d=$worst_d)")
        return max_err_iel
    end

    println("REFINED ELEMENT:")
    report(target_iel)
    println("SAME-RANK NEIGHBORS:")
    for iel in same_rank_neighbors
        report(iel)
    end
    println("FAR ELEMENTS (control):")
    for iel in far_elements
        report(iel)
    end

    return x_dg, n_dirs, same_rank_neighbors, far_elements
end

function test_angular_ncf_upwind_unified_isolated(mesh, target_iel::Int, target_e_ext::Int)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    em_iel = mesh.extra_mesh[target_iel]
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    jel = 0; test_fid = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); test_fid = fid; break
    end
    em_jel = mesh.extra_mesh[jel]

    ngl = mesh.ngl
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly=mesh.ymax-mesh.ymin; Lz=mesh.zmax-mesh.zmin
    Pspatial(t) = sum(0.3^k*t^k for k=0:ngl-1)

    overall_max_err = 0.0
    n_tested_fine = 0
    n_tested_coarse = 0

    # ── Test BOTH directions: fine (iel=target_iel) depending on coarse
    # (jel), AND coarse (iel=jel) depending on fine (target_iel) — using
    # angular_deps directly, exactly as run_angular_ncf_sweep! does. ────────
    for (test_iel, test_jel, test_fid_local) in [(target_iel, jel, test_fid),
                                                    (jel, target_iel, Int(neighbor_fid[target_iel][test_fid]))]
        em_test_iel = mesh.extra_mesh[test_iel]
        em_test_jel = mesh.extra_mesh[test_jel]

        for d = 1:em_test_iel.extra_npoin
            key = (test_iel, d)
            haskey(angular_deps, key) || continue
            found = false
            local dep_nodes
            for (dj, dn) in angular_deps[key]
                dj == test_jel || continue
                found = true; dep_nodes = dn; break
            end
            found || continue
            length(dep_nodes) == 1 && continue   # skip singletons here, focus on GROUP cases

            src_e_ext = find_containing_angular_element(em_test_jel, Set(dep_nodes))
            src_e_ext == 0 && continue

            θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(
                em_test_jel.extra_coords, em_test_jel.extra_connijk, em_test_jel.extra_nop, src_e_ext)
            near_seam = isapprox(ϕmax, 2π; atol=1e-6)
            fix_ϕ(ϕ) = (near_seam && ϕ < 1e-8) ? 2π : ϕ
            nsrc = em_test_jel.extra_nop[src_e_ext] + 1
            Pang(t) = sum(0.3^k*t^k for k=0:nsrc-1)
            field(x,y,z,θ,ϕ) = Pspatial((x-x0)/Lx) + Pspatial((y-y0)/Ly) + Pspatial((z-z0)/Lz) +
                                Pang((θ-θmin)/(θmax-θmin)) + Pang((fix_ϕ(ϕ)-ϕmin)/(ϕmax-ϕmin))

            # Fabricate x_dg for EVERY node in dep_nodes (the source element)
            x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
            for ip in dep_nodes
                θip = em_test_jel.extra_coords[1,ip]; ϕip = em_test_jel.extra_coords[2,ip]
                vals = zeros(ngl^3)
                for k=1:ngl,j=1:ngl,i=1:ngl
                    gip = mesh.connijk[test_jel,i,j,k]
                    vals[lidx(i,j,k,ngl)] = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θip, ϕip)
                end
                x_dg[(test_jel,ip)] = vals
            end

            closure = build_angular_ncf_upwind_unified(test_iel, test_fid_local, test_jel, d,
                                                        src_e_ext, mesh, em_test_jel, x_dg,
                                                        face_correspondence, ngl)

            θd = em_test_iel.extra_coords[1,d]; ϕd = em_test_iel.extra_coords[2,d]
            irange,jrange,krange = face_ip_range(test_fid_local, ngl)
            for (i,j,k) in Iterators.product(irange,jrange,krange)
                computed = closure(i,j,k)
                gip = mesh.connijk[test_iel,i,j,k]
                exact = field(mesh.x[gip], mesh.y[gip], mesh.z[gip], θd, ϕd)
                err = abs(computed-exact)
                overall_max_err = max(overall_max_err, err)
            end
            test_iel == target_iel ? (n_tested_fine += 1) : (n_tested_coarse += 1)
        end
    end

    println("Tested $n_tested_fine fine-depends-on-coarse directions, $n_tested_coarse coarse-depends-on-fine")
    println("Overall max abs error: $overall_max_err")
    return overall_max_err
end

function test_angular_dependency_graph_coverage(mesh, neighbor_elem, angular_deps; tol::Float64=1e-6)
    n_checked = 0
    n_mismatch = 0
    mismatches = Tuple{Int,Int,Int}[]   # (iel, d, jel)

    for iel = 1:mesh.nelem
        em_iel = mesh.extra_mesh[iel]
        for fid = 1:6
            jel_raw = neighbor_elem[iel][fid]
            jel_raw == zero(TInt) && continue
            jel = Int(jel_raw)
            em_jel = mesh.extra_mesh[jel]

            for d = 1:em_iel.extra_npoin
                key = (iel, d)

                # Does angular_deps have an entry for THIS key, specifically
                # pointing at THIS jel?
                has_entry_for_this_jel = false
                if haskey(angular_deps, key)
                    for (dep_jel, _) in angular_deps[key]
                        dep_jel == jel && (has_entry_for_this_jel = true; break)
                    end
                end
                has_entry_for_this_jel && continue   # covered — skip

                # NO explicit entry — the sweep would fall through to
                # "ordinary same-index" here. Verify that's ACTUALLY correct:
                # does d even exist at jel, and does it have the SAME (θ,ϕ)?
                n_checked += 1
                if d > em_jel.extra_npoin
                    n_mismatch += 1
                    length(mismatches) < 20 && push!(mismatches, (iel,d,jel))
                    continue
                end
                θi = em_iel.extra_coords[1,d]; ϕi = em_iel.extra_coords[2,d]
                θj = em_jel.extra_coords[1,d]; ϕj = em_jel.extra_coords[2,d]
                if !(isapprox(θi,θj;atol=tol) && isapprox(ϕi,ϕj;atol=tol))
                    n_mismatch += 1
                    length(mismatches) < 20 && push!(mismatches, (iel,d,jel))
                end
            end
        end
    end

    println("Checked $n_checked (iel,d,jel) triples with NO angular_deps entry (fallback path)")
    println("Mismatches (fallback assumption WRONG): $n_mismatch  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))
    if n_mismatch > 0
        println("First mismatches (iel, d, jel):")
        for (iel,d,jel) in mismatches
            em_iel = mesh.extra_mesh[iel]
            println("  iel=$iel d=$d jel=$jel: iel's (θ,ϕ)=($(em_iel.extra_coords[1,d]),$(em_iel.extra_coords[2,d]))" *
                    (d <= mesh.extra_mesh[jel].extra_npoin ?
                     ", jel's (θ,ϕ)=($(mesh.extra_mesh[jel].extra_coords[1,d]),$(mesh.extra_mesh[jel].extra_coords[2,d]))" :
                     ", jel has no direction $d at all"))
        end
    end
    return n_mismatch
end

function test_angular_dependency_graph_coverage_full(mesh, metrics; target_iel::Int=1,
                                                       target_e_ext::Int=4, tol::Float64=1e-6)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    n_checked = 0
    n_mismatch = 0
    mismatches = Tuple{Int,Int,Int}[]

    for iel = 1:mesh.nelem
        em_iel = mesh.extra_mesh[iel]
        for fid = 1:6
            jel_raw = neighbor_elem[iel][fid]
            jel_raw == zero(TInt) && continue
            jel = Int(jel_raw)
            em_jel = mesh.extra_mesh[jel]

            for d = 1:em_iel.extra_npoin
                key = (iel, d)

                has_entry_for_this_jel = false
                if haskey(angular_deps, key)
                    for (dep_jel, _) in angular_deps[key]
                        dep_jel == jel && (has_entry_for_this_jel = true; break)
                    end
                end
                has_entry_for_this_jel && continue

                n_checked += 1
                if d > em_jel.extra_npoin
                    n_mismatch += 1
                    length(mismatches) < 20 && push!(mismatches, (iel,d,jel))
                    continue
                end
                θi = em_iel.extra_coords[1,d]; ϕi = em_iel.extra_coords[2,d]
                θj = em_jel.extra_coords[1,d]; ϕj = em_jel.extra_coords[2,d]
                if !(isapprox(θi,θj;atol=tol) && isapprox(ϕi,ϕj;atol=tol))
                    n_mismatch += 1
                    length(mismatches) < 20 && push!(mismatches, (iel,d,jel))
                end
            end
        end
    end

    println("Checked $n_checked (iel,d,jel) triples with NO angular_deps entry (fallback path)")
    println("Mismatches (fallback assumption WRONG): $n_mismatch  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))
    if n_mismatch > 0
        println("First mismatches (iel, d, jel):")
        for (iel,d,jel) in mismatches
            em_iel = mesh.extra_mesh[iel]
            em_jel = mesh.extra_mesh[jel]
            if d <= em_jel.extra_npoin
                println("  iel=$iel d=$d jel=$jel: iel(θ,ϕ)=($(em_iel.extra_coords[1,d]),$(em_iel.extra_coords[2,d])) " *
                        "jel(θ,ϕ)=($(em_jel.extra_coords[1,d]),$(em_jel.extra_coords[2,d]))")
            else
                println("  iel=$iel d=$d jel=$jel: jel has no direction $d at all " *
                        "(jel only has $(em_jel.extra_npoin) directions)")
            end
        end
    end

    return n_mismatch, mismatches
end

function test_angular_ncf_two_element_trace(mesh, metrics, basis, ω_gll; target_iel=1, target_e_ext=4,
                                             κ_val=0.0001, I0=4.2)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    ψ, dψ = basis.ψ, basis.dψ
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    jel = 0; test_fid = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); test_fid = fid; break
    end
    println("Isolating: iel=$target_iel <-> jel=$jel via fid=$test_fid")

    # ── Restrict: zero out EVERY OTHER face on BOTH elements, forcing a
    # genuinely closed two-element system (every other face -> true
    # boundary, using our own already-verified BC mechanism). ───────────────
    neighbor_elem_restricted = deepcopy(neighbor_elem)
    for fid = 1:6
        fid == test_fid || (neighbor_elem_restricted[target_iel] =
            Base.setindex(neighbor_elem_restricted[target_iel], zero(TInt), fid))
    end
    jfid = Int(neighbor_fid[target_iel][test_fid])
    for fid = 1:6
        fid == jfid || (neighbor_elem_restricted[jel] =
            Base.setindex(neighbor_elem_restricted[jel], zero(TInt), fid))
    end

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
        build_connijk_dg(mesh, neighbor_elem_restricted, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem_restricted)

    x_dg, n_dirs = run_angular_ncf_sweep_restricted!(mesh, metrics, ψ, dψ, ω_gll,
                                                      neighbor_elem_restricted, neighbor_fid,
                                                      face_correspondence, angular_deps,
                                                      κ_val, 0.0, I0, [target_iel, jel])

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    ngl = mesh.ngl
    for iel in [target_iel, jel]
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                gip = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                computed = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                err = abs(computed-exact)
                if err > 1e-6
                    println("MISMATCH iel=$iel d=$d (i,j,k)=($i,$j,$k): computed=$computed exact=$exact err=$err")
                end
            end
        end
    end
    println("Two-element trace complete.")
    return x_dg, n_dirs
end

function run_angular_ncf_sweep_restricted!(mesh, metrics, ψ, dψ, ω,
                                            neighbor_elem, neighbor_fid, face_correspondence,
                                            angular_deps, κ_val::Float64, σ_val::Float64, I0::Float64,
                                            elements::Vector{Int}; max_rounds::Int = 500)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    Ω_of = Dict{Int,Vector{NTuple{3,Float64}}}()
    for iel in elements
        em = mesh.extra_mesh[iel]
        Ω_of[iel] = Vector{NTuple{3,Float64}}(undef, em.extra_npoin)
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ω_of[iel][d] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        end
    end

    solved = Dict(iel => falses(n_dirs[iel]) for iel in elements)
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = fill(σ_val, npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    total_pairs = sum(n_dirs[iel] for iel in elements)
    n_solved_total = 0
    round = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("STALLED at round $max_rounds")

        n_newly_solved = 0
        for iel in elements
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                Ωx = Ω_of[iel][d][1]; Ωy = Ω_of[iel][d][2]; Ωz = Ω_of[iel][d][3]

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                upwind_fns = Dict{Int, Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]
                    if jel_raw == zero(TInt)
                        Ωx_c,Ωy_c,Ωz_c = Ωx,Ωy,Ωz
                        upwind_fns[fid] = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            I0*exp(-κ_val*(Ωx_c*(mesh.x[gip]-r_bc[1]) +
                                           Ωy_c*(mesh.y[gip]-r_bc[2]) +
                                           Ωz_c*(mesh.z[gip]-r_bc[3])))
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]
                                    ready = false
                                else
                                    dj_c = dj; fid_c = fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved
                                    ready = false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c = d; fid_c = fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σz, intΦ_eff, mesh.connijk)
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns, fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
                println("  ROUND $round: solved (iel=$iel, d=$d), Ω=($Ωx,$Ωy,$Ωz)")
            end
        end
        n_solved_total += n_newly_solved
        n_newly_solved == 0 && n_solved_total < total_pairs && error("STALLED")
    end
    println("Converged in $round rounds")
    return x_dg, n_dirs
end

function test_angular_ncf_targeted_trace(mesh, metrics, basis, ω_gll; target_iel=1, target_e_ext=4,
                                          κ_val=0.0001, I0=4.2)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    ψ, dψ = basis.ψ, basis.dψ
    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    jel = 0; test_fid = 0
    for fid = 1:6
        cand = neighbor_elem[target_iel][fid]
        cand == zero(TInt) && continue
        jel = Int(cand); test_fid = fid; break
    end
    jfid = Int(neighbor_fid[target_iel][test_fid])
    println("iel=$target_iel <-> jel=$jel via fid=$test_fid/jfid=$jfid")

    # Restrict to the two-element system, exactly as before.
    neighbor_elem_restricted = deepcopy(neighbor_elem)
    for fid = 1:6
        fid == test_fid || (neighbor_elem_restricted[target_iel] =
            Base.setindex(neighbor_elem_restricted[target_iel], zero(TInt), fid))
    end
    for fid = 1:6
        fid == jfid || (neighbor_elem_restricted[jel] =
            Base.setindex(neighbor_elem_restricted[jel], zero(TInt), fid))
    end

    connijk_dg, npoin_dg, dg_to_cg, face_correspondence =
        build_connijk_dg(mesh, neighbor_elem_restricted, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem_restricted)

    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    exact_at(iel,d,i,j,k) = begin
        em = mesh.extra_mesh[iel]
        θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        gip = mesh.connijk[iel,i,j,k]
        I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
    end

    targets = [(target_iel,1), (target_iel,72), (jel,1), (jel,28)]

    println("─"^70)
    println("STEP 1: classify each target's dependency TYPE")
    for (iel,d) in targets
        em_iel = mesh.extra_mesh[iel]
        θ = em_iel.extra_coords[1,d]; ϕ = em_iel.extra_coords[2,d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)

        println("(iel=$iel, d=$d): Ω=($Ωx,$Ωy,$Ωz)")
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            other = neighbor_elem_restricted[iel][fid]
            if other == zero(TInt)
                println("  fid=$fid: INFLOW, target=BOUNDARY")
                continue
            end
            other = Int(other)
            key = (iel,d)
            found = false
            if haskey(angular_deps, key)
                for (dj,dn) in angular_deps[key]
                    dj == other || continue
                    found = true
                    kind = length(dn)==1 ? "SINGLETON" : "GROUP($(length(dn)) nodes)"
                    println("  fid=$fid: INFLOW, target=jel=$other, angular_deps entry: $kind, nodes=$dn")
                end
            end
            if !found
                println("  fid=$fid: INFLOW, target=jel=$other, NO angular_deps entry -> FALLBACK same-index (d=$d)")
            end
        end
    end

    println("─"^70)
    println("STEP 2: run the restricted sweep, tracing EACH target's solve + comparison")
    x_dg, n_dirs, solved_order = run_angular_ncf_sweep_traced!(mesh, metrics, ψ, dψ, ω_gll,
                                                                 neighbor_elem_restricted, neighbor_fid,
                                                                 face_correspondence, angular_deps,
                                                                 κ_val, 0.0, I0, [target_iel, jel], targets)

    println("─"^70)
    println("STEP 3: final error at EACH target, sampled over several (i,j,k)")
    for (iel,d) in targets
        vals = x_dg[(iel,d)]
        max_err = 0.0
        for k=1:ngl,j=1:ngl,i=1:ngl
            err = abs(vals[lidx(i,j,k,ngl)] - exact_at(iel,d,i,j,k))
            max_err = max(max_err, err)
        end
        println("(iel=$iel, d=$d): max abs error over element = $max_err")
    end

    return x_dg
end

function run_angular_ncf_sweep_traced!(mesh, metrics, ψ, dψ, ω,
                                        neighbor_elem, neighbor_fid, face_correspondence,
                                        angular_deps, κ_val::Float64, σ_val::Float64, I0::Float64,
                                        elements::Vector{Int}, targets::Vector{Tuple{Int,Int}};
                                        max_rounds::Int = 500)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    Ω_of = Dict{Int,Vector{NTuple{3,Float64}}}()
    for iel in elements
        em = mesh.extra_mesh[iel]
        Ω_of[iel] = Vector{NTuple{3,Float64}}(undef, em.extra_npoin)
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ω_of[iel][d] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        end
    end

    solved = Dict(iel => falses(n_dirs[iel]) for iel in elements)
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = fill(σ_val, npoin_cg); intΦ_eff = zeros(npoin_cg)

    total_pairs = sum(n_dirs[iel] for iel in elements)
    n_solved_total = 0
    round = 0
    solved_order = Tuple{Int,Int,Int}[]  # (iel,d,round)

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("STALLED at round $max_rounds")
        n_newly_solved = 0

        for iel in elements
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                Ωx = Ω_of[iel][d][1]; Ωy = Ω_of[iel][d][2]; Ωz = Ω_of[iel][d][3]
                is_trace = (iel,d) in targets

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                upwind_fns = Dict{Int, Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]
                    if jel_raw == zero(TInt)
                        Ωx_c,Ωy_c,Ωz_c = Ωx,Ωy,Ωz
                        upwind_fns[fid] = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            I0*exp(-κ_val*(Ωx_c*(mesh.x[gip]-r_bc[1]) +
                                           Ωy_c*(mesh.y[gip]-r_bc[2]) +
                                           Ωz_c*(mesh.z[gip]-r_bc[3])))
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]
                                    ready = false
                                else
                                    dj_c = dj; fid_c = fid
                                    if is_trace
                                        println("    [TRACE (iel=$iel,d=$d) fid=$fid] SINGLETON from (jel=$jel,dj=$dj_c), " *
                                                "sample x_dg value at (1,1,1)=$(x_dg[(jel,dj_c)][lidx(1,1,1,ngl)])")
                                    end
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved
                                    ready = false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    if is_trace
                                        println("    [TRACE (iel=$iel,d=$d) fid=$fid] GROUP from jel=$jel, " *
                                                "src_e_ext=$src_e_ext, dep_nodes=$dep_nodes")
                                    end
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c = d; fid_c = fid
                            if is_trace
                                println("    [TRACE (iel=$iel,d=$d) fid=$fid] FALLBACK same-index from (jel=$jel,d=$d_c)")
                            end
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σz, intΦ_eff, mesh.connijk)
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns, fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    if is_trace
                        irange, jrange, krange = face_ip_range(fid, ngl)
			sample_ijk = (first(irange), first(jrange), first(krange))
			println("    [TRACE (iel=$iel,d=$d) fid=$fid] is_inflow=$(is_inflow[fid]) " *
        		"is_active=$(is_active[fid]) sample I_upwind$sample_ijk=" *
        		"$(is_active[fid] ? I_upwind(sample_ijk...) : NaN)")
		    end
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
                push!(solved_order, (iel,d,round))
                is_trace && println("  >>> SOLVED (iel=$iel,d=$d) in round $round")
            end
        end
        n_solved_total += n_newly_solved
        n_newly_solved == 0 && n_solved_total < total_pairs && error("STALLED")
    end

    return x_dg, n_dirs, solved_order
end

function test_angular_ncf_check_dependency_sources(mesh, x_dg, κ_val, I0)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    ngl = mesh.ngl

    groups = [(2, [3, 12, 16, 4, 19, 21, 23, 25, 20, 22, 24, 26, 1, 9, 13, 2]),
              (1, [13, 28, 32, 36, 25, 29, 33, 37, 26, 30, 34, 38, 27, 31, 35, 39])]

    for (jel, dep_nodes) in groups
        em = mesh.extra_mesh[jel]
        println("─"^50)
        println("Checking jel=$jel's dependency-source directions:")
        for dn in dep_nodes
            haskey(x_dg, (jel,dn)) || (println("  d=$dn: NOT YET SOLVED (unexpected)"); continue)
            θ = em.extra_coords[1,dn]; ϕ = em.extra_coords[2,dn]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            vals = x_dg[(jel,dn)]
            max_err = 0.0
            for k=1:ngl,j=1:ngl,i=1:ngl
                gip = mesh.connijk[jel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                max_err = max(max_err, abs(vals[lidx(i,j,k,ngl)] - exact))
            end
            println("  d=$dn (Ω=($Ωx,$Ωy,$Ωz)): max error = $max_err")
        end
    end
end

function test_angular_ncf_inspect_group_inflow(mesh, metrics, x_dg, κ_val, I0)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    # ── Case: (iel=1, d=72), GROUP from jel=2, src_e_ext=2, fid=1 ───────────
    iel, d, jel, fid, src_e_ext = 1, 72, 2, 1, 2
    em_iel = mesh.extra_mesh[iel]
    em_jel = mesh.extra_mesh[jel]

    θd = em_iel.extra_coords[1,d]; ϕd = em_iel.extra_coords[2,d]
    Ωx = sin(θd)*cos(ϕd); Ωy = sin(θd)*sin(ϕd); Ωz = cos(θd)
    println("Testing (iel=$iel, d=$d), Ω=($Ωx,$Ωy,$Ωz), via fid=$fid, src_e_ext=$src_e_ext")

    closure = build_angular_ncf_upwind_unified_debug(iel, fid, jel, d, src_e_ext, mesh, em_jel,
                                                x_dg, face_correspondence, ngl)

    irange, jrange, krange = face_ip_range(fid, ngl)
    max_err = 0.0
    for (i,j,k) in Iterators.product(irange,jrange,krange)
        computed = closure(i,j,k)
        gip = mesh.connijk[iel,i,j,k]
        exact = I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
        err = abs(computed - exact)
        max_err = max(max_err, err)
        err > 1e-6 && println("  (i,j,k)=($i,$j,$k): computed=$computed exact=$exact err=$err")
    end
    println("  max inflow error over face: $max_err")
    return max_err
end

function build_angular_ncf_upwind_unified_debug(iel::Int, fid::Int, jel::Int, d::Int,
                                                 src_e_ext::Int, mesh, em_jel,
                                                 x_dg, face_correspondence, ngl::Int)
    em_iel = mesh.extra_mesh[iel]

    # ── Step 1: find (θd, ϕd) — the TARGET direction's own coordinates ──────
    θd = 0.0; ϕd = 0.0
    for e = 1:em_iel.extra_nelem
        n = em_iel.extra_nop[e] + 1
        for a = 1:n
            for b = 1:n
                if em_iel.extra_connijk[e,a,b] == d
                    θd = em_iel.extra_coords[1, d]
                    ϕd = em_iel.extra_coords[2, d]
                end
            end
        end
    end
    println("STEP 1: target direction d=$d has (θd,ϕd)=($θd,$ϕd)")

    # ── Step 2: source element's node grid, bounds, seam check ──────────────
    nsrc = em_jel.extra_nop[src_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(
        em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd
    println("STEP 2: src_e_ext=$src_e_ext bounds: θ∈[$θmin,$θmax] ϕ∈[$ϕmin,$ϕmax], near_seam=$near_seam, ϕd_fix=$ϕd_fix")

    θ_src = Vector{Float64}(undef, nsrc)
    for a = 1:nsrc
        θ_src[a] = em_jel.extra_coords[1, em_jel.extra_connijk[src_e_ext,a,1]]
    end
    ϕ_src = Vector{Float64}(undef, nsrc)
    for b = 1:nsrc
        raw = em_jel.extra_coords[2, em_jel.extra_connijk[src_e_ext,1,b]]
        ϕ_src[b] = (near_seam && raw < 1e-8) ? 2π : raw
    end
    println("STEP 2: θ_src=$θ_src")
    println("STEP 2: ϕ_src=$ϕ_src")

    # ── Step 3: 1D Lagrange weights ──────────────────────────────────────────
    Lθ = build_1d_lagrange_matrix(θ_src, [θd])
    Lϕ = build_1d_lagrange_matrix(ϕ_src, [ϕd_fix])
    println("STEP 3: Lθ=$Lθ")
    println("STEP 3: Lϕ=$Lϕ")

    src_ips = Matrix{Int}(undef, nsrc, nsrc)
    for a = 1:nsrc
        for b = 1:nsrc
            src_ips[a,b] = em_jel.extra_connijk[src_e_ext,a,b]
        end
    end
    println("STEP 3: src_ips grid=$src_ips")

    return (i,j,k) -> begin
        corr = face_correspondence[(iel,fid)][(i,j,k)]
        ii2 = corr[1]; jj2 = corr[2]; kk2 = corr[3]

        val = 0.0
        for a = 1:nsrc
            for b = 1:nsrc
                ip = src_ips[a,b]
                raw_val = x_dg[(jel,ip)][lidx(ii2,jj2,kk2,ngl)]
                contribution = Lθ[1,a]*Lϕ[1,b]*raw_val
                println("    a=$a b=$b ip=$ip Lθ[1,$a]=$(Lθ[1,a]) Lϕ[1,$b]=$(Lϕ[1,b]) " *
                        "x_dg[($jel,$ip)][($ii2,$jj2,$kk2)]=$raw_val contribution=$contribution")
                val += contribution
            end
        end
        val
    end
end

function debug_angular_scheduler_stall(mesh, metrics, neighbor_elem, angular_deps; max_rounds::Int=200)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    Ω_of = Vector{Vector{NTuple{3,Float64}}}(undef, mesh.nelem)
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        Ω_of[iel] = Vector{NTuple{3,Float64}}(undef, em.extra_npoin)
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ω_of[iel][d] = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
        end
    end

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    total_pairs = sum(n_dirs)
    n_solved_total = 0
    round = 0

    while n_solved_total < total_pairs && round < max_rounds
        round += 1
        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                Ωx,Ωy,Ωz = Ω_of[iel][d]
                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel = neighbor_elem[iel][fid]
                    jel == zero(TInt) && continue
                    jel = Int(jel)
                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            for dn in dep_nodes
                                solved[jel][dn] || (ready = false)
                            end
                        end
                    end
                    if !handled
                        (d <= n_dirs[jel] && !solved[jel][d]) || continue
                        ready = false
                    end
                    ready || break
                end
                if ready
                    solved[iel][d] = true
                    n_newly_solved += 1
                end
            end
        end
        n_solved_total += n_newly_solved
        n_newly_solved == 0 && break
    end

    println("Stopped after round $round, $n_solved_total/$total_pairs solved")
    println("─"^60)
    println("STUCK pairs and their blocking dependencies:")
    n_shown = 0
    for iel = 1:mesh.nelem
        for d = 1:n_dirs[iel]
            solved[iel][d] && continue
            n_shown >= 15 && break
            Ωx,Ωy,Ωz = Ω_of[iel][d]
            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)
            println("STUCK: (iel=$iel, d=$d), Ω=($Ωx,$Ωy,$Ωz)")
            for fid = 1:6
                (is_inflow[fid] && is_active[fid]) || continue
                jel = neighbor_elem[iel][fid]
                jel == zero(TInt) && continue
                jel = Int(jel)
                key = (iel,d)
                if haskey(angular_deps, key)
                    for (dep_jel, dep_nodes) in angular_deps[key]
                        dep_jel == jel || continue
                        unsolved = [dn for dn in dep_nodes if !solved[jel][dn]]
                        println("  fid=$fid -> jel=$jel: angular_deps entry, UNSOLVED deps=$unsolved (of $(length(dep_nodes)) total)")
                    end
                else
                    println("  fid=$fid -> jel=$jel: FALLBACK same-index, jel has direction $d? " *
                            "$(d<=n_dirs[jel]), solved? $(d<=n_dirs[jel] && solved[jel][d])")
                end
            end
            n_shown += 1
        end
    end
end

function trace_single_stuck_node(mesh, metrics, neighbor_elem, angular_deps, target_iel::Int, target_d::Int)
    ngl = mesh.ngl
    em = mesh.extra_mesh[target_iel]
    θ = em.extra_coords[1,target_d]; ϕ = em.extra_coords[2,target_d]
    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
    println("Tracing (iel=$target_iel, d=$target_d), Ω=($Ωx,$Ωy,$Ωz)")

    is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
    classify_element_faces!(is_inflow, is_active, target_iel, Ωx, Ωy, Ωz, ngl,
                             metrics.dξdx,metrics.dξdy,metrics.dξdz,
                             metrics.dηdx,metrics.dηdy,metrics.dηdz,
                             metrics.dζdx,metrics.dζdy,metrics.dζdz)

    for fid = 1:6
        (is_inflow[fid] && is_active[fid]) || continue
        jel = neighbor_elem[target_iel][fid]
        if jel == zero(TInt)
            println("  fid=$fid: TRUE BOUNDARY — always ready")
            continue
        end
        jel = Int(jel)
        key = (target_iel, target_d)
        if haskey(angular_deps, key)
            found = false
            for (dep_jel, dep_nodes) in angular_deps[key]
                dep_jel == jel || continue
                found = true
                println("  fid=$fid -> jel=$jel: angular_deps entry, depends on $(length(dep_nodes)) nodes: $dep_nodes")
            end
            found || println("  fid=$fid -> jel=$jel: has angular_deps key but NO entry for THIS jel — falls to FALLBACK")
        else
            println("  fid=$fid -> jel=$jel: NO angular_deps entry at all — FALLBACK same-index (d=$target_d)")
            if target_d > mesh.extra_mesh[jel].extra_npoin
                println("    ⚠ jel only has $(mesh.extra_mesh[jel].extra_npoin) directions — d=$target_d DOES NOT EXIST THERE")
            end
        end
    end
end

function trace_multiple_stuck_nodes(mesh, metrics, neighbor_elem, angular_deps, targets::Vector{Tuple{Int,Int}})
    ngl = mesh.ngl
    for (target_iel, target_d) in targets
        em = mesh.extra_mesh[target_iel]
        θ = em.extra_coords[1,target_d]; ϕ = em.extra_coords[2,target_d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        println("Tracing (iel=$target_iel, d=$target_d), Ω=($Ωx,$Ωy,$Ωz)")

        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, target_iel, Ωx, Ωy, Ωz, ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            jel = neighbor_elem[target_iel][fid]
            if jel == zero(TInt)
                println("  fid=$fid: TRUE BOUNDARY")
                continue
            end
            jel = Int(jel)
            key = (target_iel, target_d)
            if haskey(angular_deps, key)
                found = false
                for (dep_jel, dep_nodes) in angular_deps[key]
                    dep_jel == jel || continue
                    found = true
                    println("  fid=$fid -> jel=$jel: depends on $(length(dep_nodes)) nodes: $dep_nodes")
                end
                found || println("  fid=$fid -> jel=$jel: no entry for this jel, FALLBACK (d=$target_d)")
            else
                println("  fid=$fid -> jel=$jel: NO entry, FALLBACK (d=$target_d)")
            end
        end
        println()
    end
end

function test_dependency_group_single_element(mesh, iel::Int, node_list::Vector{Int})
    em = mesh.extra_mesh[iel]

    # For each angular element at iel, build its node set, then check
    # whether node_list is a SUBSET of (or exactly equal to) any single
    # element's node set.
    matches = Int[]
    for e = 1:em.extra_nelem
        n = em.extra_nop[e] + 1
        elem_nodes = Set{Int}()
        for a = 1:n
            for b = 1:n
                push!(elem_nodes, em.extra_connijk[e,a,b])
            end
        end
        if Set(node_list) == elem_nodes
            push!(matches, e)
        end
    end

    if length(matches) == 1
        println("✓ All $(length(node_list)) nodes belong EXACTLY to single angular element e_ext=$(matches[1])")
    elseif length(matches) > 1
        println("⚠ node_list exactly matches MULTIPLE elements: $matches (unexpected)")
    else
        # No exact match — check whether the nodes span more than one element
        contributing_elements = Set{Int}()
        for e = 1:em.extra_nelem
            n = em.extra_nop[e] + 1
            elem_nodes = Set{Int}()
            for a = 1:n
                for b = 1:n
                    push!(elem_nodes, em.extra_connijk[e,a,b])
                end
            end
            !isempty(intersect(Set(node_list), elem_nodes)) && push!(contributing_elements, e)
        end
        println("✗ No single element exactly matches node_list. " *
                "Nodes are drawn from $(length(contributing_elements)) different elements: " *
                "$(sort(collect(contributing_elements)))")
    end

    return matches
end

using Graphs

"""
At a scheduler stall (n_newly_solved==0, pairs remain), find every
strongly-connected component of size >1 among the CURRENTLY UNSOLVED
(iel,d) pairs, restricted to dependency edges pointing at OTHER unsolved
pairs (already-solved dependencies are irrelevant to the cycle itself).
"""
function find_cyclic_sccs(mesh, metrics, neighbor_elem, angular_deps, solved, n_dirs)
    unsolved_pairs = Tuple{Int,Int}[]
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        solved[iel][d] || push!(unsolved_pairs, (iel,d))
    end
    pair_to_idx = Dict(p => i for (i,p) in enumerate(unsolved_pairs))

    g = SimpleDiGraph(length(unsolved_pairs))
    for (i,(iel,d)) in enumerate(unsolved_pairs)
        em = mesh.extra_mesh[iel]
        θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
        classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, mesh.ngl,
                                 metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                 metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                 metrics.dζdx,metrics.dζdy,metrics.dζdz)
        for fid = 1:6
            (is_inflow[fid] && is_active[fid]) || continue
            jel = neighbor_elem[iel][fid]
            jel == zero(TInt) && continue
            jel = Int(jel)
            key = (iel,d)
            if haskey(angular_deps, key)
                for (dep_jel, dep_nodes) in angular_deps[key]
                    dep_jel == jel || continue
                    for dn in dep_nodes
                        haskey(pair_to_idx, (jel,dn)) || continue
                        add_edge!(g, i, pair_to_idx[(jel,dn)])
                    end
                end
            elseif d <= n_dirs[jel] && haskey(pair_to_idx, (jel,d))
                add_edge!(g, i, pair_to_idx[(jel,d)])
            end
        end
    end

    sccs = strongly_connected_components(g)
    cyclic_sccs = [ [unsolved_pairs[i] for i in scc] for scc in sccs if length(scc) > 1 ]
    return cyclic_sccs
end

"""
Locally resolve ONE cyclic SCC via Gauss-Seidel iteration: every member
solved using CURRENT guesses for other SCC members' values (treating
already-solved, non-cyclic neighbors as exact fixed data), updated each
pass, until converged.
"""
#=function resolve_angular_cycle!(mesh, metrics, ψ, dψ, ω, neighbor_elem, neighbor_fid,
                                 face_correspondence, angular_deps, x_dg, solved, n_dirs,
                                 scc::Vector{Tuple{Int,Int}}, κ_val, σ_val, I0;
                                 max_iter::Int=100, tol::Float64=1e-9)
    ngl = mesh.ngl
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = fill(σ_val, npoin_cg); intΦ_eff = zeros(npoin_cg)
    scc_set = Set(scc)

    # Initialize every SCC member's storage to zero (starting guess).
    for (iel,d) in scc
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    for iter = 1:max_iter
        max_delta = 0.0
        for (iel,d) in scc
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σz, intΦ_eff, mesh.connijk)

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
                I_upwind = (i,j,k) -> 0.0

                if is_inflow[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        I_upwind = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            I0*exp(-κ_val*(Ωx*(mesh.x[gip]-r_bc[1])+Ωy*(mesh.y[gip]-r_bc[2])+Ωz*(mesh.z[gip]-r_bc[3])))
                        end
                    else
                        jel = Int(jel)
                        key = (iel,d)
                        handled = false
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                handled = true
                                if length(dep_nodes) == 1
                                    dj = dep_nodes[1]; fid_c = fid
                                    if haskey(x_dg, (jel,dj))
                                        I_upwind = (i,j,k) -> begin
                                            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                            x_dg[(jel,dj)][lidx(corr[1],corr[2],corr[3],ngl)]
                                        end
                                    end
                                else
                                    if all(haskey(x_dg,(jel,dn)) for dn in dep_nodes)
                                        em_jel = mesh.extra_mesh[jel]
                                        src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                        I_upwind = build_angular_ncf_upwind_unified(
                                            iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    end
                                end
                            end
                        end
                        if !handled && haskey(x_dg, (jel,d))
                            d_c = d; fid_c = fid
                            I_upwind = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta = max(max_delta, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end
        max_delta < tol && (println("  cycle of size $(length(scc)) converged in $iter iterations"); break)
        iter == max_iter && @warn "cycle of size $(length(scc)) did not converge (max_delta=$max_delta)"
    end

    for (iel,d) in scc
        solved[iel][d] = true
    end
end=#

function resolve_angular_cycle!(mesh, metrics, ψ, dψ, ω, neighbor_elem, neighbor_fid,
                                 face_correspondence, angular_deps, x_dg, solved, n_dirs,
                                 scc::Vector{Tuple{Int,Int}}, κ, σ, intΦ_eff, bc_fn_all;
                                 extra_source_all = nothing, max_iter::Int=100, tol::Float64=1e-9)
    ngl = mesh.ngl

    for (iel,d) in scc
        x_dg[(iel,d)] = zeros(ngl^3)
    end

    for iter = 1:max_iter
        max_delta = 0.0
        for (iel,d) in scc
            em = mesh.extra_mesh[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

            is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
            classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz)

            A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
            build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                     metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                     metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                     metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                     metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
            if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
            end

            for fid = 1:6
                is_active[fid] || continue
                face_normal = compute_face_normal_avg(iel, fid, ngl,
                    metrics.dξdx,metrics.dξdy,metrics.dξdz,
                    metrics.dηdx,metrics.dηdy,metrics.dηdz,
                    metrics.dζdx,metrics.dζdy,metrics.dζdz)
                I_upwind = (i,j,k) -> 0.0

                if is_inflow[fid]
                    jel = neighbor_elem[iel][fid]
                    if jel == zero(TInt)
                        bc = bc_fn_all[(iel,d)]
                        I_upwind = (i,j,k) -> begin
                            gip = mesh.connijk[iel,i,j,k]
                            bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                        end
                    else
                        jel = Int(jel)
                        key = (iel,d)
                        handled = false
                        if haskey(angular_deps, key)
                            for (dep_jel, dep_nodes) in angular_deps[key]
                                dep_jel == jel || continue
                                handled = true
                                if length(dep_nodes) == 1
                                    dj = dep_nodes[1]; fid_c = fid
                                    if haskey(x_dg, (jel,dj))
                                        I_upwind = (i,j,k) -> begin
                                            corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                            x_dg[(jel,dj)][lidx(corr[1],corr[2],corr[3],ngl)]
                                        end
                                    end
                                else
                                    if all(haskey(x_dg,(jel,dn)) for dn in dep_nodes)
                                        em_jel = mesh.extra_mesh[jel]
                                        src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                        I_upwind = build_angular_ncf_upwind_unified(
                                            iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                    end
                                end
                            end
                        end
                        if !handled && haskey(x_dg, (jel,d))
                            d_c = d; fid_c = fid
                            I_upwind = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                end
                build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                      is_inflow[fid], is_active[fid], face_normal, ω,
                                      metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                      metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                      metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                      metrics.Je, iel, I_upwind)
            end

            x_new = A_vol \ RHS_local
            max_delta = max(max_delta, maximum(abs.(x_new .- x_dg[(iel,d)])))
            x_dg[(iel,d)] = x_new
        end
        max_delta < tol && break
    end

    for (iel,d) in scc
        solved[iel][d] = true
    end
end

function test_find_cyclic_sccs(mesh, metrics, target_iel::Int, target_e_ext::Int)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    # ── Run the round-based scheduler (placeholder-only, no physics) UNTIL
    # it stalls — same shape as test_angular_scheduling — so `solved`
    # reflects the ACTUAL stall state we traced by hand earlier. ────────────
    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    round = 0
    while true
        round += 1
        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)
                ready = true
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel = neighbor_elem[iel][fid]
                    jel == zero(TInt) && continue
                    jel = Int(jel)
                    key = (iel,d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            for dn in dep_nodes
                                solved[jel][dn] || (ready = false)
                            end
                        end
                    end
                    if !handled
                        (d <= n_dirs[jel] && !solved[jel][d]) && (ready = false)
                    end
                    ready || break
                end
                if ready
                    solved[iel][d] = true
                    n_newly_solved += 1
                end
            end
        end
        n_newly_solved == 0 && break
        round > 100 && break
    end

    n_unsolved = sum(count(!,solved[iel]) for iel = 1:mesh.nelem)
    println("Scheduler stalled after $round rounds, $n_unsolved pairs unsolved")

    # ── Now run cycle detection on the stalled state ─────────────────────────
    cyclic_sccs = find_cyclic_sccs(mesh, metrics, neighbor_elem, angular_deps, solved, n_dirs)

    println("Found $(length(cyclic_sccs)) cyclic SCC(s):")
    for (i, scc) in enumerate(cyclic_sccs)
        println("  SCC $i (size $(length(scc))): $scc")
    end

    # ── Cross-check against the specific cycle we already traced by hand:
    # (1,73) and (11,45) should be in the SAME SCC. ──────────────────────────
    target_pair_a = (1, 73)
    target_pair_b = (11, 45)
    found_together = false
    for scc in cyclic_sccs
        if target_pair_a in scc && target_pair_b in scc
            found_together = true
        end
    end
    println("(1,73) and (11,45) in the SAME cyclic SCC? $found_together  " *
            (found_together ? "✓ MATCHES hand-traced cycle" : "✗ does NOT match — investigate"))

    return cyclic_sccs, solved
end

function test_manufactured_isotropic_angular_ncf(mesh, metrics, basis, ω_gll;
                                                  target_iel=1, target_e_ext=4,
                                                  I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                  κ_val=0.001, σ_val=0.0005,
                                                  max_iter=500, tol=1e-9)
    adapt_single_angular_element_conforming!(mesh, target_iel, target_e_ext)
    verify_angular_phi_seam_conformity(mesh, target_iel)

    ψ, dψ = basis.ψ, basis.dψ
    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    # ── Group elements by angular-mesh signature — the actual fix. ──────────
    groups, representative_of = group_elements_by_angular_mesh(mesh)
    println("Grouped $(mesh.nelem) elements into $(length(groups)) distinct angular-mesh signatures")

    # ── Angular quadrature weights, C, Ω_m — computed ONCE PER GROUP, not
    # once per element. ───────────────────────────────────────────────────────
    w_ext_all = Dict{Int, Vector{Float64}}()
    W_all = Dict{Int, Float64}()
    C_all = Dict{Int, Matrix{Float64}}()
    Ωm_all = Dict{Int, Vector{NTuple{3,Float64}}}()

    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        W = sum(w_ext)
        w_ext_all[rep] = w_ext
        W_all[rep] = W
        C_all[rep] = [σ_val*w_ext[m]/W for kk=1:em.extra_npoin, m=1:em.extra_npoin]
        Ωm_all[rep] = [(sin(em.extra_coords[1,m])*cos(em.extra_coords[2,m]),
                        sin(em.extra_coords[1,m])*sin(em.extra_coords[2,m]),
                        cos(em.extra_coords[1,m])) for m = 1:em.extra_npoin]
    end

    # ── Manufactured source (q_ext) and BC — still per-ELEMENT (spatial
    # positions genuinely differ), but reusing the GROUP's angular data. ────
    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()

    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]
        C_iel = C_all[rep]
        Ω_m = Ωm_all[rep]

        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)

            q_vec = zeros(npoin_cg)
            for k=1:ngl, j=1:ngl, i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = sum(C_iel[dk,m]*I_manuf(x_,y_,z_, Ω_m[m]...) for m = 1:em.extra_npoin)
                q_vec[g] = Ω_dot_grad + κ_val*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
        end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        I_cg_all[(iel,d)] = zeros(ngl^3)
    end

    coeff_all = Dict{Int, Vector{Float64}}()
    for rep in keys(groups)
    	w_ext = w_ext_all[rep]; W = W_all[rep]
    	coeff_all[rep] = σ_val .* w_ext ./ W
    end

    for it = 1:max_iter
	extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
	for iel = 1:mesh.nelem
    	    em = mesh.extra_mesh[iel]
    	    rep = representative_of[iel]
    	    coeff = coeff_all[rep]

    	    # Precompute ONCE per element: the ngl^3 global ids in the SAME order
    	    # as lidx(i,j,k,ngl), so LOCAL (I_cg_all) and a VIEW into q_k at these
    	    # positions line up index-for-index.
    	    local_gids = Vector{Int}(undef, ngl^3)
    	    for k=1:ngl, j=1:ngl, i=1:ngl
        	local_gids[lidx(i,j,k,ngl)] = mesh.connijk[iel,i,j,k]
    	    end

    	    for dk = 1:em.extra_npoin
        	q_k = copy(q_ext_all[(iel,dk)])
        	q_k_local = @view q_k[local_gids]   # LOCAL-sized VIEW, same ordering as I_cg_all

        	for m = 1:em.extra_npoin
            	    coeff[m] == 0.0 && continue
            	    q_k_local .+= coeff[m] .* I_cg_all[(iel,m)]
        	end
        	extra_source_all[(iel,dk)] = q_k
    	    end
	end

        x_dg, _ = run_angular_ncf_sweep!(mesh, metrics, ψ, dψ, ω_gll,
                                         neighbor_elem, neighbor_fid, face_correspondence,
                                         angular_deps, κ, σz, intΦ_eff, bc_fn_all;
                                         extra_source_all=extra_source_all)

        max_delta = 0.0
        for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
            max_delta = max(max_delta, maximum(abs.(x_dg[(iel,d)] .- I_cg_all[(iel,d)])))
            I_cg_all[(iel,d)] = copy(x_dg[(iel,d)])
        end
        println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (println("Converged after $it iterations"); break)
    end

    max_err = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
                computed = I_cg_all[(iel,d)][lidx(i,j,k,ngl)]
                max_err = max(max_err, abs(computed - I_true))
            end
        end
    end
    println("Manufactured (isotropic, angular NCF) — max abs error: $max_err")
    return max_err
end
 
function adapt_single_angular_element_with_metrics!(mesh, target_iel::Int, target_e_ext::Int)
    iel = target_iel
    e_ext = target_e_ext

    em = mesh.extra_mesh[iel]
    ψ  = em.ψ
    dψ = em.dψ
    nop_ang     = em.extra_nop
    nelem_ang   = em.extra_nelem
    connijk_ang = em.extra_connijk
    coords_ang  = em.extra_coords
    Je_ang      = em.extra_metrics.Je
    dξdx_ang    = em.extra_metrics.dξdx
    dxdξ_ang    = em.extra_metrics.dxdξ
    dξdy_ang    = em.extra_metrics.dξdy
    dydξ_ang    = em.extra_metrics.dydξ
    dηdy_ang    = em.extra_metrics.dηdy
    dydη_ang    = em.extra_metrics.dydη
    dηdx_ang    = em.extra_metrics.dηdx
    dxdη_ang    = em.extra_metrics.dxdη

    θmax = coords_ang[1, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    θmin = coords_ang[1, connijk_ang[e_ext, 1, 1]]
    ϕmax = coords_ang[2, connijk_ang[e_ext, nop_ang[e_ext]+1, nop_ang[e_ext]+1]]
    ϕmin = coords_ang[2, connijk_ang[e_ext, 1, 1]]
    ϕmax == 0 && (ϕmax = 2π)
    θ12 = (θmax + θmin) / 2
    ϕ12 = (ϕmax + ϕmin) / 2

    nelem_ang_new = nelem_ang + 3
    npoin_ang_new_upper = em.extra_npoin + (4*(nop_ang[e_ext]+1)^2 - 4*(nop_ang[e_ext]+1))

    connijk_ang_new = zeros(Int, nelem_ang_new, nop_ang[1]+1, nop_ang[1]+1)
    coords_new      = zeros(Float64, 2, npoin_ang_new_upper)
    metrics_new     = allocate_metrics(NSD_2D(), nelem_ang_new, 0, nop_ang[e_ext], TFloat, CPU())
    nop_ang_new     = zeros(Int, nelem_ang_new)
    nop_ang_new[1:nelem_ang_new-1] .= nop_ang[e_ext]
    nop_ang_new[nelem_ang_new]      = nop_ang[1]

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

    lgl = basis_structs_ξ_ω!(LGL(), nop_ang[1], CPU())

    # ── Copy elements before e_ext, INCLUDING their metrics ──────────────────
    if e_ext > 1
        for e_ext1 = 1:e_ext-1
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                metrics_new.dxdξ[e_ext1,i,j] = dxdξ_ang[e_ext1,i,j]
                metrics_new.Je[e_ext1,i,j]   = Je_ang[e_ext1,i,j]
                metrics_new.dξdx[e_ext1,i,j] = dξdx_ang[e_ext1,i,j]
                metrics_new.dxdη[e_ext1,i,j] = dxdη_ang[e_ext1,i,j]
                metrics_new.dηdx[e_ext1,i,j] = dηdx_ang[e_ext1,i,j]
                metrics_new.dηdy[e_ext1,i,j] = dηdy_ang[e_ext1,i,j]
                metrics_new.dξdy[e_ext1,i,j] = dξdy_ang[e_ext1,i,j]
                metrics_new.dydη[e_ext1,i,j] = dydη_ang[e_ext1,i,j]
                metrics_new.dydξ[e_ext1,i,j] = dydξ_ang[e_ext1,i,j]
                insert_point!(θ, ϕ, connijk_ang_new, e_ext1, i, j)
            end
        end
    end

    child_ranges = [
        (θmin, θ12, ϕmin, ϕ12, e_ext),
        (θ12, θmax, ϕmin, ϕ12, e_ext+1),
        (θmin, θ12, ϕ12, ϕmax, e_ext+2),
        (θ12, θmax, ϕ12, ϕmax, e_ext+3),
    ]
    for (θlo, θhi, ϕlo, ϕhi, e_child) in child_ranges
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

    # ── Copy elements after e_ext, INCLUDING their metrics ───────────────────
    if e_ext < nelem_ang
        for e_ext1 = e_ext+1:nelem_ang
            dst = e_ext1 + 3
            for i = 1:nop_ang[e_ext1]+1, j = 1:nop_ang[e_ext1]+1
                θ = coords_ang[1, connijk_ang[e_ext1,i,j]]
                ϕ = coords_ang[2, connijk_ang[e_ext1,i,j]]
                metrics_new.dxdξ[dst,i,j] = dxdξ_ang[e_ext1,i,j]
                metrics_new.Je[dst,i,j]   = Je_ang[e_ext1,i,j]
                metrics_new.dξdx[dst,i,j] = dξdx_ang[e_ext1,i,j]
                metrics_new.dxdη[dst,i,j] = dxdη_ang[e_ext1,i,j]
                metrics_new.dηdx[dst,i,j] = dηdx_ang[e_ext1,i,j]
                metrics_new.dηdy[dst,i,j] = dηdy_ang[e_ext1,i,j]
                metrics_new.dξdy[dst,i,j] = dξdy_ang[e_ext1,i,j]
                metrics_new.dydη[dst,i,j] = dydη_ang[e_ext1,i,j]
                metrics_new.dydξ[dst,i,j] = dydξ_ang[e_ext1,i,j]
                insert_point!(θ, ϕ, connijk_ang_new, dst, i, j)
            end
        end
    end

    npoin_ang_new = iter - 1
    coords_new = coords_new[:, 1:npoin_ang_new]

    # ── Build metrics for the 4 NEW child elements (EXACT logic from
    # adapt_angular_grid_3Dby2D!, restored — this is the piece that was
    # missing). ────────────────────────────────────────────────────────────
    for e_ext1 = e_ext:e_ext+3
        for i = 1:nop_ang_new[e_ext1]+1, j = 1:nop_ang_new[e_ext1]+1
            ip  = connijk_ang_new[e_ext1,i,j]
            θij = coords_new[1,ip]; ϕij = coords_new[2,ip]
            for l = 1:nop_ang_new[e_ext1]+1
                for k = 1:nop_ang_new[e_ext1]+1
                    a = dψ[i,k]*ψ[j,l]; b = ψ[i,k]*dψ[j,l]
                    metrics_new.dxdξ[e_ext1,k,l] += a*θij
                    metrics_new.dxdη[e_ext1,k,l] += b*θij
                    metrics_new.dydξ[e_ext1,k,l] += a*ϕij
                    metrics_new.dydη[e_ext1,k,l] += b*ϕij
                end
            end
        end
        for l = 1:nop_ang_new[e_ext1]+1, k = 1:nop_ang_new[e_ext1]+1
            dxdξ_v = metrics_new.dxdξ[e_ext1,k,l]; dydη_v = metrics_new.dydη[e_ext1,k,l]
            dydξ_v = metrics_new.dydξ[e_ext1,k,l]; dxdη_v = metrics_new.dxdη[e_ext1,k,l]
            ip  = connijk_ang_new[e_ext1,k,l]
            θ = coords_new[1,ip]
            Je_v   = dxdξ_v*dydη_v - dydξ_v*dxdη_v
            Jinv   = 1.0 / Je_v
            metrics_new.Je[e_ext1,k,l]   = sin(θ) * Je_v
            metrics_new.dξdx[e_ext1,k,l] =  dydη_v*Jinv
            metrics_new.dξdy[e_ext1,k,l] = -dxdη_v*Jinv
            metrics_new.dηdx[e_ext1,k,l] = -dydξ_v*Jinv
            metrics_new.dηdy[e_ext1,k,l] =  dxdξ_v*Jinv
        end
    end

    # ── ϕ-periodicity merge (unchanged from before, but now ALSO needs to
    # keep metrics arrays in sync — metrics are PER-ELEMENT, indexed by
    # (e,i,j), not by global point id, so NO remapping needed here at all;
    # only coords_new/connijk_ang_new (point-indexed) were ever touched by
    # the merge). ─────────────────────────────────────────────────────────
    zero_ϕ_dict = Dict{Float64, Int}()
    for iper = 1:npoin_ang_new
        if coords_new[2, iper] < 1e-10
            zero_ϕ_dict[coords_new[1, iper]] = iper
        end
    end
    merge_pairs = Dict{Int, Int}()
    for iper = 1:npoin_ang_new
        ϕ = coords_new[2, iper]
        if abs(ϕ / π - 2.0) < 1e-10
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
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] == ip_old && (connijk_ang_new[e,i,j] = ip_new)
                end
            end
            for e = 1:nelem_ang_new
                for i = 1:nop_ang_new[e]+1, j = 1:nop_ang_new[e]+1
                    connijk_ang_new[e,i,j] > ip_old && (connijk_ang_new[e,i,j] -= 1)
                end
            end
            merge_pairs = Dict(
                (k > ip_old ? k - 1 : k) => (v > ip_old ? v - 1 : v)
                    for (k, v) in merge_pairs if k != ip_old
                )
            for i = ip_old+1:npoin_ang_new
                coords_new[1,i-1] = coords_new[1,i]
                coords_new[2,i-1] = coords_new[2,i]
            end
            npoin_ang_new -= 1
        end
    end
    coords_new = coords_new[:, 1:npoin_ang_new]

    em.extra_connijk = connijk_ang_new
    em.extra_coords   = coords_new
    em.extra_nop      = nop_ang_new
    em.extra_nelem    = nelem_ang_new
    em.extra_npoin    = npoin_ang_new
    em.extra_metrics  = metrics_new

    return nothing
end

function adapt_single_angular_element_conforming!(mesh, target_iel::Int, target_e_ext::Int)
    em = mesh.extra_mesh[target_iel]
    partner_e_ext = find_phi_seam_partner(em.extra_coords, em.extra_connijk, em.extra_nop,
                                           em.extra_nelem, target_e_ext)
    if partner_e_ext == 0
        adapt_single_angular_element_with_metrics!(mesh, target_iel, target_e_ext)
        return nothing
    end

    nop_p = em.extra_nop[partner_e_ext]
    θ_lo_p = em.extra_coords[1, em.extra_connijk[partner_e_ext,1,1]]
    θ_hi_p = em.extra_coords[1, em.extra_connijk[partner_e_ext,nop_p+1,1]]
    ϕ_first_p = em.extra_coords[2, em.extra_connijk[partner_e_ext,1,1]]
    ϕ_last_p  = em.extra_coords[2, em.extra_connijk[partner_e_ext,1,nop_p+1]]

    adapt_single_angular_element_with_metrics!(mesh, target_iel, target_e_ext)

    em2 = mesh.extra_mesh[target_iel]
    partner_e_ext_now = 0
    for e = 1:em2.extra_nelem
        nop_e = em2.extra_nop[e]
        θ_lo_e = em2.extra_coords[1, em2.extra_connijk[e,1,1]]
        θ_hi_e = em2.extra_coords[1, em2.extra_connijk[e,nop_e+1,1]]
        ϕf_e = em2.extra_coords[2, em2.extra_connijk[e,1,1]]
        ϕl_e = em2.extra_coords[2, em2.extra_connijk[e,1,nop_e+1]]
        if abs(θ_lo_e-θ_lo_p)<1e-10 && abs(θ_hi_e-θ_hi_p)<1e-10 &&
           isapprox(ϕf_e,ϕ_first_p;atol=1e-10) && isapprox(ϕl_e,ϕ_last_p;atol=1e-10)
            partner_e_ext_now = e; break
        end
    end
    partner_e_ext_now == 0 && error("could not re-locate φ-seam partner element after refining target")

    adapt_single_angular_element_with_metrics!(mesh, target_iel, partner_e_ext_now)
    return nothing
end

"""
Groups spatial elements by their angular mesh's SIGNATURE — elements
sharing the identical extra_npoin, extra_nelem, and extra_coords (within
tolerance) can share ONE computation of w_ext/W/C_iel/Ω_m, rather than
recomputing per element. Refined elements (each geometrically distinct)
naturally end up in their own singleton groups.
"""
function group_elements_by_angular_mesh(mesh; tol::Float64=1e-8)
    groups = Dict{Int, Vector{Int}}()   # representative iel -> [iel,...] sharing its mesh
    representative_of = zeros(Int, mesh.nelem)

    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        found_group = 0
        for (rep, _) in groups
            em_rep = mesh.extra_mesh[rep]
            em_rep.extra_npoin == em.extra_npoin || continue
            em_rep.extra_nelem == em.extra_nelem || continue
            isapprox(em_rep.extra_coords, em.extra_coords; atol=tol) || continue
            found_group = rep
            break
        end
        if found_group == 0
            groups[iel] = [iel]
            representative_of[iel] = iel
        else
            push!(groups[found_group], iel)
            representative_of[iel] = found_group
        end
    end
    return groups, representative_of
end

function test_group_elements_by_angular_mesh(mesh)
    groups, representative_of = group_elements_by_angular_mesh(mesh)

    println("Total spatial elements: $(mesh.nelem)")
    println("Distinct angular-mesh signatures (groups): $(length(groups))")
    println("─"^50)
    for (rep, members) in groups
        em = mesh.extra_mesh[rep]
        println("  Group rep=$rep: $(length(members)) member(s), " *
                "extra_npoin=$(em.extra_npoin), extra_nelem=$(em.extra_nelem)")
    end

    # ── Correctness check: every member of a group must ACTUALLY share the
    # representative's angular mesh exactly (not just have been grouped by
    # a bug that silently accepts mismatches). ──────────────────────────────
    n_checked = 0
    n_mismatch = 0
    for (rep, members) in groups
        em_rep = mesh.extra_mesh[rep]
        for m in members
            n_checked += 1
            em_m = mesh.extra_mesh[m]
            ok = (em_m.extra_npoin == em_rep.extra_npoin) &&
                 (em_m.extra_nelem == em_rep.extra_nelem) &&
                 isapprox(em_m.extra_coords, em_rep.extra_coords; atol=1e-8)
            ok || (n_mismatch += 1; println("  ✗ MISMATCH: element $m grouped with rep=$rep but coords differ"))
        end
    end
    println("─"^50)
    println("Checked $n_checked (element, representative) pairs, $n_mismatch mismatches  " *
            (n_mismatch == 0 ? "✓ PASS" : "✗ FAIL"))

    # ── Sanity: every element should belong to EXACTLY one group ────────────
    all_members = vcat(values(groups)...)
    n_total_members = length(all_members)
    n_unique_members = length(unique(all_members))
    println("Total group memberships: $n_total_members, unique elements covered: $n_unique_members, " *
            "mesh.nelem=$(mesh.nelem)  " *
            ((n_total_members == n_unique_members == mesh.nelem) ? "✓ PASS" : "✗ FAIL"))

    return groups, representative_of
end

function test_combined_ncf_scheduling(inputs, nparts, distribute, sem, partitioned_model, spatial_target_iel::Int,
                                       angular_target_iel::Int, angular_target_e_ext::Int;
                                       max_rounds::Int = 500)
    # ── Apply BOTH refinements — spatial far from angular, per the agreed test design ──
    println("Applying spatial refinement at iel=$spatial_target_iel...")
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    ψ, dψ, ω = sem.basis.ψ, sem.basis.dψ, sem.ω

    # NOTE: exact call signature for spatial refinement depends on your
    # sem_setup-based workflow — using the same pattern as earlier spatial
    # NCF tests; adjust if your actual call site differs.

    println("Applying angular refinement at iel=$angular_target_iel...")
    adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
    verify_angular_phi_seam_conformity(mesh, angular_target_iel)

    ngl = mesh.ngl
    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    # ── Spatial NCF tables (from the spatial-NCF effort) ─────────────────────
    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
    ncf_child_to_parent_elemid = Dict{Tuple{Int,Int}, Int}()
    ncf_parent_to_children = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ncf_idx = 1:mesh.num_ncf
        child_iel  = Int(mesh.cip[ncf_idx])
        parent_iel = Int(mesh.pip[ncf_idx])
        child_fid  = Int(mesh.lfid[ncf_idx])
        parent_fid = OPP_FACE[child_fid]
        ncf_as_child[(child_iel, child_fid)] = ncf_idx
        ncf_child_to_parent_elemid[(child_iel,child_fid)] = parent_iel
        push!(get!(ncf_as_parent, (parent_iel,parent_fid), Int[]), ncf_idx)
        push!(get!(ncf_parent_to_children, (parent_iel,parent_fid), Int[]), child_iel)
    end
    println("Spatial NCF: $(mesh.num_ncf) entries")

    # ── Angular NCF dependency graph ──────────────────────────────────────────
    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    println("Angular NCF: $(length(angular_deps)) (iel,d) entries with explicit deps")

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    total_pairs = sum(n_dirs)
    n_solved_total = 0
    round = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("STALLED after $max_rounds rounds")

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]

                    if jel_raw == zero(TInt)
                        # ── Check SPATIAL NCF first (jel_raw==0 could mean
                        # true boundary OR a spatial-NCF face) ────────────────
                        if haskey(ncf_child_to_parent_elemid, (iel,fid))
                            # iel is the CHILD (spatial) — needs parent solved,
                            # SAME direction d (spatial NCF doesn't change
                            # angular meaning, only spatial coupling).
                            p_iel = ncf_child_to_parent_elemid[(iel,fid)]
                            (d <= n_dirs[p_iel] && solved[p_iel][d]) || (ready = false)
                        elseif haskey(ncf_as_parent, (iel,fid))
                            # iel is the PARENT (spatial) — needs ALL children solved.
                            for ncf_idx in ncf_as_parent[(iel,fid)]
                                c_iel = Int(mesh.cip[ncf_idx])
                                (d <= n_dirs[c_iel] && solved[c_iel][d]) || (ready = false)
                            end
                        end
                        # else: TRUE physical boundary — always ready, no check.
                        continue
                    end
                    jel = Int(jel_raw)

                    # ── Ordinary same-rank face — check ANGULAR NCF ──────────
                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            for dn in dep_nodes
                                solved[jel][dn] || (ready = false)
                            end
                        end
                    end
                    if !handled
                        (d <= n_dirs[jel] && !solved[jel][d]) && (ready = false)
                    end
                    ready || break
                end

                if ready
                    solved[iel][d] = true
                    n_newly_solved += 1
                end
            end
        end

        n_solved_total += n_newly_solved
        println("round $round: $n_newly_solved solved, $n_solved_total/$total_pairs")

        if n_newly_solved == 0 && n_solved_total < total_pairs
            cyclic_sccs = find_cyclic_sccs(mesh, neighbor_elem, angular_deps, solved, n_dirs)
            if isempty(cyclic_sccs)
                error("STALLED with $(total_pairs-n_solved_total) unsolved and NO angular cycle " *
                      "detected — could be a genuine bug, OR a cycle involving SPATIAL NCF edges " *
                      "that find_cyclic_sccs doesn't currently know about (flagged as a real gap).")
            end
            println("  Found $(length(cyclic_sccs)) cyclic SCC(s) (angular-only detection)")
            for scc in cyclic_sccs
                for (iel,d) in scc
                    solved[iel][d] = true   # placeholder mark for THIS scheduling-only test
                end
                n_solved_total += length(scc)
            end
        end
    end

    println("Converged in $round rounds")
    return solved, n_dirs
end

function run_combined_ncf_sweep!(mesh, metrics, ψ, dψ, ω,
                                  neighbor_elem, neighbor_fid, face_correspondence,
                                  ncf_as_child, ncf_L1L2, ncf_as_parent,
                                  angular_deps,
                                  κ, σ, intΦ_eff, bc_fn_all;
                                  extra_source_all = nothing,
                                  max_rounds::Int = 500, max_cycle_iter::Int = 100, cycle_tol::Float64 = 1e-9)
    ngl = mesh.ngl
    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    solved = [falses(n_dirs[iel]) for iel = 1:mesh.nelem]
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()

    total_pairs = sum(n_dirs)
    n_solved_total = 0
    round = 0

    while n_solved_total < total_pairs
        round += 1
        round > max_rounds && error("run_combined_ncf_sweep!: exceeded max_rounds=$max_rounds — " *
                                     "$(total_pairs-n_solved_total) unresolved")

        n_newly_solved = 0
        for iel = 1:mesh.nelem
            for d = 1:n_dirs[iel]
                solved[iel][d] && continue
                em = mesh.extra_mesh[iel]
                θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
                Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)

                is_inflow = Vector{Bool}(undef,6); is_active = Vector{Bool}(undef,6)
                classify_element_faces!(is_inflow, is_active, iel, Ωx, Ωy, Ωz, ngl,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz)

                ready = true
                upwind_fns = Dict{Int, Function}()

                for fid = 1:6
                    (is_inflow[fid] && is_active[fid]) || continue
                    jel_raw = neighbor_elem[iel][fid]

                    if jel_raw == zero(TInt)
                        # ── SPATIAL NCF check first (a jel_raw==0 face is
                        # either a true boundary, a spatial-NCF child face,
                        # or a spatial-NCF parent face). ────────────────────
                        if haskey(ncf_as_child, (iel,fid))
                            ncf_idx = ncf_as_child[(iel,fid)]
                            parent_iel = Int(mesh.pip[ncf_idx])
                            if d > n_dirs[parent_iel] || !solved[parent_iel][d]
                                ready = false
                            else
                                L1,L2 = ncf_L1L2[ncf_idx]
                                IPc_nodes = mesh.IPc_list[:,ncf_idx]; IPp_nodes = mesh.IPp_list[:,ncf_idx]
                                # NOTE: build_dg_ncf_parent_to_child needs a
                                # GLOBAL-indexed x_dg (per its original
                                # spatial-NCF design) — but our combined x_dg
                                # is now (iel,d)-keyed, LOCAL ngl^3 storage.
                                # Needs an adapter — see caveat below.
                                child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, iel, parent_iel, fid, d, IPc_nodes, IPp_nodes, L1, L2, ngl)
				upwind_fns[fid] = (i,j,k) -> child_vals[face_flat_index(fid, ngl, i, j, k)]
                            end
                        elseif haskey(ncf_as_parent, (iel,fid))
                            ncf_group = ncf_as_parent[(iel,fid)]
                            children = [Int(mesh.cip[idx]) for idx in ncf_group]
                            if any(c -> d > n_dirs[c] || !solved[c][d], children)
                                ready = false
                            else
                                entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx])
                                           for idx in ncf_group]
                                    parent_grid = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, iel, fid, d, entries, ngl)
				    upwind_fns[fid] = (i,j,k) -> begin
                                    node_idx = face_flat_index(fid, ngl, i, j, k)
                                    p1 = ((node_idx-1) % ngl) + 1
                                    p2 = ((node_idx-1) ÷ ngl) + 1
                                    parent_grid[p1,p2]
                                end
                            end
                        else
                            # TRUE physical boundary
                            bc = bc_fn_all[(iel,d)]
                            upwind_fns[fid] = (i,j,k) -> begin
                                gip = mesh.connijk[iel,i,j,k]
                                bc(mesh.x[gip], mesh.y[gip], mesh.z[gip])
                            end
                        end
                        continue
                    end
                    jel = Int(jel_raw)

                    # ── ANGULAR NCF check (ordinary same-rank spatial face) ──
                    key = (iel, d)
                    handled = false
                    if haskey(angular_deps, key)
                        for (dep_jel, dep_nodes) in angular_deps[key]
                            dep_jel == jel || continue
                            handled = true
                            if length(dep_nodes) == 1
                                dj = dep_nodes[1]
                                if !solved[jel][dj]
                                    ready = false
                                else
                                    dj_c = dj; fid_c = fid
                                    upwind_fns[fid] = (i,j,k) -> begin
                                        corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                        x_dg[(jel,dj_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                                    end
                                end
                            else
                                all_solved = all(solved[jel][dn] for dn in dep_nodes)
                                if !all_solved
                                    ready = false
                                else
                                    em_jel = mesh.extra_mesh[jel]
                                    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
                                    upwind_fns[fid] = build_angular_ncf_upwind_unified(
                                        iel, fid, jel, d, src_e_ext, mesh, em_jel, x_dg, face_correspondence, ngl)
                                end
                            end
                            break
                        end
                    end
                    if !handled
                        if d > n_dirs[jel] || !solved[jel][d]
                            ready = false
                        else
                            d_c = d; fid_c = fid
                            upwind_fns[fid] = (i,j,k) -> begin
                                corr = face_correspondence[(iel,fid_c)][(i,j,k)]
                                x_dg[(jel,d_c)][lidx(corr[1],corr[2],corr[3],ngl)]
                            end
                        end
                    end
                    ready || break
                end
                ready || continue

                A_vol = zeros(ngl^3,ngl^3); RHS_local = zeros(ngl^3)
                build_dg_volume_matrix!(A_vol, iel, Ωx,Ωy,Ωz, ngl, ψ,dψ,ω,
                                         metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                         metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                         metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                         metrics.Je, κ, σ, intΦ_eff, mesh.connijk)
                if extra_source_all !== nothing && haskey(extra_source_all, (iel,d))
                    build_dg_source_term!(RHS_local, iel, ngl, ω, metrics.Je, extra_source_all[(iel,d)], mesh.connijk)
                end
                for fid = 1:6
                    is_active[fid] || continue
                    face_normal = compute_face_normal_avg(iel, fid, ngl,
                        metrics.dξdx,metrics.dξdy,metrics.dξdz,
                        metrics.dηdx,metrics.dηdy,metrics.dηdz,
                        metrics.dζdx,metrics.dζdy,metrics.dζdz)
                    I_upwind = haskey(upwind_fns, fid) ? upwind_fns[fid] : ((i,j,k)->0.0)
                    build_dg_flux_terms!(A_vol, RHS_local, fid, ngl, Ωx,Ωy,Ωz,
                                          is_inflow[fid], is_active[fid], face_normal, ω,
                                          metrics.dξdx,metrics.dξdy,metrics.dξdz,
                                          metrics.dηdx,metrics.dηdy,metrics.dηdz,
                                          metrics.dζdx,metrics.dζdy,metrics.dζdz,
                                          metrics.Je, iel, I_upwind)
                end
                x_dg[(iel,d)] = A_vol \ RHS_local
                solved[iel][d] = true
                n_newly_solved += 1
            end
        end

        n_solved_total += n_newly_solved
        println("round $round: $n_newly_solved solved, $n_solved_total/$total_pairs")

        if n_newly_solved == 0 && n_solved_total < total_pairs
            cyclic_sccs = find_cyclic_sccs(mesh, metrics, neighbor_elem, angular_deps, solved, n_dirs)
            isempty(cyclic_sccs) && error("STALLED with $(total_pairs-n_solved_total) unsolved, no cycle found")
            for scc in cyclic_sccs
                touches_spatial_ncf = any(iel -> any(fid -> haskey(ncf_as_child,(iel,fid)) || haskey(ncf_as_parent,(iel,fid)), 1:6),
                               [p[1] for p in scc])
    	 	println("  SCC of size $(length(scc)) touches a spatial-NCF face? $touches_spatial_ncf") 
	        resolve_angular_cycle!(mesh, metrics, ψ, dψ, ω, neighbor_elem, neighbor_fid,
                                       face_correspondence, angular_deps, x_dg, solved, n_dirs,
                                       scc, κ, σ, intΦ_eff, bc_fn_all;
                                       extra_source_all=extra_source_all, max_iter=max_cycle_iter, tol=cycle_tol)
                n_solved_total += length(scc)
            end
        end
    end

    return x_dg, n_dirs
end

function build_dg_ncf_parent_to_child_local(mesh, x_dg, child_iel::Int, parent_iel::Int, fid::Int, d::Int,
                                             IPc_nodes, IPp_nodes, L1::Matrix{Float64}, L2::Matrix{Float64}, ngl::Int)
    parent_vals = x_dg[(parent_iel, d)]   # LOCAL, ngl^3

    parent_ijk = Dict{Int,NTuple{3,Int}}()
    for k=1:ngl, j=1:ngl, i=1:ngl
        parent_ijk[mesh.connijk[parent_iel,i,j,k]] = (i,j,k)
    end

    parent_grid = zeros(ngl, ngl)
    for idx = 1:length(IPp_nodes)
        node = Int(IPp_nodes[idx])
        node <= 0 && continue
        p1 = ((idx-1) % ngl) + 1
        p2 = ((idx-1) ÷ ngl) + 1
        (i,j,k) = parent_ijk[node]
        parent_grid[p1,p2] = parent_vals[lidx(i,j,k,ngl)]   # LOCAL lookup, was x_dg[connijk_dg[...]]
    end

    child_vals = zeros(length(IPc_nodes))
    for node_idx = 1:length(IPc_nodes)
        node = Int(IPc_nodes[node_idx])
        node <= 0 && continue
        ci = div(node_idx-1, ngl) + 1
        cj = mod(node_idx-1, ngl) + 1
        (ci < 1 || ci > size(L1,1) || cj < 1 || cj > size(L2,1)) && continue
        w_1 = L1[cj, :]
        w_2 = L2[ci, :]
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += w_1[p1]*w_2[p2]*parent_grid[p1,p2]
        end
        child_vals[node_idx] = val
    end
    return child_vals
end

function build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, parent_iel::Int, fid::Int, d::Int,
                                                ncf_entries, ngl::Int)
    parent_result = fill(NaN, ngl, ngl)

    for (ncf_idx, child_iel, IPc_nodes, IPp_nodes) in ncf_entries
        child_vals = x_dg[(child_iel, d)]   # LOCAL

        child_ijk = Dict{Int,NTuple{3,Int}}()
        for k=1:ngl, j=1:ngl, i=1:ngl
            child_ijk[mesh.connijk[child_iel,i,j,k]] = (i,j,k)
        end

        child_line_1 = Int[]; child_line_2 = Int[]
        for i = 1:ngl
            idx = (i-1)*ngl + 1; idx <= length(IPc_nodes) || continue
            push!(child_line_1, Int(IPc_nodes[idx]))
        end
        for j = 1:ngl
            idx = j; idx <= length(IPc_nodes) || continue
            push!(child_line_2, Int(IPc_nodes[idx]))
        end

        xs2=[mesh.x[n] for n in child_line_2]; ys2=[mesh.y[n] for n in child_line_2]; zs2=[mesh.z[n] for n in child_line_2]
        axis_c1 = axis_of_max_spread((xs2,ys2,zs2))
        xs1=[mesh.x[n] for n in child_line_1]; ys1=[mesh.y[n] for n in child_line_1]; zs1=[mesh.z[n] for n in child_line_1]
        axis_c2 = axis_of_max_spread((xs1,ys1,zs1))
        coord(node,axis) = axis==1 ? mesh.x[node] : (axis==2 ? mesh.y[node] : mesh.z[node])

        c1_coords = [coord(n,axis_c1) for n in child_line_2]
        c2_coords = [coord(n,axis_c2) for n in child_line_1]
        c1_min,c1_max = extrema(c1_coords); c2_min,c2_max = extrema(c2_coords)
        tol1 = 1e-6*max(c1_max-c1_min,1.0); tol2 = 1e-6*max(c2_max-c2_min,1.0)

        child_grid = zeros(ngl,ngl)
        for idx = 1:length(IPc_nodes)
            node = Int(IPc_nodes[idx]); node<=0 && continue
            c1 = ((idx-1)%ngl)+1; c2 = ((idx-1)÷ngl)+1
            child_grid[c1,c2] = child_vals[lidx(child_ijk[node]...,ngl)]   # LOCAL lookup
        end

        for idx = 1:length(IPp_nodes)
            node = Int(IPp_nodes[idx]); node<=0 && continue
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            pc1 = coord(node,axis_c1); pc2 = coord(node,axis_c2)
            (pc1 < c1_min-tol1 || pc1 > c1_max+tol1 || pc2 < c2_min-tol2 || pc2 > c2_max+tol2) && continue

            L_c1 = build_1d_lagrange_matrix(c1_coords, [pc1])
            L_c2 = build_1d_lagrange_matrix(c2_coords, [pc2])
            val = 0.0
            for a=1:ngl, b=1:ngl
                val += L_c1[1,a]*L_c2[1,b]*child_grid[a,b]
            end
            parent_result[p1,p2] = val
        end
    end
    return parent_result
end

function test_dg_ncf_local_storage(inputs, nparts, distribute, sem, partitioned_model)
    rank = MPI.Comm_rank(get_mpi_comm())
    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[1] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh = sem.mesh; ngl = mesh.ngl

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)

    ncf_idx = 1
    child_iel  = Int(mesh.cip[ncf_idx])
    parent_iel = Int(mesh.pip[ncf_idx])
    local_facet_id = Int(mesh.lfid[ncf_idx])
    IPc_nodes = mesh.IPc_list[:, ncf_idx]
    IPp_nodes = mesh.IPp_list[:, ncf_idx]
    L1, L2 = build_spatial_interpolation_matrices_from_facet_nodes(
        child_iel, parent_iel, local_facet_id, IPc_nodes, IPp_nodes, mesh, ngl)

    x0,y0,z0 = mesh.xmin, mesh.ymin, mesh.zmin
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    P(t) = sum(0.3^k * t^k for k=0:ngl-1)
    field(x,y,z) = P((x-x0)/Lx) + P((y-y0)/Ly) + P((z-z0)/Lz)

    d_fake = 1
    x_dg = Dict{Tuple{Int,Int}, Vector{Float64}}()
    vals = zeros(ngl^3)
    for k=1:ngl,j=1:ngl,i=1:ngl
        ip = mesh.connijk[parent_iel,i,j,k]
        vals[lidx(i,j,k,ngl)] = field(mesh.x[ip], mesh.y[ip], mesh.z[ip])
    end
    x_dg[(parent_iel,d_fake)] = vals

    child_vals = build_dg_ncf_parent_to_child_local(mesh, x_dg, child_iel, parent_iel, local_facet_id, d_fake,
                                                     IPc_nodes, IPp_nodes, L1, L2, ngl)

    max_err_p2c = 0.0
    for idx = 1:length(IPc_nodes)
        node = Int(IPc_nodes[idx]); node<=0 && continue
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        max_err_p2c = max(max_err_p2c, abs(child_vals[idx] - exact))
    end
    println("Parent->child LOCAL storage: max abs error = $max_err_p2c")

    # ── Child -> parent, same style ────────────────────────────────────────
    groups = Dict{Tuple{Int,Int}, Vector{Int}}()
    for idx = 1:mesh.num_ncf
        key = (Int(mesh.pip[idx]), Int(mesh.lfid[idx]))
        push!(get!(groups, key, Int[]), idx)
    end
    valid_group = nothing
    for (key, members) in groups
        length(members) == 4 && (valid_group = (key, members); break)
    end
    (parent_iel2, lfid2), member_idxs = valid_group

    x_dg2 = Dict{Tuple{Int,Int}, Vector{Float64}}()
    ncf_entries = Tuple{Int,Int,Vector,Vector}[]
    for idx in member_idxs
        c_iel = Int(mesh.cip[idx])
        IPc = mesh.IPc_list[:,idx]; IPp = mesh.IPp_list[:,idx]
        push!(ncf_entries, (idx, c_iel, IPc, IPp))
        vals2 = zeros(ngl^3)
        for k=1:ngl,j=1:ngl,i=1:ngl
            ip = mesh.connijk[c_iel,i,j,k]
            vals2[lidx(i,j,k,ngl)] = field(mesh.x[ip], mesh.y[ip], mesh.z[ip])
        end
        x_dg2[(c_iel,d_fake)] = vals2
    end

    parent_computed = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg2, parent_iel2, lfid2, d_fake, ncf_entries, ngl)
    IPp_ref = mesh.IPp_list[:, member_idxs[1]]
    max_err_c2p = 0.0
    for idx = 1:length(IPp_ref)
        node = Int(IPp_ref[idx]); node<=0 && continue
        p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
        exact = field(mesh.x[node], mesh.y[node], mesh.z[node])
        max_err_c2p = max(max_err_c2p, abs(parent_computed[p1,p2] - exact))
    end
    println("Child->parent LOCAL storage: max abs error = $max_err_c2p")

    return max_err_p2c, max_err_c2p
end

function test_combined_ncf_exponential(inputs, nparts, distribute, sem, partitioned_model;
                                        spatial_target_iel::Int,
                                        angular_target_iel::Int, angular_target_e_ext::Int,
                                        κ_val=0.0001, I0=4.2)

    # ── Apply BOTH refinements — spatial far from angular, per the agreed test design ──
    println("Applying spatial refinement at iel=$spatial_target_iel...")
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    ψ, dψ, ω_gll = sem.basis.ψ, sem.basis.dψ, sem.ω

    # ── Apply both refinements, far apart per the agreed test design ────────
    adapt_flags = zeros(Int, mesh.nelem)
    adapt_flags[spatial_target_iel] = 1
    # (spatial refinement call site — matching your sem_setup workflow)

    adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
    verify_angular_phi_seam_conformity(mesh, angular_target_iel)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
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
    println("Spatial NCF: $(mesh.num_ncf) entries")

    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    println("Angular NCF: $(length(angular_deps)) entries")

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]
    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)
    r_bc = (mesh.xmin, mesh.ymin, mesh.zmin)

    bc_fn_all = Dict{Tuple{Int,Int}, Function}()
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:em.extra_npoin
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            bc_fn_all[(iel,d)] = (x,y,z) -> I0*exp(-κ_val*(Ωx*(x-r_bc[1])+Ωy*(y-r_bc[2])+Ωz*(z-r_bc[3])))
        end
    end

    x_dg, n_dirs_ret = run_combined_ncf_sweep!(mesh, metrics, ψ, dψ, ω_gll,
                                               neighbor_elem, neighbor_fid, face_correspondence,
                                               ncf_as_child, ncf_L1L2, ncf_as_parent,
                                               angular_deps, κ, σz, intΦ_eff, bc_fn_all)

    max_err = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                exact = I0*exp(-κ_val*(Ωx*(mesh.x[g]-r_bc[1])+Ωy*(mesh.y[g]-r_bc[2])+Ωz*(mesh.z[g]-r_bc[3])))
                computed = x_dg[(iel,d)][lidx(i,j,k,ngl)]
                max_err = max(max_err, abs(computed - exact))
            end
        end
    end
    println("Combined spatial+angular NCF exponential test — max abs error: $max_err")
    return max_err
end

function test_manufactured_isotropic_combined_ncf(inputs, nparts, distribute, sem, partitioned_model;
                                                   spatial_target_iel::Int,
                                                   angular_target_iel::Int, angular_target_e_ext::Int,
                                                   I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
                                                   κ_val=0.001, σ_val=0.0005,
                                                   max_iter=500, tol=1e-9)
    println("Applying spatial refinement at iel=$spatial_target_iel...")
    rank = MPI.Comm_rank(get_mpi_comm())

    adapt_flags = zeros(Int, sem.mesh.nelem)
    if sem.mesh.nelem > 0 && rank == 0
        adapt_flags[spatial_target_iel] = 1
    end
    sem, partitioned_model = sem_setup(inputs, nparts, distribute,
                                        adapt_flags, partitioned_model,
                                        sem.mesh, sem.interp, nothing, nothing)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    ψ, dψ, ω_gll = sem.basis.ψ, sem.basis.dψ, sem.ω

    adapt_single_angular_element_conforming!(mesh, angular_target_iel, angular_target_e_ext)
    verify_angular_phi_seam_conformity(mesh, angular_target_iel)

    neighbor_elem, neighbor_fid = build_face_adjacency_from_conn_face_el(mesh.nelem, mesh.conn_face_el)
    connijk_dg, npoin_dg, dg_to_cg, face_correspondence = build_connijk_dg(mesh, neighbor_elem, neighbor_fid)

    ncf_as_child = Dict{Tuple{Int,Int}, Int}()
    ncf_L1L2     = Dict{Int, Tuple{Matrix{Float64},Matrix{Float64}}}()
    ncf_as_parent = Dict{Tuple{Int,Int}, Vector{Int}}()
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
    println("Spatial NCF: $(mesh.num_ncf) entries")

    angular_deps = build_angular_dependency_graph(mesh, neighbor_elem)
    println("Angular NCF: $(length(angular_deps)) entries")

    n_dirs = [mesh.extra_mesh[iel].extra_npoin for iel = 1:mesh.nelem]

    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    I_manuf(x,y,z,Ωx,Ωy,Ωz) = I0 + A*sin(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz)
    dI_manuf(x,y,z,Ωx,Ωy,Ωz) = A*cos(kx*x+ky*y+kz*z)*(1+b*Ωx+c*Ωy+d*Ωz) .* (kx,ky,kz)

    npoin_cg = length(mesh.x)
    κ = fill(κ_val, npoin_cg); σz = zeros(npoin_cg); intΦ_eff = zeros(npoin_cg)

    # ── Group elements by angular-mesh signature — same optimization as the
    # angular-only manufactured test, still valid here since spatial
    # refinement doesn't change any element's ANGULAR mesh. ─────────────────
    groups, representative_of = group_elements_by_angular_mesh(mesh)
    println("Grouped $(mesh.nelem) elements into $(length(groups)) distinct angular-mesh signatures")

    w_ext_all = Dict{Int, Vector{Float64}}()
    W_all = Dict{Int, Float64}()
    C_all = Dict{Int, Matrix{Float64}}()
    Ωm_all = Dict{Int, Vector{NTuple{3,Float64}}}()
    coeff_all = Dict{Int, Vector{Float64}}()

    for rep in keys(groups)
        em = mesh.extra_mesh[rep]
        w_ext = zeros(Float64, em.extra_npoin)
        for e_ext = 1:em.extra_nelem
            nop = em.extra_nop[e_ext]
            for iθ=1:nop+1, iϕ=1:nop+1
                ip_ext = em.extra_connijk[e_ext,iθ,iϕ]
                w_ext[ip_ext] += em.extra_metrics.Je[e_ext,iθ,iϕ] * em.ωθ[iθ] * em.ωθ[iϕ]
            end
        end
        W = sum(w_ext)
        w_ext_all[rep] = w_ext; W_all[rep] = W
        C_all[rep] = [σ_val*w_ext[m]/W for kk=1:em.extra_npoin, m=1:em.extra_npoin]
        Ωm_all[rep] = [(sin(em.extra_coords[1,m])*cos(em.extra_coords[2,m]),
                        sin(em.extra_coords[1,m])*sin(em.extra_coords[2,m]),
                        cos(em.extra_coords[1,m])) for m = 1:em.extra_npoin]
        coeff_all[rep] = σ_val .* w_ext ./ W
    end

    q_ext_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    bc_fn_all = Dict{Tuple{Int,Int}, Function}()

    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        rep = representative_of[iel]
        C_iel = C_all[rep]; Ω_m = Ωm_all[rep]

        for dk = 1:em.extra_npoin
            θk = em.extra_coords[1,dk]; ϕk = em.extra_coords[2,dk]
            Ωxk = sin(θk)*cos(ϕk); Ωyk = sin(θk)*sin(ϕk); Ωzk = cos(θk)

            q_vec = zeros(npoin_cg)
            for k=1:ngl, j=1:ngl, i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                x_,y_,z_ = mesh.x[g], mesh.y[g], mesh.z[g]
                grad = dI_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                Ω_dot_grad = Ωxk*grad[1] + Ωyk*grad[2] + Ωzk*grad[3]
                I_true = I_manuf(x_,y_,z_,Ωxk,Ωyk,Ωzk)
                scatter_true = sum(C_iel[dk,m]*I_manuf(x_,y_,z_, Ω_m[m]...) for m = 1:em.extra_npoin)
                q_vec[g] = Ω_dot_grad + κ_val*I_true - scatter_true
            end
            q_ext_all[(iel,dk)] = q_vec
            bc_fn_all[(iel,dk)] = (x,y,z) -> I_manuf(x,y,z,Ωxk,Ωyk,Ωzk)
        end
    end

    I_cg_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
        I_cg_all[(iel,d)] = zeros(ngl^3)
    end

    for it = 1:max_iter
        extra_source_all = Dict{Tuple{Int,Int}, Vector{Float64}}()
        for iel = 1:mesh.nelem
            em = mesh.extra_mesh[iel]
            rep = representative_of[iel]
            coeff = coeff_all[rep]

            local_gids = Vector{Int}(undef, ngl^3)
            for k=1:ngl, j=1:ngl, i=1:ngl
                local_gids[lidx(i,j,k,ngl)] = mesh.connijk[iel,i,j,k]
            end

            for dk = 1:em.extra_npoin
                q_k = copy(q_ext_all[(iel,dk)])
                q_k_local = @view q_k[local_gids]
                for m = 1:em.extra_npoin
                    coeff[m] == 0.0 && continue
                    q_k_local .+= coeff[m] .* I_cg_all[(iel,m)]
                end
                extra_source_all[(iel,dk)] = q_k
            end
        end

        x_dg, _ = run_combined_ncf_sweep!(mesh, metrics, ψ, dψ, ω_gll,
                                          neighbor_elem, neighbor_fid, face_correspondence,
                                          ncf_as_child, ncf_L1L2, ncf_as_parent,
                                          angular_deps, κ, σz, intΦ_eff, bc_fn_all;
                                          extra_source_all=extra_source_all)

        max_delta = 0.0
        for iel = 1:mesh.nelem, d = 1:n_dirs[iel]
            max_delta = max(max_delta, maximum(abs.(x_dg[(iel,d)] .- I_cg_all[(iel,d)])))
            I_cg_all[(iel,d)] = copy(x_dg[(iel,d)])
        end
        println("  iter $it: max_delta=$max_delta")
        max_delta < tol && (println("Converged after $it iterations"); break)
    end

    max_err = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for d = 1:n_dirs[iel]
            θ = em.extra_coords[1,d]; ϕ = em.extra_coords[2,d]
            Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
            for k=1:ngl, j=1:ngl, i=1:ngl
                g = mesh.connijk[iel,i,j,k]
                I_true = I_manuf(mesh.x[g],mesh.y[g],mesh.z[g],Ωx,Ωy,Ωz)
                computed = I_cg_all[(iel,d)][lidx(i,j,k,ngl)]
                max_err = max(max_err, abs(computed - I_true))
            end
        end
    end
    println("Manufactured (isotropic, combined spatial+angular NCF) — max abs error: $max_err")
    return max_err
end

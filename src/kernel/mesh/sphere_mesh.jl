#---------------------------------------------------------------------------------
# sphere_mesh.jl
#
# Build a HIGH-ORDER (LGL) spectral-element grid on a CLOSED SPHERICAL SHELL
# starting from a linear quadrilateral gmsh grid of the shell surface
# (e.g. meshes/gmsh_grids/cubed_sphere.msh).
#
# WHY THIS IS NOT `mod_mesh_read_gmsh!`
# ------------------------------------
# The existing 2D reader in mesh.jl treats a 2D grid as a FLAT patch: it
# stores only (x,y) for the linear nodes and interpolates only (x,y) when it
# adds the high-order points (see add_high_order_nodes_edges!(…, NSD_2D)).
# A spherical shell is a 2D MANIFOLD EMBEDDED IN 3D: every node needs
# (x,y,z), the high-order points must be placed ON the sphere (not on the
# chord of the linear element), and — most importantly — the surface is
# CLOSED: it has NO boundary. That last property is what makes the node
# numbering delicate, and it is the part this file is careful about:
#
#   1. gmsh grids of a cubed sphere are very often produced panel-by-panel.
#      When the six panels are separate surfaces that do NOT share their
#      boundary curves, gmsh writes DUPLICATE (geometrically coincident)
#      nodes along the panel seams. Read naively, the "closed shell" is then
#      six disconnected patches that merely look watertight, and every
#      CG/SEM assembly across a seam is silently wrong. We therefore MERGE
#      geometrically coincident vertices before touching the topology
#      (merge_coincident_nodes!).
#
#   2. On a closed shell EVERY edge is an interior edge shared by EXACTLY
#      two elements, and there is no boundary to anchor a global ordering
#      to. High-order edge points must therefore be created ONCE per UNIQUE
#      edge and then re-used (in the correct traversal direction) by both
#      elements that touch it. Creating them per element would duplicate
#      every seam node — the classic bug on a closed surface. The unique
#      edge table is built here from scratch (build_unique_edges) instead of
#      relying on a Gridap topology, so the numbering is explicit and
#      verifiable.
#
#   3. A closed shell is ORIENTABLE but gmsh does not guarantee a globally
#      consistent element orientation. We re-orient every quad so that its
#      normal points OUTWARD (away from the sphere centre) and then VERIFY
#      that each unique edge is traversed exactly once in each direction —
#      the standard consistency test for a closed oriented manifold.
#
#   4. Everything is verified afterwards: Euler characteristic
#      (V - E + F = 2 for a sphere), edge multiplicity, single connected
#      component, exact high-order node count, and a geometric duplicate-node
#      test on the FINAL high-order grid (check_sphere_mesh).
#
# The result is a St_mesh_sphere struct + a VTK writer
# (write_vtk_sphere_grid in src/io/write_output.jl, next to
# write_vtk_grid_only) so that the grid can be inspected in ParaView
# before any equation is solved on it.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_mesh_sphere
export mod_mesh_sphere_driver
export build_sphere_shell_mesh
export read_gmsh_quad_surface
export check_sphere_mesh
export sphere_mesh_area
export merge_coincident_nodes
export build_unique_edges

#---------------------------------------------------------------------------------
# The spherical-shell spectral-element grid.
#
# Storage conventions (2D manifold in 3D):
#
#   connijk[iel, i, j]  i,j ∈ 1:ngl   the (i,j) LGL node of element iel.
#                                      i runs along ξ (v1→v2), j along η (v1→v4).
#   conn[iel, 1:4]                     the four LINEAR corners v1..v4 of iel,
#                                      counter-clockwise SEEN FROM OUTSIDE.
#   poin_in_edge[iedge, 1:ngl]         the LGL nodes of unique edge iedge, in
#                                      the canonical direction
#                                      conn_unique_edges[iedge,1] → [iedge,2].
#   edge_in_elem[iel, 1:4]             global id of the 4 local edges of iel,
#                                      local edge ordering:
#                                        1: (v1,v2)  j = 1
#                                        2: (v2,v3)  i = ngl
#                                        3: (v4,v3)  j = ngl
#                                        4: (v1,v4)  i = 1
#   node_type[ip]                      0 = linear vertex, 1 = edge node,
#                                      2 = element-interior node.
#
# Global node numbering (contiguous, no gaps, no duplicates):
#
#   1                              … npoin_linear                : vertices
#   npoin_linear + 1               … npoin_linear + nedges*(ngl-2): edge nodes
#   npoin_linear + nedges*(ngl-2)+1… npoin                        : interiors
#
#   npoin = npoin_linear + nedges*(ngl-2) + nelem*(ngl-2)^2
#
# which is EXACT on a closed shell precisely because every edge is counted
# once in `nedges` (no boundary edges, no duplicated seam edges).
#---------------------------------------------------------------------------------
mutable struct St_mesh_sphere{TFloat, TInt}

    nop::TInt                       # polynomial order
    ngl::TInt                       # nop + 1
    radius::TFloat                  # shell radius

    nelem::TInt                     # number of spectral elements (quads)
    npoin::TInt                     # total number of high-order nodes
    npoin_linear::TInt              # number of linear (gmsh) vertices
    nedges::TInt                    # number of UNIQUE edges

    x::Array{TFloat, 1}
    y::Array{TFloat, 1}
    z::Array{TFloat, 1}

    lon::Array{TFloat, 1}           # longitude λ ∈ (-π, π]
    lat::Array{TFloat, 1}           # latitude  φ ∈ [-π/2, π/2]

    conn::Array{TInt, 2}            # nelem × 4
    connijk::Array{TInt, 3}         # nelem × ngl × ngl
    conn_unique_edges::Array{TInt, 2}  # nedges × 2
    poin_in_edge::Array{TInt, 2}       # nedges × ngl
    edge_in_elem::Array{TInt, 2}       # nelem  × 4

    elem_tag::Array{TInt, 1}        # gmsh physical/entity tag (panel id)
    node_type::Array{TInt, 1}       # 0 vertex / 1 edge / 2 interior

    SD::NSD_2D                      # a shell is a 2D manifold
end


#---------------------------------------------------------------------------------
#
# gmsh reader for a QUADRILATERAL SURFACE grid (MSH 2.2 and MSH 4.1, ASCII)
#
# Returns
#   coords   :: npoin_linear × 3 node coordinates (in the order they are used
#               by `quads`, i.e. gmsh tags already compacted to 1:npoin_linear)
#   quads    :: nelem × 4 node ids
#   tags     :: nelem   physical (2.2) or entity (4.1) tag of each quad
#
# Only 4-node quads (gmsh element type 3) are retained: point/line elements
# that gmsh writes for the underlying geometry are skipped, and triangles are
# rejected loudly because the SEM kernel is tensor-product.
#---------------------------------------------------------------------------------
function read_gmsh_quad_surface(fname::String, TFloat=Float64, TInt=Int64)

    isfile(fname) || error(string(" # ERROR sphere_mesh.jl: gmsh file not found: ", fname))

    lines   = readlines(fname)
    nlines  = length(lines)

    version = 2.2
    i = 1
    while i <= nlines
        if startswith(strip(lines[i]), "\$MeshFormat")
            version = parse(Float64, split(strip(lines[i+1]))[1])
            break
        end
        i += 1
    end

    if version >= 4.0
        node_tags, node_xyz, quads_tagged = _read_msh41(lines, TFloat, TInt)
    else
        node_tags, node_xyz, quads_tagged = _read_msh22(lines, TFloat, TInt)
    end

    isempty(quads_tagged) && error(string(" # ERROR sphere_mesh.jl: no 4-node quadrilateral (gmsh type 3) elements found in ", fname))

    #
    # Compact the gmsh node tags (which may be sparse and 1-based-but-gappy)
    # to a dense 1:npoin_linear numbering.
    #
    tag2ip = Dict{TInt, TInt}()
    for (k, t) in enumerate(node_tags)
        tag2ip[t] = TInt(k)
    end

    nelem  = length(quads_tagged)
    quads  = Array{TInt}(undef, nelem, 4)
    tags   = Array{TInt}(undef, nelem)
    for (iel, (etag, nds)) in enumerate(quads_tagged)
        for k = 1:4
            haskey(tag2ip, nds[k]) || error(string(" # ERROR sphere_mesh.jl: element ", iel,
                                                   " references undefined node tag ", nds[k]))
            quads[iel, k] = tag2ip[nds[k]]
        end
        tags[iel] = etag
    end

    return node_xyz, quads, tags
end


function _read_msh22(lines, TFloat, TInt)

    nlines    = length(lines)
    node_tags = TInt[]
    node_xyz  = Array{TFloat}(undef, 0, 3)
    quads     = Tuple{TInt, NTuple{4, TInt}}[]

    i = 1
    while i <= nlines
        header = strip(lines[i])

        if header == "\$Nodes"
            nnodes    = parse(Int, strip(lines[i+1]))
            node_tags = Array{TInt}(undef, nnodes)
            node_xyz  = Array{TFloat}(undef, nnodes, 3)
            for k = 1:nnodes
                f = split(strip(lines[i+1+k]))
                node_tags[k]   = parse(TInt, f[1])
                node_xyz[k, 1] = parse(TFloat, f[2])
                node_xyz[k, 2] = parse(TFloat, f[3])
                node_xyz[k, 3] = parse(TFloat, f[4])
            end
            i += nnodes + 2

        elseif header == "\$Elements"
            nel = parse(Int, strip(lines[i+1]))
            for k = 1:nel
                f      = split(strip(lines[i+1+k]))
                etype  = parse(Int, f[2])
                ntags  = parse(Int, f[3])
                if etype == 3                      # 4-node quadrangle
                    phys = ntags >= 1 ? parse(TInt, f[4]) : TInt(0)
                    off  = 3 + ntags
                    push!(quads, (phys, (parse(TInt, f[off+1]), parse(TInt, f[off+2]),
                                         parse(TInt, f[off+3]), parse(TInt, f[off+4]))))
                elseif etype == 2 || etype == 9
                    error(" # ERROR sphere_mesh.jl: triangular elements found. The SEM kernel is tensor-product: the shell grid must be all-quadrilateral.")
                end
            end
            i += nel + 2
        else
            i += 1
        end
    end

    return node_tags, node_xyz, quads
end


function _read_msh41(lines, TFloat, TInt)

    nlines    = length(lines)
    node_tags = TInt[]
    node_xyz  = Array{TFloat}(undef, 0, 3)
    quads     = Tuple{TInt, NTuple{4, TInt}}[]

    i = 1
    while i <= nlines
        header = strip(lines[i])

        if header == "\$Nodes"
            f          = split(strip(lines[i+1]))
            nblocks    = parse(Int, f[1])
            nnodes     = parse(Int, f[2])
            node_tags  = Array{TInt}(undef, nnodes)
            node_xyz   = Array{TFloat}(undef, nnodes, 3)
            l          = i + 2
            ip         = 1
            for _ = 1:nblocks
                bh    = split(strip(lines[l]))
                parse(Int, bh[3]) == 0 ||
                    error(" # ERROR sphere_mesh.jl: MSH 4.1 files with parametric nodes are not supported. Re-export with `Mesh.SaveParametric = 0` (or in MSH 2.2 format).")
                nb    = parse(Int, bh[4])       # numNodesInBlock
                l    += 1
                # tags block, then coordinates block
                for k = 0:nb-1
                    node_tags[ip+k] = parse(TInt, strip(lines[l+k]))
                end
                l += nb
                for k = 0:nb-1
                    c = split(strip(lines[l+k]))
                    node_xyz[ip+k, 1] = parse(TFloat, c[1])
                    node_xyz[ip+k, 2] = parse(TFloat, c[2])
                    node_xyz[ip+k, 3] = parse(TFloat, c[3])
                end
                l  += nb
                ip += nb
            end
            i = l

        elseif header == "\$Elements"
            f       = split(strip(lines[i+1]))
            nblocks = parse(Int, f[1])
            l       = i + 2
            for _ = 1:nblocks
                bh     = split(strip(lines[l]))
                etag   = parse(TInt, bh[2])     # entityTag → used as panel id
                etype  = parse(Int,  bh[3])
                nb     = parse(Int,  bh[4])
                l     += 1
                if etype == 3
                    for k = 0:nb-1
                        c = split(strip(lines[l+k]))
                        push!(quads, (etag, (parse(TInt, c[2]), parse(TInt, c[3]),
                                             parse(TInt, c[4]), parse(TInt, c[5]))))
                    end
                elseif etype == 2 || etype == 9
                    error(" # ERROR sphere_mesh.jl: triangular elements found. The SEM kernel is tensor-product: the shell grid must be all-quadrilateral.")
                end
                l += nb
            end
            i = l
        else
            i += 1
        end
    end

    return node_tags, node_xyz, quads
end


#---------------------------------------------------------------------------------
# Merge geometrically coincident vertices.
#
# THE step that makes a panel-by-panel cubed sphere actually closed. Nodes are
# bucketed on a uniform spatial hash of cell size `tol`; each node is compared
# against the 27 neighbouring buckets so that a pair straddling a bucket
# boundary is still found.
#
# Returns (coords_merged, old2new, nmerged).
#---------------------------------------------------------------------------------
function merge_coincident_nodes(coords::Array{TF,2}, tol::TF) where {TF<:AbstractFloat}

    npoin   = size(coords, 1)
    old2new = zeros(Int64, npoin)
    buckets = Dict{NTuple{3,Int64}, Vector{Int64}}()
    keep    = Int64[]
    tol2    = tol*tol

    for ip = 1:npoin
        p  = (coords[ip,1], coords[ip,2], coords[ip,3])
        c  = (floor(Int64, p[1]/tol), floor(Int64, p[2]/tol), floor(Int64, p[3]/tol))
        found = 0
        for di = -1:1, dj = -1:1, dk = -1:1
            key = (c[1]+di, c[2]+dj, c[3]+dk)
            haskey(buckets, key) || continue
            for jp in buckets[key]
                d2 = (coords[jp,1]-p[1])^2 + (coords[jp,2]-p[2])^2 + (coords[jp,3]-p[3])^2
                if d2 <= tol2
                    found = old2new[jp]
                    break
                end
            end
            found != 0 && break
        end

        if found == 0
            push!(keep, ip)
            old2new[ip] = length(keep)
            push!(get!(buckets, c, Int64[]), ip)
        else
            old2new[ip] = found
        end
    end

    nmerged = npoin - length(keep)
    return coords[keep, :], old2new, nmerged
end


#---------------------------------------------------------------------------------
# Unique-edge table of a quad surface grid.
#
# An edge is keyed on the SORTED pair of its endpoints, so the two elements
# that share it (always two on a closed shell) land on the same entry. The
# stored direction conn_unique_edges[ie,:] = (a,b) is the canonical one and
# is the direction in which the high-order edge nodes are laid down; the
# element loop reverses its traversal when it walks the edge the other way.
#
# Returns (conn_unique_edges, edge_in_elem, edge_count) where edge_count[ie]
# is how many elements touch edge ie (must be 2 everywhere on a closed shell).
#---------------------------------------------------------------------------------
function build_unique_edges(conn::Array{TI,2}) where {TI<:Integer}

    nelem     = size(conn, 1)
    edge_id   = Dict{Tuple{TI,TI}, TI}()
    e1        = TI[]
    e2        = TI[]
    counts    = TI[]
    edge_in_elem = zeros(TI, nelem, 4)

    #   local edge ordering: 1:(v1,v2)  2:(v2,v3)  3:(v4,v3)  4:(v1,v4)
    local_edges = ((1,2), (2,3), (4,3), (1,4))

    for iel = 1:nelem
        for (le, (la, lb)) in enumerate(local_edges)
            a = conn[iel, la]
            b = conn[iel, lb]
            key = a < b ? (a, b) : (b, a)
            ie  = get(edge_id, key, TI(0))
            if ie == 0
                push!(e1, a); push!(e2, b); push!(counts, TI(1))
                ie = TI(length(e1))
                edge_id[key] = ie
            else
                counts[ie] += TI(1)
            end
            edge_in_elem[iel, le] = ie
        end
    end

    nedges = length(e1)
    conn_unique_edges = Array{TI}(undef, nedges, 2)
    for ie = 1:nedges
        conn_unique_edges[ie, 1] = e1[ie]
        conn_unique_edges[ie, 2] = e2[ie]
    end

    return conn_unique_edges, edge_in_elem, counts
end


#---------------------------------------------------------------------------------
# Re-orient every quad so its normal points OUTWARD (away from the centre).
#
# gmsh does not guarantee a globally consistent orientation on a closed
# surface built from several panels. An outward-consistent orientation is
# what makes (i,j) → (ξ,η) → outward normal a right-handed frame in EVERY
# element, which is what the metric terms and the tangent-plane velocity
# basis of the shallow-water equations will need.
#
# The quad normal is taken from the diagonals, n = (P3-P1) × (P4-P2), which is
# the exact normal for a planar quad and robust for a slightly warped one.
# Reversing means swapping v2 ↔ v4 (keeps v1 as the origin of the element).
#
# Returns the number of elements that had to be flipped.
#---------------------------------------------------------------------------------
function orient_elements_outward!(conn::Array{TI,2}, coords::Array{TF,2}) where {TI<:Integer, TF<:AbstractFloat}

    nflip = 0
    for iel = 1:size(conn, 1)
        i1, i2, i3, i4 = conn[iel,1], conn[iel,2], conn[iel,3], conn[iel,4]

        d1 = (coords[i3,1]-coords[i1,1], coords[i3,2]-coords[i1,2], coords[i3,3]-coords[i1,3])
        d2 = (coords[i4,1]-coords[i2,1], coords[i4,2]-coords[i2,2], coords[i4,3]-coords[i2,3])

        nx = d1[2]*d2[3] - d1[3]*d2[2]
        ny = d1[3]*d2[1] - d1[1]*d2[3]
        nz = d1[1]*d2[2] - d1[2]*d2[1]

        cx = (coords[i1,1]+coords[i2,1]+coords[i3,1]+coords[i4,1])/4
        cy = (coords[i1,2]+coords[i2,2]+coords[i3,2]+coords[i4,2])/4
        cz = (coords[i1,3]+coords[i2,3]+coords[i3,3]+coords[i4,3])/4

        if nx*cx + ny*cy + nz*cz < 0
            conn[iel,2], conn[iel,4] = i4, i2
            nflip += 1
        end
    end

    return nflip
end


#---------------------------------------------------------------------------------
# Great-circle (slerp) interpolation between two points of the SAME radius.
# t ∈ [0,1]. Falls back to the chord when the two points nearly coincide.
#---------------------------------------------------------------------------------
@inline function _slerp(p1::NTuple{3,TF}, p2::NTuple{3,TF}, t::TF, R::TF) where {TF<:AbstractFloat}

    d  = (p1[1]*p2[1] + p1[2]*p2[2] + p1[3]*p2[3]) / (R*R)
    d  = d >  1 ?  one(TF) : d
    d  = d < -1 ? -one(TF) : d
    ω  = acos(d)
    s  = sin(ω)

    if s < TF(1.0e-12)
        a = one(TF) - t
        b = t
    else
        a = sin((one(TF) - t)*ω) / s
        b = sin(t*ω) / s
    end

    return (a*p1[1] + b*p2[1], a*p1[2] + b*p2[2], a*p1[3] + b*p2[3])
end


@inline function _project_to_sphere(p::NTuple{3,TF}, R::TF) where {TF<:AbstractFloat}
    n = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3])
    n = n < TF(1.0e-300) ? one(TF) : n
    return (R*p[1]/n, R*p[2]/n, R*p[3]/n)
end


#---------------------------------------------------------------------------------
# MAIN ENTRY POINT
#
#   build_sphere_shell_mesh(gmsh_filename, nop; …) -> St_mesh_sphere
#
# Reads the linear quad shell grid, cleans it up (merge + orient), builds the
# unique-edge table, and populates the shell with the high-order LGL points:
#
#   * corners      : the (projected) gmsh vertices
#   * edge nodes   : great-circle (slerp) interpolation at the LGL abscissae,
#                    computed ONCE per unique edge and shared by both
#                    neighbouring elements
#   * interiors    : a Coons (transfinite) patch built from the four ALREADY
#                    COMPUTED great-circle edges, then projected radially onto
#                    the sphere
#
# Because the edge nodes are shared rather than recomputed, the grid is
# watertight BY CONSTRUCTION: continuity across a seam is exact node identity,
# not a floating-point coincidence.
#---------------------------------------------------------------------------------
function build_sphere_shell_mesh(gmsh_filename::String,
                                 nop::Int;
                                 interpolation_nodes = LGL(),
                                 radius              = nothing,
                                 merge_tol_rel       = 1.0e-8,
                                 lmerge_nodes        = true,
                                 lproject_to_sphere  = true,
                                 backend             = CPU(),
                                 verbose             = true,
                                 TF                  = TFloat,
                                 TI                  = TInt)

    ngl = nop + 1
    ngl >= 2 || error(" # ERROR sphere_mesh.jl: :nop must be >= 1")

    verbose && println(" # ")
    verbose && println(" # SPHERICAL SHELL GRID ........................................")
    verbose && println(" #   gmsh file: ", gmsh_filename)

    #-----------------------------------------------------------------------------
    # 1. Read the linear quad surface grid
    #-----------------------------------------------------------------------------
    coords, conn, elem_tag = read_gmsh_quad_surface(gmsh_filename, TF, TI)
    nelem = size(conn, 1)
    verbose && println(" #   linear grid as read: ", size(coords,1), " nodes, ", nelem, " quads")

    #-----------------------------------------------------------------------------
    # 2. Shell radius + optional radial projection.
    #
    #    The radius is the mean |x| over the linear nodes. A large spread
    #    means the file is not a spherical shell (or is a thick shell), which
    #    this preliminary kernel does not handle: say so instead of silently
    #    flattening it onto a sphere.
    #-----------------------------------------------------------------------------
    rnode = [sqrt(coords[i,1]^2 + coords[i,2]^2 + coords[i,3]^2) for i = 1:size(coords,1)]
    R     = radius === nothing ? TF(sum(rnode)/length(rnode)) : TF(radius)
    rdev  = maximum(abs.(rnode .- R)) / R

    verbose && @printf(" #   shell radius R = %.8e  (max relative deviation of the linear nodes: %.3e)\n", R, rdev)
    if rdev > 0.05
        @warn string(" # sphere_mesh.jl: the gmsh nodes deviate from a sphere of radius ", R,
                     " by ", rdev*100, "%. Is ", gmsh_filename, " really a spherical shell?")
    end

    if lproject_to_sphere
        for i = 1:size(coords,1)
            p = _project_to_sphere((coords[i,1], coords[i,2], coords[i,3]), R)
            coords[i,1], coords[i,2], coords[i,3] = p[1], p[2], p[3]
        end
    end

    #-----------------------------------------------------------------------------
    # 3. Merge coincident vertices (panel seams!) — see the header note.
    #-----------------------------------------------------------------------------
    if lmerge_nodes
        tol = TF(merge_tol_rel) * R
        coords_m, old2new, nmerged = merge_coincident_nodes(coords, tol)
        if nmerged > 0
            verbose && println(" #   merged ", nmerged, " duplicated (panel-seam) vertices — the shell was NOT watertight as read")
            for iel = 1:nelem, k = 1:4
                conn[iel,k] = old2new[conn[iel,k]]
            end
            coords = coords_m
        else
            verbose && println(" #   no duplicated vertices found: the gmsh shell is already watertight")
        end
    end
    npoin_linear = size(coords, 1)

    #-----------------------------------------------------------------------------
    # 4. Consistent OUTWARD orientation of every quad
    #-----------------------------------------------------------------------------
    nflip = orient_elements_outward!(conn, coords)
    verbose && println(" #   re-oriented ", nflip, " of ", nelem, " quads to an outward-pointing normal")

    #-----------------------------------------------------------------------------
    # 5. Unique edges of the closed shell
    #-----------------------------------------------------------------------------
    conn_unique_edges, edge_in_elem, edge_count = build_unique_edges(conn)
    nedges = size(conn_unique_edges, 1)

    nbdy_edges = count(c -> c == 1, edge_count)
    nbad_edges = count(c -> c > 2, edge_count)
    verbose && println(" #   unique edges: ", nedges,
                       "  (edges used by 1 element: ", nbdy_edges,
                       ", by >2 elements: ", nbad_edges, ")")
    if nbad_edges > 0
        error(string(" # ERROR sphere_mesh.jl: ", nbad_edges,
                     " edges are shared by more than two elements: the grid is not a manifold."))
    end
    if nbdy_edges > 0
        @warn string(" # sphere_mesh.jl: ", nbdy_edges, " edges belong to a single element. ",
                     "The surface is NOT closed — either the grid is a spherical PATCH, or panel-seam ",
                     "nodes were not merged (try a larger :node_merge_tol).")
    end

    #-----------------------------------------------------------------------------
    # 6. LGL points
    #-----------------------------------------------------------------------------
    lgl = basis_structs_ξ_ω!(interpolation_nodes, nop, backend)
    ξ   = TF.(collect(lgl.ξ))

    #-----------------------------------------------------------------------------
    # 7. Global node numbering + coordinates
    #-----------------------------------------------------------------------------
    n_edge_int = nedges * (ngl - 2)
    n_elem_int = nelem  * (ngl - 2) * (ngl - 2)
    npoin      = npoin_linear + n_edge_int + n_elem_int

    x = zeros(TF, npoin); y = zeros(TF, npoin); z = zeros(TF, npoin)
    node_type = zeros(TI, npoin)

    @inbounds for ip = 1:npoin_linear
        x[ip], y[ip], z[ip] = coords[ip,1], coords[ip,2], coords[ip,3]
        node_type[ip] = TI(0)
    end

    connijk      = zeros(TI, nelem, ngl, ngl)
    poin_in_edge = zeros(TI, nedges, ngl)

    #--- 7a. edge nodes: created ONCE per unique edge, in the canonical direction
    for ie = 1:nedges
        ip1 = conn_unique_edges[ie, 1]
        ip2 = conn_unique_edges[ie, 2]

        poin_in_edge[ie, 1]   = ip1
        poin_in_edge[ie, ngl] = ip2

        p1 = (x[ip1], y[ip1], z[ip1])
        p2 = (x[ip2], y[ip2], z[ip2])

        for l = 2:ngl-1
            ip = npoin_linear + (ie - 1)*(ngl - 2) + (l - 1)
            t  = (one(TF) + ξ[l]) / 2
            q  = _slerp(p1, p2, t, R)
            x[ip], y[ip], z[ip] = q[1], q[2], q[3]
            poin_in_edge[ie, l] = TI(ip)
            node_type[ip]       = TI(1)
        end
    end

    #--- 7b. element corners and edges → connijk (direction-aware)
    #
    #    local edge → (i,j) walk, in the direction of the local endpoint pair:
    #      le=1 (v1→v2): (l, 1)
    #      le=2 (v2→v3): (ngl, l)
    #      le=3 (v4→v3): (l, ngl)
    #      le=4 (v1→v4): (1, l)
    #    If the element walks the edge opposite to the canonical direction,
    #    the LGL index is mirrored (l → ngl+1-l). THIS is what keeps two
    #    neighbouring elements pointing at the SAME node on a shared seam.
    #
    local_edges = ((1,2), (2,3), (4,3), (1,4))
    for iel = 1:nelem
        v1, v2, v3, v4 = conn[iel,1], conn[iel,2], conn[iel,3], conn[iel,4]

        connijk[iel, 1,   1  ] = v1
        connijk[iel, ngl, 1  ] = v2
        connijk[iel, ngl, ngl] = v3
        connijk[iel, 1,   ngl] = v4

        for (le, (la, lb)) in enumerate(local_edges)
            ie      = edge_in_elem[iel, le]
            forward = (conn_unique_edges[ie, 1] == conn[iel, la])
            for l = 2:ngl-1
                lcan = forward ? l : ngl + 1 - l
                ip   = poin_in_edge[ie, lcan]
                if le == 1
                    connijk[iel, l,   1  ] = ip
                elseif le == 2
                    connijk[iel, ngl, l  ] = ip
                elseif le == 3
                    connijk[iel, l,   ngl] = ip
                else
                    connijk[iel, 1,   l  ] = ip
                end
            end
        end
    end

    #--- 7c. element interiors: Coons patch on the four great-circle edges,
    #        then radial projection back onto the shell.
    off_int = npoin_linear + n_edge_int
    for iel = 1:nelem
        p11 = (x[connijk[iel,1,1]],     y[connijk[iel,1,1]],     z[connijk[iel,1,1]])
        pn1 = (x[connijk[iel,ngl,1]],   y[connijk[iel,ngl,1]],   z[connijk[iel,ngl,1]])
        pnn = (x[connijk[iel,ngl,ngl]], y[connijk[iel,ngl,ngl]], z[connijk[iel,ngl,ngl]])
        p1n = (x[connijk[iel,1,ngl]],   y[connijk[iel,1,ngl]],   z[connijk[iel,1,ngl]])

        k = 0
        for j = 2:ngl-1
            v  = (one(TF) + ξ[j]) / 2
            ipL = connijk[iel, 1,   j]
            ipR = connijk[iel, ngl, j]
            for i = 2:ngl-1
                u   = (one(TF) + ξ[i]) / 2
                ipB = connijk[iel, i, 1]
                ipT = connijk[iel, i, ngl]

                cx = (one(TF)-v)*x[ipB] + v*x[ipT] + (one(TF)-u)*x[ipL] + u*x[ipR] -
                     ((one(TF)-u)*(one(TF)-v)*p11[1] + u*(one(TF)-v)*pn1[1] + u*v*pnn[1] + (one(TF)-u)*v*p1n[1])
                cy = (one(TF)-v)*y[ipB] + v*y[ipT] + (one(TF)-u)*y[ipL] + u*y[ipR] -
                     ((one(TF)-u)*(one(TF)-v)*p11[2] + u*(one(TF)-v)*pn1[2] + u*v*pnn[2] + (one(TF)-u)*v*p1n[2])
                cz = (one(TF)-v)*z[ipB] + v*z[ipT] + (one(TF)-u)*z[ipL] + u*z[ipR] -
                     ((one(TF)-u)*(one(TF)-v)*p11[3] + u*(one(TF)-v)*pn1[3] + u*v*pnn[3] + (one(TF)-u)*v*p1n[3])

                q  = _project_to_sphere((cx, cy, cz), R)

                ip = off_int + (iel - 1)*(ngl-2)*(ngl-2) + k + 1
                x[ip], y[ip], z[ip] = q[1], q[2], q[3]
                connijk[iel, i, j]  = TI(ip)
                node_type[ip]       = TI(2)
                k += 1
            end
        end
    end

    #-----------------------------------------------------------------------------
    # 8. Spherical coordinates of every node (handy for the ICs to come)
    #-----------------------------------------------------------------------------
    lon = zeros(TF, npoin)
    lat = zeros(TF, npoin)
    @inbounds for ip = 1:npoin
        lon[ip] = atan(y[ip], x[ip])
        lat[ip] = asin(clamp(z[ip]/R, -one(TF), one(TF)))
    end

    mesh = St_mesh_sphere{TF, TI}(TI(nop), TI(ngl), R,
                                  TI(nelem), TI(npoin), TI(npoin_linear), TI(nedges),
                                  x, y, z, lon, lat,
                                  conn, connijk, conn_unique_edges, poin_in_edge, edge_in_elem,
                                  elem_tag, node_type, NSD_2D())

    verbose && println(" #   high-order grid: nop = ", nop, ", ngl = ", ngl,
                       ", npoin = ", npoin, " (", npoin_linear, " vertices + ",
                       n_edge_int, " edge + ", n_elem_int, " interior)")
    verbose && println(" # SPHERICAL SHELL GRID ........................................ DONE")

    return mesh
end


#---------------------------------------------------------------------------------
# check_sphere_mesh(mesh; …)
#
# Verification of everything that can go wrong on a CLOSED shell. Runs after
# the grid is built and prints a PASS/FAIL line per test. Returns true when
# all tests pass.
#
#   T1  edge multiplicity ....... every unique edge is used by exactly 2 quads
#   T2  Euler characteristic .... V - E + F = 2  (genus-0 closed surface)
#   T3  orientation ............. every unique edge traversed once in each
#                                 direction  ⇒ globally consistent orientation
#   T4  connectivity ............ the element graph is a single component
#   T5  node accounting ......... every node index in 1:npoin is used, and the
#                                 count matches npoin_linear + nedges*(ngl-2)
#                                 + nelem*(ngl-2)^2
#   T6  no duplicate nodes ...... no two DISTINCT node indices sit at the same
#                                 point (this is the test that catches a seam
#                                 that was numbered twice)
#   T7  on-sphere ............... every node lies on the shell to round-off
#   T8  positive Jacobian ....... every sub-cell of every element has an
#                                 outward-consistent orientation (no folds)
#---------------------------------------------------------------------------------
function check_sphere_mesh(mesh::St_mesh_sphere; verbose = true, tol_rel = 1.0e-8)

    allok = true
    line(name, ok, extra="") = begin
        verbose && println("     ", ok ? "PASS" : "FAIL", "  ", name, extra == "" ? "" : string("  [", extra, "]"))
        ok
    end

    verbose && println(" # CHECK SPHERICAL SHELL GRID ..................................")

    ngl    = mesh.ngl
    nelem  = mesh.nelem
    nedges = mesh.nedges

    #--- T1 edge multiplicity
    used = zeros(Int, nedges)
    for iel = 1:nelem, le = 1:4
        used[mesh.edge_in_elem[iel, le]] += 1
    end
    nbad = count(c -> c != 2, used)
    allok &= line("T1 edge multiplicity (all edges shared by 2 elements)", nbad == 0,
                  nbad == 0 ? "" : string(nbad, " edges with count != 2"))

    #--- T2 Euler characteristic
    χ = mesh.npoin_linear - nedges + nelem
    allok &= line("T2 Euler characteristic V-E+F = 2", χ == 2, string("χ = ", χ))

    #--- T3 orientation consistency
    fwd = zeros(Int, nedges)
    bwd = zeros(Int, nedges)
    local_edges = ((1,2), (2,3), (4,3), (1,4))
    for iel = 1:nelem
        for (le, (la, lb)) in enumerate(local_edges)
            ie = mesh.edge_in_elem[iel, le]
            #
            # local edges 1 and 2 run "with" the CCW loop v1→v2→v3→v4, local
            # edges 3 and 4 are stored reversed (v4→v3, v1→v4), so the CCW
            # traversal direction is the stored one flipped.
            #
            a = mesh.conn[iel, la]
            ccw_start = (le <= 2) ? a : mesh.conn[iel, lb]
            if mesh.conn_unique_edges[ie, 1] == ccw_start
                fwd[ie] += 1
            else
                bwd[ie] += 1
            end
        end
    end
    nbad = count(ie -> !(fwd[ie] == 1 && bwd[ie] == 1), 1:nedges)
    allok &= line("T3 globally consistent orientation", nbad == 0,
                  nbad == 0 ? "" : string(nbad, " inconsistently oriented edges"))

    #--- T4 single connected component
    edge2elem = [Int[] for _ = 1:nedges]
    for iel = 1:nelem, le = 1:4
        push!(edge2elem[mesh.edge_in_elem[iel, le]], iel)
    end
    seen  = falses(nelem)
    stack = [1]
    seen[1] = true
    ncomp_elems = 0
    while !isempty(stack)
        iel = pop!(stack)
        ncomp_elems += 1
        for le = 1:4, jel in edge2elem[mesh.edge_in_elem[iel, le]]
            if !seen[jel]
                seen[jel] = true
                push!(stack, jel)
            end
        end
    end
    allok &= line("T4 single connected shell", ncomp_elems == nelem,
                  string(ncomp_elems, "/", nelem, " elements reachable"))

    #--- T5 node accounting
    touched = falses(mesh.npoin)
    inrange = true
    for iel = 1:nelem, j = 1:ngl, i = 1:ngl
        ip = mesh.connijk[iel, i, j]
        if ip < 1 || ip > mesh.npoin
            inrange = false
        else
            touched[ip] = true
        end
    end
    expected = mesh.npoin_linear + nedges*(ngl-2) + nelem*(ngl-2)^2
    ok5 = inrange && all(touched) && (mesh.npoin == expected)
    allok &= line("T5 node accounting (all of 1:npoin used, count exact)", ok5,
                  string("npoin = ", mesh.npoin, ", expected = ", expected,
                         ", unused = ", count(!, touched)))

    #--- T6 geometric duplicates on the FINAL high-order grid
    tol  = tol_rel * mesh.radius
    tol2 = tol*tol
    buckets = Dict{NTuple{3,Int64}, Vector{Int64}}()
    ndup = 0
    for ip = 1:mesh.npoin
        c = (floor(Int64, mesh.x[ip]/tol), floor(Int64, mesh.y[ip]/tol), floor(Int64, mesh.z[ip]/tol))
        for di = -1:1, dj = -1:1, dk = -1:1
            key = (c[1]+di, c[2]+dj, c[3]+dk)
            haskey(buckets, key) || continue
            for jp in buckets[key]
                d2 = (mesh.x[jp]-mesh.x[ip])^2 + (mesh.y[jp]-mesh.y[ip])^2 + (mesh.z[jp]-mesh.z[ip])^2
                d2 <= tol2 && (ndup += 1)
            end
        end
        push!(get!(buckets, c, Int64[]), ip)
    end
    allok &= line("T6 no duplicated nodes on the high-order grid", ndup == 0,
                  ndup == 0 ? "" : string(ndup, " coincident node pairs"))

    #--- T7 all nodes on the shell
    dmax = 0.0
    for ip = 1:mesh.npoin
        r = sqrt(mesh.x[ip]^2 + mesh.y[ip]^2 + mesh.z[ip]^2)
        dmax = max(dmax, abs(r - mesh.radius)/mesh.radius)
    end
    allok &= line("T7 every node on the sphere", dmax < 1.0e-10, @sprintf("max |r-R|/R = %.3e", dmax))

    #--- T8 no folded sub-cells (outward-consistent sub-quad normals)
    nneg = 0
    for iel = 1:nelem, j = 1:ngl-1, i = 1:ngl-1
        a = mesh.connijk[iel, i,   j  ]
        b = mesh.connijk[iel, i+1, j  ]
        c = mesh.connijk[iel, i+1, j+1]
        d = mesh.connijk[iel, i,   j+1]

        d1 = (mesh.x[c]-mesh.x[a], mesh.y[c]-mesh.y[a], mesh.z[c]-mesh.z[a])
        d2 = (mesh.x[d]-mesh.x[b], mesh.y[d]-mesh.y[b], mesh.z[d]-mesh.z[b])
        nx = d1[2]*d2[3] - d1[3]*d2[2]
        ny = d1[3]*d2[1] - d1[1]*d2[3]
        nz = d1[1]*d2[2] - d1[2]*d2[1]
        cx = (mesh.x[a]+mesh.x[b]+mesh.x[c]+mesh.x[d])/4
        cy = (mesh.y[a]+mesh.y[b]+mesh.y[c]+mesh.y[d])/4
        cz = (mesh.z[a]+mesh.z[b]+mesh.z[c]+mesh.z[d])/4
        (nx*cx + ny*cy + nz*cz <= 0) && (nneg += 1)
    end
    allok &= line("T8 no folded/inverted sub-cells", nneg == 0,
                  nneg == 0 ? "" : string(nneg, " sub-cells with inward normal"))

    #--- T9 the decisive one: two elements sharing an edge must list the SAME
    #        ngl node indices along it (same order or reversed). This is C0
    #        continuity expressed as node IDENTITY — if a panel seam were
    #        numbered twice, T9 fails even when every count above is right.
    edge_sides = [Tuple{Int,Int}[] for _ = 1:nedges]
    for iel = 1:nelem, le = 1:4
        push!(edge_sides[mesh.edge_in_elem[iel, le]], (iel, le))
    end
    nbad = 0
    edge_nodes(iel, le) = le == 1 ? [mesh.connijk[iel, l, 1]   for l = 1:ngl] :
                          le == 2 ? [mesh.connijk[iel, ngl, l] for l = 1:ngl] :
                          le == 3 ? [mesh.connijk[iel, l, ngl] for l = 1:ngl] :
                                    [mesh.connijk[iel, 1, l]   for l = 1:ngl]
    for ie = 1:nedges
        length(edge_sides[ie]) == 2 || continue
        (ia, la), (ib, lb) = edge_sides[ie][1], edge_sides[ie][2]
        na, nb = edge_nodes(ia, la), edge_nodes(ib, lb)
        (na == nb || na == reverse(nb)) || (nbad += 1)
    end
    allok &= line("T9 neighbouring elements share the SAME edge nodes", nbad == 0,
                  nbad == 0 ? "" : string(nbad, " torn seams"))

    verbose && println(" # CHECK SPHERICAL SHELL GRID .................................. ",
                       allok ? "ALL TESTS PASSED" : "FAILED")

    return allok
end


#---------------------------------------------------------------------------------
# Surface area of the shell obtained by summing the sub-cell areas of the
# high-order grid. Compared against 4πR² it is a global sanity check on the
# element mapping (and a first, crude, convergence indicator for `nop`).
#---------------------------------------------------------------------------------
function sphere_mesh_area(mesh::St_mesh_sphere)
    A = 0.0
    for iel = 1:mesh.nelem, j = 1:mesh.ngl-1, i = 1:mesh.ngl-1
        a = mesh.connijk[iel, i,   j  ]
        b = mesh.connijk[iel, i+1, j  ]
        c = mesh.connijk[iel, i+1, j+1]
        d = mesh.connijk[iel, i,   j+1]
        for (p, q, r) in ((a,b,c), (a,c,d))
            ux = mesh.x[q]-mesh.x[p]; uy = mesh.y[q]-mesh.y[p]; uz = mesh.z[q]-mesh.z[p]
            vx = mesh.x[r]-mesh.x[p]; vy = mesh.y[r]-mesh.y[p]; vz = mesh.z[r]-mesh.z[p]
            cx = uy*vz - uz*vy; cy = uz*vx - ux*vz; cz = ux*vy - uy*vx
            A += 0.5*sqrt(cx*cx + cy*cy + cz*cz)
        end
    end
    return A
end


#---------------------------------------------------------------------------------
# Driver called from problems/drivers.jl when :lspherical_shell => true.
#---------------------------------------------------------------------------------
function mod_mesh_sphere_driver(inputs, TFloat_in=TFloat)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    if MPI.Comm_size(comm) > 1
        error(" # ERROR sphere_mesh.jl: the spherical-shell grid builder is serial for now. Run it on a single MPI rank.")
    end

    haskey(inputs, :gmsh_filename) ||
        error(" # ERROR sphere_mesh.jl: :lspherical_shell => true requires :gmsh_filename in user_inputs.jl")

    mesh = build_sphere_shell_mesh(string(inputs[:gmsh_filename]),
                                   Int(inputs[:nop]);
                                   interpolation_nodes = inputs[:interpolation_nodes],
                                   radius              = get(inputs, :sphere_radius,  nothing),
                                   merge_tol_rel       = get(inputs, :node_merge_tol, 1.0e-8),
                                   lmerge_nodes        = get(inputs, :lmerge_coincident_nodes, true),
                                   lproject_to_sphere  = get(inputs, :lproject_to_sphere, true),
                                   backend             = inputs[:backend],
                                   verbose             = (rank == 0),
                                   TF                  = TFloat_in,
                                   TI                  = TInt)

    if get(inputs, :lcheck_grid, true) == true
        ok = check_sphere_mesh(mesh; verbose = (rank == 0))
        if !ok && get(inputs, :lstop_on_bad_grid, true) == true
            error(" # ERROR sphere_mesh.jl: the spherical shell grid failed its consistency checks (see the FAIL lines above).")
        end
    end

    if rank == 0
        A = sphere_mesh_area(mesh)
        @printf(" #   shell area (sub-cell sum) = %.10e ; 4πR² = %.10e ; relative error = %.3e\n",
                A, 4π*mesh.radius^2, abs(A - 4π*mesh.radius^2)/(4π*mesh.radius^2))
    end

    return mesh
end

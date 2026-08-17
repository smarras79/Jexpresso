#!/usr/bin/env julia
#---------------------------------------------------------------------------------
# tools/generate_cubed_sphere.jl — write a closed, all-quadrilateral cubed-sphere
# surface grid directly, without gmsh.
#
#   julia tools/generate_cubed_sphere.jl [n] [R] [out.msh] [map]
#
#   n     elements per panel edge        (default 10  -> 6n² = 600 quads)
#   R     sphere radius in metres        (default 6.371e6, Earth)
#   out   output file                    (default ./cubed_sphere.msh)
#   map   gnomonic | equiangular         (default equiangular)
#
#   julia tools/generate_cubed_sphere.jl 10 6.371e6 cubed_sphere.msh equiangular
#
# Runs in a bare Julia — no `using Jexpresso`, no gmsh, no packages. It includes
# src/kernel/mesh/cubed_sphere_maps.jl for the cube-face → sphere maps, which is
# dependency-free on purpose, so the map used here is the SAME tested code as
# :cubed_sphere_map applies at read time.
#
# WHY THIS EXISTS: no sphere-centre node.
#
# cubed_sphere.geo builds the panels as `Ruled Surface ... In Sphere {1}`, and
# that construction needs `Point(1) = {0,0,0}` as the sphere's centre reference.
# With gmsh's "Save all elements" the centre is emitted as a node of the mesh —
# a node at the origin that belongs to no element, sits at r = 0, and is
# meaningless to every downstream consumer. It also perturbs the node count that
# test/test_sphere_mesh.jl checks. Generating the grid arithmetically removes
# the whole question: every node produced here lies ON the sphere and belongs to
# at least one quad.
#
# WATERTIGHTNESS IS EXACT, NOT TOLERANCED.
#
# The usual way to build a cubed sphere is per panel, then merge coincident
# nodes with a distance tolerance. That is fragile: pick the tolerance too tight
# and the seams stay open (six disconnected patches that still LOOK like a
# sphere), too loose and distinct nodes get welded.
#
# Instead every node is identified by an INTEGER triple on the surface of the
# cube [-n, n]³: panel p, grid index (i, j) maps to
#
#     T = (2i-n)·êu + (2j-n)·êv + n·êw
#
# with (êu, êv, êw) the panel's frame from cubed_sphere_maps.jl. Two panels
# meeting along a cube edge generate the IDENTICAL triple for the same physical
# node — integers, so the dictionary lookup is exact and no tolerance appears
# anywhere. The node is positioned once, by whichever panel reaches it first;
# the seam-consistency property proved in test/test_cubed_sphere_maps.jl is what
# makes "whichever" safe.
#
# The counts that follow, for free: V = 6n²+2, F = 6n², E = 12n², so
# V - E + F = 2 exactly. The script asserts all of it before writing.
#
# NOT OFFERED: :conformal. It is singular at the eight cube corners — the local
# map scale goes as d^(1/3) — and a structured n×n panel grid necessarily puts a
# node on every corner, so the surface Jacobian there is zero and
# build_sphere_metrics rejects the element. See the note in
# remap_cubed_sphere_nodes! (src/kernel/mesh/mesh.jl).
#---------------------------------------------------------------------------------

# Printf first: @printf expands at PARSE time, so the `using` has to precede
# every macro use below, not merely every call.
using Printf

include(joinpath(@__DIR__, "..", "src", "kernel", "mesh", "cubed_sphere_maps.jl"))

const SUPPORTED_MAPS = (:gnomonic, :equidistant, :equiangular)

"""
    build_cubed_sphere(n, R, map) -> (nodes, quads)

`nodes` is a Vector of (x,y,z) on the sphere of radius `R`; `quads` is a Vector
of (panel, n1, n2, n3, n4) with 1-based node indices wound counter-clockwise as
seen from OUTSIDE the sphere.
"""
function build_cubed_sphere(n::Int, R::Float64, mapname::Symbol)
    n >= 1 || error("n must be >= 1, got $n")
    R > 0   || error("R must be > 0, got $R")
    mapname in SUPPORTED_MAPS ||
        error("map must be one of $SUPPORTED_MAPS, got :$mapname" *
              (mapname === :conformal ?
               "\n  :conformal is singular at the eight cube corners and a structured" *
               "\n  panel grid puts a node on every one of them — see the file header." : ""))

    idof   = Dict{NTuple{3,Int},Int}()          # exact integer identity -> node id
    nodes  = NTuple{3,Float64}[]
    quads  = NTuple{5,Int}[]

    # id of the node at panel-p grid index (i, j), creating it on first sight
    function nodeid(p::Int, i::Int, j::Int)
        eu, ev, ew = _PANEL_FRAMES[p]
        si, sj = 2i - n, 2j - n
        T = (round(Int, si*eu[1] + sj*ev[1] + n*ew[1]),
             round(Int, si*eu[2] + sj*ev[2] + n*ew[2]),
             round(Int, si*eu[3] + sj*ev[3] + n*ew[3]))
        haskey(idof, T) && return idof[T]
        u, v = si/n, sj/n
        U, V, W = _map_forward(mapname, u, v)
        x, y, z = _panel_to_cartesian(U, V, W, p)
        push!(nodes, (R*x, R*y, R*z))
        idof[T] = length(nodes)
        return idof[T]
    end

    for p = 1:6, j = 0:n-1, i = 0:n-1
        # counter-clockwise in (u,v); êu × êv = êw points OUT of the sphere, so
        # this is counter-clockwise seen from outside and the VTK normals face
        # outward without any later re-orientation.
        push!(quads, (p, nodeid(p,i,j), nodeid(p,i+1,j), nodeid(p,i+1,j+1), nodeid(p,i,j+1)))
    end
    return nodes, quads
end

"""
    verify(nodes, quads, n, R)

Everything that must hold for the grid to be usable, checked before it is
written rather than discovered inside a solver.
"""
function verify(nodes, quads, n::Int, R::Float64)
    nv, nf = length(nodes), length(quads)
    @assert nv == 6n^2 + 2   "node count $nv != 6n^2+2 = $(6n^2+2) — seams did not merge"
    @assert nf == 6n^2       "quad count $nf != 6n^2 = $(6n^2)"

    # every node ON the sphere, and none at the origin
    for (k, p) in enumerate(nodes)
        r = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
        @assert abs(r - R) < 1e-9*R "node $k is off the sphere (r/R = $(r/R))"
    end

    # every edge shared by EXACTLY two quads -> closed, watertight, no seam left open
    edges = Dict{Tuple{Int,Int},Int}()
    for q in quads, (a,b) in ((q[2],q[3]),(q[3],q[4]),(q[4],q[5]),(q[5],q[2]))
        k = a < b ? (a,b) : (b,a)
        edges[k] = get(edges, k, 0) + 1
    end
    bad = count(v -> v != 2, values(edges))
    @assert bad == 0 "$bad edge(s) not shared by exactly two quads — the shell is torn"
    ne = length(edges)
    @assert ne == 12n^2 "edge count $ne != 12n^2 = $(12n^2)"
    @assert nv - ne + nf == 2 "Euler characteristic $(nv-ne+nf) != 2 — not a closed sphere"

    # no two nodes coincident
    mind = Inf
    for i = 1:nv, j = i+1:nv
        d = (nodes[i][1]-nodes[j][1])^2 + (nodes[i][2]-nodes[j][2])^2 + (nodes[i][3]-nodes[j][3])^2
        d < mind && (mind = d)
    end
    mind = sqrt(mind)
    @assert mind > 1e-6*R "two nodes are coincident (separation $(mind))"

    # winding: the quad normal must point AWAY from the centre for every quad
    nin = 0
    for q in quads
        a, b, c, d = nodes[q[2]], nodes[q[3]], nodes[q[4]], nodes[q[5]]
        d1 = (c[1]-a[1], c[2]-a[2], c[3]-a[3])
        d2 = (d[1]-b[1], d[2]-b[2], d[3]-b[3])
        nrm = (d1[2]*d2[3]-d1[3]*d2[2], d1[3]*d2[1]-d1[1]*d2[3], d1[1]*d2[2]-d1[2]*d2[1])
        ctr = ((a[1]+b[1]+c[1]+d[1])/4, (a[2]+b[2]+c[2]+d[2])/4, (a[3]+b[3]+c[3]+d[3])/4)
        sum(nrm .* ctr) > 0 || (nin += 1)
    end
    @assert nin == 0 "$nin quad(s) wound inward"

    return (nv = nv, ne = ne, nf = nf, mind = mind)
end

"""
    write_msh(path, nodes, quads)

MSH 2.2 ASCII, one physical surface per panel — the format and tagging of the
grid that ships with ShallowWater/SWsphere, so it drops straight in.
"""
function write_msh(path::String, nodes, quads)
    open(path, "w") do io
        println(io, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat")
        println(io, "\$PhysicalNames\n6")
        for p = 1:6; println(io, "2 $p \"panel$p\""); end
        println(io, "\$EndPhysicalNames")
        println(io, "\$Nodes\n", length(nodes))
        for (k, p) in enumerate(nodes)
            @printf(io, "%d %.17g %.17g %.17g\n", k, p[1], p[2], p[3])
        end
        println(io, "\$EndNodes")
        println(io, "\$Elements\n", length(quads))
        for (k, q) in enumerate(quads)
            # id type=3(quad) ntags=2 physical geometrical n1 n2 n3 n4
            println(io, "$k 3 2 $(q[1]) $(q[1]) $(q[2]) $(q[3]) $(q[4]) $(q[5])")
        end
        println(io, "\$EndElements")
    end
    return path
end

function main(args)
    n    = length(args) >= 1 ? parse(Int, args[1])     : 10
    R    = length(args) >= 2 ? parse(Float64, args[2]) : 6.371e6
    out  = length(args) >= 3 ? args[3]                 : "cubed_sphere.msh"
    mapn = length(args) >= 4 ? Symbol(args[4])         : :equiangular

    @printf("# cubed sphere: n = %d per panel edge, R = %.6g m, map = :%s\n", n, R, mapn)
    nodes, quads = build_cubed_sphere(n, R, mapn)
    s = verify(nodes, quads, n, R)
    @printf("#   %d nodes, %d quads, %d edges ; V-E+F = %d ; closest node pair %.1f km\n",
            s.nv, s.nf, s.ne, s.nv - s.ne + s.nf, s.mind/1000)
    println("#   no sphere-centre node; every node lies on the sphere and in >= 1 quad")
    write_msh(out, nodes, quads)
    println("#   wrote ", abspath(out))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)

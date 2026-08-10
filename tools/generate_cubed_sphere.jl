#!/usr/bin/env julia
#---------------------------------------------------------------------------------
# generate_cubed_sphere.jl
#
# Write an all-quadrilateral CUBED-SPHERE surface grid to a gmsh MSH 2.2 ASCII
# file, without needing gmsh installed.
#
# Usage
#   julia tools/generate_cubed_sphere.jl [n] [R] [outfile] [--split-panels]
#
#     n              elements per panel edge (default 10 → 6n² = 600 quads)
#     R              sphere radius in metres (default 6.371e6, Earth)
#     outfile        default ./meshes/gmsh_grids/cubed_sphere.msh
#     --split-panels write the six panels with DUPLICATED seam nodes, i.e. a
#                    grid that only LOOKS closed. Useful to exercise the
#                    node-merging path of src/kernel/mesh/sphere_mesh.jl.
#
# The panels use the EQUIANGULAR gnomonic mapping (Ronchi et al. 1996):
# α, β uniformly spaced in [-π/4, π/4] and the cube-face coordinates are
# tan α, tan β. That is the standard cubed sphere used for atmospheric
# dynamical cores — element sizes vary by only ~√2 across a panel, against ~2
# for the equidistant (tangent-plane) variant.
#
# Physical tags 1..6 identify the six panels, so the panel decomposition is
# visible in ParaView through the `panel` cell field written by
# write_vtk_sphere_grid().
#---------------------------------------------------------------------------------

using Printf

#
# Direction on the cube of face `p` at face coordinates (X, Y) = (tan α, tan β).
# The six frames are chosen so that (i, j) → outward normal is right-handed on
# every panel, i.e. the quads come out already counter-clockwise seen from
# outside. (sphere_mesh.jl re-orients anyway, but a generator that gets this
# right makes `nflip = 0` a meaningful sanity signal.)
#
function panel_dir(p::Int, X::Float64, Y::Float64)
    p == 1 && return ( 1.0,  X,   Y  )   # +x
    p == 2 && return (-X,    1.0, Y  )   # +y
    p == 3 && return (-1.0, -X,   Y  )   # -x
    p == 4 && return ( X,   -1.0, Y  )   # -y
    p == 5 && return (-Y,    X,   1.0)   # +z  (north)
    p == 6 && return ( Y,    X,  -1.0)   # -z  (south)
    error("panel index must be in 1:6")
end


#
# Merge geometrically coincident nodes (same spatial-hash algorithm as
# merge_coincident_nodes in src/kernel/mesh/sphere_mesh.jl). This is what
# stitches the six panels into ONE closed shell at generation time.
#
function merge_nodes(coords::Vector{NTuple{3,Float64}}, tol::Float64)
    buckets = Dict{NTuple{3,Int64}, Vector{Int}}()
    old2new = zeros(Int, length(coords))
    keep    = Int[]
    tol2    = tol*tol
    for ip in eachindex(coords)
        p = coords[ip]
        c = (floor(Int64, p[1]/tol), floor(Int64, p[2]/tol), floor(Int64, p[3]/tol))
        found = 0
        for di = -1:1, dj = -1:1, dk = -1:1
            key = (c[1]+di, c[2]+dj, c[3]+dk)
            haskey(buckets, key) || continue
            for jp in buckets[key]
                q = coords[jp]
                if (q[1]-p[1])^2 + (q[2]-p[2])^2 + (q[3]-p[3])^2 <= tol2
                    found = old2new[jp]; break
                end
            end
            found != 0 && break
        end
        if found == 0
            push!(keep, ip); old2new[ip] = length(keep)
            push!(get!(buckets, c, Int[]), ip)
        else
            old2new[ip] = found
        end
    end
    return coords[keep], old2new
end


function generate_cubed_sphere(n::Int, R::Float64; split_panels::Bool = false)

    α = [-π/4 + i*(π/2)/n for i = 0:n]
    t = tan.(α)

    coords = NTuple{3,Float64}[]
    quads  = NTuple{4,Int}[]
    tags   = Int[]

    for p = 1:6
        base = length(coords)
        for j = 0:n, i = 0:n
            d  = panel_dir(p, t[i+1], t[j+1])
            nn = sqrt(d[1]^2 + d[2]^2 + d[3]^2)
            push!(coords, (R*d[1]/nn, R*d[2]/nn, R*d[3]/nn))
        end
        for j = 0:n-1, i = 0:n-1
            v1 = base + j*(n+1) + i + 1
            push!(quads, (v1, v1+1, v1+(n+1)+1, v1+(n+1)))
            push!(tags, p)
        end
    end

    if !split_panels
        coords, old2new = merge_nodes(coords, 1.0e-8*R)
        quads = [(old2new[q[1]], old2new[q[2]], old2new[q[3]], old2new[q[4]]) for q in quads]
    end

    return coords, quads, tags
end


function write_msh22(fname::String, coords, quads, tags)
    mkpath(dirname(fname))
    open(fname, "w") do f
        println(f, "\$MeshFormat")
        println(f, "2.2 0 8")
        println(f, "\$EndMeshFormat")
        println(f, "\$PhysicalNames")
        println(f, 6)
        for p = 1:6
            println(f, "2 ", p, " \"panel", p, "\"")
        end
        println(f, "\$EndPhysicalNames")
        println(f, "\$Nodes")
        println(f, length(coords))
        for (k, c) in enumerate(coords)
            @printf(f, "%d %.17g %.17g %.17g\n", k, c[1], c[2], c[3])
        end
        println(f, "\$EndNodes")
        println(f, "\$Elements")
        println(f, length(quads))
        for (k, q) in enumerate(quads)
            @printf(f, "%d 3 2 %d %d %d %d %d %d\n", k, tags[k], tags[k], q[1], q[2], q[3], q[4])
        end
        println(f, "\$EndElements")
    end
end


if abspath(PROGRAM_FILE) == @__FILE__

    args   = filter(a -> !startswith(a, "--"), ARGS)
    split  = "--split-panels" in ARGS

    n      = length(args) >= 1 ? parse(Int,     args[1]) : 10
    R      = length(args) >= 2 ? parse(Float64, args[2]) : 6.371e6
    out    = length(args) >= 3 ? args[3] : joinpath(@__DIR__, "..", "meshes", "gmsh_grids", "cubed_sphere.msh")

    coords, quads, tags = generate_cubed_sphere(n, R; split_panels = split)
    write_msh22(out, coords, quads, tags)

    #  6n² quads, 6n² + 2 vertices on the merged (closed) shell
    @printf(" # Wrote %s\n", out)
    @printf(" #   %d panels of %dx%d  ->  %d quads, %d nodes  (R = %.6e, %s)\n",
            6, n, n, length(quads), length(coords), R,
            split ? "panels SPLIT: seam nodes duplicated" : "watertight: seam nodes shared")
end

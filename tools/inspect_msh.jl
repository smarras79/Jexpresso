#!/usr/bin/env julia
#
# inspect_msh.jl -- report a structured .msh's element counts and extents.
#
#   julia tools/inspect_msh.jl <mesh.msh> [--machine]
#
# Exists so the setup script can take a user-supplied mesh's dimensions from
# the FILE rather than from a variable the user also has to keep in sync. Every
# mismatch of that kind in this project has cost a queue cycle: a mesh whose
# real shape differed from what the deck claimed produced empty ranks, a wrong
# resolution, or an out-of-memory kill, none of which name the mismatch.
#
# Assumes a structured box mesh (the LESICP family): element centroids fall on
# a regular lattice, so the number of DISTINCT centroid coordinates along each
# axis is the element count along that axis.

function parse_msh(path)
    lines = split(read(path, String), '\n')
    i = findfirst(==("\$Nodes"), lines)
    i === nothing && error("no \$Nodes section in $path (is this a .msh 4.1 file?)")
    nb = parse(Int, split(lines[i+1])[1])
    coords = Dict{Int,NTuple{3,Float64}}()
    p = i + 2
    for _ in 1:nb
        _, _, _, n = parse.(Int, split(lines[p])); p += 1
        tags = [parse(Int, lines[p+k]) for k in 0:n-1]; p += n
        for k in 0:n-1
            v = split(lines[p+k])          # may carry trailing parametric coords
            coords[tags[k+1]] = (parse(Float64, v[1]), parse(Float64, v[2]), parse(Float64, v[3]))
        end
        p += n
    end
    j = findfirst(==("\$Elements"), lines)
    j === nothing && error("no \$Elements section in $path")
    nb = parse(Int, split(lines[j+1])[1])
    p = j + 2
    cells = Vector{Vector{Int}}()
    for _ in 1:nb
        dim, _, _, n = parse.(Int, split(lines[p])); p += 1
        if dim == 3
            for k in 0:n-1
                push!(cells, [parse(Int, t) for t in split(lines[p+k])[2:end]])
            end
        end
        p += n
    end
    isempty(cells) && error("no 3D elements in $path")
    return coords, cells
end

function main()
    length(ARGS) < 1 && (println("usage: julia tools/inspect_msh.jl <mesh.msh> [--machine]"); exit(1))
    path = ARGS[1]
    isfile(path) || (println(stderr, "ERROR: no such file: $path"); exit(1))
    coords, cells = parse_msh(path)

    ndist(f) = begin
        s = Set{Float64}()
        for c in cells
            push!(s, round(sum(coords[t][f] for t in c) / length(c), digits=6))
        end
        length(s)
    end
    nx, ny, nz = ndist(1), ndist(2), ndist(3)

    xs = (minimum(v[1] for v in values(coords)), maximum(v[1] for v in values(coords)))
    ys = (minimum(v[2] for v in values(coords)), maximum(v[2] for v in values(coords)))
    zs = (minimum(v[3] for v in values(coords)), maximum(v[3] for v in values(coords)))
    lx, ly, lz = xs[2]-xs[1], ys[2]-ys[1], zs[2]-zs[1]
    structured = (nx * ny * nz == length(cells))

    if "--machine" in ARGS
        println("MESH_NX=$nx"); println("MESH_NY=$ny"); println("MESH_NZ=$nz")
        println("MESH_LX=$(round(Int, lx))"); println("MESH_LY=$(round(Int, ly))")
        println("MESH_LZ=$(round(Int, lz))")
        println("MESH_NELEM=$(length(cells))")
        println("MESH_STRUCTURED=$(structured ? 1 : 0)")
    else
        println("mesh      : $path")
        println("elements  : $(length(cells))  ->  $(nx) x $(ny) x $(nz)")
        println("extents   : $(round(lx,digits=1)) x $(round(ly,digits=1)) x $(round(lz,digits=1)) m")
        println("columns   : $(nx*ny)   (upper bound on ranks that can each own ground)")
        structured || println("WARNING: $(nx)x$(ny)x$(nz) = $(nx*ny*nz) != $(length(cells)) elements — " *
                              "this mesh is not a regular box lattice, so the rank planning " *
                              "in plan_run.jl does not apply to it.")
    end
end

main()

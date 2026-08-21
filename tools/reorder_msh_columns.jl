#!/usr/bin/env julia
#
# reorder_msh_columns.jl -- rewrite a .msh so its 3D elements are ordered
#                           COLUMN-MAJOR (all cells of one vertical column,
#                           bottom to top, then the next column).
#
#   julia tools/reorder_msh_columns.jl <in.msh> <out.msh> [nranks]
#
# WHY
# ---
# With :linitial_refine => true the mesh goes through p4est, and p4est gives
# each rank a CONTIGUOUS range of the space-filling curve. That curve follows
# coarse-tree order, which follows the order elements appear in the .msh.
#
# gmsh writes a box mesh in horizontal LAYERS. A contiguous range of a
# layer-ordered mesh is therefore a horizontal SLAB: the ranks holding the
# bottom layer own the entire ground surface and every other rank owns none.
# On the 8x2x60 case that produced
#
#     # Bottom faces (min/avg/max) : 0 / 1.0 / 60   [imbalance: 60.0x]
#
# i.e. 62 of 64 ranks doing no wall-model work at all, and one rank doing 94%
# of it, every timestep.
#
# Order the same elements column-major and a contiguous range becomes a
# VERTICAL COLUMN instead, so the ground is spread over as many ranks as there
# are columns.
#
# WHAT THIS DOES NOT FIX
# ----------------------
# If you run more ranks than the mesh has horizontal columns, some ranks still
# get no ground -- a column can only be split vertically, and only the bottom
# piece touches the surface. That is inherent to column decomposition, not to
# the ordering. To give EVERY rank ground you need nranks <= ncolumns; pass
# `nranks` below and this script reports exactly what you would get.
#
# The rewrite is deliberately minimal: node data, entity blocks, physical
# groups and all lower-dimensional (boundary) elements are copied verbatim.
# Only the order of the 3D element lines changes, with their tags renumbered
# ascending to match so that a reader keying on either file order or tag order
# sees the same sequence.

function read_nodes(lines)
    i = findfirst(==("\$Nodes"), lines)
    nblocks = parse(Int, split(lines[i+1])[1])
    coords = Dict{Int,NTuple{3,Float64}}()
    p = i + 2
    for _ in 1:nblocks
        _, _, _, n = parse.(Int, split(lines[p])); p += 1
        tags = [parse(Int, lines[p+k]) for k in 0:n-1]; p += n
        for k in 0:n-1
            v = split(lines[p+k])          # may carry extra parametric coords
            coords[tags[k+1]] = (parse(Float64, v[1]), parse(Float64, v[2]), parse(Float64, v[3]))
        end
        p += n
    end
    return coords
end

function main()
    if length(ARGS) < 2
        println("usage: julia tools/reorder_msh_columns.jl <in.msh> <out.msh> [nranks]")
        return
    end
    infile, outfile = ARGS[1], ARGS[2]
    nranks = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

    lines  = split(read(infile, String), '\n')
    coords = read_nodes(lines)

    e = findfirst(==("\$Elements"), lines)
    e === nothing && error("no \$Elements section in $infile")
    nblocks = parse(Int, split(lines[e+1])[1])

    out = copy(lines)
    total_reordered = 0
    ncols_report = 0

    p = e + 2
    for _ in 1:nblocks
        dim, ent, etype, n = parse.(Int, split(lines[p]))
        hdr = p; p += 1
        if dim == 3 && n > 1
            body = lines[p:p+n-1]
            # centroid of each element from its node list (tag is field 1)
            cent = map(body) do ln
                v  = split(ln)
                ns = parse.(Int, v[2:end])
                cx = sum(coords[t][1] for t in ns) / length(ns)
                cy = sum(coords[t][2] for t in ns) / length(ns)
                cz = sum(coords[t][3] for t in ns) / length(ns)
                (cx, cy, cz)
            end
            key(c) = (round(c[1], digits=6), round(c[2], digits=6), round(c[3], digits=6))
            # column-major: group by (x,y), then ascending z so the GROUND cell
            # is first in each column's run -- that is what puts the surface at
            # the start of a rank's contiguous range rather than buried in it.
            perm = sortperm(1:n; by = i -> key(cent[i]))
            tags = [parse(Int, split(body[i])[1]) for i in 1:n]
            sort!(tags)
            newbody = similar(body)
            for (newpos, oldidx) in enumerate(perm)
                v = split(body[oldidx])
                newbody[newpos] = string(tags[newpos], " ", join(v[2:end], " "))
            end
            out[p:p+n-1] = newbody
            total_reordered += n
            ncols_report = length(Set((key(c)[1], key(c)[2]) for c in cent))

            if nranks > 0
                zmin = minimum(c[3] for c in cent)
                isbot = [c[3] - zmin < 1e-6 for c in cent]
                for (label, order) in (("before (file order)", 1:n), ("after  (column-major)", perm))
                    per = zeros(Int, nranks)
                    chunk = cld(n, nranks)
                    for (pos, idx) in enumerate(order)
                        isbot[idx] && (per[min(cld(pos, chunk), nranks)] += 1)
                    end
                    println("  ", rpad(label, 24),
                            " ground cells/rank min=", minimum(per),
                            " max=", maximum(per),
                            "  ranks with none: ", count(iszero, per), "/", nranks)
                end
            end
        end
        p += n
    end

    open(outfile, "w") do io
        write(io, join(out, '\n'))
    end
    println("reordered $(total_reordered) 3D elements over $(ncols_report) columns")
    println("wrote $(outfile)")
    if nranks > 0
        println("\nnranks=$(nranks) vs $(ncols_report) columns: ",
                nranks <= ncols_report ?
                "OK -- every rank can own ground." :
                "NOTE only $(ncols_report) ranks can own ground; a column splits vertically.")
    end
end

main()

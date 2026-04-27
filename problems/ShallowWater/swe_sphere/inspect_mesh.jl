#
# Cubed-sphere mesh diagnostics for the swe_sphere case.
#
# Usage:
#   julia --project=. problems/ShallowWater/swe_sphere/inspect_mesh.jl \
#         [path/to/cubed_sphere.msh]
#
# Defaults to ./cubed_sphere.msh in this directory.
#
# Reads a Gmsh .msh format 4.1 file, verifies that the surface elements are
# 4-node quads embedded in ℝ³ with all node radii close to a common R,
# tags each quad with a panel id (0..5) by the dominant axis at the element
# centroid, writes mesh_diagnostics.txt, and writes a Paraview-readable
# legacy ASCII VTK with cell scalar `panel_id`.
#
# Self-contained: stdlib only, no GridapGmsh / Jexpresso dependency.
#

const PANEL_NAMES = ("+x", "-x", "+y", "-y", "+z", "-z")

# Centroid -> panel id 0..5 by dominant axis (ties go to lower index).
function panel_id(cx::Float64, cy::Float64, cz::Float64)
    ax, ay, az = abs(cx), abs(cy), abs(cz)
    if ax >= ay && ax >= az
        return cx >= 0 ? 0 : 1
    elseif ay >= az
        return cy >= 0 ? 2 : 3
    else
        return cz >= 0 ? 4 : 5
    end
end

# Minimal parser for Gmsh .msh format 4.1. Only the data we need:
#   - Node coordinates (x, y, z) keyed by node tag.
#   - 4-node quad connectivity (element type 3); points/lines are counted
#     and ignored.
function parse_msh4(path::AbstractString)
    nodes  = Dict{Int, NTuple{3, Float64}}()
    quads  = Vector{NTuple{4, Int}}()
    n_pts  = 0
    n_lines = 0
    n_other = 0

    open(path, "r") do io
        while !eof(io)
            line = strip(readline(io))
            if line == "\$Nodes"
                hdr = split(strip(readline(io)))
                nblocks = parse(Int, hdr[1])
                for _ in 1:nblocks
                    bhdr = split(strip(readline(io)))
                    nIn  = parse(Int, bhdr[4])
                    tags = Vector{Int}(undef, nIn)
                    for k in 1:nIn
                        tags[k] = parse(Int, strip(readline(io)))
                    end
                    for k in 1:nIn
                        toks = split(strip(readline(io)))
                        # Parametric blocks append (u) or (u,v); first 3 are x,y,z.
                        nodes[tags[k]] = (parse(Float64, toks[1]),
                                          parse(Float64, toks[2]),
                                          parse(Float64, toks[3]))
                    end
                end
            elseif line == "\$Elements"
                hdr = split(strip(readline(io)))
                nblocks = parse(Int, hdr[1])
                for _ in 1:nblocks
                    bhdr  = split(strip(readline(io)))
                    etype = parse(Int, bhdr[3])
                    nIn   = parse(Int, bhdr[4])
                    for _ in 1:nIn
                        toks = split(strip(readline(io)))
                        if etype == 3            # 4-node quadrilateral
                            push!(quads, (parse(Int, toks[2]),
                                          parse(Int, toks[3]),
                                          parse(Int, toks[4]),
                                          parse(Int, toks[5])))
                        elseif etype == 15
                            n_pts += 1
                        elseif etype == 1
                            n_lines += 1
                        else
                            n_other += 1
                        end
                    end
                end
            end
        end
    end

    return nodes, quads, n_pts, n_lines, n_other
end

function write_vtk_ascii(path::AbstractString,
                         nodes::Dict{Int, NTuple{3, Float64}},
                         quads::Vector{NTuple{4, Int}},
                         panel_ids::Vector{Int})
    # Compact node remapping: only nodes referenced by quads are exported.
    used = Set{Int}()
    for q in quads
        push!(used, q[1]); push!(used, q[2]); push!(used, q[3]); push!(used, q[4])
    end
    sorted_tags = sort!(collect(used))
    idx = Dict(t => i - 1 for (i, t) in enumerate(sorted_tags))
    nq = length(quads)

    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "cubed_sphere panel coloring")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        println(io, "POINTS $(length(sorted_tags)) double")
        for t in sorted_tags
            c = nodes[t]
            println(io, "$(c[1]) $(c[2]) $(c[3])")
        end
        println(io, "CELLS $(nq) $(5 * nq)")
        for q in quads
            println(io, "4 $(idx[q[1]]) $(idx[q[2]]) $(idx[q[3]]) $(idx[q[4]])")
        end
        println(io, "CELL_TYPES $(nq)")
        for _ in 1:nq
            println(io, "9")    # VTK_QUAD
        end
        println(io, "CELL_DATA $(nq)")
        println(io, "SCALARS panel_id int 1")
        println(io, "LOOKUP_TABLE default")
        for p in panel_ids
            println(io, p)
        end
    end
end

function main(args::Vector{String})
    msh_path = length(args) >= 1 ? args[1] : joinpath(@__DIR__, "cubed_sphere.msh")
    isfile(msh_path) || error("Mesh file not found: $(msh_path)")

    out_dir   = dirname(abspath(msh_path))
    diag_path = joinpath(out_dir, "mesh_diagnostics.txt")
    vtk_path  = joinpath(out_dir, "cubed_sphere_panels.vtk")

    nodes, quads, n_pts, n_lines, n_other =
        parse_msh4(msh_path)

    isempty(quads) && error("No 4-node quad elements (gmsh type 3) in $(msh_path).")

    # Sphere radius / deviation across *all* node tags, not just the ones
    # referenced by quads -- catches stray nodes (origin point, etc.).
    radii = Float64[sqrt(c[1]^2 + c[2]^2 + c[3]^2) for c in values(nodes)]
    nz_radii = filter(r -> r > 1e-12, radii)
    R   = sum(nz_radii) / length(nz_radii)
    dev = maximum(abs.(nz_radii .- R))
    n_offsphere = count(r -> r <= 1e-12, radii)

    # Panel tagging + edge length sweep over quads.
    counts    = zeros(Int, 6)
    panel_ids = Vector{Int}(undef, length(quads))
    minedge   =  Inf
    maxedge   = -Inf

    for (i, q) in enumerate(quads)
        p1 = nodes[q[1]]; p2 = nodes[q[2]]; p3 = nodes[q[3]]; p4 = nodes[q[4]]
        cx = (p1[1] + p2[1] + p3[1] + p4[1]) / 4
        cy = (p1[2] + p2[2] + p3[2] + p4[2]) / 4
        cz = (p1[3] + p2[3] + p3[3] + p4[3]) / 4
        pid = panel_id(cx, cy, cz)
        panel_ids[i] = pid
        counts[pid + 1] += 1

        @inbounds for (a, b) in ((p1, p2), (p2, p3), (p3, p4), (p4, p1))
            d = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)
            d < minedge && (minedge = d)
            d > maxedge && (maxedge = d)
        end
    end

    # mesh_diagnostics.txt
    open(diag_path, "w") do io
        println(io, "cubed_sphere mesh diagnostics")
        println(io, "  source              : $(msh_path)")
        println(io, "  nodes (total)       : $(length(nodes))")
        println(io, "  nodes at origin     : $(n_offsphere)")
        println(io, "  quad elements       : $(length(quads))")
        println(io, "  point elements      : $(n_pts)")
        println(io, "  line  elements      : $(n_lines)")
        println(io, "  other elements      : $(n_other)")
        println(io, "  R (mean node norm)  : $(R)")
        println(io, "  max |‖x‖ - R|       : $(dev)")
        println(io, "  edge length min     : $(minedge)")
        println(io, "  edge length max     : $(maxedge)")
        println(io, "")
        println(io, "  nelem per panel     :")
        for pid in 0:5
            println(io, "    panel $(pid) ($(PANEL_NAMES[pid + 1])): $(counts[pid + 1])")
        end
    end

    write_vtk_ascii(vtk_path, nodes, quads, panel_ids)

    # Console summary
    println("nodes      : $(length(nodes))   (origin nodes: $(n_offsphere))")
    println("quads      : $(length(quads))")
    println("R          : $(R)")
    println("max|r - R| : $(dev)")
    println("edge min   : $(minedge)")
    println("edge max   : $(maxedge)")
    for pid in 0:5
        println("panel $(pid) ($(PANEL_NAMES[pid + 1])): $(counts[pid + 1])")
    end
    println("wrote      : $(diag_path)")
    println("wrote      : $(vtk_path)")
end

main(ARGS)

#---------------------------------------------------------------------------------
# check_mesh.jl — is naca64A210.msh a grid this case can actually be curved on?
#
#     julia problems/CompEuler/naca64A210/check_mesh.jl
#
# Bare Julia, no package: it reads the .msh and user_exactGeo.jl and nothing
# else. Run it after changing the resolution in make_mesh.jl, because the two
# ways a curved-boundary grid goes wrong are both invisible in a mesh viewer:
#
#   1. THE GRID IS NOT THIS AEROFOIL. If the .msh and the deck's placement have
#      drifted apart, snapping to the section drags the whole wall somewhere
#      nobody asked for. user_exactGeo_setup refuses when that is gross (> 1% of
#      chord); this prints the actual number, which on a mesh made by
#      make_mesh.jl is round-off.
#
#   2. THE ELEMENTS FOLD. Curving pushes the mid-edge nodes of a boundary edge
#      out to the true section, by up to the SAGITTA of that edge's arc. If an
#      element is thinner, normal to the wall, than its own sagitta, it turns
#      inside out and the Jacobian changes sign. exact_geometry.jl stops the run
#      when that happens (_check_curved_elements), but at that point the grid is
#      already built; this says it beforehand, and says by how much margin.
#
# It also checks the two things Jexpresso's 2D gmsh reader assumes: every
# element is a quadrilateral, and the physical groups the deck and user_bc.jl
# name are all present.
#---------------------------------------------------------------------------------

using Printf

const CASE_DIR = @__DIR__
include(joinpath(CASE_DIR, "user_exactGeo.jl"))

const MSH = joinpath(CASE_DIR, "naca64A210.msh")

#---------------------------------------------------------------------------------
# Minimal gmsh 4.1 reader: physical names, node coordinates, and the elements of
# each entity. Enough for the checks above and no more — the real reader is
# src/kernel/mesh/mesh.jl.
#---------------------------------------------------------------------------------
function read_msh4(path)
    lines = readlines(path)
    i = 1
    names   = Dict{Int,String}()                  # physical tag  -> name
    curvephys = Dict{Int,Vector{Int}}()           # curve entity  -> physical tags
    coord   = Dict{Int,Tuple{Float64,Float64}}()  # node tag      -> (x, y)
    edges   = Tuple{Int,Int,Int}[]                # (curve entity, n1, n2)
    quads   = NTuple{4,Int}[]

    while i <= length(lines)
        sec = strip(lines[i]); i += 1

        if sec == "\$PhysicalNames"
            n = parse(Int, strip(lines[i])); i += 1
            for _ = 1:n
                m = match(r"^\s*(\d+)\s+(\d+)\s+\"(.*)\"", lines[i]); i += 1
                names[parse(Int, m.captures[2])] = m.captures[3]
            end

        elseif sec == "\$Entities"
            np, nc, ns, nv = parse.(Int, split(strip(lines[i]))); i += 1
            i += np                                            # points: skipped
            for _ = 1:nc
                f   = split(strip(lines[i])); i += 1
                tag = parse(Int, f[1])
                npt = parse(Int, f[8])                         # numPhysicalTags
                curvephys[tag] = [parse(Int, f[8+k]) for k = 1:npt]
            end
            i += ns + nv

        elseif sec == "\$Nodes"
            nb, _, _, _ = parse.(Int, split(strip(lines[i]))); i += 1
            for _ = 1:nb
                f = split(strip(lines[i])); i += 1
                nn = parse(Int, f[4])
                tags = [parse(Int, strip(lines[i + k - 1])) for k = 1:nn]; i += nn
                for k = 1:nn
                    c = parse.(Float64, split(strip(lines[i]))); i += 1
                    coord[tags[k]] = (c[1], c[2])
                end
            end

        elseif sec == "\$Elements"
            nb, _, _, _ = parse.(Int, split(strip(lines[i]))); i += 1
            for _ = 1:nb
                f = split(strip(lines[i])); i += 1
                edim, etag, etype, ne = parse(Int, f[1]), parse(Int, f[2]),
                                        parse(Int, f[3]), parse(Int, f[4])
                for _ = 1:ne
                    g = parse.(Int, split(strip(lines[i]))); i += 1
                    if etype == 1                              # 2-node line
                        push!(edges, (etag, g[2], g[3]))
                    elseif etype == 3                          # 4-node quad
                        push!(quads, (g[2], g[3], g[4], g[5]))
                    elseif edim == 2
                        error("check_mesh.jl: element type $etype in a surface block — " *
                              "the 2D reader wants quadrilaterals only (type 3). " *
                              "Is Mesh.RecombinationAlgorithm = 3 still set in the .geo?")
                    end
                end
            end
        end
    end
    return names, curvephys, coord, edges, quads
end

#---------------------------------------------------------------------------------
# The checks. In a function, not at top level: Julia's soft scope would make
# every accumulator inside a `for` a fresh local and the totals would come out
# undefined.
#---------------------------------------------------------------------------------
function main()
    names, curvephys, coord, edges, quads = read_msh4(MSH)

    println("naca64A210.msh")
    @printf("  %d nodes, %d quads, %d boundary line elements\n",
            length(coord), length(quads), length(edges))

    #---------------------------------------------------------------------------------
    # The physical groups the deck and user_bc.jl dispatch on.
    #---------------------------------------------------------------------------------
    have = Set(values(names))
    for want in ("bottom", "outflow", "top", "inflow", "airfoil", "domain")
        want in have || error("check_mesh.jl: the mesh has no physical group \"$want\". " *
                              "Present: " * join(sort(collect(have)), ", "))
    end
    println("  physical groups: ", join(sort(collect(have)), ", "))

    #---------------------------------------------------------------------------------
    # The aerofoil boundary, and the section the deck places.
    #---------------------------------------------------------------------------------
    airfoil_phys = only([t for (t, n) in names if n == "airfoil"])
    airfoil_ents = Set(t for (t, ps) in curvephys if airfoil_phys in ps)
    wall = [(e[2], e[3]) for e in edges if e[1] in airfoil_ents]
    isempty(wall) && error("check_mesh.jl: no boundary edge carries the \"airfoil\" tag.")

    function deck_placement(file)
        for line in eachline(file)
            code = lstrip(line)
            startswith(code, "#") && continue
            m = match(r":naca64A210((?:\s*,\s*[-\d.eE+]+){5})\s*\)", code)
            m === nothing && continue
            return ntuple(i -> parse(Float64, strip(split(m.captures[1], ",")[i+1])), 5)
        end
        error("check_mesh.jl: no (:naca64A210, ...) spec in " * file)
    end
    xle, yle, chord, alpha, rte = deck_placement(joinpath(CASE_DIR, "user_inputs.jl"))
    sec = naca64A210_section(xle, yle, chord, alpha, rte)
    @printf("  section: NACA 64A210, chord %.4g, LE (%.4g, %.4g), α %.4g deg ; %d wall edges\n",
            chord, xle, yle, alpha, length(wall))

    #---------------------------------------------------------------------------------
    # 1. Are the grid's own vertices on the section?
    #---------------------------------------------------------------------------------
    wallnodes = Set{Int}()
    for (a, b) in wall; push!(wallnodes, a); push!(wallnodes, b); end
    off = 0.0
    for n in wallnodes
        x, y = coord[n]
        px, py = user_exactGeo("airfoil", sec, x, y)
        off = max(off, hypot(px - x, py - y))
    end
    @printf("  1. vertices off the section : %.3e  (%.2e %% of chord)%s\n",
            off, 100*off/chord, off > 0.01*chord ? "   <-- REFUSED by user_exactGeo_setup" : "")

    #---------------------------------------------------------------------------------
    # 2. Sagitta of each wall edge against the thickness of the element behind it.
    #
    # The sagitta is how far the curving moves the mid-edge nodes; the thickness is
    # how much room the element has to give. Their ratio is the fold margin, and it
    # has to stay comfortably above 1.
    #---------------------------------------------------------------------------------
    # For each wall edge, the quad that owns it.
    edgequad = Dict{Tuple{Int,Int},NTuple{4,Int}}()
    for q in quads, k = 1:4
        a, b = q[k], q[k == 4 ? 1 : k+1]
        edgequad[minmax(a, b)] = q
    end

    mid_of(a, b) = (0.5*(a[1] + b[1]), 0.5*(a[2] + b[2]))
    L2sq(a, b)   = (b[1]-a[1])^2 + (b[2]-a[2])^2

    # Distance from p to the segment (a, b).
    function seg_dist(p, a, b)
        vx, vy = b[1] - a[1], b[2] - a[2]
        wx, wy = p[1] - a[1], p[2] - a[2]
        L2 = vx*vx + vy*vy
        s  = L2 > 0 ? clamp((wx*vx + wy*vy)/L2, 0.0, 1.0) : 0.0
        return hypot(wx - s*vx, wy - s*vy)
    end

    # Parameter on the section nearest to a point: user_exactGeo's projection,
    # but returning the PARAMETER rather than the point, because the sagitta
    # needs the arc between two nodes and not just their images.
    #
    # Sweep AND Newton, Newton started from EVERY local minimum of the sweep —
    # the same construction, and for the same reason, as user_exactGeo itself
    # (see the note above _NACA_SWEEP_PER_SEGMENT there). A single-start search
    # near the sharp trailing edge locks onto whichever surface happened to hold
    # the nearest sample, which is not always the one the node is on, and the
    # arc it then reports spans the wrong part of the section. That shows up
    # here as a sagitta of a millimetre where the true one is 1.7e-4 —
    # measurement noise dressed as geometry.
    function param_of(sec, x, y)
        # Only the spline BETWEEN the tangency points is boundary; the rest was
        # replaced by the trailing-edge arc.
        tlo = sec.te === nothing ? sec.t[1]   : sec.te.tu
        thi = sec.te === nothing ? sec.t[end] : sec.te.tl
        n  = (length(sec.t) - 1)*32
        dt = (thi - tlo)/n
        function refine(t0)
            lo, hi = max(tlo, t0 - dt), min(thi, t0 + dt)
            t = t0
            for _ = 1:40
                px, py, dx, dy, d2x, d2y = naca64A210_point(sec, t)
                g  = (px - x)*dx + (py - y)*dy
                gp = dx*dx + dy*dy + (px - x)*d2x + (py - y)*d2y
                abs(gp) > 0.0 || break
                tn = clamp(t - g/gp, lo, hi)
                abs(tn - t) <= 1.0e-15*max(1.0, abs(t)) && (t = tn; break)
                t = tn
            end
            px, py, _, _, _, _ = naca64A210_point(sec, t)
            return t, (px - x)^2 + (py - y)^2
        end
        tbest, fbest = tlo, Inf
        consider(t0) = ((tc, fc) = refine(t0); fc < fbest && (fbest = fc; tbest = tc))
        fm2, fm1, tm1 = Inf, Inf, tlo
        for k = 0:n
            t = tlo + (thi - tlo)*k/n
            px, py, _, _, _, _ = naca64A210_point(sec, t)
            f = (px - x)^2 + (py - y)^2
            fm1 <= fm2 && fm1 <= f && consider(tm1)
            fm2, fm1, tm1 = fm1, f, t
        end
        consider(tlo); consider(thi)
        fm1 <= fm2 && consider(tm1)
        return tbest
    end

    worst_ratio = Inf
    worst_sag   = 0.0
    hmin        = Inf
    at_ratio    = (0.0, 0.0); at_sag = (0.0, 0.0)
    hr_ratio    = 0.0; hr_sag = 0.0
    for (a, b) in wall
        pa, pb = coord[a], coord[b]
        hmin = min(hmin, hypot(pb[1] - pa[1], pb[2] - pa[2]))

        # Sagitta. An edge lying on the ROUNDED TRAILING EDGE is a chord of a
        # circle and its sagitta is r - sqrt(r^2 - (L/2)^2) exactly; only the
        # spline part needs the parameter search below. (Before the trailing
        # edge was rounded every edge was on the spline, and an arc edge sent
        # through param_of comes back with a parameter span across the whole
        # aerofoil, because its nearest spline point is a tangency point.)
        onarc(q) = sec.te !== nothing &&
                   abs(hypot(q[1] - sec.te.cx, q[2] - sec.te.cy) - sec.te.r) <
                       1.0e-8*sec.chord
        if onarc(pa) && onarc(pb)
            r   = sec.te.r
            sag = r - sqrt(max(r*r - 0.25*L2sq(pa, pb), 0.0))
            worst_sag < sag && (worst_sag = sag; at_sag = mid_of(pa, pb); hr_sag = sqrt(L2sq(pa,pb)))
            q = get(edgequad, minmax(a, b), nothing)
            if q !== nothing
                others = [n for n in q if n != a && n != b]
                if length(others) == 2 && sag > 0
                    thick = seg_dist(mid_of(pa,pb), coord[others[1]], coord[others[2]])
                    if thick/sag < worst_ratio
                        worst_ratio = thick/sag; at_ratio = mid_of(pa,pb); hr_ratio = sqrt(L2sq(pa,pb))
                    end
                end
            end
            continue
        end
        ta, tb = param_of(sec, pa...), param_of(sec, pb...)
        # The trailing edge is the one point of the section with TWO parameters,
        # t = 0 and t = t_end, and the sweep above always reports the first. So
        # the edge that runs into the trailing edge along the LOWER surface comes
        # back as a parameter span across the entire aerofoil. Whenever the span
        # is absurd for the edge's length, that is what happened: move the
        # endpoint that sits at a parameter end to the other one.
        L = hypot(pb[1] - pa[1], pb[2] - pa[2])
        if abs(tb - ta) > 5*L
            ta < tb ? (ta = sec.t[end]) : (tb = sec.t[end])
        end
        abs(tb - ta) > 5*L && error(
            "check_mesh.jl: wall edge of length " * string(L) * " spans a section " *
            "parameter range of " * string(abs(tb - ta)) * ". Either the mesh is not " *
            "this aerofoil (see check 1) or an element straddles the trailing edge.")
        sag = 0.0
        for k = 1:19
            t = ta + (tb - ta)*k/20
            px, py, _, _, _, _ = naca64A210_point(sec, t)
            sag = max(sag, seg_dist((px, py), pa, pb))
        end
        sag > worst_sag && (worst_sag = sag; at_sag = mid_of(pa, pb); hr_sag = L)

        # Thickness: how far the opposite edge of the owning quad is from this one.
        q = get(edgequad, minmax(a, b), nothing)
        q === nothing && continue
        others = [n for n in q if n != a && n != b]
        length(others) == 2 || continue
        mid  = (0.5*(pa[1] + pb[1]), 0.5*(pa[2] + pb[2]))
        thick = seg_dist(mid, coord[others[1]], coord[others[2]])
        if sag > 0 && thick/sag < worst_ratio
            worst_ratio = thick/sag
            at_ratio    = mid
            hr_ratio    = L
        end
    end

    @printf("  2. largest wall-edge sagitta: %.3e at (%.4f, %.4f), edge %.3e  (x/c = %.3f)\n",
            worst_sag, at_sag[1], at_sag[2], hr_sag, (at_sag[1] - xle)/chord)
    @printf("     smallest wall edge       : %.3e\n", hmin)
    @printf("     worst thickness/sagitta  : %.1f at (%.4f, %.4f), edge %.3e (x/c = %.3f)%s\n",
            worst_ratio, at_ratio[1], at_ratio[2], hr_ratio, (at_ratio[1] - xle)/chord,
            worst_ratio < 3 ? "   <-- THIN. Elements may fold; refine normal to the wall." : "   (fold margin, wants >> 1)")

    println(worst_ratio > 3 && off <= 0.01*chord ?
            "  OK — this grid can be curved onto the section." :
            "  PROBLEM — see the markers above.")
end

main()

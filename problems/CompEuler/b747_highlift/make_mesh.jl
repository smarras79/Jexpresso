#---------------------------------------------------------------------------------
# make_mesh.jl — write b747_highlift.geo, a FOUR-ELEMENT high-lift section.
#
#     julia problems/CompEuler/b747_highlift/make_mesh.jl
#     gmsh -2 problems/CompEuler/b747_highlift/b747_highlift.geo \
#          -format msh4 -o problems/CompEuler/b747_highlift/b747_highlift.msh
#
# Both outputs are committed, so nothing here has to be run to USE the case.
#
# THE CONFIGURATION. A leading-edge slat, the main element, and two trailing-edge
# flaps, deployed in the arrangement of a transport wing in the landing
# configuration — a 747-style layout, not a reproduction of any published 747
# section (those ordinates are not public). Every element is a NACA 64A210 of its
# own chord, deflection and trailing-edge radius, which is what makes the whole
# thing exactly definable: see ./user_exactGeo.jl.
#
# EVERY EDGE IS ROUND. The 64A210 has a rounded leading edge already (radius
# 0.0056c); the trailing edge of each element is truncated and capped with a
# circle tangent to both surfaces. So there is no point anywhere on any of the
# four boundaries at which the normal is undefined, which on a single element was
# what made the free-slip condition misbehave at the cusp. With four elements and
# eight edges the problem would have been four times over.
#
# THE SLOTS ARE THE MESH. Element sizing here is not about the aerofoils, it is
# about the gaps between them: the slot flow is the whole point of a high-lift
# section, and a slot resolved by two cells is not resolved. As shipped the three
# gaps are 0.032, 0.049 and 0.044 m against H_SURF, i.e. 6-10 elements across the
# narrowest one. check_mesh.jl prints them.
#
# The placement lives in the deck's :exact_geometry — one entry per element, and
# this script reads all four out of it, so the mesh and the analytic geometry
# cannot drift apart.
#---------------------------------------------------------------------------------

using Printf

const CASE_DIR = @__DIR__

include(joinpath(CASE_DIR, "user_exactGeo.jl"))

#---------------------------------------------------------------------------------
# The elements, read out of the deck.
#
# Read TEXTUALLY rather than by evaluating user_inputs(), which would need
# CarpenterKennedy2N54, DSGS, TOTAL and the rest of the package: this script has
# to run on a bare Julia, before anything is instantiated, because its whole job
# is to produce the mesh the package then reads.
#
# Returns [(name, xle, yle, chord, α_deg, r_TE/chord), ...] in deck order, which
# is also the order the .geo lays them out.
#---------------------------------------------------------------------------------
function deck_elements(file)
    out = Tuple{String,Float64,Float64,Float64,Float64,Float64}[]
    for line in eachline(file)
        code = lstrip(line)
        startswith(code, "#") && continue
        m = match(r"\"([A-Za-z0-9_]+)\"\s*=>\s*\(:naca64A210((?:\s*,\s*[-\d.eE+]+){5})\s*\)", code)
        m === nothing && continue
        v = [parse(Float64, strip(s)) for s in split(m.captures[2], ",")[2:end]]
        push!(out, (m.captures[1], v[1], v[2], v[3], v[4], v[5]))
    end
    isempty(out) && error("make_mesh.jl: no \"<name>\" => (:naca64A210, x_LE, y_LE, chord, " *
                          "α, r_TE) entries found in " * file * ".\n  The deck's " *
                          ":exact_geometry is what places the elements; this script reads it.")
    return out
end

const ELEMS = deck_elements(joinpath(CASE_DIR, "user_inputs.jl"))

# The tunnel. It lives here rather than in the deck because the deck never needs
# it: the grid IS the domain, and mod_mesh_read_gmsh! reads the bounds off it.
const XMIN, XMAX =  0.0, 3.0
const YMIN, YMAX = -1.0, 1.0

#---------------------------------------------------------------------------------
# Resolution.
#
# NSURF is how finely the .geo DRAWS each element's two surfaces (the arcs are
# exact gmsh Circles, not sampled). gmsh's `Spline` is Catmull-Rom and the
# section is the interpolating cubic spline, so the two agree at the sampled
# points and differ a little in between; check_mesh.jl measures the difference.
#
#   H_TIP   at every leading edge and every trailing-edge apex. The binding one
#           is the smallest TE arc: at r = 0.003 m the arc is 0.0088 m long, so
#           H_TIP must be well under that or the arc gets two elements and gmsh
#           returns inverted quads.
#   H_SURF  along the rest of every surface, and therefore across the slots.
#   H_FAR   in the far field.
#
# H_TIP IS THE TIME STEP. The DynSGS diffusive limit is Δt ≲ ½Δx²/(μ⃗·C2·Δ·(|u|+c))
# with Δx ∝ Δ ∝ the smallest element, so Δt ∝ H_TIP. Four elements each needing
# a resolved round tip is what makes this case cost what it does; see the :Δt
# note in user_inputs.jl.
#---------------------------------------------------------------------------------
const NSURF    = 201
const H_TIP    = 0.0012
const H_SURF   = 0.0050
const H_FAR    = 0.1000
const D_NEAR   = 0.030             # distance out to which H_SURF holds
const D_FAR    = 0.900             # distance at which H_FAR is reached
const D_TIP    = 0.010             # radius of the H_TIP patch
const D_TIPFAR = 0.120             # ...and how far it relaxes

# Wake band behind the last flap. X0 starts where flap2 ends (x = 1.869), so
# the band is wake and not aerofoil.
const WAKE_H    = 0.010
const WAKE_X0   = 1.87
const WAKE_X1   = 2.70
const WAKE_YMIN = -0.42
const WAKE_YMAX =  0.10
const WAKE_THK  = 0.10

#---------------------------------------------------------------------------------
# Sample one element's two surfaces between the trailing-edge tangency points.
# Cosine clustering on each surface separately, so both the leading edge and the
# two tangency points get the dense end of the distribution.
#---------------------------------------------------------------------------------
function sample_surfaces(sec, nper)
    tle = sec.t[argmin(sec.x)]
    tlo = sec.te === nothing ? sec.t[1]   : sec.te.tu
    thi = sec.te === nothing ? sec.t[end] : sec.te.tl
    pts = Tuple{Float64,Float64}[]
    for (a, b, skipfirst) in ((tlo, tle, false), (tle, thi, true))
        for k = 0:nper
            skipfirst && k == 0 && continue
            s = 0.5*(1.0 - cos(π*k/nper))
            x, y, _, _, _, _ = naca64A210_point(sec, a + (b - a)*s)
            push!(pts, (x, y))
        end
    end
    return pts
end

const SECS = [naca64A210_section(e[2], e[3], e[4], e[5], e[6]) for e in ELEMS]
const PTS  = [sample_surfaces(s, (NSURF - 1) ÷ 2) for s in SECS]

#---------------------------------------------------------------------------------
# Write the .geo.
#---------------------------------------------------------------------------------
open(joinpath(CASE_DIR, "b747_highlift.geo"), "w") do io

    println(io, "// b747_highlift.geo — GENERATED by make_mesh.jl. Do not edit by hand:")
    println(io, "// re-run  julia problems/CompEuler/b747_highlift/make_mesh.jl  instead.")
    println(io, "//")
    println(io, "// A four-element high-lift section in a Mach-3 wind tunnel: leading-edge")
    println(io, "// slat, main element, and two trailing-edge flaps. Each element is a")
    println(io, "// NACA 64A210 of its own chord and deflection with a ROUNDED trailing")
    println(io, "// edge — a circle tangent to both surfaces — so no point of any of the")
    println(io, "// four boundaries is a cusp. user_exactGeo.jl holds the definition and")
    println(io, "// the solver curves the high-order nodes onto it.")
    println(io, "//")
    for (k, e) in enumerate(ELEMS)
        @printf(io, "//   %-6s chord %.4g, LE (%.4g, %.4g), deflection %+.4g deg, TE r = %.4g m\n",
                e[1], e[4], e[2], e[3], -e[5], SECS[k].te.r)
    end
    @printf(io, "// tunnel  : [%.6g, %.6g] x [%.6g, %.6g]\n", XMIN, XMAX, YMIN, YMAX)
    println(io)

    @printf(io, "hfar = %.10g;\n\n", H_FAR)

    println(io, "// ---- tunnel ----------------------------------------------------------")
    @printf(io, "Point(1) = {%.10g, %.10g, 0, hfar};\n", XMIN, YMIN)
    @printf(io, "Point(2) = {%.10g, %.10g, 0, hfar};\n", XMAX, YMIN)
    @printf(io, "Point(3) = {%.10g, %.10g, 0, hfar};\n", XMAX, YMAX)
    @printf(io, "Point(4) = {%.10g, %.10g, 0, hfar};\n", XMIN, YMAX)
    println(io, "Line(1) = {1, 2};   // bottom")
    println(io, "Line(2) = {2, 3};   // outflow")
    println(io, "Line(3) = {3, 4};   // top")
    println(io, "Line(4) = {4, 1};   // inflow")
    println(io)

    # Per element: points 1000e+.., curves 10e+1..10e+4, Curve Loop 1+e.
    allcurves = Int[]
    tips      = Int[]
    loops     = Int[]
    for (k, e) in enumerate(ELEMS)
        P   = PTS[k]
        te  = SECS[k].te
        b   = 1000*k
        ile = b + argmin([p[1] for p in P])
        ic  = b + length(P) + 1
        ia  = b + length(P) + 2
        c1, c2, c3, c4 = 10k + 1, 10k + 2, 10k + 3, 10k + 4

        @printf(io, "// ---- %s ----\n", e[1])
        for (i, p) in enumerate(P)
            @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};\n", b + i, p[1], p[2])
        end
        @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};   // TE arc centre\n", ic, te.cx, te.cy)
        @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};   // TE apex\n",       ia, te.xa, te.ya)
        @printf(io, "Spline(%d) = {%d:%d};   // upper surface, TE tangency -> LE\n", c1, b + 1, ile)
        @printf(io, "Spline(%d) = {%d:%d};   // lower surface, LE -> TE tangency\n", c2, ile, b + length(P))
        @printf(io, "Circle(%d) = {%d, %d, %d};   // lower tangency -> apex\n", c3, b + length(P), ic, ia)
        @printf(io, "Circle(%d) = {%d, %d, %d};   // apex -> upper tangency\n", c4, ia, ic, b + 1)
        @printf(io, "Curve Loop(%d) = {%d, %d, %d, %d};\n", 1 + k, c1, c2, c3, c4)
        println(io)

        append!(allcurves, (c1, c2, c3, c4))
        append!(tips, (ile, ia))
        push!(loops, 1 + k)
    end

    println(io, "Curve Loop(1) = {1, 2, 3, 4};      // tunnel, counter-clockwise")
    @printf(io, "Plane Surface(1) = {1, %s};        // the tunnel with four holes\n",
            join(loops, ", "))
    println(io)

    println(io, "// ---- element size ----------------------------------------------------")
    println(io, "// Distance/Threshold rather than a background mesh so the grading follows")
    println(io, "// the geometry. The slots between the elements sit inside the near band,")
    println(io, "// so they are meshed at H_SURF.")
    println(io, "Field[1] = Distance;")
    @printf(io, "Field[1].CurvesList = {%s};\n", join(allcurves, ", "))
    println(io, "Field[1].Sampling = 400;")
    println(io, "Field[2] = Threshold;")
    println(io, "Field[2].InField = 1;")
    @printf(io, "Field[2].SizeMin = %.10g;\n", H_SURF)
    @printf(io, "Field[2].SizeMax = %.10g;\n", H_FAR)
    @printf(io, "Field[2].DistMin = %.10g;\n", D_NEAR)
    @printf(io, "Field[2].DistMax = %.10g;\n", D_FAR)
    println(io)
    println(io, "// Every leading edge and every trailing-edge apex. The smallest TE arc")
    println(io, "// must carry several elements or gmsh returns inverted quads.")
    println(io, "Field[3] = Distance;")
    @printf(io, "Field[3].PointsList = {%s};\n", join(tips, ", "))
    println(io, "Field[4] = Threshold;")
    println(io, "Field[4].InField = 3;")
    @printf(io, "Field[4].SizeMin = %.10g;\n", H_TIP)
    @printf(io, "Field[4].SizeMax = %.10g;\n", H_FAR)
    @printf(io, "Field[4].DistMin = %.10g;\n", D_TIP)
    @printf(io, "Field[4].DistMax = %.10g;\n", D_TIPFAR)
    println(io)
    println(io, "// The wake. Everything the four elements shed goes through here.")
    println(io, "Field[5] = Box;")
    @printf(io, "Field[5].VIn  = %.10g;\n", WAKE_H)
    @printf(io, "Field[5].VOut = %.10g;\n", H_FAR)
    @printf(io, "Field[5].XMin = %.10g;\n", WAKE_X0)
    @printf(io, "Field[5].XMax = %.10g;\n", WAKE_X1)
    @printf(io, "Field[5].YMin = %.10g;\n", WAKE_YMIN)
    @printf(io, "Field[5].YMax = %.10g;\n", WAKE_YMAX)
    @printf(io, "Field[5].Thickness = %.10g;\n", WAKE_THK)
    println(io)
    println(io, "Field[6] = Min;")
    println(io, "Field[6].FieldsList = {2, 4, 5};")
    println(io, "Background Field = 6;")
    println(io, "Mesh.MeshSizeExtendFromBoundary = 0;")
    println(io, "Mesh.MeshSizeFromPoints = 0;")
    println(io, "Mesh.MeshSizeFromCurvature = 0;")
    println(io)

    println(io, "// ---- quads -----------------------------------------------------------")
    println(io, "// Jexpresso's 2D reader wants quadrilaterals, and ALL of them: a")
    println(io, "// quad-dominant mesh with leftover triangles is not one it can read.")
    println(io, "Recombine Surface {1};")
    println(io, "Mesh.RecombinationAlgorithm = 3;   // blossom, full-quad")
    println(io, "Mesh.Algorithm = 8;                // frontal-Delaunay for quads")
    println(io)

    println(io, "// ---- physical groups -------------------------------------------------")
    println(io, "// user_bc.jl treats anything that is not \"inflow\" or \"outflow\" as a")
    println(io, "// free-slip wall, so all four elements and the two tunnel walls are walls.")
    println(io, "// These names are also the keys of :exact_geometry in user_inputs.jl.")
    println(io, "Physical Curve(\"bottom\",  1) = {1};")
    println(io, "Physical Curve(\"outflow\", 2) = {2};")
    println(io, "Physical Curve(\"top\",     3) = {3};")
    println(io, "Physical Curve(\"inflow\",  4) = {4};")
    for (k, e) in enumerate(ELEMS)
        @printf(io, "Physical Curve(\"%s\", %d) = {%d, %d, %d, %d};\n",
                e[1], 4 + k, 10k + 1, 10k + 2, 10k + 3, 10k + 4)
    end
    println(io, "Physical Surface(\"domain\") = {1};")
end

println("wrote ", joinpath(CASE_DIR, "b747_highlift.geo"), "  (", length(ELEMS), " elements, ",
        sum(length, PTS), " surface points)")
for (k, e) in enumerate(ELEMS)
    @printf("   %-6s chord %.4g  deflection %+6.2f deg  TE r = %.4g m  (%d pts)\n",
            e[1], e[4], -e[5], SECS[k].te.r, length(PTS[k]))
end
println("now:  gmsh -2 ", joinpath(CASE_DIR, "b747_highlift.geo"),
        " -format msh4 -o ", joinpath(CASE_DIR, "b747_highlift.msh"))

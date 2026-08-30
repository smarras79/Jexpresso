#---------------------------------------------------------------------------------
# make_mesh.jl — write naca64A210.geo from the section user_exactGeo.jl defines.
#
#     julia problems/CompEuler/naca64A210/make_mesh.jl
#     gmsh -2 problems/CompEuler/naca64A210/naca64A210.geo \
#          -format msh4 -o problems/CompEuler/naca64A210/naca64A210.msh
#
# Both outputs are committed, so nothing here has to be run to USE the case.
# Run it when you change the placement (chord, incidence, where the aerofoil
# sits in the tunnel) or the resolution.
#
# WHY THIS SCRIPT EXISTS AT ALL, rather than a hand-written .geo: a NACA
# 6A-series section is a table of ordinates, and gmsh's own `Spline` through
# those 51 points is a Catmull-Rom curve, which is NOT the interpolating cubic
# spline that everyone means by "the NACA 64A210". Near the nose, where the
# published points are 0.004c apart across a 0.0056c nose radius, the two differ by
# far more than the curving is trying to fix. So the .geo is written from the
# ACTUAL section: this script includes user_exactGeo.jl and samples
# naca64A210_section, and the geometry the mesh generator sees is the geometry
# the solver snaps to, by construction rather than by coincidence.
#
# It also reads the placement out of user_inputs.jl — the very tuple the deck
# hands to user_exactGeo — so the mesh and the exact geometry cannot drift
# apart. Edit :exact_geometry, re-run this, and both follow.
#---------------------------------------------------------------------------------

using Printf

const CASE_DIR = @__DIR__

include(joinpath(CASE_DIR, "user_exactGeo.jl"))

#---------------------------------------------------------------------------------
# The placement, read out of the deck.
#
# Read TEXTUALLY rather than by evaluating user_inputs(), which would need
# CarpenterKennedy2N54, DSGS, TOTAL and the rest of the package: this script has
# to run on a bare Julia, before anything is instantiated, because its whole job
# is to produce the mesh the package then reads. Same trick, and the same
# reason, as deck_mesh() in test/ci_cases.jl.
#---------------------------------------------------------------------------------
function deck_placement(file)
    for line in eachline(file)
        code = lstrip(line)
        startswith(code, "#") && continue
        m = match(r":naca64A210((?:\s*,\s*[-\d.eE+]+){5})\s*\)", code)
        m === nothing && continue
        return ntuple(i -> parse(Float64, strip(split(m.captures[1], ",")[i+1])), 5)
    end
    error("make_mesh.jl: no (:naca64A210, x_LE, y_LE, chord, α, r_TE) spec found in " * file *
          ".\n  The deck's :exact_geometry entry is what places the aerofoil; this " *
          "script reads it\n  so the mesh and the exact geometry cannot disagree.")
end

const XLE, YLE, CHORD, ALPHA, RTE = deck_placement(joinpath(CASE_DIR, "user_inputs.jl"))

# The tunnel. It lives here rather than in the deck because the deck never needs
# it: the grid IS the domain, and mod_mesh_read_gmsh! reads the bounds off it.
# These are shock_circle's tunnel exactly — same box, different obstacle.
const XMIN, XMAX =  0.0, 3.0
const YMIN, YMAX = -1.0, 1.0

#---------------------------------------------------------------------------------
# Resolution.
#
# NSURF is the number of points the aerofoil boundary is drawn with. They are
# NOT the mesh — gmsh meshes the spline through them under the size fields below
# — they are how finely the .geo states the curve. gmsh's `Spline` is Catmull-Rom
# and the section is the interpolating cubic spline, so the two agree at the
# sampled points and differ a little in between; the only requirement is that
# the difference be well below the correction the curving is there to make. At
# 241 cosine-clustered points check_mesh.jl measures the mesh vertices as 3.2e-6
# (5.4e-6 of chord) off the true section, against a wall-edge sagitta of ~1e-4
# at the nose — a thirtieth of it. Raise NSURF if you refine the surface a lot.
#
# The size fields are what set the actual element count:
#
#   H_LE    at the leading and trailing edges. The nose radius is 0.0056c, so
#           this is what decides whether the stagnation region is a curve or a
#           polygon, and it is the number to lower if the bow shock looks
#           faceted.
#   H_SURF  along the rest of the aerofoil.
#   H_FAR   in the far field, reached over a distance D_FAR.
#
# As shipped this is 2468 quads — the same order as shock_circle, i.e. a case
# that runs on a laptop. It is a demonstration, not a grid-converged aerofoil
# calculation. check_mesh.jl measures what it is worth: the wall-edge sagitta
# peaks at 3.2e-4 m at the nose (that is the error the curving removes) and the
# thinnest element there is still 9.4 times its own sagitta, so nothing folds.
#
# H_LE IS ALSO THE TIME STEP. The DynSGS diffusive stability limit is
# Δt ≲ ½Δx²/(μ⃗·C2·Δ·(|u|+c)) with Δx ∝ Δ ∝ H_LE, so Δt ∝ H_LE — halving the
# nose size halves the time step and doubles the cost of the run. See the :Δt
# note in user_inputs.jl, which carries the measurements. Refine the nose only
# as far as the fold margin actually needs.
#---------------------------------------------------------------------------------
const NSURF  = 241
# H_LE also sets the size at the trailing edge, and the trailing edge is now a
# circle of radius RTE·c. Keeping H_LE = RTE·c puts TWO elements across the base,
# which is the point of rounding it: a base one element wide is a base the
# boundary condition cannot see.
const H_LE   = 0.005*CHORD         # = RTE·c; nose radius 0.0056c gets ~3 elements
# 0.008c, not 0.020c. Aft of ~0.85c the aerofoil is thinner than 0.020c, so an
# element there was larger than the scale over which the upper and lower streams
# differ — the recompression and the shear layer leaving the trailing edge. A
# high-order scheme cannot resolve a feature smaller than its element and answers
# with oscillations, which is what drove p (and hence T) negative in the wake.
#
# This costs NOTHING in time step: Δt is set by the SMALLEST element, which is
# the nose at H_LE, and 0.008c is above it. Free resolution exactly where the
# case was failing.
const H_SURF = 0.008*CHORD
const H_FAR  = 0.100
const D_NEAR = 0.05*CHORD          # distance out to which H_SURF holds
const D_FAR  = 1.50*CHORD          # distance at which H_FAR is reached
const D_TIP  = 0.03*CHORD          # radius of the H_LE patch at LE and TE
# ...and the distance over which the tip patch relaxes back to the far-field
# size, where Field[2] takes over. Short, so the refinement stays a patch on the
# nose and the trailing edge rather than a halo around the whole aerofoil: at
# 1.5c it is the same grid with 50% more elements and the same fold margin.
const D_TIPFAR = 0.25*CHORD

# Wake band behind the trailing edge. The shear layer and the recompression
# shock system live here, and the graded field alone was letting the elements
# grow to 0.024 m by x = 1.75 — which is exactly where the first unphysical
# states appeared. Also free in Δt for the same reason as H_SURF.
const WAKE_H    = 0.010
const WAKE_X0   = XLE + CHORD                 # the trailing edge
const WAKE_X1   = XLE + CHORD + 1.4*CHORD
const WAKE_HALF = 0.20*CHORD
const WAKE_THK  = 0.15*CHORD                  # smooth transition out of the band

#---------------------------------------------------------------------------------
# Sample the section.
#
# Cosine clustering in the parameter, applied separately to the two surfaces so
# that both the nose and the trailing edge get the dense end of the distribution
# (the parameter is arclength from the upper trailing edge, so the leading edge
# sits in the middle of the range).
#---------------------------------------------------------------------------------
function sample_section(sec, nper)
    tle = sec.t[argmin(sec.x)]           # parameter at the leading edge
    # With a rounded trailing edge the SPLINE part of the section runs only
    # between the two tangency points; the rest is the arc, emitted separately
    # below as two gmsh Circle segments so the .geo carries the true circle.
    tlo = sec.te === nothing ? sec.t[1]   : sec.te.tu
    thi = sec.te === nothing ? sec.t[end] : sec.te.tl

    pts = Tuple{Float64,Float64}[]
    for (a, b, skipfirst) in ((tlo, tle, false), (tle, thi, true))
        for k = 0:nper
            skipfirst && k == 0 && continue
            # cosine: dense at both ends of each surface
            s = 0.5*(1.0 - cos(π*k/nper))
            t = a + (b - a)*s
            x, y, _, _, _, _ = naca64A210_point(sec, t)
            push!(pts, (x, y))
        end
    end
    return pts
end

const SEC = naca64A210_section(XLE, YLE, CHORD, ALPHA, RTE)
# First and last sample are the two tangency points, which are distinct: the
# arc joins them round the back.
const PTS = sample_section(SEC, (NSURF - 1) ÷ 2)

#---------------------------------------------------------------------------------
# Write the .geo.
#---------------------------------------------------------------------------------
open(joinpath(CASE_DIR, "naca64A210.geo"), "w") do io

    println(io, "// naca64A210.geo — GENERATED by make_mesh.jl. Do not edit by hand:")
    println(io, "// re-run  julia problems/CompEuler/naca64A210/make_mesh.jl  instead.")
    println(io, "//")
    println(io, "// A NACA 64A210 in a Mach-3 wind tunnel. The aerofoil is a hole in the")
    println(io, "// tunnel, drawn as two splines through points sampled from the section")
    println(io, "// defined in user_exactGeo.jl — which is what the solver curves the")
    println(io, "// high-order nodes onto (:exact_geometry in user_inputs.jl).")
    println(io, "//")
    @printf(io, "// section : NACA 64A210, chord %.6g, LE at (%.6g, %.6g), incidence %.6g deg\n",
            CHORD, XLE, YLE, ALPHA)
    @printf(io, "// TE      : ROUNDED, radius %.6g c = %.6g m (base %.6g m over h = %.6g)\n",
            RTE, SEC.te.r, 2*SEC.te.r, H_LE)
    @printf(io, "// tunnel  : [%.6g, %.6g] x [%.6g, %.6g]\n", XMIN, XMAX, YMIN, YMAX)
    @printf(io, "// aerofoil: %d points, splined; h = %.4g at the tips, %.4g on the surface\n",
            length(PTS), H_LE, H_SURF)
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

    println(io, "// ---- aerofoil --------------------------------------------------------")
    println(io, "// Points 100.. run from the UPPER trailing-edge tangency point forward")
    println(io, "// to the leading edge, then aft over the LOWER surface to the lower")
    println(io, "// tangency point. The two splines meet at the leading edge; the two")
    println(io, "// Circle arcs below close the loop round the rounded base. Every point")
    println(io, "// of this boundary has a well-defined normal, which a cusp does not.")
    for (k, p) in enumerate(PTS)
        @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};\n", 99 + k, p[1], p[2])
    end
    println(io)

    nle   = argmin([p[1] for p in PTS])         # index of the leading-edge point
    ile   = 99 + nle
    ilast = 99 + length(PTS)
    te    = SEC.te
    ic, ia = ilast + 1, ilast + 2               # arc centre, arc apex
    @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};   // trailing-edge arc centre\n",
            ic, te.cx, te.cy)
    @printf(io, "Point(%d) = {%.12g, %.12g, 0, hfar};   // trailing-edge apex\n",
            ia, te.xa, te.ya)
    println(io)
    @printf(io, "Spline(5) = {%d:%d};   // upper surface, tangency -> LE\n", 100, ile)
    @printf(io, "Spline(6) = {%d:%d};   // lower surface, LE -> tangency\n", ile, ilast)
    println(io, "// The trailing edge is a CIRCLE tangent to both surfaces, not a cusp.")
    println(io, "// Two arcs, not one: gmsh's built-in Circle spans less than pi, and")
    println(io, "// this one sweeps about 168 deg. They meet at the apex.")
    @printf(io, "Circle(7) = {%d, %d, %d};   // lower tangency -> apex\n", ilast, ic, ia)
    @printf(io, "Circle(8) = {%d, %d, %d};   // apex -> upper tangency\n", ia, ic, 100)
    println(io)

    println(io, "Curve Loop(1) = {1, 2, 3, 4};      // tunnel, counter-clockwise")
    println(io, "Curve Loop(2) = {5, 6, 7, 8};      // aerofoil, the hole")
    println(io, "Plane Surface(1) = {1, 2};")
    println(io)

    println(io, "// ---- element size ----------------------------------------------------")
    println(io, "// Fine on the aerofoil, finer still at the two tips, coarse far away.")
    println(io, "// Distance/Threshold rather than a background mesh so the grading follows")
    println(io, "// the geometry rather than the box.")
    println(io, "Field[1] = Distance;")
    println(io, "Field[1].CurvesList = {5, 6, 7, 8};")
    println(io, "Field[1].Sampling = 400;")
    println(io, "Field[2] = Threshold;")
    println(io, "Field[2].InField = 1;")
    @printf(io, "Field[2].SizeMin = %.10g;\n", H_SURF)
    @printf(io, "Field[2].SizeMax = %.10g;\n", H_FAR)
    @printf(io, "Field[2].DistMin = %.10g;\n", D_NEAR)
    @printf(io, "Field[2].DistMax = %.10g;\n", D_FAR)
    println(io)
    println(io, "// The leading edge (radius 0.0056c) and the trailing-edge arc, whose")
    println(io, "// base must carry at least a couple of elements across it.")
    println(io, "Field[3] = Distance;")
    @printf(io, "Field[3].PointsList = {%d, %d};\n", ile, ia)
    println(io, "Field[4] = Threshold;")
    println(io, "Field[4].InField = 3;")
    @printf(io, "Field[4].SizeMin = %.10g;\n", H_LE)
    @printf(io, "Field[4].SizeMax = %.10g;\n", H_FAR)
    @printf(io, "Field[4].DistMin = %.10g;\n", D_TIP)
    @printf(io, "Field[4].DistMax = %.10g;\n", D_TIPFAR)
    println(io)
    println(io)
    println(io, "// The wake. Everything the trailing edge sheds goes through here.")
    println(io, "Field[5] = Box;")
    @printf(io, "Field[5].VIn  = %.10g;\n", WAKE_H)
    @printf(io, "Field[5].VOut = %.10g;\n", H_FAR)
    @printf(io, "Field[5].XMin = %.10g;\n", WAKE_X0)
    @printf(io, "Field[5].XMax = %.10g;\n", WAKE_X1)
    @printf(io, "Field[5].YMin = %.10g;\n", -WAKE_HALF)
    @printf(io, "Field[5].YMax = %.10g;\n", WAKE_HALF)
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
    println(io, "// Jexpresso's 2D reader wants quadrilaterals, and it wants ALL of them:")
    println(io, "// a quad-dominant mesh with a few leftover triangles is not a mesh it can")
    println(io, "// read. Plain blossom (RecombinationAlgorithm = 1) leaves a few hundred")
    println(io, "// triangles behind on this geometry, so use the full-quad variant, which")
    println(io, "// splits whatever it cannot pair.")
    println(io, "Recombine Surface {1};")
    println(io, "Mesh.RecombinationAlgorithm = 3;   // blossom, full-quad")
    println(io, "Mesh.Algorithm = 8;                // frontal-Delaunay for quads")
    println(io)

    println(io, "// ---- physical groups -------------------------------------------------")
    println(io, "// These names are what user_bc.jl dispatches on, and \"airfoil\" is what")
    println(io, "// :exact_geometry names in user_inputs.jl.")
    println(io, "Physical Curve(\"bottom\",  1) = {1};")
    println(io, "Physical Curve(\"outflow\", 2) = {2};")
    println(io, "Physical Curve(\"top\",     3) = {3};")
    println(io, "Physical Curve(\"inflow\",  4) = {4};")
    println(io, "Physical Curve(\"airfoil\", 5) = {5, 6, 7, 8};")
    println(io, "Physical Surface(\"domain\") = {1};")
end

println("wrote ", joinpath(CASE_DIR, "naca64A210.geo"),
        "  (", length(PTS), " aerofoil points)")
println("now:  gmsh -2 ", joinpath(CASE_DIR, "naca64A210.geo"),
        " -format msh4 -o ", joinpath(CASE_DIR, "naca64A210.msh"))

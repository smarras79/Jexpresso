#---------------------------------------------------------------------------------
# user_exactGeo.jl — the ANALYTIC definition of this case's geometry.
#
# THIS IS THE FILE TO COPY when you want a curved boundary in your own case. It
# is deliberately the only place a circle appears anywhere in Jexpresso: the
# kernel (src/kernel/mesh/exact_geometry.jl) owns the curving — putting the
# boundary nodes on the curve, blending the correction into the touching
# elements with the Gordon-Hall transfinite map, checking that the grid still
# conforms and that no element folded — and owns no geometry at all. What shape
# a given wall IS cannot be known there; it is a property of the case, the way
# the flux, the source and the boundary conditions are.
#
# THE CONTRACT
# ------------
# The deck (user_inputs.jl) names the boundaries:
#
#     :exact_geometry => Dict("circle_boundary" => (:circle, 1.0, 0.0, 0.2)),
#
# — a gmsh `Physical Curve` name, and a SPEC, which is any value at all: the
# kernel never looks inside it, it only hands it back here. (A bare list,
# `:exact_geometry => ["circle_boundary"]`, is the same thing with spec =
# `nothing`, for a case that prefers to keep its numbers in this file.)
#
# This file then answers two questions:
#
#   user_exactGeo(tag, spec, x, y) -> (x, y)             REQUIRED
#       Given a point at or near boundary `tag`, where does it belong ON the
#       exact geometry? For a closed shape, the nearest point. Return (x, y)
#       unchanged to leave a node alone — which is also how you decline a tag
#       this case does not handle. Called for every node of every named
#       boundary, so keep it cheap and keep it PURE: two ranks holding the same
#       node must get bit-for-bit the same answer, or the grid tears along the
#       partition boundary.
#
#   user_exactGeo_setup(tag, spec, xs, ys) -> spec, or nothing   OPTIONAL
#       Called ONCE per boundary, before any node moves, with the LINEAR (gmsh
#       vertex) coordinates of that boundary gathered over all ranks — the
#       geometry the mesh file itself asserts. Return the spec to use (possibly
#       a different one, fitted from the grid), or `nothing` to leave the
#       boundary alone. Every rank calls it, so a collective is allowed here and
#       nowhere else. Omit it entirely and the deck's spec is used as written.
#
# WHY A NEAREST-POINT PROJECTION IS ENOUGH. Kopriva (J. Sci. Comput. 26(3):301,
# 2006) Theorem 3 asks only that the element mapping be in P^N. The mapping is
# the degree-N interpolant through the nodes, so it is in P^N however the nodes
# are placed along the curve — the theorem does not care WHICH parametrization
# of the arc they end up at, only that they lie on it. Free-stream preservation
# survives the radial snap below exactly.
#---------------------------------------------------------------------------------

#
# The geometry of this case. plate_hole_circle_unit.geo cuts a circular hole of
# radius 0.2 centred at (1, 0) out of the unit plate, and tags its boundary
# "circle_boundary". Those three numbers are what the deck passes in as
# (:circle, 1.0, 0.0, 0.2); they are read here rather than hardcoded so the same
# file serves a second circle without editing.
#
# Radial projection onto the circle: the 2D twin of the radial snap that
# project_nodes_to_shell! does on a spherical shell.
#
function user_exactGeo(tag, spec, x, y)

    if spec isa Tuple && length(spec) == 4 && spec[1] === :circle
        xc, yc, r = spec[2], spec[3], spec[4]
        dx = x - xc
        dy = y - yc
        d  = sqrt(dx*dx + dy*dy)
        d > 0 || return (x, y)              # node sits on the centre: nothing to do
        s = r/d
        return (xc + s*dx, yc + s*dy)
    end

    # A tag, or a spec, this case does not know: leave the node where it is.
    return (x, y)
end

#
# Resolve the deck's spec before any node is moved.
#
#   (:circle, xc, yc, r)   taken as written. This is what the deck uses, because
#                          :linitial_refine is on — see below.
#
#   :circle                centre and radius FITTED from the boundary vertices
#                          the mesh file already carries, and REFUSED if they do
#                          not actually lie on a circle. Same "trust the grid,
#                          but check it" policy as project_nodes_to_shell!.
#
# je_fit_circle (src/kernel/mesh/exact_geometry.jl) is the least-squares fit;
# it is a kernel utility rather than a geometry because it is an MPI collective,
# and it returns the largest radial residual relative to r precisely so that the
# caller can refuse a bad fit instead of deforming the wall onto a shape nobody
# asked for.
#
# WHY THE FIT CANNOT BE USED ON A REFINED GRID. Refinement happens on the
# straight-sided Gridap model and the curving runs on the result, so every new
# boundary vertex arrives at the MIDPOINT OF A CHORD — a sagitta (3.8e-3 here at
# :init_refine_lvl => 1) inside the circle. Least squares then lands exactly
# between the two clusters, the on-circle vertices and the pulled-in ones, at
# the same residual from both, so there is no residual-based way to tell them
# apart and guessing would be worse than refusing. Refusing is what happens
# below; the deck states the circle instead, and the snap fixes the pulled-in
# vertices along with everything else (it moves vertices too), so the wall comes
# out at ~3e-16 from the circle, refined or not.
#
function user_exactGeo_setup(tag, spec, xs, ys)

    spec === :circle || return spec

    fit = je_fit_circle(xs, ys)
    if fit === nothing
        @warn string(":exact_geometry \"", tag, "\" => :circle ignored: its vertices are ",
                     "collinear or too few to define a circle.")
        return nothing
    end

    xc, yc, r, rresid = fit
    if rresid > 1.0e-6
        @warn string(":exact_geometry \"", tag, "\" => :circle REFUSED: the boundary vertices ",
                     "are not on a circle (best fit centre = (", xc, ", ", yc, "), r = ", r,
                     " leaves a relative residual of ", rresid, "). Nothing was curved. ",
                     "If this grid is h-refined (:linitial_refine / :ladapt), that is expected — ",
                     "refinement puts the new boundary vertices on the chords, so the grid no ",
                     "longer states the circle. Give it outright instead: ",
                     ":exact_geometry => Dict(\"", tag, "\" => (:circle, xc, yc, r)), ",
                     "which works on a refined grid because the vertices are snapped too.")
        return nothing
    end

    return (:circle, xc, yc, r)
end

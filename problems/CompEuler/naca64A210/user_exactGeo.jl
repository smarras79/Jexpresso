#---------------------------------------------------------------------------------
# user_exactGeo.jl — the NACA 64A210 section, analytically.
#
# This is the second worked example of the :exact_geometry hook (the first,
# problems/CompEuler/shock_circle/user_exactGeo.jl, is a circle) and it is the
# one that shows why the geometry cannot live in the kernel.
#
# A NACA 6A-series section has NO closed form. Unlike the 4- and 5-digit
# families, whose thickness distribution is a polynomial in sqrt(x), the
# 6-series sections were derived by conformal mapping from a prescribed pressure
# distribution and are DEFINED by their published table of ordinates (Abbott &
# von Doenhoff, "Theory of Wing Sections", Dover 1959, Appendix IV). The exact
# geometry is therefore "the smooth curve through those points", and every
# aerodynamicist means the same thing by that: the interpolating cubic spline.
#
# So this file is:
#
#   1. the published ordinates,                       NACA64A210_ORDINATES
#   2. a natural cubic spline through them, in arclength,   _spline / _section
#   3. the nearest point on that spline,               user_exactGeo
#
# and nothing about it could have been anticipated by src/kernel/mesh/
# exact_geometry.jl, which owns the curving and owns no geometry: the
# isoparametric boundary, the Gordon-Hall blend, the conformity guard and the
# fold check are all shape-independent, and this is the shape.
#
# WHAT THE CURVING BUYS HERE. The grid is straight-sided: naca64A210.msh puts
# its vertices on the section and every high-order node in between on the CHORD.
# At the leading edge the radius of curvature is 0.0056c, so a boundary edge
# there is a chord across a strongly curved arc and its midpoint sits ~1e-4 c
# INSIDE the nose — on an aerofoil at Mach 3 that is a row of small forward-
# facing steps around the stagnation point, each shedding its own entropy. The
# snap puts every boundary node on the true section and blends the correction
# into the touching elements, so the wall the solver integrates over is the
# degree-:nop interpolant of the real aerofoil.
#
# THE CONTRACT (see src/kernel/mesh/exact_geometry.jl for the full statement):
#
#   user_exactGeo(tag, spec, x, y) -> (x, y)                       REQUIRED
#       the analytic definition: where does this point belong on the section?
#
#   user_exactGeo_setup(tag, spec, xs, ys) -> spec, or false        OPTIONAL
#       called ONCE per boundary, before any node moves, with that boundary's
#       linear (gmsh vertex) cloud. Here it does the two things the hook exists
#       for: it BUILDS the spline once, so the per-node call is cheap, and it
#       CHECKS that the grid's own vertices really are on the section before
#       anything is moved onto it. Returning `false` declines the boundary —
#       `false` rather than `nothing`, since `nothing` is itself a valid spec.
#
# The deck (user_inputs.jl) names the boundary and states the placement:
#
#   :exact_geometry => Dict("airfoil" => (:naca64A210, 1.0, 0.0, 0.6, 0.0)),
#                                        kind        x_LE y_LE chord α [deg]
#
# make_mesh.jl reads the SAME tuple to write the .geo, so the mesh and the exact
# geometry cannot drift apart: change the chord or the incidence in the deck,
# re-run make_mesh.jl, and both follow.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# THE SECTION, as published.
#
# NACA 64A210: 6-series, minimum pressure at 0.4c on the basic symmetrical
# section, design lift coefficient 0.2 on the a = 0.8 (modified) mean line, 10%
# thick, "A" = the 6A-series thickness modification (straight surfaces aft of
# 0.8c, a trailing-edge angle small enough to build).
#
# 51 points, in the standard aerofoil ordering: trailing edge, forward along the
# UPPER surface to the leading edge at (0,0), then aft along the LOWER surface
# back to the trailing edge. x and y are fractions of chord.
#
# Source: the NACA 64A210 ordinates as distributed in the UIUC Airfoil
# Coordinate Database (naca64a210.dat), which reproduces the Abbott & von
# Doenhoff table. Ordinates are measurements of a standard section, not
# something to be recomputed — copy them, do not "improve" them.
#
# WHAT THIS TABLE DOES AND DOES NOT PIN DOWN. Splined, it reproduces the section
# where the section is smooth: maximum thickness comes out at t/c = 0.0998 at
# x/c = 0.39, which is the "10" of 64A2-10 and the 6-series' 0.4c thickness
# station. The NOSE is the exception. Abbott & von Doenhoff quote a leading-edge
# radius of 0.00687c, and they quote it separately precisely BECAUSE the table
# cannot carry it: the three points nearest the leading edge here are at
# x/c = 0.0042, 0, 0.0058, and a cubic spline through points that far apart
# realises a nose radius of 0.0056c — 18% under the published figure. That is a
# property of the standard 51-point table, not of this file, and it is what any
# mesh generator handed the same .dat produces. It is the geometry used
# throughout: check_mesh.jl's sagitta, the element sizing in make_mesh.jl and
# the estimates in the deck are all quoted against 0.0056c, not the published
# value. Fairing the nose to the published radius means adding points to the
# table above — a deliberate change to the geometry, not a tweak.
#---------------------------------------------------------------------------------
const NACA64A210_ORDINATES = [
    1.00000   0.00021
    0.95027   0.00785
    0.90052   0.01551
    0.85074   0.02301
    0.80076   0.03037
    0.75063   0.03702
    0.70054   0.04310
    0.65042   0.04852
    0.60028   0.05323
    0.55012   0.05714
    0.49994   0.06014
    0.44975   0.06208
    0.39955   0.06274
    0.34935   0.06192
    0.29917   0.05984
    0.24900   0.05656
    0.19885   0.05200
    0.14874   0.04592
    0.09868   0.03792
    0.07369   0.03288
    0.04874   0.02685
    0.02387   0.01895
    0.01153   0.01342
    0.00665   0.01044
    0.00424   0.00856
    0.00000   0.00000
    0.00576  -0.00744
    0.00835  -0.00886
    0.01347  -0.01100
    0.02613  -0.01473
    0.05126  -0.01963
    0.07631  -0.02316
    0.10132  -0.02600
    0.15126  -0.03030
    0.20115  -0.03340
    0.25100  -0.03554
    0.30083  -0.03688
    0.35065  -0.03744
    0.40045  -0.03716
    0.45025  -0.03580
    0.50006  -0.03354
    0.54988  -0.03062
    0.59972  -0.02719
    0.64958  -0.02342
    0.69946  -0.01944
    0.74937  -0.01542
    0.79924  -0.01167
    0.84926  -0.00859
    0.89948  -0.00571
    0.94974  -0.00295
    1.00000  -0.00021
]

# The published section has a trailing edge 0.042% of chord thick. Closed here
# by the usual linear taper — each surface is pulled toward the camber line in
# proportion to x, so the leading edge is untouched and the two surfaces meet at
# (1, 0). It removes 2.1e-4 c of thickness at the trailing edge and less
# everywhere else, i.e. nothing an inviscid solution can see, and it buys a
# geometry that is a single closed curve with one corner instead of a blunt base
# whose two sharp corners would each want their own resolution.
const NACA64A210_TE_HALF_GAP = 0.00021


#---------------------------------------------------------------------------------
# Natural cubic spline through (t[i], v[i]), returned as the vector of second
# derivatives M. Natural = M vanishes at both ends.
#
# Thomas algorithm on the standard tridiagonal system
#
#   h[i-1]/6 M[i-1] + (h[i-1]+h[i])/3 M[i] + h[i]/6 M[i+1]
#       = (v[i+1]-v[i])/h[i] - (v[i]-v[i-1])/h[i-1]
#
# (Press et al., Numerical Recipes, §3.3). Written out rather than pulled from a
# package because it has to give the SAME answer here and in make_mesh.jl, and
# because a spline is four lines of arithmetic.
#---------------------------------------------------------------------------------
function _spline_second_derivatives(t::Vector{Float64}, v::Vector{Float64})
    n = length(t)
    M = zeros(Float64, n)
    u = zeros(Float64, n)
    for i = 2:n-1
        hm = t[i] - t[i-1]
        hp = t[i+1] - t[i]
        sig = hm/(hm + hp)
        p   = sig*M[i-1] + 2.0
        M[i] = (sig - 1.0)/p
        u[i] = ((v[i+1] - v[i])/hp - (v[i] - v[i-1])/hm)
        u[i] = (6.0*u[i]/(hm + hp) - sig*u[i-1])/p
    end
    for i = n-1:-1:1
        M[i] = M[i]*M[i+1] + u[i]
    end
    return M
end

# Value and first two derivatives of the spline on the interval that contains t.
# `i` is the knot index, found by the caller, so the search is not repeated for
# x and y.
@inline function _spline_eval(t, ts, v, M, i)
    h = ts[i+1] - ts[i]
    A = (ts[i+1] - t)/h
    B = (t - ts[i])/h
    f   = A*v[i] + B*v[i+1] + ((A^3 - A)*M[i] + (B^3 - B)*M[i+1])*h*h/6.0
    df  = (v[i+1] - v[i])/h + ((3.0*B*B - 1.0)*M[i+1] - (3.0*A*A - 1.0)*M[i])*h/6.0
    d2f = A*M[i] + B*M[i+1]
    return f, df, d2f
end

@inline function _spline_interval(t, ts)
    n = length(ts)
    t <= ts[1]   && return 1
    t >= ts[n]   && return n - 1
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        ts[mid] > t ? (hi = mid) : (lo = mid)
    end
    return lo
end


#---------------------------------------------------------------------------------
# naca64A210_section(xle, yle, chord, α_deg) -> the curve, as a NamedTuple
#
# The published ordinates, trailing edge closed, scaled by `chord`, rotated by
# `α_deg` NOSE UP about the leading edge, and placed with the leading edge at
# (xle, yle). Splined against cumulative chord length, so the parameter is very
# nearly arclength and the round leading edge — where dy/dx is infinite and a
# y(x) spline would fail — is no different from anywhere else. The two ends of
# the parameter range are the two sides of the sharp trailing edge, which is
# exactly why the curve is left OPEN rather than made periodic: a periodic
# spline would round off the one corner the section actually has.
#
# A NamedTuple rather than a struct so that re-including this file (which
# src/run.jl does whenever the case is edited) can never collide with a type it
# defined on a previous include.
#---------------------------------------------------------------------------------
function naca64A210_section(xle::Real, yle::Real, chord::Real, α_deg::Real)

    xr = NACA64A210_ORDINATES[:, 1]
    yr = NACA64A210_ORDINATES[:, 2]
    n  = length(xr)

    # Close the trailing edge: the upper surface (the first half of the table,
    # down to the leading edge at index ile) comes down, the lower surface goes
    # up, both in proportion to x.
    ile = argmin(xr)
    y   = similar(yr)
    for i = 1:n
        y[i] = yr[i] + (i <= ile ? -1.0 : 1.0)*NACA64A210_TE_HALF_GAP*xr[i]
    end

    # Scale, rotate nose-up by α, translate. Nose-up means the SECTION turns,
    # which for a wind tunnel with a horizontal free stream is the same thing as
    # incidence, and keeps the free stream aligned with the tunnel walls.
    α  = α_deg*π/180.0
    ca = cos(α); sa = sin(α)
    px = Vector{Float64}(undef, n)
    py = Vector{Float64}(undef, n)
    for i = 1:n
        xs = chord*xr[i]
        ys = chord*y[i]
        px[i] = xle + ca*xs + sa*ys
        py[i] = yle - sa*xs + ca*ys
    end

    # Cumulative chord length as the spline parameter.
    ts = Vector{Float64}(undef, n)
    ts[1] = 0.0
    for i = 2:n
        ts[i] = ts[i-1] + hypot(px[i] - px[i-1], py[i] - py[i-1])
    end

    return (t = ts, x = px, y = py,
            Mx = _spline_second_derivatives(ts, px),
            My = _spline_second_derivatives(ts, py),
            chord = Float64(chord))
end

# A point on the section, and the tangent there.
@inline function naca64A210_point(sec, t)
    i = _spline_interval(t, sec.t)
    x, dx, d2x = _spline_eval(t, sec.t, sec.x, sec.Mx, i)
    y, dy, d2y = _spline_eval(t, sec.t, sec.y, sec.My, i)
    return x, y, dx, dy, d2x, d2y
end


#---------------------------------------------------------------------------------
# THE ANALYTIC DEFINITION. Nearest point of the section to (x, y).
#
# Minimise f(t) = |P(t) - q|² over the parameter range, i.e. solve the
# orthogonality condition
#
#     g(t) = (P(t) - q) · P'(t) = 0,      g'(t) = |P'|² + (P - q) · P''.
#
# A coarse sweep brackets the global minimum — necessary, because f has a local
# minimum against every part of the section and Newton from a bad start would
# happily project a node on the upper surface onto the lower one — and Newton
# then converges on it. The sweep is a fixed sampling, so this is deterministic
# to the last bit: two MPI ranks that own the same node get the same answer, and
# the grid cannot tear along a partition boundary.
#
# Clamped to the parameter range, so a point beyond the trailing edge projects
# onto the trailing edge itself rather than running off the end of the spline.
#
# Kopriva (J. Sci. Comput. 26(3):301, 2006) Theorem 3 asks only that the element
# mapping be in P^N. It is the degree-N interpolant through the nodes, so it is
# in P^N wherever ON the curve they land — the theorem does not care about the
# parametrization, only that the nodes lie on the section. A nearest-point
# projection is simply the choice that moves them least.
#---------------------------------------------------------------------------------
const _NACA_SWEEP_PER_SEGMENT = 6        # sweep samples per published ordinate
const _NACA_NEWTON_ITERS      = 40

function user_exactGeo(tag, spec, x, y)

    spec isa NamedTuple && haskey(spec, :Mx) || return (x, y)   # not ours

    sec  = spec
    ts   = sec.t
    tlo  = ts[1]
    thi  = ts[end]

    # Coarse sweep for the global minimum.
    nseg = (length(ts) - 1)*_NACA_SWEEP_PER_SEGMENT
    tbest = tlo
    fbest = Inf
    for k = 0:nseg
        t = tlo + (thi - tlo)*k/nseg
        px, py, _, _, _, _ = naca64A210_point(sec, t)
        f = (px - x)^2 + (py - y)^2
        if f < fbest
            fbest = f; tbest = t
        end
    end

    # Newton on g(t) = 0, confined to one sweep cell either side of the best
    # sample: that is where the minimum is, and clamping there stops a large
    # Newton step from jumping to a different branch of the surface.
    dt   = (thi - tlo)/nseg
    lo   = max(tlo, tbest - dt)
    hi   = min(thi, tbest + dt)
    t    = tbest
    for _ = 1:_NACA_NEWTON_ITERS
        px, py, dx, dy, d2x, d2y = naca64A210_point(sec, t)
        g  = (px - x)*dx + (py - y)*dy
        gp = dx*dx + dy*dy + (px - x)*d2x + (py - y)*d2y
        abs(gp) > 0.0 || break
        tn = t - g/gp
        tn = tn < lo ? lo : (tn > hi ? hi : tn)
        abs(tn - t) <= 1.0e-15*max(1.0, abs(t)) && (t = tn; break)
        t = tn
    end

    px, py, _, _, _, _ = naca64A210_point(sec, t)
    return (px, py)
end


#---------------------------------------------------------------------------------
# Build the spline once, and check the grid before anything is moved onto it.
#
# The check is the aerofoil's version of the circle fit's residual test in
# shock_circle: if the mesh vertices are NOT on the section this deck describes
# — the wrong .msh, the wrong chord, the wrong incidence, a section edited in
# one place and not the other — then snapping to it would drag the whole
# boundary somewhere nobody asked for, and every symptom of that appears much
# later and looks like a physics bug. So measure it here, where the numbers can
# still be read, and refuse rather than deform.
#
# The tolerance is generous on purpose (1% of chord): the vertices are supposed
# to be ON the section, so this is a "did you hand me the right mesh" test, not
# an accuracy test. The number it prints is the accuracy test, and on
# naca64A210.msh it is round-off.
#---------------------------------------------------------------------------------
function user_exactGeo_setup(tag, spec, xs, ys)

    (spec isa Tuple && length(spec) == 5 && spec[1] === :naca64A210) || return spec

    sec = naca64A210_section(spec[2], spec[3], spec[4], spec[5])

    # How far the grid's own vertices are from the section. Local only: the
    # kernel hands each rank its own share, and a rank that owns none of this
    # boundary simply contributes nothing.
    worst = 0.0
    for k in eachindex(xs)
        px, py = user_exactGeo(tag, sec, xs[k], ys[k])
        worst  = max(worst, hypot(px - Float64(xs[k]), py - Float64(ys[k])))
    end

    if worst > 0.01*sec.chord
        @warn string(":exact_geometry \"", tag, "\" REFUSED: the boundary vertices of this ",
                     "grid are up to ", worst, " (", 100*worst/sec.chord, "% of chord) away ",
                     "from the NACA 64A210 section at x_LE = ", spec[2], ", y_LE = ", spec[3],
                     ", chord = ", spec[4], ", α = ", spec[5], " deg. That is not this ",
                     "aerofoil, so nothing was curved. Regenerate the mesh from the SAME ",
                     "placement with make_mesh.jl, which reads this very tuple out of ",
                     "user_inputs.jl.")
        return false
    end

    return sec
end

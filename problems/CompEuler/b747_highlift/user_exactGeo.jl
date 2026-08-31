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
#   :exact_geometry => Dict("airfoil" => (:naca64A210, 1.0, 0.0, 0.6, 0.0, 0.005)),
#                                        kind        x_LE y_LE chord α  r_TE
#
# r_TE is the trailing-edge radius as a fraction of chord; 0 leaves the published
# cusp. See _naca_te_round below for why this case does not use 0.
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
# the parameter range are the two sides of the trailing edge, which is
# exactly why the curve is left OPEN rather than made periodic: a periodic
# spline would round off the one corner the section actually has.
#
# A NamedTuple rather than a struct so that re-including this file (which
# src/run.jl does whenever the case is edited) can never collide with a type it
# defined on a previous include.
#---------------------------------------------------------------------------------
function naca64A210_section(xle::Real, yle::Real, chord::Real, α_deg::Real,
                            te_radius_c::Real = 0.0)

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

    sec = (t = ts, x = px, y = py,
           Mx = _spline_second_derivatives(ts, px),
           My = _spline_second_derivatives(ts, py),
           chord = Float64(chord), te = nothing)

    te_radius_c > 0 || return sec
    return (t = sec.t, x = sec.x, y = sec.y, Mx = sec.Mx, My = sec.My,
            chord = sec.chord,
            te = _naca_te_round(sec, Float64(te_radius_c)*Float64(chord)))
end


#---------------------------------------------------------------------------------
# ROUNDING THE TRAILING EDGE.
#
# WHY. The published section ends in a cusp: the two surfaces meet at one point
# with an included angle of about 12°. A single point is a place where the
# boundary has NO normal — the two edges meeting there carry normals almost
# 170° apart — and the free-slip condition is applied edge by edge, so that node
# gets both projections in succession and neither is the wall it is on. It is
# the origin of the oscillations that spread along the aft surface and into the
# wake, and no amount of refinement removes it: the thinner the elements, the
# sharper the wedge they have to represent.
#
# WHAT. Truncate the section where a circle of radius `r` is tangent to BOTH
# surfaces, and cap it with that circle. The result is G1-continuous: the normal
# rotates smoothly through the whole 180° around the base, every boundary node
# has one well-defined normal, and there is no corner anywhere on the aerofoil.
#
# WHY ROUND AND NOT SQUARE. A blunt (squared-off) base is the other way to give
# the trailing edge a finite thickness, and it is worse here: it replaces one
# 12° cusp with TWO corners of about 90°, and at a 90° corner the two successive
# slip projections remove BOTH momentum components — a spurious stagnation point
# planted in the middle of the recompression. That is the failure documented at
# the convex step corner in problems/CompEuler/ffs_step/user_bc.jl. A circle has
# no corner to get wrong.
#
# HOW. Find the centre C that is exactly `r` from each surface. Newton on the
# pair of distance equations: with N_u, N_l the nearest points on the upper and
# lower surface and n_u, n_l the unit vectors from them to C,
#
#     n_u·δ = r - |C - N_u| ,    n_l·δ = r - |C - N_l|
#
# is the linearisation, because moving C by δ changes its distance to a surface
# by n·δ to first order. Two equations, two unknowns, quadratic convergence.
# The tangency points N_u, N_l are then where the spline hands over to the arc.
#---------------------------------------------------------------------------------
function _naca_te_round(sec, r::Float64)

    ts   = sec.t
    tlo, thi = ts[1], ts[end]
    tle  = ts[argmin(sec.x)]                 # parameter at the leading edge

    # Start from the point midway between the surfaces where the section is
    # about 2r thick. Bisect on the parameter offset from the trailing edge.
    lo, hi = 0.0, 0.45*(tle - tlo)
    cx = cy = 0.0
    for _ = 1:60
        d = 0.5*(lo + hi)
        xu, yu, = naca64A210_point(sec, tlo + d)
        xl, yl, = naca64A210_point(sec, thi - d)
        hypot(xu - xl, yu - yl) < 2r ? (lo = d) : (hi = d)
        cx = 0.5*(xu + xl); cy = 0.5*(yu + yl)
    end

    tu = tlo; tl = thi
    for _ = 1:60
        tu, xu, yu, du = _nearest_on_spline(sec, cx, cy, tlo, tle)
        tl, xl, yl, dl = _nearest_on_spline(sec, cx, cy, tle, thi)
        du = sqrt(du); dl = sqrt(dl)
        (du > 0 && dl > 0) || error(" # ERROR user_exactGeo.jl: the trailing-edge " *
                                    "rounding centre landed on the section.")
        nux, nuy = (cx - xu)/du, (cy - yu)/du
        nlx, nly = (cx - xl)/dl, (cy - yl)/dl
        det = nux*nly - nuy*nlx
        abs(det) > 1.0e-12 ||
            error(" # ERROR user_exactGeo.jl: :te_radius = " * string(r) *
                  " is too large for this section — the two surfaces are parallel " *
                  "where the tangent circle would have to sit. Use a smaller radius.")
        bu, bl = r - du, r - dl
        δx = ( bu*nly - bl*nuy)/det
        δy = (-bu*nlx + bl*nux)/det
        cx += δx; cy += δy
        hypot(δx, δy) <= 1.0e-14*max(1.0, r) && break
    end

    xu, yu, = naca64A210_point(sec, tu)
    xl, yl, = naca64A210_point(sec, tl)
    θu = atan(yu - cy, xu - cx)
    θl = atan(yl - cy, xl - cx)

    # The arc runs from the LOWER tangency point round the BACK to the upper
    # one. "Back" is the direction from the leading edge toward the trailing
    # edge, so pick the sweep whose midpoint lies on that side.
    xle_, yle_ = naca64A210_point(sec, tle)[1], naca64A210_point(sec, tle)[2]
    dirx, diry = cx - xle_, cy - yle_
    dn = hypot(dirx, diry); dirx /= dn; diry /= dn
    Δ⁺ = mod(θu - θl, 2π)                       # counter-clockwise sweep
    for (Δ, s) in ((Δ⁺, 1.0), (Δ⁺ - 2π, -1.0))
        θm = θl + 0.5*Δ
        if cos(θm)*dirx + sin(θm)*diry > 0      # midpoint is downstream: this one
            return (cx = cx, cy = cy, r = r, tu = tu, tl = tl,
                    θl = θl, Δθ = Δ,
                    xu = xu, yu = yu, xl = xl, yl = yl,
                    xa = cx + r*dirx, ya = cy + r*diry)   # the apex, for the .geo
        end
    end
    error(" # ERROR user_exactGeo.jl: could not orient the trailing-edge arc.")
end


# Nearest point on the arc, or on whichever end of it is closer when the query
# falls outside its angular span.
@inline function _naca_nearest_on_arc(te, x, y)
    dx = x - te.cx; dy = y - te.cy
    d  = hypot(dx, dy)
    if d > 0
        θ = atan(dy, dx)
        s = mod((θ - te.θl)*sign(te.Δθ), 2π)
        if s <= abs(te.Δθ)
            px = te.cx + te.r*dx/d; py = te.cy + te.r*dy/d
            return px, py, (px - x)^2 + (py - y)^2
        end
    end
    fu = (te.xu - x)^2 + (te.yu - y)^2
    fl = (te.xl - x)^2 + (te.yl - y)^2
    return fu <= fl ? (te.xu, te.yu, fu) : (te.xl, te.yl, fl)
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
# Sweep resolution, and why it is not "a few samples per ordinate".
#
# THE TRAILING EDGE IS WHY. The sweep exists to pick the branch of the section
# that Newton then converges on, and Newton, once started, stays on its branch.
# Near a sharp trailing edge the two surfaces come within a few 1e-5 m of each
# other while the sweep spacing was 4e-3 m — so the sample nearest a node on the
# LOWER surface was routinely a point on the UPPER one, and the projection
# happily moved that node 4.1e-4 m ACROSS the trailing edge. That is bigger than
# the 1.7e-4 m sagitta the whole exercise is there to remove, it crosses the
# elements at the trailing edge, and a crossed element is a negative Jacobian,
# which is a NaN on the first time step.
#
# Two changes, and BOTH are needed:
#
#   * a finer sweep, so the true branch is sampled at all near the trailing
#     edge, and
#   * Newton from EVERY local minimum of the swept distance, keeping the
#     closest result — not from the single best sample. A point on the lower
#     surface has a local minimum on its own branch with distance ~0, which
#     beats anything on the other branch no matter which sample happened to be
#     nearest. This is what makes the answer independent of the sampling, and it
#     costs a handful of extra Newton solves on a few hundred boundary nodes,
#     once, at grid-build time.
#
# The two parameter ends are always candidates as well: they are the two sides
# of the trailing edge itself, and a node just downstream of it projects there.
const _NACA_SWEEP_PER_SEGMENT = 100      # sweep samples per published ordinate
const _NACA_NEWTON_ITERS      = 40

# Newton on g(t) = (P(t) - q)·P'(t) = 0, confined to [lo, hi] so a large step
# cannot jump to another branch. Returns the point and its squared distance.
@inline function _naca_refine(sec, x, y, t0, lo, hi)
    t = t0
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
    return t, px, py, (px - x)^2 + (py - y)^2
end

# Nearest point on the spline, restricted to the parameter window [ta, tb], and
# the parameter it sits at. Newton from EVERY local minimum of the sweep — see
# the note above _NACA_SWEEP_PER_SEGMENT for why the single-best-sample version
# is not good enough near a thin trailing edge.
function _nearest_on_spline(sec, x, y, ta, tb)

    ts   = sec.t
    span = ts[end] - ts[1]
    nseg = max(16, ceil(Int, (tb - ta)/span*(length(ts) - 1)*_NACA_SWEEP_PER_SEGMENT))
    dt   = (tb - ta)/nseg

    bt = ta; bx = 0.0; by = 0.0; bf = Inf
    @inline function consider(t0)
        lo = t0 - dt < ta ? ta : t0 - dt
        hi = t0 + dt > tb ? tb : t0 + dt
        t, px, py, f = _naca_refine(sec, x, y, t0, lo, hi)
        if f < bf
            bf = f; bx = px; by = py; bt = t
        end
        return nothing
    end

    fm2 = Inf; fm1 = Inf; tm1 = ta
    for k = 0:nseg
        t = ta + (tb - ta)*k/nseg
        px, py, _, _, _, _ = naca64A210_point(sec, t)
        f = (px - x)^2 + (py - y)^2
        fm1 <= fm2 && fm1 <= f && consider(tm1)
        fm2 = fm1; fm1 = f; tm1 = t
    end
    consider(ta); consider(tb)
    fm1 <= fm2 && consider(tm1)

    return bt, bx, by, bf
end


function user_exactGeo(tag, spec, x, y)

    spec isa NamedTuple && haskey(spec, :Mx) || return (x, y)   # not ours

    sec = spec
    te  = sec.te

    # With a rounded trailing edge the section is the spline BETWEEN the two
    # tangency points plus the arc that caps it; the spline outside that window
    # is the cusp the rounding removed and must not be projected onto.
    ta, tb = te === nothing ? (sec.t[1], sec.t[end]) : (te.tu, te.tl)
    _, bx, by, bf = _nearest_on_spline(sec, x, y, ta, tb)

    if te !== nothing
        ax, ay, af = _naca_nearest_on_arc(te, x, y)
        af < bf && ((bx, by) = (ax, ay))
    end

    return (bx, by)
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

    (spec isa Tuple && length(spec) == 6 && spec[1] === :naca64A210) || return spec

    sec = naca64A210_section(spec[2], spec[3], spec[4], spec[5], spec[6])

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
                     ", chord = ", spec[4], ", α = ", spec[5], " deg, r_TE = ", spec[6],
                     ". That is not this ",
                     "aerofoil, so nothing was curved. Regenerate the mesh from the SAME ",
                     "placement with make_mesh.jl, which reads this very tuple out of ",
                     "user_inputs.jl.")
        return false
    end

    return sec
end

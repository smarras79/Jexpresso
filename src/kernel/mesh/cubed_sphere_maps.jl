#---------------------------------------------------------------------------------
# cubed_sphere_maps.jl
#
# Re-position the nodes of a cubed-sphere shell by changing the map that carries
# the face of the inscribed cube onto the sphere, WITHOUT touching connectivity.
#
# The grid Jexpresso reads (problems/ShallowWater/SWsphere/cubed_sphere.geo, or
# any equivalent .msh) is a classical GNOMONIC cubed sphere: the six panels are
# the central projection of the faces of an inscribed cube (Sadourny 1972,
# Mon. Wea. Rev. 100, 136). That map is the simplest one but the worst
# conditioned — the panel corners carry the whole of its distortion.
#
# What this file adds is a POST-READ REMAP. Every node is
#
#   1. assigned to the panel whose outward face normal its position points at;
#   2. read back to its logical coordinate (u, v) ∈ [-1,1]² on that cube face,
#      which is the gnomonic inverse — see remap_direction on why the source
#      map is fixed rather than selectable;
#   3. pushed forward through a DIFFERENT map to give the new position.
#
# Step 2/3 leave (u, v) — the node's LOGICAL identity — alone, so the element
# topology, the panel decomposition and the panel boundaries are all unchanged.
# Only where the nodes sit inside each panel changes. That is what makes this
# safe to do after mod_mesh_read_gmsh! has already built the connectivity.
#
# THE MAPS
#
#   :gnomonic     (default source) equidistant central projection. The cube face
#                 point (u, v) ∈ [-1,1]² goes to (u, v, 1)/‖(u, v, 1)‖ in the
#                 panel frame. Sadourny (1972).
#
#   :equiangular  the same projection but with the face coordinate measured as
#                 an ANGLE: u = tan(α), α ∈ [-π/4, π/4]. Ronchi, Iacono &
#                 Paolucci (1996), JCP 124, 93. Equal angular increments give a
#                 far more even node spacing along a panel edge than equal
#                 distances on the cube face do; this is the variant most
#                 modern cubed-sphere codes actually run.
#
#   :conformal    the conformal cube-to-sphere map of Rančić, Purser & Mesinger
#                 (1996), QJRMS 122, 959-982. Angles are preserved everywhere
#                 except at the eight cube corners, so the grid is locally
#                 orthogonal — the property that makes the metric terms
#                 diagonal. The price is a strongly non-uniform Jacobian: cells
#                 near the corners are much smaller than at the panel centre,
#                 which is the opposite trade-off from :equiangular.
#
# THE CONFORMAL MAP, and where its numbers come from
#
# The forward map has no closed form. RPM96 give it as a truncated Taylor
# series whose 30 coefficients are their Table B1, reproduced verbatim in
# A_RANCIC below; the algorithm that wraps the series (fold to the first
# quadrant, Z = ((xᶜ + i yᶜ)/2)⁴, the cube-root/Möbius pair, stereographic
# lift) is theirs as well, and reached this file through the lineage
# Purser & Rančić Fortran 77 → MITgcm matlab-grid-generator/map_xy2xyz.m →
# CubedSphere.jl. The inverse series B is NOT tabulated anywhere: it is the
# formal reversion of A, computed here by _reverted_series (a Lagrange
# inversion), which is checked against A ∘ B = identity in the unit tests.
#
# All of this operates on the UNIT sphere; the caller rescales by mesh.radius.
#
# NOTHING HERE DEPENDS ON THE REST OF JEXPRESSO — no St_mesh, no MPI, no
# KernelAbstractions — so test/test_cubed_sphere_maps.jl can `include` this file
# on its own and check the maps without building the package. The driver that
# walks a St_mesh and applies a remap lives with its sibling
# project_nodes_to_shell! in mesh.jl. Keep it that way.
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# Panel frames.
#
# Panel p is the set of directions whose largest-magnitude Cartesian component
# is the one p names: 1..6 = +x, -x, +y, -y, +z, -z. This is the SAME rule as
# _cubed_sphere_panel in src/io/write_output.jl, so the `panel` field in the VTK
# output and the panel used here agree.
#
# Each panel carries a right-handed orthonormal triple (êu, êv, êw), êw the
# outward face normal and êu × êv = êw. A point of the panel is
#
#     r̂ = U êu + V êv + W êw ,   (U, V, W) = map(u, v)
#
# and the face coordinate that comes back out is (u, v) = (r̂·êu, r̂·êv)/(r̂·êw).
#
# WHY THE HANDEDNESS MATTERS. Two panels share a cube edge, and a node ON that
# edge is a SINGLE node of the watertight grid — it gets exactly one panel
# assignment, decided by a floating-point tie-break between two equal
# components. The remap is only well posed if both panels would have moved that
# node to the same place. They do, for every map here, because
#
#   * the maps are symmetric, map(v, u) = swap₁₂(map(u, v)), and
#   * an edge of the square maps into the plane of the cube edge, U = W on
#     u = 1,
#
# and with right-handed frames throughout those two facts make the two panels'
# images of a shared edge identical point by point. test/test_cubed_sphere_maps.jl
# measures it on all twelve cube edges rather than trusting the argument.
#---------------------------------------------------------------------------------
const _PANEL_FRAMES = (
    # (êu, êv, êw)
    ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0), ( 1.0, 0.0, 0.0)),   # 1: +x
    ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)),   # 2: -x
    ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), ( 0.0, 1.0, 0.0)),   # 3: +y
    ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), ( 0.0,-1.0, 0.0)),   # 4: -y
    ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), ( 0.0, 0.0, 1.0)),   # 5: +z
    ((0.0, 1.0, 0.0), (1.0, 0.0, 0.0), ( 0.0, 0.0,-1.0)),   # 6: -z
)

# Which panel a direction belongs to. Ties (a point exactly on a cube edge or
# corner) are broken in a fixed order; the maps agree there, so the choice does
# not matter — see the header note above _PANEL_FRAMES.
@inline function cubed_sphere_panel_of(x::Real, y::Real, z::Real)
    ax, ay, az = abs(x), abs(y), abs(z)
    if ax >= ay && ax >= az
        return x >= 0 ? 1 : 2
    elseif ay >= az
        return y >= 0 ? 3 : 4
    else
        return z >= 0 ? 5 : 6
    end
end

# Cartesian direction -> (u, v) face coordinate of the EQUIDISTANT gnomonic map.
@inline function _panel_face_coords(x::Real, y::Real, z::Real, panel::Int)
    eu, ev, ew = _PANEL_FRAMES[panel]
    U = x*eu[1] + y*eu[2] + z*eu[3]
    V = x*ev[1] + y*ev[2] + z*ev[3]
    W = x*ew[1] + y*ew[2] + z*ew[3]
    return U/W, V/W
end

# (U, V, W) in the panel frame -> Cartesian direction.
@inline function _panel_to_cartesian(U::Real, V::Real, W::Real, panel::Int)
    eu, ev, ew = _PANEL_FRAMES[panel]
    return (U*eu[1] + V*ev[1] + W*ew[1],
            U*eu[2] + V*ev[2] + W*ew[2],
            U*eu[3] + V*ev[3] + W*ew[3])
end


#---------------------------------------------------------------------------------
# The gnomonic and equiangular maps, both closed form.
#
# Forward: face coordinate -> point of the unit sphere, in the panel frame.
# Inverse: the face coordinate of a point of the unit sphere, in the panel frame.
#
# `_gnomonic_forward` is the identity of this file's face coordinate: the (u, v)
# that _panel_face_coords returns is BY DEFINITION the equidistant gnomonic one.
#---------------------------------------------------------------------------------
@inline function _gnomonic_forward(u::Real, v::Real)
    s = 1.0/sqrt(1.0 + u*u + v*v)
    return u*s, v*s, s
end

@inline _gnomonic_inverse(u::Real, v::Real) = (u, v)

# Equiangular: the face coordinate is the tangent of an angle that runs
# uniformly over [-π/4, π/4], so a uniform grid in (α, β) has u = tan α.
@inline function _equiangular_forward(u::Real, v::Real)
    a = tan(0.25*π*u)
    b = tan(0.25*π*v)
    s = 1.0/sqrt(1.0 + a*a + b*b)
    return a*s, b*s, s
end

@inline function _equiangular_inverse(u::Real, v::Real)
    return (4.0/π)*atan(u), (4.0/π)*atan(v)
end


#---------------------------------------------------------------------------------
# The conformal map of Rančić, Purser & Mesinger (1996).
#
# A_RANCIC is their Table B1, index 1 holding the constant term (zero) so that
# A_RANCIC[k] multiplies Z^(k-1).
#---------------------------------------------------------------------------------
const A_RANCIC = (
    +0.00000000000000,
    +1.47713062600964,
    -0.38183510510174,
    -0.05573058001191,
    -0.00895883606818,
    -0.00791315785221,
    -0.00486625437708,
    -0.00329251751279,
    -0.00235481488325,
    -0.00175870527475,
    -0.00135681133278,
    -0.00107459847699,
    -0.00086944475948,
    -0.00071607115121,
    -0.00059867100093,
    -0.00050699063239,
    -0.00043415191279,
    -0.00037541003286,
    -0.00032741060100,
    -0.00028773091482,
    -0.00025458777519,
    -0.00022664642371,
    -0.00020289261022,
    -0.00018254510830,
    -0.00016499474461,
    -0.00014976117168,
    -0.00013646173946,
    -0.00012478875823,
    -0.00011449267279,
    -0.00010536946150,
    -0.00009725109376,
)

#
# Formal reversion of a power series.
#
# Given W(Z) = Σ_{n≥1} a[n+1] Zⁿ with a[1] = 0 and a[2] ≠ 0, return b with
# Z(W) = Σ_{n≥1} b[n+1] Wⁿ and Z(W(Z)) = Z + O(Z^{N+1}). Both vectors are
# 1-based on the CONSTANT term, matching A_RANCIC's layout.
#
# Solving order by order: [Zⁿ](Σ_k b_k Wᵏ) = δ_{n1}, and [Zⁿ]Wᵏ vanishes for
# k > n, so b_n falls out of the orders already known. `pow` carries Wᵏ as it
# grows, truncated at N, which is what keeps this O(N³) rather than symbolic.
#
function _reverted_series(a::NTuple{M,Float64}) where {M}
    N = M - 1                                  # highest power retained
    A = collect(a)
    b = zeros(Float64, M)
    a1 = A[2]
    a1 == 0.0 && error(" # ERROR cubed_sphere_maps.jl: series is not invertible (a₁ = 0).")

    # Wᵏ for k = 1, built up in place. wk[n+1] = [Zⁿ] Wᵏ.
    wk = zeros(Float64, M)
    copyto!(wk, A)
    # coefficient of Zⁿ in Wᵏ, for every k, stored as we go
    powcoef = Vector{Vector{Float64}}(undef, N)
    powcoef[1] = copy(wk)
    for k = 2:N
        nxt = zeros(Float64, M)
        prev = powcoef[k-1]
        for i = 0:N, j = 0:N-i
            (prev[i+1] == 0.0 || A[j+1] == 0.0) && continue
            nxt[i+j+1] += prev[i+1]*A[j+1]
        end
        powcoef[k] = nxt
    end

    for n = 1:N
        acc = (n == 1) ? 1.0 : 0.0
        for k = 1:n-1
            acc -= b[k+1]*powcoef[k][n+1]
        end
        b[n+1] = acc/powcoef[n][n+1]
    end
    return Tuple(b)
end

const B_RANCIC = _reverted_series(A_RANCIC)

@inline function _series_eval(c::NTuple{M,Float64}, Z::Complex) where {M}
    # Horner, from the highest power down.
    acc = complex(c[M], 0.0)
    @inbounds for k = M-1:-1:1
        acc = acc*Z + c[k]
    end
    return acc
end

@inline _W_rancic(Z::Complex) = _series_eval(A_RANCIC, Z)
@inline _Z_rancic(W::Complex) = _series_eval(B_RANCIC, W)

#
# Face coordinate (u, v) ∈ [-1,1]² -> unit-sphere point in the panel frame.
#
# The eight-fold symmetry of the square is used first: the series is only valid
# on one octant, so |u|, |v| and the u/v ordering are folded away and restored
# at the end. The interior of the algorithm is RPM96's.
#
function _conformal_forward(u::Real, v::Real)
    (abs(u) > 1.0 || abs(v) > 1.0) &&
        error(" # ERROR cubed_sphere_maps.jl: conformal map needs (u,v) ∈ [-1,1]², got ($u, $v).")

    au, av = abs(u), abs(v)
    swapped = av > au                       # fold onto the octant au ≥ av

    uc = 1.0 - (swapped ? av : au)
    vc = 1.0 - (swapped ? au : av)

    Z = ((uc + im*vc)/2)^4
    W = _W_rancic(Z)

    im3 = im^(1/3)
    ra  = sqrt(3.0) - 1.0
    cb  = -1.0 + im
    cc  = ra*cb/2

    W = im3*(W*im)^(1/3)
    W = (W - ra)/(cb + cc*W)
    U, V = reim(W)

    # inverse stereographic lift onto the unit sphere
    H = 2.0/(1.0 + U*U + V*V)
    U *= H
    V *= H
    W3 = H - 1.0

    swapped && ((U, V) = (V, U))
    u < 0 && (U = -U)
    v < 0 && (V = -V)
    # the folds above leave a ±0 rather than an exact zero on the axes
    u == 0 && (U = 0.0)
    v == 0 && (V = 0.0)

    return U, V, W3
end

#
# The inverse: a unit-sphere point in the panel frame -> face coordinate.
# Only the first quadrant is handled by the series; the rest is symmetry.
#
function _conformal_inverse_UVW(U::Real, V::Real, W::Real)
    su, sv = sign(U), sign(V)
    aU, aV = abs(U), abs(V)
    swapped = aV > aU
    swapped && ((aU, aV) = (aV, aU))

    H  = W + 1.0
    ωs = (aU/H) + im*(aV/H)

    ra = sqrt(3.0) - 1.0
    cb = -1.0 + im
    cc = ra*cb/2
    ω0 = (ωs*cb + ra)/(1.0 - ωs*cc)
    W0 = im*ω0^3*im
    Z  = _Z_rancic(W0)
    z  = 2*Z^(1/4)
    uc, vc = reim(z)

    u = 1.0 - abs(vc)
    v = 1.0 - abs(uc)
    swapped && ((u, v) = (v, u))

    return (su == 0 ? 0.0 : su*u), (sv == 0 ? 0.0 : sv*v)
end

@inline function _conformal_inverse(u::Real, v::Real)
    # (u, v) here is the EQUIDISTANT gnomonic face coordinate of the point, so
    # rebuild the panel-frame direction first and then invert the conformal map.
    U, V, W = _gnomonic_forward(u, v)
    return _conformal_inverse_UVW(U, V, W)
end


#---------------------------------------------------------------------------------
# Dispatch table. `forward` takes a face coordinate to the panel-frame point,
# `inverse` takes the EQUIDISTANT gnomonic face coordinate of a point (which is
# what _panel_face_coords hands out) back to that map's own face coordinate.
#---------------------------------------------------------------------------------
const CUBED_SPHERE_MAPS = (:gnomonic, :equidistant, :equiangular, :conformal)

@inline function _map_forward(name::Symbol, u::Real, v::Real)
    name === :gnomonic || name === :equidistant ? _gnomonic_forward(u, v)    :
    name === :equiangular                       ? _equiangular_forward(u, v) :
    name === :conformal                         ? _conformal_forward(u, v)   :
    error(" # ERROR cubed_sphere_maps.jl: unknown cubed-sphere map ", name)
end

@inline function _map_inverse(name::Symbol, u::Real, v::Real)
    name === :gnomonic || name === :equidistant ? _gnomonic_inverse(u, v)    :
    name === :equiangular                       ? _equiangular_inverse(u, v) :
    name === :conformal                         ? _conformal_inverse(u, v)   :
    error(" # ERROR cubed_sphere_maps.jl: unknown cubed-sphere map ", name)
end

#
# The composite the remap actually applies: a Cartesian direction in, a
# Cartesian direction out, both on the unit sphere.
#
# The SOURCE map is always the gnomonic one. That is not a limitation dressed
# up as a default: reading a node's face coordinate back off the sphere,
# (u, v) = (r̂·êu, r̂·êv)/(r̂·êw), IS the gnomonic inverse, and it is what
# cubed_sphere.geo produces. Offering a "source map" input instead would be a
# claim about the .msh that nothing in the code can verify, and getting it
# wrong yields a grid that is neither map — with no error.
#
function remap_direction(x::Real, y::Real, z::Real, to::Symbol)
    r = sqrt(x*x + y*y + z*z)
    r > 0.0 || return (x, y, z)
    xn, yn, zn = x/r, y/r, z/r

    panel = cubed_sphere_panel_of(xn, yn, zn)
    u, v  = _panel_face_coords(xn, yn, zn, panel)
    # clamp: a node exactly on a cube edge can land at 1+ε through round-off,
    # and the conformal series rejects anything outside the closed square.
    U, V, W = _map_forward(to, clamp(u, -1.0, 1.0), clamp(v, -1.0, 1.0))
    return _panel_to_cartesian(U, V, W, panel)
end

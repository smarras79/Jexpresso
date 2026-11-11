#=================================================================================
 Shape-memory alloy material for the quasi-static solid path
 (src/kernel/solvers/quasistatic_solid.jl): the Souza–Auricchio–Petrini model
 as written by Artioli, Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60
 (2017), Sec. 3.3, in plane stress or plane strain.

 THE MODEL (infinitesimal strain, temperature T a given function of time):

     ε = εᵉ + eᵗʳ,   eᵗʳ deviatoric (the transformation is isochoric),
     ‖eᵗʳ‖ ≤ εL                                  saturation
     ψ = ½K (tr ε)² + G‖dev ε − eᵗʳ‖²            elastic
         + τ(T) ‖eᵗʳ‖,   τ = β⟨T − Mf⟩          chemical
         + ½h ‖eᵗʳ‖²                             transformation hardening
     σ = K tr ε 1 + s,   s = 2G (dev ε − eᵗʳ)
     X = s − α,   α = τ ∂‖eᵗʳ‖ + h eᵗʳ + ∂I_S(eᵗʳ)      (paper eqs. 3.32–3.36)
     f(X) = ‖X‖ − R,   R = √(2/3) σy0,   ėᵗʳ = ζ̇ X/‖X‖,  ζ̇ ≥ 0, ζ̇ f = 0

 At T ≤ Mf (τ = 0) a transformation strain left by a load stays after
 unloading as long as h εL < R — the shape-memory effect; heating raises τ
 and the chemical term drives eᵗʳ back to zero — the shape recovery. Above Mf
 the forward transformation needs ‖s‖ > R + τ and reverses on unloading —
 superelasticity. The paper's device (Sec. 4.4) runs the first, then heats.

 THE STEP. Backward Euler of the flow rule is the incremental minimisation
 (the dissipation is R‖ėᵗʳ‖, the support function of the elastic domain)

     eᵗʳₙ₊₁ = argmin over ‖x‖ ≤ εL of
              Φ(x) = G‖e − x‖² + ½h‖x‖² + τ‖x‖ + R‖x − xₙ‖,     e = dev εₙ₊₁,

 a strictly convex problem with one minimiser — the state update of Artioli &
 Bisegna (2016) that the paper uses, solved here exactly. Its optimality
 conditions are the discrete flow rule and Kuhn–Tucker conditions. Every term
 depends on x only through ‖x‖, x·e, x·xₙ, so the minimiser lies in the plane
 of e and xₙ, and the cases are:

   elastic    x = xₙ             0 ∈ ∂Φ(xₙ) + N_S(xₙ): closed-form tests
   reverse    x = 0 ≠ xₙ         ‖2G e + R xₙ/‖xₙ‖‖ ≤ τ (full reverse
                                 transformation, e.g. on heating)
   flow       0 < ‖x‖ < εL       X = R (x − xₙ)/‖x − xₙ‖; with c = 2G + h + τ/‖x‖
                                 x = (1 − R/D) a/c + (R/D) xₙ,
                                 a = 2G e,  D = ‖a − c xₙ‖ > R
                                 — closed form for τ = 0 (c fixed: the radial
                                 return of linear kinematic hardening), one
                                 scalar root in ‖x‖ otherwise
   saturated  ‖x‖ = εL           the same with c = 2G + h + τ/εL + γ, γ ≥ 0 the
                                 multiplier of the constraint: one scalar root
                                 in c (‖x(c)‖ decreases with c).

 THE TANGENT. Elastic and reverse: the elastic one. Flow with τ = 0: the
 analytic consistent tangent of the radial return. Saturated or τ > 0: central
 differences of the update — exact to ~1e-8, enough for quadratic Newton.

 PLANE STRESS (the paper's device). σzz = 0 is enforced at every quadrature
 point by Newton on εzz with the 3D update as a black box ("nested
 iterations", paper Sec. 4.4); εzz is kept with the history, and the in-plane
 tangent is the 3D one with the zz row and column condensed out,
 D − D:z Dz:/Dzz. PLANE STRAIN: εzz = 0.

 History per quadrature point: [eᵗʳ_xx, eᵗʳ_yy, eᵗʳ_zz, eᵗʳ_xy, εzz].
 Plottable: etr_xx … etr_xy, etr_norm (‖eᵗʳ‖), etr_sat (‖eᵗʳ‖/εL, 0 → 1:
 saturated), etr_eq (√(2/3)‖eᵗʳ‖, the uniaxial-equivalent transformation
 strain) and ezz.
=================================================================================#

struct SolidSMA <: AbstractSolidModel
    E::Float64; ν::Float64; K::Float64; G::Float64
    h::Float64                       # transformation hardening       [stress]
    β::Float64                       # dτ/dT                          [stress/temperature]
    Mf::Float64                      # martensite finish temperature
    εL::Float64                      # saturation: ‖eᵗʳ‖ ≤ εL
    σy0::Float64                     # f(X) = ‖X‖ − √(2/3) σy0
    plane_stress::Bool
    temperature::Any                 # t ↦ T(t), from the case
    T::Base.RefValue{Float64}        # T of the increment being solved (solid_set_time!)
end

nhistory(::SolidSMA)        = 5
islinear(::SolidSMA)        = false
solid_hypothesis(m::SolidSMA) = m.plane_stress ? "plane stress" : "plane strain"

sma_tau(m::SolidSMA, T) = m.β*max(T - m.Mf, 0.0)

function solid_set_time!(m::SolidSMA, t)
    m.T[] = Float64(m.temperature(t))
    return nothing
end

function solid_history_fields(m::SolidSMA)
    nrm = h -> sqrt(h[1]^2 + h[2]^2 + h[3]^2 + 2.0*h[4]^2)
    f = Pair{String,Function}[]
    for (c, cn) in enumerate(("xx", "yy", "zz", "xy"))
        push!(f, string("etr_", cn) => (h -> h[c]))
    end
    push!(f, "etr_norm" => nrm)
    push!(f, "etr_sat"  => (h -> nrm(h)/m.εL))
    push!(f, "etr_eq"   => (h -> sqrt(2.0/3.0)*nrm(h)))
    push!(f, "ezz"      => (h -> h[5]))
    return f
end

#
# Deviatoric tensors with xz = yz = 0 as 4-tuples (xx, yy, zz, xy); the inner
# product counts the shear twice.
#
const _SMA4 = NTuple{4,Float64}
@inline _dd(a, b)         = a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + 2.0*a[4]*b[4]
@inline _nn(a)            = sqrt(max(_dd(a, a), 0.0))
@inline _lin(α, a, β, b)  = (α*a[1] + β*b[1], α*a[2] + β*b[2], α*a[3] + β*b[3], α*a[4] + β*b[4])
@inline _sc(α, a)         = (α*a[1], α*a[2], α*a[3], α*a[4])
const _SMA0 = (0.0, 0.0, 0.0, 0.0)

# The transformation strain x(c) that minimises ½c‖x‖² − a·x + R‖x − xₙ‖.
@inline function _sma_xc(a, xn, c, R)
    u = _lin(1.0, a, -c, xn)
    D = _nn(u)
    D <= R && return xn
    return _lin((1.0 - R/D)/c, a, R/D, xn)
end

"""
    sma_return_map(m, a, xn, τ) -> (x, case)

One backward-Euler step of the transformation strain: the minimiser of
`½(2G + h)‖x‖² − a·x + τ‖x‖ + R‖x − xn‖` over `‖x‖ ≤ εL`, with `a = 2G dev ε`
(see the header). `case` is `:elastic`, `:reverse`, `:flow` or `:saturated`.
"""
function sma_return_map(m::SolidSMA, a::_SMA4, xn::_SMA4, τ::Float64)
    c0  = 2.0*m.G + m.h
    R   = sqrt(2.0/3.0)*m.σy0
    εL  = m.εL
    rn  = _nn(xn)
    tol = 1.0e-12*εL

    # x = xₙ?  0 ∈ −(a − c0 xₙ) + τ∂‖xₙ‖ + R B + N_S(xₙ)
    y = _lin(1.0, a, -c0, xn)
    if rn <= tol
        _nn(y) <= R + τ && return (xn, :elastic)
    else
        nx = _sc(1.0/rn, xn)
        z  = _lin(1.0, y, -τ, nx)
        if rn < εL - tol
            _nn(z) <= R && return (xn, :elastic)
        else                                           # saturated: the normal cone takes z·n > 0
            p  = _dd(z, nx)
            zt = p > 0 ? _lin(1.0, z, -p, nx) : z
            _nn(zt) <= R && return (xn, :elastic)
        end
        # x = 0 (full reverse transformation)?  a + R xₙ/‖xₙ‖ ∈ τ B
        _nn(_lin(1.0, a, R, nx)) <= τ && return (_SMA0, :reverse)
    end

    # Transformation. Unsaturated root first.
    xun = _SMA0
    sat = false
    if τ == 0.0
        xun = _sma_xc(a, xn, c0, R)                    # c fixed: closed form
        sat = _nn(xun) > εL
    elseif rn <= tol                                   # from xₙ = 0: ‖x‖ = (‖a‖ − R − τ)/c0
        na  = _nn(a)
        r   = (na - R - τ)/c0
        sat = r > εL
        sat || (xun = _sc(r/na, a))
    else                                               # g(r) = ‖x(c0 + τ/r)‖ − r = 0 on (0, εL)
        g = r -> _nn(_sma_xc(a, xn, c0 + τ/r, R)) - r
        if g(εL) >= 0
            sat = true
        else
            lo = εL; glo = g(lo)
            while glo <= 0 && lo > 1.0e-300           # g > 0 near 0 (the reverse case is excluded)
                lo *= 0.5; glo = g(lo)
            end
            hi = 2.0*lo > εL ? εL : 2.0*lo
            hi == lo && (hi = εL)
            for _ = 1:200                              # bisection: g(lo) > 0 ≥ g(hi)
                mid = 0.5*(lo + hi)
                (hi - lo) <= 4.0*eps(hi) && break
                g(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            r   = 0.5*(lo + hi)
            xun = _sma_xc(a, xn, c0 + τ/r, R)
        end
    end
    sat || return (xun, :flow)

    # Saturated: c ≥ cs = c0 + τ/εL with ‖x(c)‖ = εL (‖x(c)‖ decreases with c).
    cs = c0 + τ/εL
    lo = cs; hi = 2.0*cs
    while _nn(_sma_xc(a, xn, hi, R)) > εL
        lo = hi; hi *= 2.0
    end
    for _ = 1:200
        mid = 0.5*(lo + hi)
        (hi - lo) <= 4.0*eps(hi) && break
        _nn(_sma_xc(a, xn, mid, R)) > εL ? (lo = mid) : (hi = mid)
    end
    x = _sma_xc(a, xn, 0.5*(lo + hi), R)
    return (_sc(εL/_nn(x), x), :saturated)             # exactly on the constraint
end

#
# The 3D point: strain (xx, yy, zz, xy — tensor shear) ↦ stress (xx, yy, zz, xy).
#
@inline function _sma_stress(m::SolidSMA, ε::_SMA4, xn::_SMA4, τ::Float64)
    tr = ε[1] + ε[2] + ε[3]
    e  = (ε[1] - tr/3, ε[2] - tr/3, ε[3] - tr/3, ε[4])
    a  = _sc(2.0*m.G, e)
    x, case = sma_return_map(m, a, xn, τ)
    s  = _lin(1.0, a, -2.0*m.G, x)
    p  = m.K*tr
    return (p + s[1], p + s[2], p + s[3], s[4]), x, case, a
end

"""
    sma_tangent3d!(D, m, ε, xn, τ, case, a)

The 4×4 consistent tangent ∂σ/∂ε (rows σxx, σyy, σzz, σxy; columns εxx, εyy,
εzz and the ENGINEERING shear γxy = 2εxy) of the 3D update at `ε`.
"""
function sma_tangent3d!(D, m::SolidSMA, ε::_SMA4, xn::_SMA4, τ::Float64, case::Symbol, a::_SMA4)
    K, G = m.K, m.G
    if case === :elastic || case === :reverse || (case === :flow && τ == 0.0)
        θ = 1.0; θb = 0.0; N = _SMA0
        if case === :flow                            # radial return, c = 2G + h fixed
            c  = 2.0*G + m.h
            R  = sqrt(2.0/3.0)*m.σy0
            u  = _lin(1.0, a, -c, xn)
            Dn = _nn(u)
            θ  = 1.0 - (2.0*G/c)*(1.0 - R/Dn)
            θb = (2.0*G/c)*(R/Dn)
            N  = _sc(1.0/Dn, u)
        end
        @inbounds for j = 1:3, i = 1:3
            D[i,j] = K + 2.0*G*θ*((i == j ? 1.0 : 0.0) - 1.0/3.0) - 2.0*G*θb*N[i]*N[j]
        end
        @inbounds for i = 1:3
            D[i,4] = -2.0*G*θb*N[i]*N[4]
            D[4,i] = -2.0*G*θb*N[4]*N[i]
        end
        D[4,4] = G*θ - 2.0*G*θb*N[4]*N[4]
        return D
    end
    # saturated, or τ > 0: central differences of the update
    δ = 1.0e-8*max(m.εL, abs(ε[1]), abs(ε[2]), abs(ε[3]), abs(ε[4]))
    @inbounds for j = 1:4
        hj = j == 4 ? 0.5*δ : δ                      # column 4: dγ = 2 dεxy
        εp = ntuple(i -> ε[i] + (i == j ? hj : 0.0), 4)
        εm = ntuple(i -> ε[i] - (i == j ? hj : 0.0), 4)
        σp = first(_sma_stress(m, εp, xn, τ))
        σm = first(_sma_stress(m, εm, xn, τ))
        for i = 1:4
            D[i,j] = (σp[i] - σm[i])/(2.0*δ)
        end
    end
    return D
end

"""
    sma_point!(σ, D3, m, εxx, εyy, εxy, xn, εzz0, τ) -> (x, εzz, case)

The quadrature point in the 2D hypothesis of the model: the in-plane strain,
the converged transformation strain `xn` and temperature parameter `τ` give
the stress `σ = (σxx, σyy, σzz, σxy)`, the in-plane tangent `D3` (xx, yy, xy;
engineering shear), the new transformation strain and εzz. Plane stress: εzz
from σzz = 0 by Newton, starting at `εzz0`; plane strain: εzz = 0.
"""
function sma_point!(σ, D3, m::SolidSMA, εxx, εyy, εxy, xn::_SMA4, εzz0, τ)
    D4 = zeros(4, 4)
    εzz = m.plane_stress ? Float64(εzz0) : 0.0
    ε  = (Float64(εxx), Float64(εyy), εzz, Float64(εxy))
    σ4, x, case, a = _sma_stress(m, ε, xn, τ)
    if m.plane_stress
        # σzz to ~1e-14 of the stress of a saturated strain: the in-plane stress
        # (and the global residual) inherit this error, so it must sit well
        # below the global Newton tolerance
        scale = m.E*m.εL
        for it = 1:60
            abs(σ4[3]) <= 1.0e-14*scale && break
            # ∂σzz/∂εzz: analytic where the tangent is, else a central difference
            if case === :elastic || case === :reverse || (case === :flow && τ == 0.0)
                sma_tangent3d!(D4, m, ε, xn, τ, case, a)
                dzz = D4[3,3]
            else
                δ  = 1.0e-8*max(m.εL, abs(εzz))
                dzz = (_sma_stress(m, (ε[1], ε[2], εzz + δ, ε[4]), xn, τ)[1][3] -
                       _sma_stress(m, (ε[1], ε[2], εzz - δ, ε[4]), xn, τ)[1][3])/(2.0*δ)
            end
            εzz -= σ4[3]/dzz
            ε = (ε[1], ε[2], εzz, ε[4])
            σ4, x, case, a = _sma_stress(m, ε, xn, τ)
            it == 60 && error(@sprintf(" # ERROR quasistatic_solid_sma.jl: plane stress: σzz = %.3e after 60 iterations on εzz.", σ4[3]))
        end
    end
    sma_tangent3d!(D4, m, ε, xn, τ, case, a)
    idx = (1, 2, 4)
    @inbounds for (jj, j) in enumerate(idx), (ii, i) in enumerate(idx)
        D3[ii,jj] = m.plane_stress ? D4[i,j] - D4[i,3]*D4[3,j]/D4[3,3] : D4[i,j]
    end
    σ[1] = σ4[1]; σ[2] = σ4[2]; σ[3] = σ4[3]; σ[4] = σ4[4]
    return x, εzz, case
end

function solid_constitutive_update!(σ, D, hnew, m::SolidSMA, εxx, εyy, εxy, hold, δt)
    xn = (hold[1], hold[2], hold[3], hold[4])
    x, εzz, _ = sma_point!(σ, D, m, εxx, εyy, εxy, xn, hold[5], sma_tau(m, m.T[]))
    hnew[1] = x[1]; hnew[2] = x[2]; hnew[3] = x[3]; hnew[4] = x[4]; hnew[5] = εzz
    return nothing
end

#
# The case's material NamedTuple → SolidSMA (called by solid_model_from_case).
#
function solid_sma_from_case(m, E, ν, K, G)
    h   = Float64(m.h);   β  = Float64(m.β);  Mf = Float64(m.Mf)
    εL  = Float64(m.εL);  σy0 = Float64(m.σy0)
    hyp = Symbol(get(m, :hypothesis, :plane_stress))
    hyp in (:plane_stress, :plane_strain) ||
        error(" # ERROR quasistatic_solid_sma.jl: hypothesis = :", hyp, " — use :plane_stress or :plane_strain.")
    (h > 0 && β >= 0 && εL > 0 && σy0 > 0) ||
        error(" # ERROR quasistatic_solid_sma.jl: need h > 0, β ≥ 0, εL > 0, σy0 > 0 (got h = $h, β = $β, εL = $εL, σy0 = $σy0).")
    temp = get(m, :temperature, t -> Mf)
    temp isa Number && (temp = let T0 = Float64(temp); t -> T0; end)
    return SolidSMA(E, ν, K, G, h, β, Mf, εL, σy0, hyp === :plane_stress, temp, Ref(Float64(temp(0.0))))
end

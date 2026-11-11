#---------------------------------------------------------------------------------
# References for the SMA arch.
#
# 1. The device has no closed form: user_solid_exact returns `nothing` and the
#    run is compared with the paper's Fig. 5.6 (tools/artioli2017_data/fig5_6.csv,
#    read from the vector graphics of the arXiv PDF).
#
# 2. VERIFICATION MODE (JEXPRESSO_SMA_VERIFY=elastic): the same arch, linear
#    elastic, bent by a force P along +x at its free end — the curved bar of
#    Timoshenko & Goodier (Theory of Elasticity, Sec. 31), exact for any opening
#    angle. With the Airy function φ = (A r³ + B/r + C r + D r ln r) sin θ,
#
#      σr = (2Ar − 2B/r³ + D/r) sin θ,   σθ = (6Ar + 2B/r³ + D/r) sin θ,
#      τrθ = −(2Ar − 2B/r³ + D/r) cos θ,
#      A = P/2N,  B = −P a²b²/2N,  D = −P (a² + b²)/N,
#      N = b² − a² − (a² + b²) ln(b/a),
#
#    traction-free on r = a and r = b, and on the end θ = 0 a radial (here: x)
#    shear traction (2Ar − 2B/r³ + D/r) whose resultant is P. Integrating the
#    plane-stress strains (derived for this case; checked by test/test_sma_arch.jl
#    against finite differences of the strains):
#
#      u_r = (sin θ/E)[A(1−3ν) r² + B(1+ν)/r² + D(1−ν) ln r] − (2D/E)(sin θ + θ cos θ)
#      u_θ = −(cos θ/E)[A(5+ν) r² + B(1+ν)/r² + D(1−ν)(1 − ln r)] + (2D/E) θ sin θ
#
#    plus a rigid motion, chosen so that the centroid of the clamped end
#    (r = (a+b)/2, θ = π) neither moves nor rotates. Plane strain: E → E/(1−ν²),
#    ν → ν/(1−ν). The verification run prescribes this displacement on the clamp
#    and this traction on the loaded end, so the SEM solution must reproduce the
#    field everywhere: the solver reports the error (solid_history.csv).
#---------------------------------------------------------------------------------

sma_verify_force() = 1.0           # the force of the verification run [N/mm]

function curved_bar_constants(a, b, P)
    N = b^2 - a^2 - (a^2 + b^2)*log(b/a)
    return (P/(2N), -P*a^2*b^2/(2N), -P*(a^2 + b^2)/N)
end

# polar displacements and stresses, without the rigid motion
function curved_bar_polar(r, θ, a, b, P, E, ν)
    A, B, D = curved_bar_constants(a, b, P)
    ur = (sin(θ)/E)*(A*(1 - 3ν)*r^2 + B*(1 + ν)/r^2 + D*(1 - ν)*log(r)) - (2D/E)*(sin(θ) + θ*cos(θ))
    ut = -(cos(θ)/E)*(A*(5 + ν)*r^2 + B*(1 + ν)/r^2 + D*(1 - ν)*(1 - log(r))) + (2D/E)*θ*sin(θ)
    f  = 2A*r - 2B/r^3 + D/r
    q  = 6A*r + 2B/r^3 + D/r
    return ur, ut, f*sin(θ), q*sin(θ), -f*cos(θ)
end

_arch_theta(x, y) = atan(max(y, 0.0), x)                   # the arch lies in y ≥ 0

function _curved_bar_cartesian(x, y, a, b, P, E, ν)
    r = hypot(x, y); θ = _arch_theta(x, y)
    ur, ut, σr, σθ, τ = curved_bar_polar(r, θ, a, b, P, E, ν)
    c, s = cos(θ), sin(θ)
    return (ur*c - ut*s, ur*s + ut*c,
            σr*c^2 + σθ*s^2 - 2τ*s*c, σr*s^2 + σθ*c^2 + 2τ*s*c, (σr - σθ)*s*c + τ*(c^2 - s^2))
end

"""
    curved_bar_exact(x, y; P, hypothesis) -> (; ux, uy, sxx, syy, sxy, szz)

The exact elastic field of the arch bent by a force P along +x at its free
end, the centroid of the clamped end fixed (no displacement, no rotation).
"""
function curved_bar_exact(x, y; P = sma_verify_force(), hypothesis = sma_arch_hypothesis())
    M = sma_arch_properties()
    a, b = M.a, M.b
    E, ν = hypothesis === :plane_strain ? (M.E/(1 - M.ν^2), M.ν/(1 - M.ν)) : (M.E, M.ν)
    ux, uy, sxx, syy, sxy = _curved_bar_cartesian(x, y, a, b, P, E, ν)
    # rigid motion (Ux, Uy, ω): u += (Ux − ω y, Uy + ω x) with u = 0 and
    # ½(∂uy/∂x − ∂ux/∂y) = 0 at the clamped end's centroid (−r0, 0)
    r0 = (a + b)/2; δ = 1.0e-5*r0
    u0 = _curved_bar_cartesian(-r0, 0.0, a, b, P, E, ν)
    ω0 = 0.5*((_curved_bar_cartesian(-r0 + δ, 0.0, a, b, P, E, ν)[2] - _curved_bar_cartesian(-r0 - δ, 0.0, a, b, P, E, ν)[2]) -
              (_curved_bar_cartesian(-r0, δ, a, b, P, E, ν)[1] - _curved_bar_cartesian(-r0, -δ, a, b, P, E, ν)[1]))/(2δ)
    ω  = -ω0
    Ux = -u0[1]; Uy = -u0[2] + ω*r0
    szz = hypothesis === :plane_strain ? M.ν*(sxx + syy) : 0.0
    return (; ux = ux + Ux - ω*y, uy = uy + Uy + ω*x, sxx, syy, sxy, szz)
end

# the end traction of the curved-bar solution: on θ = 0 (x = r), σ·n = (f(r), 0)
function curved_bar_end_traction(x; P = sma_verify_force())
    M = sma_arch_properties()
    A, B, D = curved_bar_constants(M.a, M.b, P)
    return 2A*x - 2B/x^3 + D/x
end

user_solid_exact(x, y, t) = sma_verify_mode() ? curved_bar_exact(x, y) : nothing

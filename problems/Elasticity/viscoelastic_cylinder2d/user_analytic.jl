#---------------------------------------------------------------------------------
# Reference solutions for the viscoelastic thick-walled cylinder.
#
# 1. THE CLOSED FORM — pressure p switched on as a Heaviside step at t = 0.
#
# The elastic plane-strain (Lamé) solution of a cylinder a ≤ r ≤ b with
# pressure p inside and a free outer surface is
#
#     u_r = P [ r/M + b²/(2G r) ],     P = p a²/(b² − a²),   M = 2K + 2G/3,
#
# and the correspondence principle turns it into the viscoelastic one: the step
# load becomes p/s, G becomes sĜ(s), K (elastic in volume) stays K, and the two
# moduli invert into creep compliances,
#
#     u_r(r,t) = P [ r C(t) + (b²/2r) J(t) ]
#
#     J(t) = (1/G) [ 1/μ₀ − (1/μ₀ − 1) e^{−μ₀ t/λ} ]           shear creep compliance
#     C(t) = 1/M∞ − (1/M∞ − 1/M₀) e^{−M∞ t/(M₀ λ)},
#            M₀ = 2K + 2G/3,  M∞ = 2K + 2μ₀G/3,
#
# (μ₀ + μ₁ = 1). At t = 0, J = 1/G and C = 1/M₀: the Lamé solution,
# u_r = P (1+ν)/E [(1−2ν) r + b²/r].
#
# THE STRESSES. σ_rr and σ_θθ of this problem do not depend on the material at
# all — both surfaces carry prescribed tractions and the problem is statically
# determinate in the plane — so they are the Lamé stresses at EVERY t:
#
#     σ_rr = P (1 − b²/r²),    σ_θθ = P (1 + b²/r²).
#
# Only σ_zz, which plane strain needs to keep ε_zz = 0, evolves. Elastically
# σ_zz = ν(σ_rr + σ_θθ) = 2νP; by the same correspondence argument
#
#     σ_zz(t) = 2P [ ν∞ + (ν − ν∞) e^{−t/τ_z} ],
#     ν∞ = (3K − 2μ₀G)/(2(3K + μ₀G)),   τ_z = λ (3K + G)/(3K + μ₀G).
#
# So in the free-outer-surface case the stresses do NOT redistribute in the
# plane: the displacement (and σ_zz) are what carries the viscoelasticity, and
# the displacement is what the case checks. (In the case-bonded variant the
# in-plane stresses do redistribute.)
#
# How this was checked before being trusted as a regression reference:
# re-derived independently; its t = 0 limit is the Lamé solution to round-off;
# the exact-in-space recursion below converges to it at exactly second order
# for both exponential integrators and first order for backward Euler (both
# Prony sets; u_r for all three, σ_zz for the exponential one); and that
# recursion, run the way FEAP runs it,
# reproduces the paper's Fig. 5.2 to the precision of the figure (see the
# README of this case and test/test_viscoelastic_cylinder.jl).
#
# NOTE λ is the RELAXATION TIME (τ in the code), not the Lamé constant.
#---------------------------------------------------------------------------------

@inline function cylinder_creep_compliances(t)
    P  = cylinder_properties()
    μ0 = P.g0
    J  = (1.0/P.G)*(1.0/μ0 - (1.0/μ0 - 1.0)*exp(-μ0*t/P.τ))
    M0 = 2.0*P.K + 2.0*P.G/3.0
    Mi = 2.0*P.K + 2.0*μ0*P.G/3.0
    C  = 1.0/Mi - (1.0/Mi - 1.0/M0)*exp(-Mi*t/(M0*P.τ))
    return (; J, C)
end

"""
    cylinder_ur_exact(r, t)

Radial displacement of the cylinder at radius `r`, time `t ≥ 0`, under a
Heaviside internal pressure applied at `t = 0` (closed form).
"""
@inline function cylinder_ur_exact(r, t)
    P  = cylinder_properties()
    Pf = P.p*P.a^2/(P.b^2 - P.a^2)
    c  = cylinder_creep_compliances(t)
    return Pf*(r*c.C + P.b^2/(2.0*r)*c.J)
end

@inline function cylinder_szz_exact(t)
    P   = cylinder_properties()
    Pf  = P.p*P.a^2/(P.b^2 - P.a^2)
    μ0  = P.g0
    νi  = (3.0*P.K - 2.0*μ0*P.G)/(2.0*(3.0*P.K + μ0*P.G))
    τz  = P.τ*(3.0*P.K + P.G)/(3.0*P.K + μ0*P.G)
    return 2.0*Pf*(νi + (P.ν - νi)*exp(-t/τz))
end

#
# Hook read by the quasi-static solver: Cartesian displacement and stress.
#
function user_solid_exact(x, y, t)
    P   = cylinder_properties()
    r   = hypot(x, y)
    c   = x/r; s = y/r
    ur  = cylinder_ur_exact(r, t)
    Pf  = P.p*P.a^2/(P.b^2 - P.a^2)
    srr = Pf*(1.0 - P.b^2/r^2)
    stt = Pf*(1.0 + P.b^2/r^2)
    return (; ux  = ur*c,
              uy  = ur*s,
              sxx = srr*c*c + stt*s*s,
              syy = srr*s*s + stt*c*c,
              sxy = (srr - stt)*s*c,
              szz = cylinder_szz_exact(t))
end

#---------------------------------------------------------------------------------
# 2. THE TIME-DISCRETE REFERENCE — exact in space, discrete in time.
#
# For this geometry every displacement the solver can reach is of Lamé form,
# u_r = A r + B/r: its deviatoric strain is (A/3 − B/r², A/3 + B/r², −2A/3)
# in (rr, θθ, zz), every Maxwell internal variable driven by it stays of the
# same form (α/3 − β/r², ...), and equilibrium dσ_rr/dr + (σ_rr − σ_θθ)/r = 0
# then holds IDENTICALLY. Each time step therefore reduces to the two traction
# conditions σ_rr(a) = −p, σ_rr(b) = 0 for (A, B):
#
#     B = [ p a²b²/(b² − a²) − 2Gμ₁ β* ] / (2G̃)
#     A = [ P − 2Gμ₁ α*/3 ] / (2K + 2G̃/3)
#
# with G̃ = G(μ₀ + μ₁γ), α* = d α − γ A, β* = d β − γ B (the previous step), and
# afterwards α ← α* + γA, β ← β* + γB — the very recursion the SEM solver runs
# at every quadrature point, with the SAME (d, γ) (maxwell_step_coefficients).
#
# What it is for:
#   * SEM minus this = the SPATIAL error alone, at any Δt;
#   * this minus the closed form = the TIME-integration error alone;
#   * with initial_solve = false it is what FEAP computes in the paper (the load
#     entering during the first increment), i.e. the paper's Fig. 5.2.
#
# Returns the step times and, per time, (A, B) and σ_zz; u_r(r) = A r + B/r.
#---------------------------------------------------------------------------------
function cylinder_time_discrete_reference(times::AbstractVector;
                                          integrator::Symbol = :exponential,
                                          initial_solve::Bool = true)
    P  = cylinder_properties()
    μ0 = P.g0; μ1 = P.g1
    Pf = P.p*P.a^2/(P.b^2 - P.a^2)

    t  = Float64[times[1]]
    As = [0.0]; Bs = [0.0]; szz = [0.0]
    α = β = A = B = 0.0
    function step!(δt)
        d, γ = maxwell_step_coefficients(δt, P.τ, integrator)
        Gt = P.G*(μ0 + μ1*γ)
        αs = d*α - γ*A
        βs = d*β - γ*B
        B  = (P.p*P.a^2*P.b^2/(P.b^2 - P.a^2) - 2.0*P.G*μ1*βs)/(2.0*Gt)
        A  = (Pf - 2.0*P.G*μ1*αs/3.0)/(2.0*P.K + 2.0*Gt/3.0)
        α  = αs + γ*A
        β  = βs + γ*B
        return nothing
    end
    stress_zz() = 2.0*P.K*A + 2.0*P.G*(μ0*(-2.0*A/3.0) + μ1*(-2.0*α/3.0))

    if initial_solve
        step!(0.0)
        As[1] = A; Bs[1] = B; szz[1] = stress_zz()
    end
    for n = 2:length(times)
        step!(times[n] - times[n-1])
        push!(t, times[n]); push!(As, A); push!(Bs, B); push!(szz, stress_zz())
    end
    return (; t, A = As, B = Bs, szz)
end

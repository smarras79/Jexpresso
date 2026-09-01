#---------------------------------------------------------------------------------
# 2D linear elastodynamics, velocity-stress form, with displacement recovery.
#
# The first five equations are the plane-strain elastodynamic system, exactly
# as in problems/Elasticity/plane_wave:
#
#     U₁₋₅ = (ρu, ρv, σxx, σyy, σxy)ᵀ
#     F₁₋₅ = -(σxx, σxy, (λ+2μ)u, λu,      μv)ᵀ
#     G₁₋₅ = -(σxy, σyy, λv,      (λ+2μ)v, μu)ᵀ
#
# — balance of linear momentum plus the rate form of Hooke's law, with
# eigenvalues ±cp = ±√((λ+2μ)/ρ), ±cs = ±√(μ/ρ) and 0 in every direction, and
# the mechanical energy ½ρ|v|² + ½σ:C⁻¹:σ as convex entropy.
#
# Two more are carried here so the run produces the DEFORMED BEAM and not only
# the velocity and stress fields:
#
#     U₆ = uₓ     ∂ₜuₓ + ∂ₓ(0) + ∂_y(0) = u
#     U₇ = u_y    ∂ₜu_y + ∂ₓ(0) + ∂_y(0) = v
#
# Zero flux, so they add two zero eigenvalues and change neither the
# characteristic structure nor the stable time step. They are reconstructions,
# not conservation laws — the physics is in rows 1-5 — but they are what lets
# the output be read as a beam bending. Warping the mesh by (uₓ, u_y) in
# ParaView is the picture.
#
# Unlike the 1D Timoshenko cases, the two recovery rows are the ONLY source in
# this system: 2D elastodynamics itself is pure divergence.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Material and beam geometry, in one place; every other user_*.jl of this case
# reads them from here.
#
# A METHOD and not a set of `const`s on purpose: src/run.jl loads a case by
# `include`ing its user_*.jl files into the Jexpresso module, so a `const`
# would collide with the sibling case (Elasticity/plane_wave) the first time
# the two are run in the same session.
#
# Same material and same beam as the 1D Timoshenko cases — ν = 0.3, L/h = 10 —
# so the frequency this case produces can be read directly against them.
#---------------------------------------------------------------------------------
@inline function elastic_properties()

    ρ = 1.0        # density
    E = 1.0        # Young's modulus
    ν = 0.3        # Poisson's ratio

    L = 1.0        # beam length      (must match the mesh)
    h = 0.1        # beam depth       (must match the mesh: y ∈ [-h/2, h/2])

    λ  = E*ν/((1.0 + ν)*(1.0 - 2.0*ν))   # Lamé λ   (plane strain)
    μs = E/(2.0*(1.0 + ν))               # Lamé μ = shear modulus. NOT `μ`:
                                         # :μ is also an inputs key.

    return (ρ = ρ, E = E, ν = ν, L = L, h = h,
            λ = λ, μ = μs, λ2μ = λ + 2.0*μs,
            cp = sqrt((λ + 2.0*μs)/ρ),   # P-wave speed
            cs = sqrt(μs/ρ),             # S-wave speed
            #
            # Beam-theory reference numbers, for comparison only — nothing in
            # the solver uses them. In PLANE STRAIN the bending modulus is
            # E/(1-ν²), not E: a plane-strain beam is stiffer than a plane-
            # stress one by 1/(1-ν²) = 1.0989 here.
            #
            Ebend = E/(1.0 - ν*ν),
            Ib    = h^3/12.0,            # second moment of area, unit thickness
            ρA    = ρ*h)
end

#---------------------------------------------------------------------------------
# First bending mode of a clamped-free Euler-Bernoulli beam, normalised to 1 at
# the tip. β₁L = 1.8751… is the first root of cos(βL)cosh(βL) + 1 = 0, and
#
#     Φ(x) = cosh βx - cos βx - σ₁(sinh βx - sin βx),
#     σ₁   = (cosh β₁L + cos β₁L)/(sinh β₁L + sin β₁L) = 0.7340955…
#
# Φ(0) = Φ'(0) = 0, which is what makes it a legal initial velocity for a
# clamped end: no motion and no rotation at the wall.
#
# It is used here only to SHAPE the initial velocity. It is a mode of the
# reduced 1D model, not of the 2D elastic beam, so the run will contain a
# little of the neighbouring modes as well — visible as a slow beat, and
# honest: it is the difference between the two models.
#---------------------------------------------------------------------------------
@inline function cantilever_mode1(x)

    P  = elastic_properties()
    βL = 1.8751040687119611
    σ1 = 0.7340955138
    β  = βL/P.L

    Φ  = (cosh(β*x) - cos(β*x)) - σ1*(sinh(β*x) - sin(β*x))

    return Φ/2.0        # Φ(L) = 2 exactly, so this is 1 at the tip
end

#
# Euler-Bernoulli estimate of that mode's frequency, in PLANE STRAIN.
#
@inline function cantilever_omega1()
    P  = elastic_properties()
    βL = 1.8751040687119611
    return (βL/P.L)^2*sqrt(P.Ebend*P.Ib/P.ρA)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=7, ip=1)

    P = elastic_properties()

    u   = q[1]/P.ρ
    v   = q[2]/P.ρ
    σxx = q[3]
    σyy = q[4]
    σxy = q[5]

    F[1] = -σxx
    F[2] = -σxy
    F[3] = -P.λ2μ*u
    F[4] = -P.λ*u
    F[5] = -P.μ*v
    F[6] = 0.0          # uₓ and u_y are advanced by their sources alone
    F[7] = 0.0

    G[1] = -σxy
    G[2] = -σyy
    G[3] = -P.λ*v
    G[4] = -P.λ2μ*v
    G[5] = -P.μ*u
    G[6] = 0.0
    G[7] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    P = elastic_properties()

    u = q[1]/P.ρ
    v = q[2]/P.ρ

    return (T(-q[3]), T(-q[5]), T(-P.λ2μ*u), T(-P.λ*u), T(-P.μ*v), T(0.0), T(0.0),
            T(-q[5]), T(-q[4]), T(-P.λ*v),   T(-P.λ2μ*v), T(-P.μ*u), T(0.0), T(0.0))
end

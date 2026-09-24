#---------------------------------------------------------------------------------
# 2D linear elastodynamics, velocity-stress form, with displacement recovery.
#
# THE PROBLEM: a beam resting on two hinged supports with a load pressed onto
# the middle of its top surface. It sags, overshoots, and rings about the
# static deflection. Textbook statics, solved as a time-dependent 2D elastic
# body with no beam theory anywhere in it.
#
# State (plane strain), rows 1-5:
#
#     U₁₋₅ = (ρvx, ρvy, σxx, σyy, σxy)ᵀ
#
#     ∂ₜ(ρvx) = ∂ₓσxx + ∂_y σxy
#     ∂ₜ(ρvy) = ∂ₓσxy + ∂_y σyy
#     ∂ₜσxx   = (λ+2μ) ∂ₓvx + λ ∂_y vy
#     ∂ₜσyy   = λ ∂ₓvx + (λ+2μ) ∂_y vy
#     ∂ₜσxy   = μ (∂_y vx + ∂ₓvy)
#
# i.e. ∂ₜU + ∂ₓF(U) + ∂_y G(U) = S with
#
#     F(U) = -(σxx, σxy, (λ+2μ)vx, λvx,      μvy)ᵀ
#     G(U) = -(σxy, σyy, λvy,      (λ+2μ)vy, μvx)ᵀ
#
# Rows 6-7 recover the displacement so the output is the DEFORMED BEAM and not
# only velocity and stress:
#
#     U₆ = ux    ∂ₜux + ∂ₓ(0) + ∂_y(0) = vx
#     U₇ = uy    ∂ₜuy + ∂ₓ(0) + ∂_y(0) = vy
#
# Zero flux, so they add two zero eigenvalues and change neither the
# characteristic structure nor the stable time step. They are reconstructions,
# not conservation laws, and they are the only source in this system: rows 1-5
# are pure divergence.
#
# NAMING, because the output does not warn you: vx/vy are the particle
# VELOCITY, ux/uy are the DISPLACEMENT. The sagging beam is ux/uy.
#
# The flux Jacobian in any direction n has eigenvalues ±cp = ±√((λ+2μ)/ρ),
# ±cs = ±√(μ/ρ) and 0, so two incoming characteristics and exactly two
# conditions at every boundary. Convex entropy is the mechanical energy
# ½ρ|v|² + ½σ:C⁻¹:σ, whose flux is -σ·v, the mechanical power.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Material, geometry and load, in one place; every other user_*.jl of this case
# reads them from here.
#
# A METHOD and not a set of `const`s on purpose: src/run.jl loads a case by
# `include`ing its user_*.jl files into the Jexpresso module, so a `const`
# would collide with the sibling case (Elasticity/plane_wave2d) the first time
# the two are run in the same session.
#
# Non-dimensional: unit density and Young's modulus, ν = 0.3. Span 1, depth
# 0.2, so L/h = 5 — stocky enough that shear deformation is 12% of the midspan
# sag, which is exactly the regime where a 2D calculation earns its keep over
# Euler-Bernoulli.
#---------------------------------------------------------------------------------
@inline function elastic_properties()

    ρ = 1.0        # density
    E = 1.0        # Young's modulus
    ν = 0.3        # Poisson's ratio

    L = 1.0        # span        (must match the mesh: x ∈ [0, L])
    h = 0.2        # depth       (must match the mesh: y ∈ [-h/2, h/2])

    λ  = E*ν/((1.0 + ν)*(1.0 - 2.0*ν))   # Lamé λ   (plane strain)
    μs = E/(2.0*(1.0 + ν))               # Lamé μ = shear modulus. NOT `μ`:
                                         # :μ is also an inputs key.

    # --- the load ------------------------------------------------------------
    Ptot = 3.0e-4   # total downward force pressed onto the top surface
    aload = 0.2     # width of the loaded patch, centred at midspan
    tramp = 4.0     # the load is raised smoothly over this long

    return (ρ = ρ, E = E, ν = ν, L = L, h = h,
            λ = λ, μ = μs, λ2μ = λ + 2.0*μs,
            cp = sqrt((λ + 2.0*μs)/ρ),       # P-wave speed
            cs = sqrt(μs/ρ),                 # S-wave speed
            Ptot = Ptot, aload = aload, tramp = tramp,
            #
            # Beam-theory reference numbers, for comparison only — nothing in
            # the solver uses them. In PLANE STRAIN the bending modulus is
            # E/(1-ν²), not E.
            #
            Ebend = E/(1.0 - ν*ν),
            Ib    = h^3/12.0,                # second moment of area
            ρA    = ρ*h,
            κGA   = (5.0/6.0)*μs*h)          # shear rigidity, rectangular κₛ
end

#---------------------------------------------------------------------------------
# The load pressed onto the top surface: a downward pressure p(x,t), so the
# prescribed normal traction there is σyy = -p.
#
#   in x   a cos² bump of width `aload` centred at midspan, which integrates to
#          Ptot and reaches the surface smoothly at both edges of the patch —
#          no corner for a non-dissipative scheme to ring on;
#
#   in t   a raised cosine from 0 to full over `tramp`. This matters as much as
#          the shape does: switching a load on instantaneously is a step in
#          time at that face, and it launches a front. `tramp` = 4 is about 4.6
#          wave transits of the span (long enough that the elastic waves stay
#          resolved) and 0.36 of a bending period (short enough that the beam
#          still swings and does not simply creep to the static answer).
#---------------------------------------------------------------------------------
@inline function beam_load_pressure(x, t)

    P = elastic_properties()

    ξ = (x - 0.5*P.L)/P.aload            # -1/2 … +1/2 across the patch
    abs(ξ) > 0.5 && return 0.0

    # ∫ p dx = Ptot: the cos² bump has mean 1/2 over its width.
    p0 = 2.0*P.Ptot/P.aload
    shape = cos(π*ξ)^2

    ramp = t >= P.tramp ? 1.0 : 0.5*(1.0 - cos(π*t/P.tramp))

    return p0*shape*ramp
end

#---------------------------------------------------------------------------------
# Static midspan deflection this load produces, in beam theory: the number the
# run should ring about. Printed by initialize.jl.
#---------------------------------------------------------------------------------
@inline function beam_static_midspan()
    P = elastic_properties()
    bending = P.Ptot*P.L^3/(48.0*P.Ebend*P.Ib)
    shear   = P.Ptot*P.L/(4.0*P.κGA)
    return (bending = bending, shear = shear, total = bending + shear)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=7, ip=1)

    P = elastic_properties()

    vx  = q[1]/P.ρ
    vy  = q[2]/P.ρ
    σxx = q[3]
    σyy = q[4]
    σxy = q[5]

    F[1] = -σxx
    F[2] = -σxy
    F[3] = -P.λ2μ*vx
    F[4] = -P.λ*vx
    F[5] = -P.μ*vy
    F[6] = 0.0          # ux and uy are advanced by their sources alone
    F[7] = 0.0

    G[1] = -σxy
    G[2] = -σyy
    G[3] = -P.λ*vy
    G[4] = -P.λ2μ*vy
    G[5] = -P.μ*vx
    G[6] = 0.0
    G[7] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    P = elastic_properties()

    vx = q[1]/P.ρ
    vy = q[2]/P.ρ

    return (T(-q[3]), T(-q[5]), T(-P.λ2μ*vx), T(-P.λ*vx), T(-P.μ*vy), T(0.0), T(0.0),
            T(-q[5]), T(-q[4]), T(-P.λ*vy),   T(-P.λ2μ*vy), T(-P.μ*vx), T(0.0), T(0.0))
end

#---------------------------------------------------------------------------------
# Timoshenko (shear-deformable) beam — flux vector, with displacement recovery.
#
# The first four equations are the Timoshenko conservation laws exactly as in
# problems/Elasticity/simply_supported:
#
#     U₁₋₄ = (ρA v, ρI ω, γ, χ)ᵀ
#     F₁₋₄ = -(κₛGA γ, EI χ, v, ω)ᵀ
#     S₁₋₄ =  (q, κₛGA γ, -ω, 0)ᵀ
#
#     v = ẇ,  ω = φ̇,  γ = w_x - φ,  χ = φ_x
#
# with fluxes minus the shear force Q = κₛGA γ and the bending moment M = EI χ,
# and rows 3-4 the compatibility conditions γ̇ = v_x - ω, χ̇ = ω_x.
#
# Two more are carried here so the run produces the DEFORMED SHAPE and not only
# the strain/velocity measures:
#
#     U₅ = w      ∂ₜw + ∂ₓ(0) = v
#     U₆ = φ      ∂ₜφ + ∂ₓ(0) = ω
#
# They are passive: zero flux, so they add two zero eigenvalues and change
# neither the characteristic structure nor the stable time step. They are
# reconstructions, not conservation laws — the physics is in rows 1-4.
#
# Eigenvalues of the 4×4 block split into two decoupled pairs,
#
#     λ = ±√(κₛG/ρ)   (shear wave)      λ = ±√(E/ρ)   (extensional/rotary wave)
#
# so it is strictly hyperbolic away from κₛG = E and symmetrizable by the
# energy norm ½(ρA v² + ρI ω² + κₛGA γ² + EI χ²).
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Beam properties and the reference load, in one place; every other user_*.jl of
# this case reads them from here.
#
# A METHOD and not a set of `const`s on purpose: src/run.jl loads a case by
# `include`ing its user_*.jl files into the Jexpresso module, so a `const` would
# collide with the sibling case (Elasticity/simply_supported) the first time the
# two are run in the same session. Redefining a method at an identical signature
# is exactly what an `include` is supposed to do.
#
# Same beam as the sibling case: non-dimensional, unit density and Young's
# modulus, rectangular section of unit width and depth h = 0.1 on a beam of
# length 1, so L/h = 10 and shear deformation is a percent-level correction to
# Euler-Bernoulli.
#---------------------------------------------------------------------------------
@inline function timoshenko_properties()

    L  = 1.0        # beam length
    ρ  = 1.0        # mass density
    E  = 1.0        # Young's modulus
    ν  = 0.3        # Poisson's ratio
    κs = 5.0/6.0    # shear correction factor (rectangular section)

    h  = 0.1        # section depth
    bw = 1.0        # section width

    q0 = 1.0e-5     # uniform transverse load that produced the initial shape

    Gm = E/(2.0*(1.0 + ν))   # shear modulus (never spelled `G`: that name is
                             # taken by the y-flux argument of user_flux!)
    Ar = bw*h                # cross-section area
    Iy = bw*h^3/12.0         # second moment of area

    return (L = L, ρ = ρ, E = E, ν = ν, κs = κs, h = h, q0 = q0,
            A = Ar, I = Iy,
            ρA = ρ*Ar, ρI = ρ*Iy, EI = E*Iy, κGA = κs*Gm*Ar,
            cs = sqrt(κs*Gm/ρ),   # shear wave speed
            ce = sqrt(E/ρ))       # extensional/rotary wave speed
end

#---------------------------------------------------------------------------------
# Static Timoshenko cantilever under a uniform transverse load q₀: clamped at
# x = 0, free at x = L. This is the shape the run is released from.
#
#     Q(x) = q₀(L - x)                        (∂ₓQ + q₀ = 0, Q(L) = 0)
#     M(x) = q₀(L - x)²/2                     (∂ₓM + Q  = 0, M(L) = 0)
#     γ    = Q/κₛGA                           shear strain
#     χ    = M/EI                             curvature
#     φ    = q₀[L³ - (L - x)³]/(6EI)          from φ_x = χ, φ(0) = 0
#     w    = q₀x(2L - x)/(2κₛGA)
#            + q₀[L³x + (L - x)⁴/4 - L⁴/4]/(6EI)     from w_x = γ + φ, w(0) = 0
#
# The two terms of w are the shear and bending contributions; at the tip they
# are q₀L²/(2κₛGA) and q₀L⁴/(8EI), the second about 100× the first at L/h = 10.
#---------------------------------------------------------------------------------
@inline function timoshenko_static_cantilever(x)

    p  = timoshenko_properties()
    q0 = p.q0
    L  = p.L
    s  = L - x

    γ = q0*s/p.κGA
    χ = q0*s*s/(2.0*p.EI)
    φ = q0*(L^3 - s^3)/(6.0*p.EI)
    w = q0*x*(2.0*L - x)/(2.0*p.κGA) +
        q0*(L^3*x + s^4/4.0 - L^4/4.0)/(6.0*p.EI)

    return (γ = γ, χ = χ, φ = φ, w = w)
end

#---------------------------------------------------------------------------------
# F(U) = -(κₛGA γ, EI χ, v, ω, 0, 0)ᵀ.   G is unused in 1D.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=6, ip=1)

    p = timoshenko_properties()

    v = q[1]/p.ρA
    ω = q[2]/p.ρI
    γ = q[3]
    χ = q[4]

    F[1] = -p.κGA*γ     # -Q, shear force
    F[2] = -p.EI*χ      # -M, bending moment
    F[3] = -v
    F[4] = -ω
    F[5] = 0.0          # w and φ are advanced by their sources alone
    F[6] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    p = timoshenko_properties()

    return T(-p.κGA*q[3]), T(-p.EI*q[4]), T(-q[1]/p.ρA), T(-q[2]/p.ρI), T(0.0), T(0.0)
end

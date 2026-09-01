#---------------------------------------------------------------------------------
# Timoshenko (shear-deformable) beam — flux vector.
#
# State (all densities per unit length of beam):
#
#     U = (ρA v, ρI ω, γ, χ)ᵀ
#
#     v = ẇ        transverse velocity        (w  = transverse displacement)
#     ω = φ̇        cross-section angular rate  (φ  = cross-section rotation)
#     γ = w_x - φ  shear strain
#     χ = φ_x      curvature
#
# Conservation form ∂ₜU + ∂ₓF(U) = S with
#
#     F(U) = -(κₛGA γ, EI χ, v, ω)ᵀ
#     S    =  (q, κₛGA γ, -ω, 0)ᵀ
#
# Rows 1-2 are balance of linear and angular momentum; their fluxes are minus
# the shear force Q = κₛGA γ and the bending moment M = EI χ. Rows 3-4 are the
# compatibility conditions γ̇ = v_x - ω and χ̇ = ω_x. κₛGA γ and -ω sit in the
# source because they are undifferentiated, not because they are nonlinear:
# the whole system is linear with a constant Jacobian A = ∂F/∂U whose spectrum
# splits into two decoupled pairs,
#
#     λ = ±√(κₛG/ρ)   (shear wave)      λ = ±√(E/ρ)   (extensional/rotary wave)
#
# so it is strictly hyperbolic away from κₛG = E and symmetrizable by the
# energy norm ½(ρA v² + ρI ω² + κₛGA γ² + EI χ²).
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Beam properties, in one place. Every other user_*.jl of this case reads them
# from here.
#
# A METHOD and not a set of `const`s on purpose: src/run.jl loads a case by
# `include`ing its user_*.jl files into the Jexpresso module, so a `const` would
# collide with the sibling case (Elasticity/cantilever) the first time the two
# are run in the same session. Redefining a method at an identical signature is
# exactly what an `include` is supposed to do, and it constant-folds into the
# per-node loops all the same.
#
# Non-dimensional but self-consistent: unit density, unit Young's modulus, a
# rectangular section of unit width and depth h = 0.1 on a beam of length 1, so
# the slenderness L/h = 10 puts the case in the regime where shear deformation
# is a percent-level correction to Euler-Bernoulli — visible, but not dominant.
#---------------------------------------------------------------------------------
@inline function timoshenko_properties()

    L  = 1.0        # beam length
    ρ  = 1.0        # mass density
    E  = 1.0        # Young's modulus
    ν  = 0.3        # Poisson's ratio
    κs = 5.0/6.0    # shear correction factor (rectangular section)

    h  = 0.1        # section depth
    bw = 1.0        # section width

    Gm = E/(2.0*(1.0 + ν))   # shear modulus (never spelled `G`: that name is
                             # taken by the y-flux argument of user_flux!)
    Ar = bw*h                # cross-section area
    Iy = bw*h^3/12.0         # second moment of area

    return (L = L, ρ = ρ, E = E, ν = ν, κs = κs, h = h,
            A = Ar, I = Iy,
            ρA = ρ*Ar, ρI = ρ*Iy, EI = E*Iy, κGA = κs*Gm*Ar,
            cs = sqrt(κs*Gm/ρ),   # shear wave speed
            ce = sqrt(E/ρ))       # extensional/rotary wave speed
end

#---------------------------------------------------------------------------------
# The standing mode this case is initialised with, and against which
# user_analytic.jl compares.
#
# For a SIMPLY SUPPORTED beam, w = W sin(kx) cos(Ωt), φ = Φ cos(kx) cos(Ωt) with
# k = nπ/L satisfies w = 0 and M = EIφ_x = 0 at both ends exactly. Substituting
# into the two momentum balances leaves the Timoshenko frequency equation
#
#     (ρIρA/κGA) Ω⁴ - [ρI k² + EI ρA k²/κGA + ρA] Ω² + EI k⁴ = 0
#
# whose two roots are the flexural branch (the lower one, taken here) and the
# much faster shear branch. Φ then follows from the linear-momentum row:
#
#     Φ = W k (1 - ρA Ω²/(κGA k²))
#---------------------------------------------------------------------------------
@inline function timoshenko_mode()

    p = timoshenko_properties()

    n = 1                   # mode number
    W = 0.01                # transverse amplitude (1% of L: linear theory)

    k = n*π/p.L

    a = p.ρI*p.ρA/p.κGA
    b = p.ρI*k^2 + p.EI*p.ρA*k^2/p.κGA + p.ρA
    c = p.EI*k^4

    Ω2 = (b - sqrt(b*b - 4.0*a*c))/(2.0*a)   # flexural (lower) branch
    Ω  = sqrt(Ω2)

    Φ = W*k*(1.0 - p.ρA*Ω2/(p.κGA*k*k))

    return (n = n, k = k, Ω = Ω, W = W, Φ = Φ)
end

#---------------------------------------------------------------------------------
# F(U) = -(κₛGA γ, EI χ, v, ω)ᵀ.   G is unused in 1D.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    p = timoshenko_properties()

    v = q[1]/p.ρA
    ω = q[2]/p.ρI
    γ = q[3]
    χ = q[4]

    F[1] = -p.κGA*γ     # -Q, shear force
    F[2] = -p.EI*χ      # -M, bending moment
    F[3] = -v
    F[4] = -ω
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    p = timoshenko_properties()

    return T(-p.κGA*q[3]), T(-p.EI*q[4]), T(-q[1]/p.ρA), T(-q[2]/p.ρI)
end

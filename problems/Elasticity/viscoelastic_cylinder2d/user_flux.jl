#---------------------------------------------------------------------------------
# Thick-walled viscoelastic cylinder under internal pressure, plane strain.
#
#   E. Artioli, L. Beirão da Veiga, C. Lovadina, E. Sacco, "Arbitrary order 2D
#   virtual elements for polygonal meshes: Part II, inelastic problem",
#   Comput. Mech. 60 (2017) 643-657, arXiv:1701.06676 — Secs. 3.1 and 4.2.
#
# There is no flux here: this case runs on the QUASI-STATIC solid path
# (:lquasistatic_solid, src/kernel/solvers/quasistatic_solid.jl), which never
# calls user_flux!. This file holds what every other file of the case reads:
# material, geometry and load, in one place.
#
# A METHOD and not a set of `const`s for the same reason as in beam2d: src/run.jl
# `include`s a case into the Jexpresso module, and a `const` would collide with
# a sibling case run earlier in the same session.
#
# MATERIAL (paper Sec. 3.1, eqs. 3.5-3.7): generalized Maxwell in the deviator,
# elastic in volume, one Prony term (M = 1 — a "standard linear solid"):
#
#     G(t) = G (μ₀ + μ₁ e^{−t/λ}),   λ = 1,   K constant,
#     E = 1000, ν = 0.3   →   G = 384.615..., K = 833.333...
#
# NOTATION. The paper's λ is the RELAXATION TIME. It is called τ in the code,
# because the Lamé constant λ = K − 2G/3 lives in the same solver and is a
# different thing altogether. Likewise μ₀, μ₁ are Prony WEIGHTS here, not shear
# moduli (they are spelled g0, g1 below; `:μ` is also an inputs key).
#
# Two Prony sets, chosen with the environment variable JEXPRESSO_VE_SET:
#
#     ve1  (μ₀, μ₁) = (0.01, 0.99)   default. K/G(0) = 2.167, K/G(∞) = 216.7:
#                                    nearly incompressible under sustained load
#                                    (ν(∞) = 0.498). Creeps for a LONG time — the
#                                    shear creep time is λ/μ₀ = 100.
#     ve2  (μ₀, μ₁) = (0.30, 0.70)   an intermediate response, settled by t ≈ 15.
#
# GEOMETRY AND LOAD (paper Sec. 4.2, Fig. 5.1): inner radius a = 2, outer
# radius b = 4, a quarter of the cross-section with symmetry conditions on the
# two straight edges and a traction-free outer surface; internal pressure
# p = 10 applied as a STEP at t = 0 and held to t = 20.
#---------------------------------------------------------------------------------
@inline function cylinder_properties()

    E = 1000.0     # Young's modulus
    ν = 0.3        # Poisson's ratio (instantaneous)

    set = lowercase(strip(get(ENV, "JEXPRESSO_VE_SET", "ve1")))
    if set == "ve1"
        g0, g1 = 0.01, 0.99
    elseif set == "ve2"
        g0, g1 = 0.30, 0.70
    else
        error(" # ERROR viscoelastic_cylinder2d: JEXPRESSO_VE_SET = \"", set,
              "\" — use ve1 (μ₀ = 0.01) or ve2 (μ₀ = 0.3).")
    end
    τ = 1.0        # relaxation time: the paper's λ. NOT the Lamé constant.

    a = 2.0        # inner radius  (must match the mesh)
    b = 4.0        # outer radius  (must match the mesh)
    p = 10.0       # internal pressure

    K = E/(3.0*(1.0 - 2.0*ν))
    G = E/(2.0*(1.0 + ν))

    return (; E, ν, K, G, g0, g1, τ, a, b, p, set)
end

#
# Hook read by the quasi-static solver: Young's modulus, Poisson's ratio and the
# normalised Prony series G(t)/G = g0 + Σ g[m] exp(−t/τ[m]).
#
function user_solid_material()
    P = cylinder_properties()
    return (; E = P.E, ν = P.ν, g0 = P.g0, g = [P.g1], τ = [P.τ])
end

#
# The internal pressure: a step, on from t = 0. Whether the solver treats it as
# already on at tinit (an instantaneous elastic solve first, the default) or as
# entering during the first increment is :linitial_elastic_solve, not this.
#
@inline cylinder_pressure(t) = cylinder_properties().p

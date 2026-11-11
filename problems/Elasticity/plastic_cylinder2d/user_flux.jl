#---------------------------------------------------------------------------------
# Thick-walled cylinder, elastic-perfectly plastic (J2), plane strain, opened
# up by a prescribed radial displacement of its bore until the whole wall has
# yielded — the quadrature-point-history test of the quasi-static solid path
# (src/kernel/solvers/quasistatic_solid.jl, model = :j2_plasticity).
#
# Same quarter cylinder as viscoelastic_cylinder2d (a = 2, b = 4, symmetry
# rollers, free outer surface). The material is the one Artioli et al. (2017)
# use for their plasticity benchmark (Sec. 4.3): E = 7000, ν = 0.3,
# σy0 = 24.3, no hardening.
#
# Why a cylinder: two numbers of the elastic-plastic solution are EXACT and
# known in closed form, independently of the discretisation —
#
#   first yield   the Lamé stresses (σ_zz = 2νP in plane strain) reach the von
#                 Mises surface at r = a:  p_e = σy0 / φ(a),  user_analytic.jl;
#   limit load    once the whole wall flows, plastic incompressibility and
#                 plane strain give s_zz = 0, σθθ − σrr = 2σy0/√3, and
#                 equilibrium dσrr/dr = (σθθ − σrr)/r integrates to
#                     p_L = (2/√3) σy0 ln(b/a)
#                 (Hill, The Mathematical Theory of Plasticity, 1950) —
#
# and in between, the problem stays axisymmetric, so a fine 1D radial model
# with the same constitutive routine is a reference for the whole curve.
#
# DISPLACEMENT CONTROL. A perfectly plastic body cannot carry a prescribed
# pressure above p_L, so the bore is driven by its radial displacement δ(t)
# instead and the pressure is read back as the reaction: p = R_inner,x / a.
# t ∈ [0, 1] is a load parameter, not a time: J2 plasticity is rate
# independent.
#
# A METHOD and not `const`s, as in the sibling cases (src/run.jl includes a
# case into the Jexpresso module).
#---------------------------------------------------------------------------------
@inline function plastic_cylinder_properties()

    E    = 7000.0     # Young's modulus           (Artioli et al. Sec. 4.3)
    ν    = 0.3        # Poisson's ratio
    σy0  = 24.3       # initial yield stress
    Hiso = 0.0        # isotropic hardening modulus   (perfectly plastic)
    Hkin = 0.0        # kinematic hardening modulus

    a    = 2.0        # inner radius  (must match the mesh)
    b    = 4.0        # outer radius  (must match the mesh)
    δmax = 0.05       # bore displacement at t = 1: ≈ 8.7x first yield, past full plastification

    K = E/(3.0*(1.0 - 2.0*ν))
    G = E/(2.0*(1.0 + ν))

    return (; E, ν, σy0, Hiso, Hkin, a, b, δmax, K, G)
end

#
# Hook read by the quasi-static solver.
#
function user_solid_material()
    P = plastic_cylinder_properties()
    return (; model = :j2_plasticity, E = P.E, ν = P.ν, σy0 = P.σy0, Hiso = P.Hiso, Hkin = P.Hkin)
end

# The radial displacement imposed on the bore, ramped linearly with t ∈ [0, 1].
@inline bore_displacement(t) = plastic_cylinder_properties().δmax*t

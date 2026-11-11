function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the SMA arch .................................. ")

    #-----------------------------------------------------------------------------
    # The quasi-static solver's fields: SOLID_FIELDS (displacement and stress;
    # szz ≈ 0 in plane stress), the von Mises stress "svm", the transformation
    # strain — its saturation ratio "etr_sat" = ‖eᵗʳ‖/εL (0: none, 1:
    # saturated), "etr_eq" = √(2/3)‖eᵗʳ‖ (the uniaxial-equivalent strain) and
    # the components — the thickness strain "ezz" of plane stress, and the
    # nodal reactions "rx", "ry" (non-zero on the clamp only).
    # Everything starts at zero: undeformed, unstressed, no transformation
    # strain (twinned martensite at T = Mf).
    #-----------------------------------------------------------------------------
    qvars = ["ux", "uy", "sxx", "syy", "sxy", "szz", "svm", "etr_sat", "etr_eq",
             "etr_xx", "etr_yy", "etr_zz", "etr_xy", "ezz", "rx", "ry"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qvars)

    P = sma_arch_properties()
    R = sqrt(2/3)*P.σy0
    println("   Artioli et al. (2017) Sec. 4.4: semicircular arch, ", P.a, " ≤ r ≤ ", P.b,
            " mm, clamped at θ = π, traction q(t) along +x on the end θ = 0")
    println("   SMA (Souza–Auricchio–Petrini), ", sma_arch_hypothesis(), ": E = ", P.E, " MPa, ν = ", P.ν,
            ", εL = ", P.εL, ", h = ", P.h, " MPa, β = ", P.β, " MPa/K, Mf = ", P.Mf, " K, σy0 = ", P.σy0, " MPa")
    @printf("   Table 3: q = 0 → %g → 0 → %g → 0 N/mm at T = %g K, then heating to %g K at q = 0\n",
            P.qmax, -P.qmax, P.Troom, P.Troom + P.ΔT)
    @printf("   at T = Mf: transformation from ‖s‖ = R = %.2f MPa (uniaxial %g MPa); h εL = %g < R:\n",
            R, P.σy0, P.h*P.εL)
    println("              the transformation strain stays after unloading (shape memory)")
    @printf("   heating: τ = β (T − Mf) reaches %g MPa > R: the strain is recovered from T = %.1f K on\n",
            P.β*P.ΔT, P.Mf + R/P.β)
    sma_verify_mode() &&
        println("   VERIFICATION MODE: elastic (σy0 → ∞), exact curved-bar displacement on the clamp and ",
                "traction on the free end (user_analytic.jl): the SEM solution must equal the exact one")

    println(" Initialize fields for the SMA arch .................................. DONE")

    return q
end

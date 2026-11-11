function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the viscoelastic thick-walled cylinder ........ ")

    #-----------------------------------------------------------------------------
    # The fields of the quasi-static solver (SOLID_FIELDS in
    # src/kernel/solvers/quasistatic_solid.jl): the DISPLACEMENT (ux, uy), which
    # is what it solves for, and the stress it derives from it, σzz included —
    # plane strain keeps ε_zz = 0, not σ_zz. Everything starts at zero:
    # undeformed, unstressed, no history. That is the only state the solver
    # accepts, since what a prior deformation did to the Maxwell internal
    # variables would not be known.
    #
    # ASCII names: they become VTK array names (see beam2d/initialize.jl).
    # qoutvars is what user_uout! (user_primitives.jl) writes: the six fields,
    # the closed-form displacement and the difference from it.
    #-----------------------------------------------------------------------------
    qvars    = ["ux", "uy", "sxx", "syy", "sxy", "szz"]
    qoutvars = ["ux", "uy", "sxx", "syy", "sxy", "szz", "ux_exact", "uy_exact", "ux_err", "uy_err"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    P  = cylinder_properties()
    tend = Float64(inputs[:tend])
    println("   Artioli et al. (2017) Sec. 4.2 — Prony set ", P.set, ": G(t)/G = ", P.g0,
            " + ", P.g1, " exp(-t/", P.τ, ")   (λ = ", P.τ, " is a relaxation time, not Lamé's λ)")
    println("   E = ", P.E, ", ν = ", P.ν, "  →  G = ", round(P.G, sigdigits = 7), ", K = ", round(P.K, sigdigits = 7))
    println("   a = ", P.a, ", b = ", P.b, ", internal pressure p = ", P.p, ", outer surface free")

    # The numbers the run is aiming at: the closed form at the two control points.
    lame(r) = P.p*P.a^2/(P.b^2 - P.a^2)*(1.0 + P.ν)/P.E*((1.0 - 2.0*P.ν)*r + P.b^2/r)
    @printf("   closed form  u_r(A = (%g,0)):  t=0⁺ %.8f (Lamé %.8f)   t=%g %.8f\n",
            P.a, cylinder_ur_exact(P.a, 0.0), lame(P.a), tend, cylinder_ur_exact(P.a, tend))
    @printf("   closed form  u_r(B = (%g,0)):  t=0⁺ %.8f (Lamé %.8f)   t=%g %.8f\n",
            P.b, cylinder_ur_exact(P.b, 0.0), lame(P.b), tend, cylinder_ur_exact(P.b, tend))
    @printf("   σ_rr, σ_θθ stay the Lamé stresses; σ_zz goes from %.6f (t=0⁺) to %.6f (t=%g)\n",
            cylinder_szz_exact(0.0), cylinder_szz_exact(tend), tend)

    println(" Initialize fields for the viscoelastic thick-walled cylinder ........ DONE")

    return q
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the elastic-plastic thick-walled cylinder ..... ")

    #-----------------------------------------------------------------------------
    # The quasi-static solver's fields: SOLID_FIELDS (displacement and stress),
    # then two extras it fills on request — the von Mises stress "svm" and the
    # accumulated plastic strain "ep_bar", a QUADRATURE-POINT history value that
    # the solver projects onto the nodes so it can be plotted like any field.
    # Any other J2 history field may be added here (ep_xx … ep_xy, alpha_xx …).
    # Everything starts at zero: undeformed, unstressed, no plastic strain.
    #-----------------------------------------------------------------------------
    qvars = ["ux", "uy", "sxx", "syy", "sxy", "szz", "svm", "ep_bar"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qvars)

    P  = plastic_cylinder_properties()
    fy = plastic_cylinder_first_yield()
    pL = plastic_cylinder_limit_pressure()
    println("   J2, perfectly plastic: E = ", P.E, ", ν = ", P.ν, ", σy0 = ", P.σy0,
            "   (Artioli et al. 2017, Sec. 4.3 material)")
    println("   a = ", P.a, ", b = ", P.b, ", bore displaced radially to δ = ", P.δmax, " at t = 1")
    @printf("   first yield (Lamé, exact):   p_e = %.6f at δ_e = %.6e\n", fy.p_e, fy.δ_e)
    @printf("   limit pressure (exact):      p_L = (2/√3) σy0 ln(b/a) = %.6f   (p_L/p_e = %.4f)\n", pL, pL/fy.p_e)
    @printf("   δ = %.3g is %.1f x δ_e: the wall yields through and the pressure levels off at p_L\n",
            P.δmax, P.δmax/fy.δ_e)

    println(" Initialize fields for the elastic-plastic thick-walled cylinder ..... DONE")

    return q
end

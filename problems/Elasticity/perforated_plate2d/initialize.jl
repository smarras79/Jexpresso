function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the perforated plastic plate ................. ")

    #-----------------------------------------------------------------------------
    # The quasi-static solver's fields: SOLID_FIELDS (displacement and stress),
    # the von Mises stress "svm", the accumulated plastic strain "ep_bar" — a
    # quadrature-point history value, projected onto the nodes for output — and
    # the nodal reactions "rx", "ry" (zero off the constrained edges; summed
    # over the top edge, ry is the force per mm of thickness behind the
    # paper's "reaction sum", Fig. 5.4, which is 10 times it).
    # Everything starts at zero: undeformed, unstressed, no plastic strain.
    #-----------------------------------------------------------------------------
    qvars = ["ux", "uy", "sxx", "syy", "sxy", "szz", "svm", "ep_bar", "rx", "ry"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qvars)

    P   = perforated_plate_properties()
    ref = perforated_plate_paper_reference()
    nst = round(Int, (inputs[:tend] - inputs[:tinit])/inputs[:Δt])
    println("   Artioli et al. (2017) Sec. 4.3: strip 2L x 2H = ", 2P.L, " x ", 2P.H,
            " mm, hole 2R = ", 2P.R, " mm, plane strain, a quarter modelled")
    println("   J2, perfectly plastic: E = ", P.E, ", ν = ", P.ν, ", σy0 = ", P.σy0, " kg/mm²")
    println("   top edge pulled up by δ = ", P.δmax, " mm in ", nst, " increments; ",
            "u_x = 0 on the lateral edge and on x = 0, u_y = 0 on y = 0")
    @printf("   paper, Table 1 (Ref):  u_A = %.3f mm,  v_B = %.3f mm     (FEM-Q9: %.3f, %.3f)\n",
            ref.uA, ref.vB, ref.quadratic.FEM_Q9...)
    @printf("   paper, Table 2: %.2f-%.2f Newton iterations per increment for quadratic elements\n",
            minimum(values(ref.newton_quadratic)), maximum(values(ref.newton_quadratic)))

    println(" Initialize fields for the perforated plastic plate ................. DONE")

    return q
end

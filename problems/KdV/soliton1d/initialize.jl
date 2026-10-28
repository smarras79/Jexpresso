function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for the KdV soliton (hyperbolized KdV) ... ")

    #-----------------------------------------------------------------------------
    # U = (u, v, w)ᵀ, in the order used by user_flux! / user_source!.
    #
    # u is the KdV soliton. v and w are set to its EXACT first and second
    # derivatives (kdv_soliton in user_flux.jl), i.e. the well-prepared data the
    # relaxation converges to, so no initial layer of width O(τ) is excited.
    #-----------------------------------------------------------------------------
    qvars = ["u", "v", "w"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]

            η, ηx, ηxx = kdv_soliton(x, 0.0)

            q.qn[ip, 1] = η
            q.qn[ip, 2] = ηx
            q.qn[ip, 3] = ηxx

            q.qe[ip, 1] = 0.0
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 0.0
        end
    end

    p = kdv_parameters()
    println("   τ = ", p.τ, ", soliton speed c = ", p.c, ", amplitude 3c = ", 3.0*p.c, ", x0 = ", p.x0)
    println("   characteristic speeds at u = 0: -1/τ = ", -1.0/p.τ, ", ±1/√τ = ±", 1.0/sqrt(p.τ))
    println(" Initialize fields for the KdV soliton (hyperbolized KdV) ... DONE ")

    return q
end

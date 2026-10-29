function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for Kuramoto-Sivashinsky (Kassam & Trefethen 2005) ... ")

    #-----------------------------------------------------------------------------
    # U = (u, w)ᵀ, in the order used by user_flux! / user_source!.
    #
    # u(x,0) = cos(x/16)(1 + sin(x/16)), and w starts from its exact second
    # derivative (ks_initial in user_flux.jl): the relaxed manifold w = ∂ₓ²u.
    #-----------------------------------------------------------------------------
    qvars = ["u", "w"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]

            u, uxx = ks_initial(x)

            q.qn[ip, 1] = u
            q.qn[ip, 2] = uxx

            q.qe[ip, 1] = 0.0
            q.qe[ip, 2] = 0.0
        end
    end

    # The space-time diagram of user_plot.jl is accumulated over the run.
    ks_reset_history!()

    p = ks_parameters()
    println("   L = 32π = ", p.L, ", τ = ", p.τ)
    println(" Initialize fields for Kuramoto-Sivashinsky (Kassam & Trefethen 2005) ... DONE ")

    return q
end

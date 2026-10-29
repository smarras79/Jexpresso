function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for Kuramoto-Sivashinsky (Kassam & Trefethen 2005) ... ")

    #-----------------------------------------------------------------------------
    # One unknown: u, the slope of the KS front, with Kassam & Trefethen's
    # u(x,0) = cos(x/16)(1 + sin(x/16)) (ks_initial in user_flux.jl).
    #-----------------------------------------------------------------------------
    qvars = ["u"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]

            q.qn[ip, 1] = ks_initial(mesh.coords[1, ip])
            q.qe[ip, 1] = 0.0
        end
    end

    # The space-time diagram of user_plot.jl is accumulated over the run.
    ks_reset_history!()

    println("   L = 32π = ", ks_parameters().L)
    println(" Initialize fields for Kuramoto-Sivashinsky (Kassam & Trefethen 2005) ... DONE ")

    return q
end

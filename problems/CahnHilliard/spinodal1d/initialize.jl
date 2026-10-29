function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for Cahn-Hilliard spinodal growth ... ")

    #-----------------------------------------------------------------------------
    # U = (c, μ)ᵀ, in the order used by user_flux! / user_source!.
    #
    # c is a small single-mode perturbation of the unstable mixture c = 0. μ is
    # its chemical potential, built from the exact second derivative of c
    # (ch_initial in user_flux.jl), so the run starts on the relaxed manifold.
    #-----------------------------------------------------------------------------
    qvars = ["c", "μ"]

    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]

            c, _, _, μ = ch_initial(x)

            q.qn[ip, 1] = c
            q.qn[ip, 2] = μ

            q.qe[ip, 1] = 0.0
            q.qe[ip, 2] = 0.0
        end
    end

    p = ch_parameters()
    println("   ε = ", p.ε, ", τ = ", p.τ, ", mode k = ", p.k, ", amplitude δ = ", p.δ)
    println("   Cahn-Hilliard growth rate σ(k) = k²(1 - ε²k²) = ", ch_growth_rate(p.k))
    println(" Initialize fields for Cahn-Hilliard spinodal growth ... DONE ")

    return q
end

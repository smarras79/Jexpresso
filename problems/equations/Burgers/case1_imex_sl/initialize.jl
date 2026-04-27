function initialize(SD, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
       Initial condition for the 1D viscous Burgers equation:

           q(x,0) = sin(2π x / L),    L = xmax - xmin

        which steepens and forms a (viscosity-smoothed) shock.

        This case uses the IMEX_ARS(2,3,2) Runge-Kutta time integrator:
        advection (q²/2)_x is treated explicitly, diffusion ν q_xx is
        treated implicitly.
    """
    @info " Initialize fields for 1D viscous Burgers (IMEX) ................ "

    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    xmin = inputs[:xmin]
    xmax = inputs[:xmax]
    L    = xmax - xmin

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl

            ip = mesh.connijk[iel_g,i,1]
            x  = mesh.x[ip]

            q.qn[ip,1] = sin(2.0*π*(x - xmin)/L)

            # Background/reference state (here zero) for perturbation plotting.
            q.qe[ip,1] = 0.0

        end
    end

    @info " Initialize fields for 1D viscous Burgers (IMEX) ................ DONE "

    return q
end

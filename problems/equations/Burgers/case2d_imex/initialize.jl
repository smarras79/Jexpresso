function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
       Initial condition for the 2D viscous Burgers equation:

           q(x,y,0) = sin(2π x / Lx) · sin(2π y / Ly),
               Lx = xmax - xmin,  Ly = ymax - ymin

       The sine tensor-product steepens along the diagonal and forms a
       viscosity-smoothed front. On a doubly-periodic mesh this is a clean
       IMEX test: advection (q²/2)_x + (q²/2)_y is treated explicitly while
       the viscous Laplacian ν (q_xx + q_yy) is solved implicitly with
       ARS(2,3,2).
    """
    @info " Initialize fields for 2D viscous Burgers (IMEX) ................ "

    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    xmin = minimum(mesh.x); xmax = maximum(mesh.x); Lx = xmax - xmin
    ymin = minimum(mesh.y); ymax = maximum(mesh.y); Ly = ymax - ymin

    for iel_g = 1:mesh.nelem
        for j = 1:mesh.ngl, i = 1:mesh.ngl

            ip = mesh.connijk[iel_g, i, j]
            x  = mesh.x[ip]
            y  = mesh.y[ip]

            q.qn[ip,1] = sin(2.0*π*(x - xmin)/Lx) * sin(2.0*π*(y - ymin)/Ly)

            # Background/reference state (here zero) for perturbation plotting.
            q.qe[ip,1] = 0.0

        end
    end

    @info " Initialize fields for 2D viscous Burgers (IMEX) ................ DONE "

    return q
end

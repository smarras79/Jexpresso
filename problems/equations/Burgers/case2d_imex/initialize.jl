function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
       Initial condition for the 2D inviscid Burgers equation in flux form

           ∂ₜ u + ∇·(½ u² v) = 0,    v = (1, 1)

       Piecewise-constant four-quadrant Riemann data (Sec. 4.1 of the
       reference paper), with discontinuities at x = x_mid and y = y_mid:

                   │  y > y_mid : -0.2   │  y > y_mid : -1.0
           x < x_m │                     │   x > x_mid
                   │  y < y_mid : +0.5   │  y < y_mid : +0.8

       The discontinuity locations use the domain midpoints so the test
       survives any rescaling of the underlying periodic gmsh mesh; on
       the unit square [0,1]² this reproduces the paper's x = y = 0.5
       cut lines exactly. The IMEX integrator treats the (weak) viscous
       regularization ν Δu implicitly — set ν=0 in user_inputs.jl to run
       the purely inviscid problem from the paper, or leave a small ν
       to get a vanishing-viscosity regularized front.
    """
    @info " Initialize fields for 2D Burgers Riemann (IMEX) ................ "

    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    xmin = minimum(mesh.x); xmax = maximum(mesh.x); x_mid = 0.5*(xmin + xmax)
    ymin = minimum(mesh.y); ymax = maximum(mesh.y); y_mid = 0.5*(ymin + ymax)

    for iel_g = 1:mesh.nelem
        for j = 1:mesh.ngl, i = 1:mesh.ngl

            ip = mesh.connijk[iel_g, i, j]
            x  = mesh.x[ip]
            y  = mesh.y[ip]

            if x < x_mid && y > y_mid
                q.qn[ip,1] = -0.2
            elseif x > x_mid && y > y_mid
                q.qn[ip,1] = -1.0
            elseif x < x_mid && y < y_mid
                q.qn[ip,1] =  0.5
            else   # x > x_mid && y < y_mid (and the two zero-measure lines)
                q.qn[ip,1] =  0.8
            end

            # Background/reference state (here zero) for perturbation plotting.
            q.qe[ip,1] = 0.0

        end
    end

    @info " Initialize fields for 2D Burgers Riemann (IMEX) ................ DONE "

    return q
end

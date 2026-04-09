function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    TC2: Inertia-Gravity Wave on f-plane (Bishnu et al. 2024, Sec. 2.3)

    Linearized rotating shallow water equations:
      ∂η/∂t + H₀(∂u/∂x + ∂v/∂y) = 0
      ∂u/∂t - f₀v + g∂η/∂x = 0
      ∂v/∂t + f₀u + g∂η/∂y = 0

    Exact solution: superposition of two inertia-gravity wave modes.
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D Shallow Water TC2 (Inertia-Gravity Wave) ..."
    end

    #---------------------------------------------------------------------------------
    # Solution variables: q = [η, u, v]
    #---------------------------------------------------------------------------------
    qvars    = ["η", "u", "v"]
    qoutvars = ["η", "u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        #---------------------------------------------------------------------------------
        # Physical parameters
        #---------------------------------------------------------------------------------
        Lx  = 1.0e7   # domain length x [m]
        Ly  = 1.0e7   # domain length y [m]
        H0  = 1.0e3   # mean depth [m]
        g   = 9.81     # gravity [m/s²]
        f0  = 1.0e-4   # Coriolis parameter [s⁻¹]

        # Mode 1
        A1  = 0.10     # amplitude [m]
        mx1 = 1; my1 = 1
        kx1 = 2π * mx1 / Lx
        ky1 = 2π * my1 / Ly
        K1sq = kx1^2 + ky1^2
        ω1  = sqrt(f0^2 + g * H0 * K1sq)

        # Mode 2
        A2  = 0.20     # amplitude [m]
        mx2 = 2; my2 = 2
        kx2 = 2π * mx2 / Lx
        ky2 = 2π * my2 / Ly
        K2sq = kx2^2 + ky2^2
        ω2  = sqrt(f0^2 + g * H0 * K2sq)

        t0 = inputs[:tinit]

        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            y = mesh.y[ip]

            # Mode 1 at t=t0
            θ1 = kx1 * x + ky1 * y - ω1 * t0
            η1 = A1 * cos(θ1)
            fac1 = A1 / (H0 * K1sq)
            u1 = fac1 * (kx1 * ω1 * cos(θ1) - f0 * ky1 * sin(θ1))
            v1 = fac1 * (ky1 * ω1 * cos(θ1) + f0 * kx1 * sin(θ1))

            # Mode 2 at t=t0
            θ2 = kx2 * x + ky2 * y - ω2 * t0
            η2 = A2 * cos(θ2)
            fac2 = A2 / (H0 * K2sq)
            u2 = fac2 * (kx2 * ω2 * cos(θ2) - f0 * ky2 * sin(θ2))
            v2 = fac2 * (ky2 * ω2 * cos(θ2) + f0 * kx2 * sin(θ2))

            # Superposition
            q.qn[ip, 1] = η1 + η2
            q.qn[ip, 2] = u1 + u2
            q.qn[ip, 3] = v1 + v2

            # Store exact solution at t=t0 for reference
            q.qe[ip, 1] = q.qn[ip, 1]
            q.qe[ip, 2] = q.qn[ip, 2]
            q.qe[ip, 3] = q.qn[ip, 3]
        end
    end

    if rank == 0
        @info " Initialize fields for 2D Shallow Water TC2 (Inertia-Gravity Wave) ... DONE"
    end

    return q
end

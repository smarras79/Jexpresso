function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    TC3: Planetary Rossby Wave on beta-plane (Bishnu et al. 2024, Sec. 2.7)

    Linearized rotating shallow water equations with f(y) = f₀ + β₀y:
      ∂η/∂t + H₀(∂u/∂x + ∂v/∂y) = -σ η
      ∂u/∂t + g∂η/∂x              = f(y)v - σ u
      ∂v/∂t + g∂η/∂y              = -f(y)u - σ v

    Initial condition: Gaussian perturbation in geostrophic balance.
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D Shallow Water TC3 (Planetary Rossby Wave) ..."
    end

    #---------------------------------------------------------------------------------
    # Solution variables: q = [η, u, v]
    #---------------------------------------------------------------------------------
    qvars    = ["η", "u", "v"]
    qoutvars = ["η", "u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        #---------------------------------------------------------------------------------
        # Physical parameters (TC3: Bishnu et al. 2024)
        #---------------------------------------------------------------------------------
        Lx   = 1.0e6     # domain length x [m]
        Ly   = 1.0e6     # domain length y [m]
        H0   = 1.0e3     # mean depth [m]
        g    = 9.81       # gravity [m/s²]
        f0   = 1.0e-4     # base Coriolis [s⁻¹]

        eta_hat  = 0.01   # perturbation amplitude [m]
        Rx_frac  = 0.10   # Gaussian radius fraction
        Ry_frac  = 0.10
        Rx = Rx_frac * Lx
        Ry = Ry_frac * Ly

        # Domain center
        comm_mpi = MPI.COMM_WORLD
        x0 = (MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm_mpi) +
               MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm_mpi)) / 2.0
        y0 = (MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm_mpi) +
               MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm_mpi)) / 2.0

        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            y = mesh.y[ip]

            # Gaussian perturbation
            gauss = exp(-((x - x0)^2 / (2.0 * Rx^2) + (y - y0)^2 / (2.0 * Ry^2)))
            η0 = eta_hat * gauss

            # Geostrophic initial velocities
            u0 =  (g / (f0 * Ry^2)) * (y - y0) * η0
            v0 = -(g / (f0 * Rx^2)) * (x - x0) * η0

            q.qn[ip, 1] = η0
            q.qn[ip, 2] = u0
            q.qn[ip, 3] = v0

            # Store initial state for reference
            q.qe[ip, 1] = η0
            q.qe[ip, 2] = u0
            q.qe[ip, 3] = v0
        end
    end

    if rank == 0
        @info " Initialize fields for 2D Shallow Water TC3 (Planetary Rossby Wave) ... DONE"
    end

    return q
end

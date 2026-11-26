function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """
    TC1: Coastal Kelvin Wave (Bishnu et al. 2024, Sec. 2.5)

    Linearized rotating shallow water equations, coast at x=0:
      ∂η/∂t + H₀(∂u/∂x + ∂v/∂y) = 0
      ∂u/∂t - f₀v + g∂η/∂x = 0
      ∂v/∂t + f₀u + g∂η/∂y = 0

    Exact solution (Cushman-Roisin & Beckers 2011, as used by Bishnu et al. 2024):
      u = 0
      v =  √(gH₀) F(y + ct) e^{-x/R}
      η = -H₀     F(y + ct) e^{-x/R}
    with c = √(gH₀), R = c/f₀ the barotropic Rossby radius of deformation,
    and F(y) a superposition of two sinusoidal meridional modes:
      F(y) = η̂₁ sin(ky₁ y) + η̂₂ sin(ky₂ y),   η̂₂ = 2η̂₁,  ky₂ = 2ky₁.
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D Shallow Water TC1 (Coastal Kelvin Wave) ...")
    end

    #---------------------------------------------------------------------------------
    # Solution variables: q = [η, u, v]
    #---------------------------------------------------------------------------------
    qvars    = ["η", "u", "v"]
    qoutvars = ["η", "u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        #---------------------------------------------------------------------------------
        # Get domain size from the actual mesh
        #---------------------------------------------------------------------------------
        Ly = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm) - MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)

        if rank == 0
            println("  Domain size: Ly = $Ly m")
        end

        #---------------------------------------------------------------------------------
        # Physical parameters (Bishnu et al. 2024, Sec. 2.5)
        #---------------------------------------------------------------------------------
        H0 = 1.0e3     # mean depth [m]
        g  = 10.0      # gravity [m/s²]
        f0 = 1.0e-4    # Coriolis parameter [s⁻¹]
        c  = sqrt(g * H0)  # gravity-wave speed [m/s]
        R  = c / f0        # barotropic Rossby radius of deformation [m]

        # Mode 1
        η̂1 = 1.0e-4           # [m]
        ky1 = 2π / Ly

        # Mode 2: twice the amplitude and wavenumber of mode 1
        η̂2 = 2.0 * η̂1
        ky2 = 2.0 * ky1

        t0 = inputs[:tinit]

        if rank == 0
            println("  c = $c m/s, R = $R m")
            println("  Mode 1: ky=$(ky1), η̂=$(η̂1)")
            println("  Mode 2: ky=$(ky2), η̂=$(η̂2)")
        end

        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            y = mesh.y[ip]

            F  = η̂1 * sin(ky1 * (y + c * t0)) + η̂2 * sin(ky2 * (y + c * t0))
            decay = exp(-x / R)

            η = -H0 * F * decay
            u = 0.0
            v = c * F * decay

            q.qn[ip, 1] = η
            q.qn[ip, 2] = u
            q.qn[ip, 3] = v

            # Store exact solution at t=t0 for reference
            q.qe[ip, 1] = q.qn[ip, 1]
            q.qe[ip, 2] = q.qn[ip, 2]
            q.qe[ip, 3] = q.qn[ip, 3]
        end
    end

    if rank == 0
        println(" Initialize fields for 2D Shallow Water TC1 (Coastal Kelvin Wave) ... DONE")
    end

    return q
end

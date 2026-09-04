function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """
    ReducedGravity: Wind-driven double-gyre circulation in a rectangular basin with flat bathymetry.

    Sec 3.1 of Choi et al. (2004).

    Conventions (see user_flux.jl): z = 0 is the flat bottom, 
    q[1] = H is the water depth measured from the top of the bathymetry, 
     
    Initial state:

      State vector q = [H, Hu, Hv] with

          H(x, y, 0) = 150
          u(x, y, 0) = 0
          v(x, y, 0) = 0

      The reference state qe = [150, 0, 0], stores the lake at rest. 
      user_flux.jl/user_source.jl advance the perturbation about
      it (pressure flux g'(H² - He²)/2, source -g'(H - He)∇Hb) so that qe is an
      exact equilibrium of the discrete operators. There is no initial dry states in this case
      but they can develop as the wind stress drives the flow and the water piles up in the corners of the domain. 
      The wet/dry threshold is set to 1e-3 m.
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D NLSWE reduced gravity double-gyre...")
    end

    qvars    = ["H", "Hu", "Hv"]
    qoutvars = ["H", "Hu", "Hv"] 
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())

        for ip = 1:mesh.npoin
            #x = mesh.x[ip]
            #y = mesh.y[ip]

            H_water = 150.0
            q.qn[ip, 1] = H_water
            q.qn[ip, 2] = 0.0
            q.qn[ip, 3] = 0.0

            # reference state used by the kernels: lake at rest (well-balanced split)            
            q.qe[ip, 1] = H_water
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 0.0
        end
    end

    if rank == 0
        println(" Initialize fields for 2D NLSWE reduced gravity double-gyre ... DONE")
    end

    return q
end

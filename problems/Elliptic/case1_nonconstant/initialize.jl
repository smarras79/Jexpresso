function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0 println(" Initialize fields for the non-constant diffusivity test .............. ") end
    #---------------------------------------------------------------------------------
    # Solution variables:
    #---------------------------------------------------------------------------------
    qvars    = ["u"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            y = mesh.y[ip]
            # zero initial guess; exact field qe = manufactured solution u_ex
            # (used for the post-solve L2 error check in standard_linsolve!).
            q.qn[ip,1] = 0.0
            q.qe[ip,1] = ncd_u(x, y)
        end
    else
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y; ndrange = mesh.npoin)
    end

    if rank == 0 println(" Initialize fields for the non-constant diffusivity test .............. DONE ") end

    return q
end

@kernel function initialize_gpu!(qn, qe, x, y)
    T = eltype(qn)

    ip = @index(Global, Linear)
    xip = x[ip]
    yip = y[ip]
    qn[ip,1] = T(0.0)
    qe[ip,1] = T(sin(xip) * cos(yip))   # mirror of ncd_u with kx=ky=1
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0 println(" Initialize fields for the 2D periodic Poisson problem ........................ ") end

    qvars = ["u"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    # qn: zero initial guess (the solve is direct). qe: the exact solution,
    # which the SEM solve's automatic L2-error check compares against.
    if (inputs[:backend] == CPU())
        for ip = 1:mesh.npoin
            x = mesh.coords[1,ip]
            y = mesh.coords[2,ip]
            q.qn[ip,1] = 0.0
            q.qe[ip,1] = user_fft_exact(x, y)
        end
    else
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, @view(mesh.coords[1,:]), @view(mesh.coords[2,:]); ndrange = mesh.npoin)
    end

    if rank == 0 println(" Initialize fields for the 2D periodic Poisson problem ........................ DONE") end

    return q
end

@kernel function initialize_gpu!(qn, qe, x, y)
    ip = @index(Global, Linear)
    xip = x[ip]
    yip = y[ip]
    qn[ip,1] = 0.0
    qe[ip,1] = user_fft_exact(xip, yip)
end

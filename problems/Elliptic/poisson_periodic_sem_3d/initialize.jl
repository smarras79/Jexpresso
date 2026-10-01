function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0 println(" Initialize fields for the 3D periodic Poisson problem ........................ ") end

    qvars = ["u"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    # qn: zero initial guess (the solve is direct). qe: the exact solution,
    # which the SEM solve's automatic error check compares against.
    inputs[:backend] == CPU() ||
        error(" # poisson_periodic_sem_3d: the periodic linear solves are CPU-only.")
    for ip = 1:mesh.npoin
        q.qn[ip,1] = 0.0
        q.qe[ip,1] = user_fft_exact(mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip])
    end

    if rank == 0 println(" Initialize fields for the 3D periodic Poisson problem ........................ DONE") end

    return q
end

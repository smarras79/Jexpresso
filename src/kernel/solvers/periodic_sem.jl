#=============================================================================
 periodic_sem.jl — direct SEM solve of -∇²u = f on a fully PERIODIC mesh

 How periodicity reaches a continuous-Galerkin SEM operator in Jexpresso: the
 local node numbering is NOT merged. restructure4periodicity_2D gives the two
 copies of a seam node (x = xmin and x = xmax, …) the same GLOBAL id in
 mesh.ip2gip, and the global DSS (assemble_mpi! through g_dss_cache) sums the
 contributions of all local nodes sharing a global id. The time-dependent
 solvers go through that DSS; the sparse Laplacian of the linear solve does not:

   • sem.matrix.L  is assembled on the LOCAL nodes (DSS_laplace_sparse), i.e.
                   the stiffness matrix of the box with natural (Neumann)
                   boundaries — the seam copies are not coupled;
   • sem.matrix.M  (lumped) HAS been summed over each periodic class by
                   DSS_global_mass!: every copy holds the class's total mass.

 So the periodic system has to be formed here. With the classes c = 1..m of
 local nodes sharing a global id and the npoin × m incidence matrix P
 (P[ip, class(ip)] = 1),

       K = Pᵀ L P          periodic stiffness (partner rows/columns summed)
       b_c = M[rep_c]·f(rep_c)   = Σ_{ip∈c} M_local[ip] f(ip) for periodic f

 THE NULL SPACE. K is symmetric positive SEMI-definite with K·1 = 0: the
 periodic problem determines u only up to a constant, and has a solution only
 if the RHS is orthogonal to the constants (Σ b = 0 ⇔ ∫f = 0). The solve
   1. projects b onto the range: b ← b - (Σb / Σw) w,  w_c = M[rep_c]
      (i.e. solves for f - mean(f); the removed mean is reported),
   2. pins the first class (u_1 = 0) and solves the remaining SPD system,
   3. shifts u to zero M-weighted mean — the same gauge as the FFT solver,
      so the two return comparable solutions.

 periodic_sem_system builds (K, b, …) on its own so that other solvers of the
 SAME SEM system (AlgebraicMultigrid) reuse it unchanged.

 Scope: 2D, serial, every boundary edge periodic. A mesh that is periodic in
 one direction only (Dirichlet on the other boundary) is refused, not solved
 wrongly.
=============================================================================#

const _PERIODIC_EDGE_TAGS = ("periodicx", "periodicz", "periodic1", "periodic3")

"""
    sem_mesh_is_periodic(mesh) -> Bool

`true` when every boundary edge of the (2D) mesh is periodic. Errors on a mesh
that mixes periodic and non-periodic boundary edges, which the linear solves
do not support yet.
"""
function sem_mesh_is_periodic(mesh)
    mesh.nsd == 2 || return false
    nb = Int(mesh.nedges_bdy)
    nb == 0 && return false
    types = @view mesh.bdy_edge_type[1:nb]
    nper  = count(t -> t in _PERIODIC_EDGE_TAGS, types)
    nper == 0  && return false
    nper == nb && return true
    error(" # standard_linsolve!: the mesh has $nper periodic and $(nb - nper) non-periodic ",
          "boundary edges. Linear solves on partially periodic meshes are not supported yet.")
end

"""
    periodic_sem_system(sem, f_nodal) -> (; K, b, cls, rep, w, fmean)

The SEM system of -∇²u = f on a fully periodic mesh, reduced to one unknown
per periodic node class (see the file header). `f_nodal[ip]` is f at local
node ip. Returns the reduced stiffness `K` (m × m, sparse, singular),
the RHS `b` projected onto the range of K, the class of every local node
`cls`, one representative local node per class `rep`, the class weights `w`
(lumped mass) and the M-weighted mean of f that the projection removed.
"""
function periodic_sem_system(sem, f_nodal::AbstractVector)
    mesh  = sem.mesh
    npoin = Int(mesh.npoin)
    MPI.Comm_size(get_mpi_comm()) == 1 ||
        error(" # periodic_sem_system: periodic linear solves are serial only.")
    M = sem.matrix.M
    M isa AbstractVector ||
        error(" # periodic_sem_system: needs the lumped (diagonal) mass matrix.")

    # periodic classes = local nodes sharing a global id
    gip = @view mesh.ip2gip[1:npoin]
    ug  = unique(gip)
    cls = Int.(indexin(gip, ug))
    m   = length(ug)
    rep = zeros(Int, m)
    @inbounds for ip = npoin:-1:1
        rep[cls[ip]] = ip
    end

    P = sparse(1:npoin, cls, ones(Float64, npoin), npoin, m)
    K = sparse(P' * (sem.matrix.L * P))

    w = Float64[M[rep[c]] for c in 1:m]              # class mass (already summed)
    b = Float64[M[rep[c]] * f_nodal[rep[c]] for c in 1:m]
    fmean = sum(b) / sum(w)
    b .-= fmean .* w                                   # Σ b = 0: in the range of K

    return (; K, b, cls, rep, w, fmean)
end

"""
    periodic_sem_direct_solve(K, b, w) -> u (per class)

Zero-mean (w-weighted) solution of the singular system K u = b (Σb = 0) by
pinning the first unknown and a sparse direct solve of the rest.
"""
function periodic_sem_direct_solve(K, b, w)
    u = zeros(Float64, length(b))
    u[2:end] = K[2:end, 2:end] \ b[2:end]
    u .-= sum(w .* u) / sum(w)
    return u
end

function periodic_sem_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    inputs[:backend] == CPU() ||
        error(" # periodic_sem_linsolve!: periodic linear solves are CPU-only.")
    mesh  = sem.mesh
    npoin = Int(mesh.npoin)

    f = Vector{Float64}(undef, npoin)
    for ip = 1:npoin
        f[ip] = user_source!(0.0,
                             params.qp.qn[ip],
                             params.qp.qe[ip],
                             mesh.npoin,
                             inputs[:CL], inputs[:SOL_VARS_TYPE];
                             neqs=1, x=mesh.coords[1,ip], y=mesh.coords[2,ip],
                             xmax=mesh.xmax, xmin=mesh.xmin,
                             ymax=mesh.ymax, ymin=mesh.ymin)
    end

    sys = periodic_sem_system(sem, f)
    println(YELLOW_FG(string(" # Periodic SEM system: ", npoin, " local nodes → ",
                             length(sys.b), " periodic unknowns (singular; zero-mean solution)")))
    if abs(sys.fmean) > 1e-10 * max(1.0, maximum(abs, f))
        println(string(" # periodic_sem_linsolve!: mean of f = ", sys.fmean,
                       " ≠ 0; solved the projected problem -∇²u = f - mean(f) ",
                       "(periodic compatibility condition)."))
    end

    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage ..............")))
    uc = jx_robust_solve("direct SEM (Ax=b, periodic)",
                         () -> periodic_sem_direct_solve(sys.K, sys.b, sys.w);
                         robust  = get(inputs, :lbenchmark_solve, true),
                         seconds = Float64(get(inputs, :EL_timing_seconds, 2.0)))
    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage .............. DONE")))

    sol = TFloat.(uc[sys.cls])                 # every seam copy gets its class value

    # Error check against qe with quadrature weights that count each periodic
    # class ONCE: sem.matrix.M holds the class total on every copy.
    ncopies = zeros(Int, length(sys.b))
    for c in sys.cls
        ncopies[c] += 1
    end
    wq = Float64[sem.matrix.M[ip] / ncopies[sys.cls[ip]] for ip in 1:npoin]
    print_solution_L2_error(sol, params.qp.qe, wq, npoin; label="direct periodic SEM solve")

    args = (params.SD, sol, params.uaux, 1, 1,
            mesh, nothing,
            nothing, nothing,
            0.0, 0.0, 0.0,
            OUTPUT_DIR, inputs,
            params.qp.qvars,
            params.qp.qoutvars,
            inputs[:outformat])
    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe, metrics=params.metrics)

    return nothing
end

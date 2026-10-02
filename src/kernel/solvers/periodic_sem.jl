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
 SAME SEM system (AlgebraicMultigrid) reuse it unchanged. The direct solve is
 split into a one-time sparse factorisation (periodic_sem_factorize, setup)
 and the triangular solves (periodic_sem_direct_solve, the timed solve step).

 Scope: 2D and 3D, serial, every boundary edge (2D) or face (3D) periodic.
 A mesh that is periodic in one direction only (Dirichlet on another
 boundary) is refused, not solved wrongly. In 3D sem.matrix.L comes from
 DSS_laplace_sparse_3D (element_matrices.jl).
=============================================================================#

const _PERIODIC_EDGE_TAGS = ("periodicx", "periodicz", "periodic1", "periodic3")
const _PERIODIC_FACE_TAGS = ("periodicx", "periodicy", "periodicz", "periodic1", "periodic2", "periodic3")

"""
    sem_mesh_is_periodic(mesh) -> Bool

`true` when every boundary edge (2D) or face (3D) of the mesh is periodic.
Errors on a mesh that mixes periodic and non-periodic boundaries, which the
linear solves do not support yet.
"""
function sem_mesh_is_periodic(mesh)
    if mesh.nsd == 2
        nb = Int(mesh.nedges_bdy); tags = _PERIODIC_EDGE_TAGS; what = "edges"
        nb == 0 && return false
        types = @view mesh.bdy_edge_type[1:nb]
    elseif mesh.nsd == 3
        nb = Int(mesh.nfaces_bdy); tags = _PERIODIC_FACE_TAGS; what = "faces"
        nb == 0 && return false
        types = @view mesh.bdy_face_type[1:nb]
    else
        return false
    end
    nper  = count(t -> t in tags, types)
    nper == 0  && return false
    nper == nb && return true
    error(" # standard_linsolve!: the mesh has $nper periodic and $(nb - nper) non-periodic ",
          "boundary $what. Linear solves on partially periodic meshes are not supported yet.")
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
    K = P' * (sem.matrix.L * P)                   # SparseMatrixCSC; exactly symmetric

    w = Float64[M[rep[c]] for c in 1:m]              # class mass (already summed)
    b = Float64[M[rep[c]] * f_nodal[rep[c]] for c in 1:m]
    fmean = sum(b) / sum(w)
    b .-= fmean .* w                                   # Σ b = 0: in the range of K

    return (; K, b, cls, rep, w, fmean)
end

"""
    periodic_sem_factorize(K) -> F

Sparse factorisation of the periodic system with its first unknown pinned
(K[2:end, 2:end], symmetric positive definite): the one-time setup of the
direct solve. `factorize` picks what `K \\ b` would use — a sparse Cholesky
factorisation for this exactly symmetric matrix, with CHOLMOD's own
fill-reducing ordering (AMD). `ordering = :metis` (deck: :sparse_ordering =>
"metis") factorises with METIS nested dissection instead, whose fill grows
like n^(4/3) in 3D against markedly faster for AMD.
"""
periodic_sem_factorize(K, ordering::Symbol = :cholmod) =
    ordering === :metis ? (Kp = K[2:end, 2:end]; cholesky(Symmetric(Kp); perm = jx_metis_perm(Kp))) :
                          factorize(K[2:end, 2:end])

"""
    periodic_sem_direct_solve(F, b, w) -> u (per class)

Zero-mean (w-weighted) solution of the singular system K u = b (Σb = 0) from
the factorisation `F = periodic_sem_factorize(K)`: the pinned first unknown is
0, the rest come from the triangular solves, then the w-weighted mean is
removed.
"""
function periodic_sem_direct_solve(F, b, w)
    u = zeros(Float64, length(b))
    u[2:end] = F \ b[2:end]
    u .-= sum(w .* u) / sum(w)
    return u
end

# f at every local node, from the case's user_source! (function barrier: all
# arguments arrive concretely typed, so the loop is compiled for them).
function _periodic_sem_rhs(coords::AbstractMatrix, qn, qe, npoin::Int, CL, SV,
                           xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64)
    f = Vector{Float64}(undef, npoin)
    for ip = 1:npoin
        f[ip] = user_source!(0.0, qn[ip], qe[ip], npoin, CL, SV;
                             neqs=1, x=coords[1,ip], y=coords[2,ip],
                             xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin)
    end
    return f
end

function _periodic_sem_rhs(coords::AbstractMatrix, qn, qe, npoin::Int, CL, SV,
                           xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64,
                           zmin::Float64, zmax::Float64)
    f = Vector{Float64}(undef, npoin)
    for ip = 1:npoin
        f[ip] = user_source!(0.0, qn[ip], qe[ip], npoin, CL, SV;
                             neqs=1, x=coords[1,ip], y=coords[2,ip], z=coords[3,ip],
                             xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin, zmax=zmax, zmin=zmin)
    end
    return f
end

function periodic_sem_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    inputs[:backend] == CPU() ||
        error(" # periodic_sem_linsolve!: periodic linear solves are CPU-only.")
    mesh  = sem.mesh
    npoin = Int(mesh.npoin)

    f = jx_phase(:rhs) do
        # St_mesh fields and the inputs Dict are untyped: hand the concrete
        # values to the loop through a function barrier (_periodic_sem_rhs).
        mesh.nsd == 3 ?
            _periodic_sem_rhs(mesh.coords, params.qp.qn, params.qp.qe, npoin,
                              inputs[:CL], inputs[:SOL_VARS_TYPE],
                              Float64(mesh.xmin), Float64(mesh.xmax),
                              Float64(mesh.ymin), Float64(mesh.ymax),
                              Float64(mesh.zmin), Float64(mesh.zmax)) :
            _periodic_sem_rhs(mesh.coords, params.qp.qn, params.qp.qe, npoin,
                              inputs[:CL], inputs[:SOL_VARS_TYPE],
                              Float64(mesh.xmin), Float64(mesh.xmax),
                              Float64(mesh.ymin), Float64(mesh.ymax))
    end

    # setup: the periodic reduction (which also forms b = M f)
    sys = jx_phase(:setup) do
        jx_phase(() -> periodic_sem_system(sem, f), :reduce)
    end
    println(YELLOW_FG(string(" # Periodic SEM system: ", npoin, " local nodes → ",
                             length(sys.b), " periodic unknowns (singular; zero-mean solution)")))
    if abs(sys.fmean) > 1e-10 * max(1.0, maximum(abs, f))
        println(string(" # periodic_sem_linsolve!: mean of f = ", sys.fmean,
                       " ≠ 0; solved the projected problem -∇²u = f - mean(f) ",
                       "(periodic compatibility condition)."))
    end

    # Three solves of the same periodic SEM system, selected by the deck:
    #   :lstatic_condensation => true   element-learning static condensation
    #                                   (T^ie from the SEM matrix); skeleton by
    #                                   :EL_skeleton_solver ("direct" | "amg")
    #   :linsolve_amg         => true   AMG-preconditioned CG on the full system
    #   default                         sparse direct (factorise + triangular solves)
    if get(inputs, :lstatic_condensation, false)
        uc = periodic_sem_sc_solve(sem, sys, inputs)
        label = string("static condensation (", el_skeleton_options(inputs).skeleton_solver, ")")
    elseif _pmg_kind(inputs) !== :none
        uc = periodic_sem_pmg_solve(sem, sys, inputs)
        label = string("p-multigrid CG periodic SEM solve (", _pmg_kind(inputs), ")")
    elseif get(inputs, :linsolve_amg, false)
        uc = periodic_sem_amg_solve(sys, inputs)
        label = "AMG-CG periodic SEM solve"
    else
        F = jx_phase(:setup) do
            jx_phase(() -> periodic_sem_factorize(sys.K, jx_sparse_ordering(inputs)), :factorize)
        end
        println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage ..............")))
        uc = jx_robust_solve("direct SEM (triangular solves, periodic)",
                             () -> periodic_sem_direct_solve(F, sys.b, sys.w);
                             robust  = get(inputs, :lbenchmark_solve, true),
                             seconds = Float64(get(inputs, :EL_timing_seconds, 2.0)))
        println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage .............. DONE")))
        label = "direct periodic SEM solve"
    end

    sol = TFloat.(uc[sys.cls])                 # every seam copy gets its class value

    # Error check against qe with quadrature weights that count each periodic
    # class ONCE: sem.matrix.M holds the class total on every copy.
    ncopies = zeros(Int, length(sys.b))
    for c in sys.cls
        ncopies[c] += 1
    end
    wq = Float64[sem.matrix.M[ip] / ncopies[sys.cls[ip]] for ip in 1:npoin]
    print_solution_L2_error(sol, params.qp.qe, wq, npoin; label=label)

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

"""
    periodic_sem_amg_solve(sys, inputs) -> u (per class)

AMG-preconditioned CG on the full periodic SEM system: the first unknown is
pinned (K[2:end,2:end] is SPD), the result shifted to zero M-weighted mean.
Records :setup (AMG hierarchy) and :solve (CG).
"""
function periodic_sem_amg_solve(sys, inputs)
    opts = jx_amg_options(inputs)
    S  = jx_phase(() -> jx_amg_setup(sys.K[2:end, 2:end]; method = opts.method), :setup)
    uc = jx_time_solve("AMG-CG on the full periodic SEM system", () -> begin
             u = zeros(Float64, length(sys.b))
             u[2:end] = jx_amg_solve(S, sys.b[2:end]; rtol = opts.rtol, itmax = opts.itmax)
             u .-= sum(sys.w .* u) / sum(sys.w)
             u
         end)
    _el_print_amg_stats()
    return uc
end

# :linsolve_pmg => "none" (default) | "amg" | "gmg"  (p-multigrid CG, pmultigrid.jl)
function _pmg_kind(inputs)
    k = Symbol(lowercase(string(get(inputs, :linsolve_pmg, "none"))))
    k in (:none, :amg, :gmg) || error(" # :linsolve_pmg => \"$k\"; expected \"none\", \"amg\" or \"gmg\".")
    return k
end

"""
    periodic_sem_pmg_solve(sem, sys, inputs) -> u (per class)

CG on the full periodic SEM system, preconditioned by the p-multigrid V-cycle
(jx_pmg_setup: SEM levels N, N/2, …, 1, then AMG ("amg") or geometric
h-multigrid ("gmg") on the p = 1 level); the result shifted to zero
M-weighted mean. Records :setup (the hierarchy) and :solve (CG).
:pmg_degree (Chebyshev degree, default 3) and :pmg_lower (lower end of the
smoothed spectrum, fraction of λ_max, default 0.25) tune the smoother.
"""
function periodic_sem_pmg_solve(sem, sys, inputs)
    opts = jx_amg_options(inputs)
    M  = jx_phase(:setup) do
             jx_pmg_setup(sem, sys.K, sys.cls; coarse = _pmg_kind(inputs),
                          degree = Int(get(inputs, :pmg_degree, 3)),
                          lower = Float64(get(inputs, :pmg_lower, 0.25)),
                          amg_method = opts.method)
         end
    uc = jx_time_solve("p-multigrid CG on the full periodic SEM system", () -> begin
             u = jx_pmg_cg(sys.K, sys.b, M; rtol = opts.rtol, itmax = opts.itmax)
             u .-= sum(sys.w .* u) / sum(sys.w)
             u
         end)
    st = JX_AMG_STATS[]
    println(GREEN_FG(string(" # p-multigrid (", st.method, ", orders ", join(M.orders, "→"),
                            ", ", st.levels, " levels): CG converged in ", st.iters,
                            " iterations, relative residual ", st.rel_resid)))
    return uc
end

# 3D element connectivity in periodic-class numbering, each element's boundary
# nodes first (any local index 1 or ngl), its (ngl-2)³ interior nodes last,
# each group in connijk's (i fastest) order: the layout elementLearning_Axb!
# expects. (2D uses mesh.conn, which is already boundary-first.)
function _sc_conn_3d(connijk::AbstractArray{<:Integer, 4}, cls::Vector{Int}, nelem::Int, ngl::Int)
    order = NTuple{3, Int}[]
    onb(i) = i == 1 || i == ngl
    for k = 1:ngl, j = 1:ngl, i = 1:ngl
        (onb(i) || onb(j) || onb(k)) && push!(order, (i, j, k))
    end
    for k = 2:ngl-1, j = 2:ngl-1, i = 2:ngl-1
        push!(order, (i, j, k))
    end
    conn = Matrix{Int}(undef, nelem, ngl^3)
    @inbounds for (a, (i, j, k)) in enumerate(order), e = 1:nelem
        conn[e, a] = cls[connijk[e, i, j, k]]
    end
    return conn
end

"""
    periodic_sem_sc_solve(sem, sys, inputs) -> u (per class)

The element-learning STATIC CONDENSATION (elementLearning_Axb! with the local
operators T^ie computed from the SEM matrix, not predicted by the network)
applied to the periodic SEM system.

elementLearning_Axb! works on a mesh-like description (2D or 3D): the element
connectivity `conn` (element-boundary nodes first, interior last), the
Dirichlet set Γ, the internal skeleton ∂O, the skeleton ∂τ = Γ ∪ ∂O and the
strictly interior nodes Io. For the periodic system all of them are expressed
in PERIODIC-CLASS numbering (one unknown per class, the numbering of sys.K):
  conn = class of every local node of mesh.conn
  Γ    = ∅                    (no Dirichlet boundary)
  ∂O   = ∂τ = classes of the element-boundary nodes
  Io   = classes of the element-interior nodes
They are set on a separate St_mesh so the real mesh is not touched. With Γ
empty the condensed system B_{∂O,∂O} has the constants in its null space;
el_skeleton_solve pins one unknown and the result is shifted here to zero
M-weighted mean — the gauge of every other periodic solve.
"""
function periodic_sem_sc_solve(sem, sys, inputs)
    mesh  = sem.mesh
    ngl   = Int(mesh.ngl);  nelem = Int(mesh.nelem)
    nsd   = Int(mesh.nsd)
    npel  = ngl^nsd
    nint  = (ngl - 2)^nsd
    nb    = npel - nint
    m     = length(sys.b)
    opts  = el_skeleton_options(inputs)

    pm, EL, wbuf = jx_phase(:sc_alloc) do
        conn  = nsd == 3 ? _sc_conn_3d(mesh.connijk, sys.cls, nelem, ngl) :
                           Int[sys.cls[mesh.conn[e, a]] for e in 1:nelem, a in 1:npel]
        skel  = sort!(unique(vec(conn[:, 1:nb])))
        inter = vec(conn[:, nb+1:end])
        # mesh.conn must list each element's boundary nodes first: every
        # interior class then belongs to exactly one element and to no skeleton
        (allunique(inter) && isempty(intersect(skel, inter)) &&
         length(skel) + length(inter) == m) ||
            error(" # periodic_sem_sc_solve: mesh.conn is not boundary-first; cannot condense.")
        pm = typeof(mesh)(; SD = mesh.SD)
        pm.conn  = conn;       pm.ngl   = ngl;   pm.nelem = nelem
        pm.Γ     = Int[];      pm.lengthΓ  = 0
        pm.∂O    = skel;       pm.length∂O = length(skel)
        pm.∂τ    = copy(skel); pm.length∂τ = length(skel)       # ∂τ = Γ ∪ ∂O
        pm.Io    = inter;      pm.lengthIo = length(inter)
        pm.O     = vcat(inter, skel);  pm.lengthO = m
        EL   = allocate_elemLearning(nelem, ngl, pm.length∂O, pm.length∂τ, 0,
                                     Float64, CPU(); Nsamp = 1, lEL_Sample = true, nsd = nsd,
                                     blocks = opts.kernel === :legacy)
        # no model (NNfile = nothing): the inference buffers, and the skeleton
        # submatrix they would copy, are not allocated
        wbuf = EL_WorkBuffers(pm, sys.K, spzeros(0, 0), npel, nint, nb, nothing)
        pm, EL, wbuf
    end

    println(YELLOW_FG(string(" # Static condensation (T^ie from the SEM matrix) on the periodic system: ",
                             m, " unknowns → ", pm.length∂O, " skeleton unknowns; skeleton solver: ",
                             opts.skeleton_solver, " ..............")))
    u = zeros(Float64, m, 1)
    elementLearning_Axb!(u, nothing, pm, sys.K, reshape(copy(sys.b), m, 1), EL,
                         zeros(Float64, 1, npel), nothing, nothing,
                         zeros(pm.length∂O), zeros(0), wbuf;
                         skeleton_solver = opts.skeleton_solver,
                         amg_method = opts.amg_method, amg_rtol = opts.amg_rtol,
                         skeleton_ordering = opts.ordering,
                         record_tensors = false, sc_kernel = opts.kernel,
                         sc_amg = opts.amg_mode, amg_itmax = opts.amg_itmax,
                         tp1d = nsd == 3 ? (sc_gll_stiffness(sem.basis.dψ, sem.ω), sem.ω) : nothing)
    uc = vec(u)
    uc .-= sum(sys.w .* uc) / sum(sys.w)
    _el_sc_record_phases!()
    opts.skeleton_solver === :amg && _el_print_amg_stats()
    println(YELLOW_FG(string(" # Static condensation ...................................... DONE")))
    return uc
end

#=============================================================================
 tools/poisson3d_benchmark/mpi/poisson3d_mpi.jl

 Distributed-memory (MPI) version of the 3D periodic Poisson solver
 comparison of ../poisson3d.jl:
     -∇²u = f  on [0,2π]³,  u periodic in x, y, z,
 the same SEM discretisation (ne³ elements of order N, GLL collocation,
 lumped mass), the same exact solution and right-hand side, solved by

   :mumps       MUMPS sparse Cholesky (LDLᵀ for SPD), distributed matrix
                input, METIS nested-dissection ordering (on the host: the
                MUMPS_jll binary has no parallel ordering)
   :boomeramg   hypre BoomerAMG (one V-cycle) preconditioned CG (hypre PCG)
   :jacobi      Jacobi (diagonal) preconditioned CG (hypre PCG + DiagScale),
                the CEED BP5 solver

 THE OPERATOR. With GLL collocation the SEM stiffness matrix is exactly
     K = Mz⊗My⊗Kx + Mz⊗Ky⊗Mx + Kz⊗My⊗Mx,    M = Mz⊗My⊗Mx (lumped),
 with the periodic 1D matrices K1 (stiffness) and m1 (lumped mass), the same
 in every direction. So the row of node (i,j,k) is known from K1 and m1
 alone: every rank builds its own rows, with no communication.

 THE PARTITION. The Ng³ nodes (Ng = ne·N) are split in boxes over a
 px×py×pz process grid (MPI.Dims_create): rank r owns the nodes of its box,
 numbered contiguously (global row range of rank r = its box). Boxes keep
 the halo of the stencil (±N nodes in each direction) small compared with the
 owned nodes, unlike slabs.

 SINGULAR SYSTEM. As in the serial benchmark, the first unknown (node (0,0,0),
 global row 1) is pinned: its row and column are replaced by the identity,
 b₁ = 0. This is the same system as the serial K[2:end,2:end] solve. The
 solution is shifted to zero M-weighted mean afterwards.

 The 1D ingredients (LGL nodes and weights by Kopriva's algorithm, Lagrange
 derivative matrix) are computed here, not taken from Jexpresso, so that the
 MPI ranks do not load Jexpresso (hundreds of MB per rank). verify_mpi.jl
 checks that the errors equal those of the serial Jexpresso-based benchmark.
=============================================================================#
module P3DMPI

using MPI, LinearAlgebra, SparseArrays, Printf, Libdl
using HYPRE
using HYPRE.LibHYPRE                        # exports every HYPRE_* function
using HYPRE.LibHYPRE: @check, HYPRE_BigInt, HYPRE_Int, HYPRE_Complex, HYPRE_jll
import MUMPS

export run_config, MPI_SOLVERS

const MPI_SOLVERS = (:mumps, :boomeramg, :jacobi)
const LBOX = 2π
const R_DEFAULT = 0.5

# ---------------------------------------------------------------------------
# Exact solution (as ../poisson3d.jl, d = 3)
#   u = A ( p(x)p(y)p(z) - (c²-1)^(-3/2) ),  p(s) = 1/(c - cos s)
# ---------------------------------------------------------------------------
_c(r) = (r + 1 / r) / 2
_p(s, c)   = 1 / (c - cos(s))
_pss(s, c) = (q = c - cos(s); -cos(s) / q^2 + 2 * sin(s)^2 / q^3)
amp(c) = 1 / ((c - 1)^(-3) - (c^2 - 1)^(-3 / 2))
umean(c) = (c^2 - 1)^(-3 / 2)

# ---------------------------------------------------------------------------
# 1D GLL ingredients (Kopriva, "Implementing Spectral Methods for PDEs",
# Algorithms 22, 24, 25, 37: the algorithms Jexpresso's basis_structs_ξ_ω!
# and LagrangeInterpolatingPolynomials_classic implement)
# ---------------------------------------------------------------------------
function _qandL(N, x)                          # q = L_{N+1} - L_{N-1}, q', L_N
    Lm2 = 1.0; Lm1 = x; dLm2 = 0.0; dLm1 = 1.0
    L = x; dL = 1.0
    for k in 2:N
        L = (2k - 1) / k * x * Lm1 - (k - 1) / k * Lm2
        dL = dLm2 + (2k - 1) * Lm1
        Lm2, Lm1 = Lm1, L; dLm2, dLm1 = dLm1, dL
    end
    k = N + 1
    Lp1 = (2k - 1) / k * x * L - (k - 1) / k * Lm2
    dLp1 = dLm2 + (2k - 1) * Lm1
    return Lp1 - Lm2, dLp1 - dLm2, L
end

function lgl(N::Int)
    ξ = zeros(N + 1); ω = zeros(N + 1)
    if N == 1
        return [-1.0, 1.0], [1.0, 1.0]
    end
    ξ[1] = -1; ω[1] = 2 / (N * (N + 1))
    ξ[N+1] = 1; ω[N+1] = ω[1]
    for j in 1:(N + 1) ÷ 2 - 1
        x = -cos((j + 0.25) * π / N - 3 / (8N * π * (j + 0.25)))
        for _ in 1:100
            q, dq, _ = _qandL(N, x)
            Δ = -q / dq; x += Δ
            abs(Δ) <= 4eps() * abs(x) && break
        end
        _, _, L = _qandL(N, x)
        ξ[j+1] = x; ω[j+1] = 2 / (N * (N + 1) * L^2)
        ξ[N+1-j] = -x; ω[N+1-j] = ω[j+1]
    end
    if iseven(N)
        _, _, L = _qandL(N, 0.0)
        ξ[N÷2+1] = 0.0; ω[N÷2+1] = 2 / (N * (N + 1) * L^2)
    end
    return ξ, ω
end

# derivative of the Lagrange basis on the nodes ξ: D[i,k] = ψ_i'(ξ_k)
function lagrange_derivative(ξ)
    n = length(ξ)
    w = [1 / prod(ξ[j] - ξ[k] for k in 1:n if k != j) for j in 1:n]   # barycentric weights
    D = zeros(n, n)                                # D[k,i] = ψ_i'(ξ_k)
    for k in 1:n, i in 1:n
        i == k && continue
        D[k, i] = (w[i] / w[k]) / (ξ[k] - ξ[i])
    end
    for k in 1:n
        D[k, k] = -sum(D[k, i] for i in 1:n if i != k)
    end
    return Matrix(transpose(D))
end

"""
    sem1d(N, ne) -> (; K1, m1, x, n1)

Periodic 1D SEM operators on ne elements of order N over [0, 2π): assembled
stiffness K1 (exactly symmetric), lumped mass m1, node coordinates x
(as sem1d in ../poisson3d.jl).
"""
function sem1d(N::Int, ne::Int)
    ξ, ω = lgl(N)
    dψ = lagrange_derivative(ξ)
    h = LBOX / ne; J = h / 2
    Ke = [sum(ω[k] * dψ[i, k] * dψ[l, k] for k in 1:N+1) / J for i in 1:N+1, l in 1:N+1]
    Ke = (Ke + Ke') ./ 2
    n1 = ne * N
    I = Int[]; Jv = Int[]; V = Float64[]
    m1 = zeros(n1); x = zeros(n1)
    for e in 0:ne-1, i in 0:N, l in 0:N
        push!(I, mod(e * N + i, n1) + 1); push!(Jv, mod(e * N + l, n1) + 1); push!(V, Ke[i+1, l+1])
    end
    for e in 0:ne-1, i in 0:N
        m1[mod(e * N + i, n1) + 1] += J * ω[i+1]
        i < N && (x[e * N + i + 1] = e * h + (ξ[i+1] + 1) * J)
    end
    return (; K1 = sparse(I, Jv, V, n1, n1), m1, x, n1, ξ, ω)
end

# ---------------------------------------------------------------------------
# Box partition of the Ng³ nodes over a px×py×pz process grid
# ---------------------------------------------------------------------------
struct Partition
    Ng::Int
    dims::NTuple{3, Int}           # px, py, pz
    coords::NTuple{3, Int}         # this rank's block (0-based)
    blk::NTuple{3, Vector{Int}}    # node index (1-based) -> block (0-based), per direction
    loc::NTuple{3, Vector{Int}}    # node index -> index within its block (0-based)
    sz::NTuple{3, Vector{Int}}     # block -> number of nodes, per direction
    st::NTuple{3, Vector{Int}}     # block -> first node (1-based), per direction
    off::Vector{Int}               # rank -> first global row - 1
    ilower::Int; iupper::Int       # this rank's global rows (1-based, inclusive)
end

rank_of(P::Partition, b) = b[1] + P.dims[1] * (b[2] + P.dims[2] * b[3])

function Partition(Ng::Int, comm)
    np = MPI.Comm_size(comm); me = MPI.Comm_rank(comm)
    dims = Tuple(MPI.Dims_create(np, [0, 0, 0]))
    any(d -> d > Ng, dims) && error("P3DMPI: $np ranks need a grid of at least $(maximum(dims)) nodes per direction (have $Ng)")
    split(p) = (s = [Ng ÷ p + (b < Ng % p ? 1 : 0) for b in 0:p-1]; (s, cumsum([1; s[1:end-1]])))
    sz = ntuple(d -> split(dims[d])[1], 3); st = ntuple(d -> split(dims[d])[2], 3)
    blk = ntuple(d -> zeros(Int, Ng), 3); loc = ntuple(d -> zeros(Int, Ng), 3)
    for d in 1:3, b in 0:dims[d]-1, t in 0:sz[d][b+1]-1
        blk[d][st[d][b+1] + t] = b; loc[d][st[d][b+1] + t] = t
    end
    off = zeros(Int, np)
    acc = 0
    for bz in 0:dims[3]-1, by in 0:dims[2]-1, bx in 0:dims[1]-1
        r = bx + dims[1] * (by + dims[2] * bz)
        off[r+1] = acc
        acc += sz[1][bx+1] * sz[2][by+1] * sz[3][bz+1]
    end
    coords = (me % dims[1], (me ÷ dims[1]) % dims[2], me ÷ (dims[1] * dims[2]))
    nloc = prod(sz[d][coords[d]+1] for d in 1:3)
    return Partition(Ng, dims, coords, blk, loc, sz, st, off, off[me+1] + 1, off[me+1] + nloc)
end

# global row (1-based) of node (i,j,k) (1-based)
@inline function gid(P::Partition, i, j, k)
    bx = P.blk[1][i]; by = P.blk[2][j]; bz = P.blk[3][k]
    r = bx + P.dims[1] * (by + P.dims[2] * bz)
    return P.off[r+1] + 1 + P.loc[1][i] + P.sz[1][bx+1] * (P.loc[2][j] + P.sz[2][by+1] * P.loc[3][k])
end

# this rank's nodes, in the order of its global rows
function owned_nodes(P::Partition)
    r = ntuple(d -> P.st[d][P.coords[d]+1]:(P.st[d][P.coords[d]+1] + P.sz[d][P.coords[d]+1] - 1), 3)
    return [(i, j, k) for k in r[3] for j in r[2] for i in r[1]]
end

# ---------------------------------------------------------------------------
# Local rows of the pinned system, in CSR form (global column indices),
# the right-hand side and the exact solution at the owned nodes
# ---------------------------------------------------------------------------
function local_system(P::Partition, o1, r::Float64, comm; upper_only::Bool = false)
    K1 = o1.K1; m1 = o1.m1; x = o1.x
    nodes = owned_nodes(P); nloc = length(nodes)
    rv = rowvals(K1); nzv = nonzeros(K1)                 # K1 symmetric: column = row
    rowptr = zeros(Int, nloc + 1); rowptr[1] = 1
    cap = nloc * (3 * maximum(diff(K1.colptr)))
    cols = Vector{HYPRE_BigInt}(undef, cap); vals = Vector{Float64}(undef, cap)
    c = _c(r); A = amp(c); um = umean(c)
    b = zeros(nloc); U = zeros(nloc); m = zeros(nloc)
    p = 0
    @inbounds for (l, (i, j, k)) in enumerate(nodes)
        g = gid(P, i, j, k)
        mi, mj, mk = m1[i], m1[j], m1[k]
        if g == 1                                        # pinned unknown
            p += 1; cols[p] = 1; vals[p] = 1.0
        else
            pd = 0; diag = 0.0
            for (dir, a) in ((1, i), (2, j), (3, k)), q in nzrange(K1, a)
                a2 = rv[q]
                s = dir == 1 ? mj * mk : dir == 2 ? mi * mk : mi * mj
                v = s * nzv[q]
                if a2 == a
                    diag += v; continue
                end
                gc = dir == 1 ? gid(P, a2, j, k) : dir == 2 ? gid(P, i, a2, k) : gid(P, i, j, a2)
                gc == 1 && continue                      # pinned column
                upper_only && gc < g && continue
                p += 1; cols[p] = gc; vals[p] = v
            end
            p += 1; cols[p] = g; vals[p] = diag
        end
        rowptr[l+1] = p + 1
        px, py, pz = _p(x[i], c), _p(x[j], c), _p(x[k], c)
        qx, qy, qz = _pss(x[i], c), _pss(x[j], c), _pss(x[k], c)
        m[l] = mi * mj * mk
        U[l] = A * (px * py * pz - um)
        b[l] = -A * (qx * py * pz + px * qy * pz + px * py * qz) * m[l]
    end
    resize!(cols, p); resize!(vals, p)
    # project b onto the range of K (zero sum), then pin
    fmean = MPI.Allreduce(sum(b), +, comm) / MPI.Allreduce(sum(m), +, comm)
    b .-= fmean .* m
    P.ilower == 1 && (b[1] = 0.0)
    return (; rowptr, cols, vals, b, U, m, nloc, nnz = p)
end

function errors(u, sys, comm)
    um = MPI.Allreduce(dot(sys.m, u), +, comm) / MPI.Allreduce(sum(sys.m), +, comm)
    e = u .- um .- sys.U
    linf = MPI.Allreduce(maximum(abs, e; init = 0.0), max, comm)
    num = MPI.Allreduce(dot(sys.m, e .^ 2), +, comm); den = MPI.Allreduce(dot(sys.m, sys.U .^ 2), +, comm)
    return (linf = linf, l2rel = sqrt(num / den))
end

# ---------------------------------------------------------------------------
# hypre: IJ matrix and vectors from the local rows; PCG with BoomerAMG or
# diagonal scaling
# ---------------------------------------------------------------------------
function hypre_system(P::Partition, sys, comm)
    A = HYPREMatrix(comm, P.ilower, P.iupper)
    nrows = HYPRE_Int(sys.nloc)
    ncols = HYPRE_Int.(diff(sys.rowptr))
    rows = collect(HYPRE_BigInt, P.ilower:P.iupper)
    @check HYPRE_IJMatrixSetValues(A, nrows, ncols, rows, sys.cols, sys.vals)
    HYPRE.Internals.assemble_matrix(A)
    b = HYPREVector(comm, sys.b, P.ilower, P.iupper)
    x = HYPREVector(comm, zeros(sys.nloc), P.ilower, P.iupper)
    return A, b, x
end

# BoomerAMG settings for 3D: HMIS coarsening, extended+i interpolation
# truncated to 4 entries per row, l1-scaled symmetric hybrid Gauss-Seidel
# smoothing (symmetric, so the V-cycle is a valid CG preconditioner)
amg_options(theta) = (; CoarsenType = 10, InterpType = 6, PMaxElmts = 4, StrongThreshold = theta,
                        RelaxType = 8, NumSweeps = 1, MaxLevels = 25, PrintLevel = 0)

function hypre_solve!(T, info, sys, P, comm, solver; rtol, theta, itmax = 100_000)
    A, b, x = phase(T, :assembly, comm) do
        hypre_system(P, sys, comm)
    end
    pcg = HYPRE.PCG(comm; Tol = rtol, MaxIter = itmax, TwoNorm = 0, PrintLevel = 0, Logging = 1)
    amg = nothing
    if solver === :boomeramg
        amg = HYPRE.BoomerAMG(; amg_options(theta)...)
        HYPRE.Internals.set_precond_defaults(amg)        # one V-cycle, Tol = 0
        HYPRE.Internals.set_precond(pcg, amg)
    else
        lib = Libdl.dlopen(HYPRE_jll.libHYPRE)
        @check HYPRE_ParCSRPCGSetPrecond(pcg, Libdl.dlsym(lib, :HYPRE_ParCSRDiagScale),
                                         Libdl.dlsym(lib, :HYPRE_ParCSRDiagScaleSetup), C_NULL)
    end
    phase(T, :setup, comm) do
        @check HYPRE_ParCSRPCGSetup(pcg, A, b, x)
    end
    phase(T, :solve, comm) do
        HYPRE_ParCSRPCGSolve(pcg, A, b, x)               # non-zero if not converged: checked below
        HYPRE_ClearAllErrors()
    end
    info[:iters] = HYPRE.GetNumIterations(pcg)
    info[:final_res] = HYPRE.GetFinalRelativeResidualNorm(pcg)
    info[:final_res] > 10 * rtol && MPI.Comm_rank(comm) == 0 &&
        @warn "P3DMPI: $solver stopped at relative residual $(info[:final_res]) after $(info[:iters]) iterations"
    u = zeros(sys.nloc)
    copy!(u, x)
    finalize(pcg); amg === nothing || finalize(amg)
    finalize(A); finalize(b); finalize(x)
    return u
end

# ---------------------------------------------------------------------------
# MUMPS: distributed assembled input (ICNTL(18) = 3), upper triangle (SPD),
# centralised right-hand side and solution on the host
# ---------------------------------------------------------------------------
# name => (ICNTL(28), ICNTL(29) or ICNTL(7)). :parmetis and :ptscotch need a
# MUMPS built with them (MUMPS_jll 5.8: INFOG(1) = -38 / failure)
const MUMPS_ORDERINGS = Dict(
    :parmetis => (2, 2), :ptscotch => (2, 1), :metis => (1, 5), :auto => (0, 7))

function mumps_solve!(T, info, sys, P, comm; ordering = :metis, mem_relax = 50)
    me = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    n = P.Ng^3
    irn, jcn, a = phase(T, :assembly, comm) do
        irn = Vector{MUMPS.MUMPS_INT}(undef, sys.nnz); jcn = Vector{MUMPS.MUMPS_INT}(undef, sys.nnz)
        for l in 1:sys.nloc, q in sys.rowptr[l]:sys.rowptr[l+1]-1
            irn[q] = P.ilower + l - 1; jcn[q] = sys.cols[q]
        end
        irn, jcn, sys.vals
    end
    mumps = MUMPS.Mumps{Float64}(MUMPS.mumps_definite, 1)        # SPD, host works too, COMM_WORLD
    seticntl(i, v) = MUMPS.set_icntl!(mumps, i, v; displaylevel = 0)
    seticntl(1, 6); seticntl(2, 0); seticntl(3, 0); seticntl(4, 1)   # errors only
    seticntl(5, 0); seticntl(18, 3)                                  # assembled, distributed
    seticntl(20, 0); seticntl(21, 0)                                 # dense rhs / solution on host
    seticntl(14, mem_relax)                                          # workspace relaxation, %
    seticntl(35, 0)                                                  # no BLR: exact factorisation
    o28, o = MUMPS_ORDERINGS[ordering]
    seticntl(28, o28)
    o28 == 2 ? seticntl(29, o) : seticntl(7, o)
    mumps.n = n
    mumps.nnz_loc = length(a); mumps.nz_loc = length(a) <= typemax(Int32) ? length(a) : 0
    mumps.irn_loc = pointer(irn); mumps.jcn_loc = pointer(jcn); mumps.a_loc = pointer(a)
    counts = MPI.Allgather(Int32(sys.nloc), comm)
    bglob = me == 0 ? zeros(n) : zeros(0)
    u = zeros(sys.nloc)
    GC.@preserve irn jcn a bglob begin
        check(job) = (mumps.infog[1] < 0 && error("MUMPS job $job failed: INFOG(1) = $(mumps.infog[1]), INFOG(2) = $(mumps.infog[2])"))
        phase(T, :setup, comm) do
            MUMPS.set_job!(mumps, 1); MUMPS.invoke_mumps!(mumps); check(1)     # analysis (ordering)
            # factorisation. INFOG(1) = -9 / -8: the work array estimated by the
            # analysis is too small (happens for small problems on many ranks):
            # enlarge it (ICNTL(14)) and refactorise, as the MUMPS guide says.
            # The retries are part of the setup time.
            relax = mem_relax; info[:mumps_retries] = 0
            while true
                MUMPS.set_job!(mumps, 2); MUMPS.invoke_mumps!(mumps)
                (mumps.infog[1] in (-8, -9) && info[:mumps_retries] < 5) || break
                relax *= 2; info[:mumps_retries] += 1
                seticntl(14, relax)
            end
            check(2)
        end
        phase(T, :solve, comm) do
            MPI.Gatherv!(sys.b, me == 0 ? MPI.VBuffer(bglob, counts) : nothing, comm; root = 0)
            if me == 0
                mumps.rhs = pointer(bglob); mumps.lrhs = n; mumps.nrhs = 1
            end
            MUMPS.set_job!(mumps, 3); MUMPS.invoke_mumps!(mumps); check(3)     # solve
            MPI.Scatterv!(me == 0 ? MPI.VBuffer(bglob, counts) : nothing, u, comm; root = 0)
        end
    end
    # INFOG(29): entries in the factors (negative: millions); INFOG(22): MB, all ranks
    f = mumps.infog[29]; info[:factor_nnz] = f < 0 ? -Int(f) * 10^6 : Int(f)
    info[:mumps_mem_gb] = mumps.infog[22] / 1024
    info[:ordering] = ordering
    MUMPS.finalize!(mumps)
    return u
end

# ---------------------------------------------------------------------------
# Timing: a phase starts and ends at a barrier, so it is the slowest rank's
# time
# ---------------------------------------------------------------------------
function phase(f, T, key, comm)
    MPI.Barrier(comm); t0 = MPI.Wtime()
    v = f()
    MPI.Barrier(comm)
    T[key] = get(T, key, 0.0) + MPI.Wtime() - t0
    return v
end

"""
    run_config(solver, ne, N; r = 0.5, rtol = 1e-12, ordering = :metis, theta = 0.5,
               comm = MPI.COMM_WORLD) -> row

Solve the 3D periodic problem on ne³ elements of order N with `solver`
(:mumps, :boomeramg, :jacobi) on all ranks of `comm`, and return (on every
rank) the row: errors, phase timings (slowest rank), iterations, sizes,
memory.
"""
function run_config(solver::Symbol, ne::Int, N::Int; r::Float64 = R_DEFAULT, rtol = 1e-12,
                    ordering::Symbol = :metis, theta = 0.5, comm = MPI.COMM_WORLD)
    solver in MPI_SOLVERS || error("unknown MPI solver $solver (have $(MPI_SOLVERS))")
    T = Dict{Symbol, Float64}(); info = Dict{Symbol, Any}()
    Ng = ne * N; n = Ng^3
    P, o1 = phase(T, :assembly, comm) do
        Partition(Ng, comm), sem1d(N, ne)
    end
    sys = phase(T, :assembly, comm) do
        local_system(P, o1, r, comm; upper_only = solver === :mumps)
    end
    nnz_glob = MPI.Allreduce(sys.nnz, +, comm)
    u = solver === :mumps ? mumps_solve!(T, info, sys, P, comm; ordering = ordering) :
                            hypre_solve!(T, info, sys, P, comm, solver; rtol = rtol, theta = theta)
    err = errors(u, sys, comm)
    rss = Sys.maxrss() / 2^30
    g(k) = get(T, k, 0.0)
    return (solver = solver, d = 3, r = r, ne = ne, nop = N, Ng = Ng, n = n, solved = n - 1,
            linf = err.linf, l2rel = err.l2rel,
            assembly = g(:assembly), rhs = 0.0, setup = g(:setup), solve = g(:solve),
            total = g(:assembly) + g(:setup) + g(:solve),
            iters = get(info, :iters, 0),
            nnz = solver === :mumps ? 2 * nnz_glob - n : nnz_glob,       # full matrix (MUMPS gets the upper half)
            factor_nnz = get(info, :factor_nnz, 0), skeleton_nnz = 0,
            mumps_retries = get(info, :mumps_retries, 0),
            ordering = solver === :mumps ? ordering : (solver === :boomeramg ? Symbol("hmis_theta", theta) : :none),
            nranks = MPI.Comm_size(comm), grid = join(P.dims, "x"),
            mumps_mem_gb = round(get(info, :mumps_mem_gb, 0.0), digits = 3),
            maxrss_gb = round(MPI.Allreduce(rss, +, comm), digits = 3),
            maxrss_rank_gb = round(MPI.Allreduce(rss, max, comm), digits = 3))
end

end # module

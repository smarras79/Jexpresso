# ═════════════════════════════════════════════════════════════════════════════
#  asm_preconditioner.jl
#
#  Parallel GMRES preconditioners for CG-SEM radiative transfer (5D, MPI).
#
#  ── Required packages (always needed) ────────────────────────────────────
#  using SparseArrays, LinearAlgebra, KLU, ILUZero, IncompleteLU
#  using Krylov, LinearOperators, AMD, MPI
#
#  ── Optional packages (load BEFORE calling solver) ───────────────────────
#
#  :global_paru
#    import ParU_jll          # must come BEFORE using LinearSolve
#    using LinearSolve
#
#  :global_pardiso
#    using Pardiso            # loads as LinearSolve extension automatically
#    using LinearSolve
#    # Requires Intel MKL — load module intel-oneapi-mkl on cluster
#    # Recommended env vars:
#    #   export MKL_NUM_THREADS=<cores_per_node>
#
#  :global_superlu
#    using SuperLUDIST        # ] add SuperLUDIST
#    # Requires system SuperLU_DIST libs:
#    #   SuperLUDIST.set_libraries!(libsuperlu_dist_Int32="/path/libsuperlu_dist.so")
#
#  :boomeramg
#    using HYPRE; HYPRE.Init()
#    ] add HYPRE
#    # Check for system HYPRE: module spider hypre
#    # Key kwargs: amg_coarsen, amg_relax, amg_threshold
#
#  :superlu_dist
#    using SuperLUDIST, SparseBase, CIndices
#    SuperLUDIST.set_libraries!(libsuperlu_dist_Int64="$HOME/local/lib/libsuperlu_dist.so")
#    ] add SuperLUDIST SparseBase MatrixMarket
#    export OMP_NUM_THREADS=<cores_per_rank>
#
#  :mumps_dist
#    using MUMPS              # standalone, no LinearSolve needed
#    # Requires system MPI — set JULIA_MPI_BINARY=system and rebuild MPI.jl
#    # Recommended env vars:
#    #   export OPENBLAS_NUM_THREADS=1   # prevent OpenBLAS/OpenMP conflict
#    #   export OMP_NUM_THREADS=<cores_per_node>
#    # NOTE: Known MPI communicator conflicts on macOS — use on cluster only
#
#  ── Preconditioner options ────────────────────────────────────────────────
#  Local ASM:      :klu :rcmsplu :rcmilu :rcmilu0 :splu :lu :ilu :ilu0 :spilu
#  Global (rank0): :global_klu :global_lu :global_ilu :global_ilu0
#                  :global_paru :global_pardiso :global_superlu
#  Distributed:    :mumps_dist :superlu_dist :boomeramg
#  Other:          :jacobi :inner_gmres :none
# ═════════════════════════════════════════════════════════════════════════════

using SparseArrays
using LinearAlgebra
using KLU
using ILUZero
using IncompleteLU
using Krylov
using LinearOperators
using AMD
using MUMPS
using Pardiso
using Random
import Libdl

# ─────────────────────────────────────────────────────────────────────────────
#  ASMPreconditioner
# ─────────────────────────────────────────────────────────────────────────────

struct ASMPreconditioner
    local_factor  :: Any
    solve_l2g     :: Vector{Int}
    zerorow_l2g   :: Vector{Int}
    zero_l2g      :: Vector{Int}
    x_sub         :: Vector{Float64}
    y_sub         :: Vector{Float64}
    gnpoin        :: Int
end

# ─────────────────────────────────────────────────────────────────────────────
#  complete_preconditioner_matrix
#
#  Completes split interface entries via Allgatherv.
#  Returns npoin × npoin completed matrix.
# ─────────────────────────────────────────────────────────────────────────────

function complete_preconditioner_matrix(
        A        :: SparseMatrixCSC{Float64},
        ip2gip   :: AbstractVector{Int},
        gip2owner:: AbstractVector{Int},
        gnpoin   :: Int,
        npoin_g  :: Int,
        g_ip2gip :: AbstractVector{Int},
        comm     :: MPI.Comm = MPI.COMM_WORLD)

    rank   = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    npoin  = length(ip2gip)

    ip2gip_g = Vector{Int}(undef, npoin_g)
    ip2gip_g[1:npoin] .= ip2gip
    if npoin_g > npoin
        ip2gip_g[npoin+1:npoin_g] .= g_ip2gip
    end

    gip_to_local = Dict{Int,Int}(ip2gip_g[ip] => ip for ip = 1:npoin_g)

    nonowned_local = Set{Int}()
    for ip = 1:npoin
        gip2owner[ip] != rank && push!(nonowned_local, ip)
    end
    for ip = npoin+1:npoin_g
        push!(nonowned_local, ip)
    end

    nonowned_gids  = Int[ip2gip_g[ip] for ip in nonowned_local]
    n_local        = Int32(length(nonowned_gids))
    gid_counts     = MPI.Allgather([n_local], comm)
    all_nonowned_gids = Set(MPI.Allgatherv(nonowned_gids, gid_counts, comm))

    rank == 0 && @info "Interface DOFs (nonowned+ghost): $(length(all_nonowned_gids)) " *
                       "out of $gnpoin total " *
                       "($(round(100*length(all_nonowned_gids)/gnpoin, digits=1))%)"

    rows_sp = rowvals(A); vals_sp = nonzeros(A)
    local_I = Int[]; local_J = Int[]; local_V = Float64[]

    for col_ip = 1:npoin_g
        col_gid         = ip2gip_g[col_ip]
        col_is_nonowned = col_gid in all_nonowned_gids
        for idx in nzrange(A, col_ip)
            row_ip = rows_sp[idx]
            row_ip > npoin && continue
            vals_sp[idx] == 0.0 && continue
            row_gid         = ip2gip_g[row_ip]
            row_is_nonowned = row_gid in all_nonowned_gids
            (col_is_nonowned || row_is_nonowned) || continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end

    n_local_entries = Int32(length(local_V))
    entry_counts    = MPI.Allgather([n_local_entries], comm)
    send_ij = Vector{Int}(undef, 2*n_local_entries)
    for k = 1:n_local_entries
        send_ij[2k-1] = local_I[k]; send_ij[2k] = local_J[k]
    end
    all_ij = MPI.Allgatherv(send_ij, Int32.(entry_counts .* 2), comm)
    all_v  = MPI.Allgatherv(local_V, entry_counts, comm)

    global_entries = Dict{Tuple{Int,Int},Float64}()
    for k = 1:length(all_v)
        key = (all_ij[2k-1], all_ij[2k])
        global_entries[key] = get(global_entries, key, 0.0) + all_v[k]
    end

    owned_gip_to_local = Dict{Int,Int}(ip2gip[ip] => ip for ip = 1:npoin)

    base   = A[1:npoin, 1:npoin]
    I_coo  = Int[];     sizehint!(I_coo,  nnz(base))
    J_coo  = Int[];     sizehint!(J_coo,  nnz(base))
    V_coo  = Float64[]; sizehint!(V_coo,  nnz(base))
    rows_b = rowvals(base); vals_b = nonzeros(base)
    for j = 1:npoin
        for k in nzrange(base, j)
            push!(I_coo, rows_b[k]); push!(J_coo, j); push!(V_coo, vals_b[k])
        end
    end

    completed_pos = Dict{Tuple{Int,Int},Float64}()
    for ((rg, cg), v) in global_entries
        ri = get(owned_gip_to_local, rg, 0); ri == 0 && continue
        ci = get(owned_gip_to_local, cg, 0); ci == 0 && continue
        completed_pos[(ri, ci)] = v
    end

    for k = 1:length(I_coo)
        key = (I_coo[k], J_coo[k])
        if haskey(completed_pos, key)
            V_coo[k] = completed_pos[key]; delete!(completed_pos, key)
        end
    end
    for ((ri, ci), v) in completed_pos
        push!(I_coo, ri); push!(J_coo, ci); push!(V_coo, v)
    end

    return sparse(I_coo, J_coo, V_coo, npoin, npoin)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Helpers: _sparsify, _reorder_for_ilu, ReorderedILU
# ─────────────────────────────────────────────────────────────────────────────

function _sparsify(A::SparseMatrixCSC{Float64}, tau::Float64)
    n = size(A,1); d = abs.(diag(A))
    rows = rowvals(A); vals = nonzeros(A)
    I_k = Int[]; J_k = Int[]; V_k = Float64[]
    for j = 1:n
        for idx in nzrange(A, j)
            r = rows[idx]; v = vals[idx]
            if r == j || abs(v) >= tau * max(d[r], d[j], 1e-30)
                push!(I_k, r); push!(J_k, j); push!(V_k, v)
            end
        end
    end
    return sparse(I_k, J_k, V_k, n, n)
end

function _reorder_for_ilu(A::SparseMatrixCSC{Float64})
    A_sym = abs.(A) + abs.(A')
    perm  = AMD.amd(A_sym)
    iperm = invperm(perm)
    return A[perm, perm], perm, iperm
end

struct ReorderedILU
    factor  :: Any
    perm    :: Vector{Int}
    iperm   :: Vector{Int}
    tmp_in  :: Vector{Float64}
    tmp_out :: Vector{Float64}
end

function LinearAlgebra.ldiv!(y::AbstractVector, F::ReorderedILU, x::AbstractVector)
    @inbounds for i in eachindex(F.perm); F.tmp_in[i] = x[F.perm[i]]; end
    LinearAlgebra.ldiv!(F.tmp_out, F.factor, F.tmp_in)
    @inbounds for i in eachindex(F.perm); y[F.perm[i]] = F.tmp_out[i]; end
    return y
end

# ─────────────────────────────────────────────────────────────────────────────
#  _factorize_matrix
#
#  Central factorization dispatcher used by both local and global preconditioners.
#
#  Solver options:
#    :klu      — KLU (exploits block structure, fastest for physics operators)
#    :splu     — UMFPACK sparse LU
#    :lu       — dense LU (small matrices only)
#    :rcmsplu  — AMD reordering + UMFPACK LU (best fill-reducing for general sparse)
#    :rcmilu   — AMD reordering + ILU(τ) (fast setup, approximate)
#    :rcmilu0  — AMD reordering + ILU(0) (zero fill-in, very fast setup)
#    :ilu      — ILU(τ) without reordering
#    :ilu0     — ILU(0) without reordering
#    :spilu    — ILU(0) on sparsified matrix
# ─────────────────────────────────────────────────────────────────────────────

# ── LinearSolve.jl wrapper for advanced solvers ──────────────────────────────
#  Used for :paru, :pardiso solvers.
#  Returns a wrapper that acts like a factorization via ldiv!.
struct LinearSolveFactor
    cache :: Any   # LinearSolve cache object
    n     :: Int
end

function LinearAlgebra.ldiv!(y::AbstractVector, F::LinearSolveFactor, x::AbstractVector)
    F.cache.b .= x
    LinearSolve.solve!(F.cache)
    copyto!(y, F.cache.u)
    return y
end

function _factorize_matrix(A::SparseMatrixCSC{Float64}, solver::Symbol, ilu_tau::Float64)
    if solver == :klu
        return klu(A)
    elseif solver == :splu
        return lu(A)
    elseif solver == :lu
        return lu(Matrix(A))
    elseif solver == :rcmsplu
        A_rcm, perm, iperm = _reorder_for_ilu(A)
        return ReorderedILU(lu(A_rcm), perm, iperm,
                            zeros(length(perm)), zeros(length(perm)))
    elseif solver == :rcmilu
        A_rcm, perm, iperm = _reorder_for_ilu(A)
        return ReorderedILU(IncompleteLU.ilu(A_rcm, τ=ilu_tau), perm, iperm,
                            zeros(length(perm)), zeros(length(perm)))
    elseif solver == :rcmilu0
        A_rcm, perm, iperm = _reorder_for_ilu(A)
        return ReorderedILU(ILUZero.ilu0(A_rcm), perm, iperm,
                            zeros(length(perm)), zeros(length(perm)))
    elseif solver == :ilu
        return IncompleteLU.ilu(A, τ=ilu_tau)
    elseif solver == :ilu0
        return ILUZero.ilu0(A)
    elseif solver == :spilu
        return ILUZero.ilu0(_sparsify(A, ilu_tau))

    # ── Advanced parallel solvers via LinearSolve.jl ─────────────────────────

    elseif solver == :paru
        # ParU: parallel multifrontal LU (OpenMP, shared memory).
        # Reuses UMFPACK symbolic analysis, parallelises numeric phase.
        # Requires: import ParU_jll; using LinearSolve
        # Install: ] add LinearSolve ParU_jll
        @isdefined(LinearSolve) ||
            error(":paru requires `import ParU_jll; using LinearSolve`")
        isdefined(LinearSolve, :ParUFactorization) ||
            error(":paru: ParUFactorization not found — " *
                  "ensure `import ParU_jll` is called BEFORE `using LinearSolve`")
        n = size(A, 1)
        b_dummy = zeros(Float64, n)
        prob  = LinearSolve.LinearProblem(A, b_dummy)
        cache = LinearSolve.init(prob, LinearSolve.ParUFactorization())
        LinearSolve.solve!(cache)
        return LinearSolveFactor(cache, n)

    elseif solver == :pardiso
        # MKL Pardiso: highly optimised parallel sparse direct solver (OpenMP).
        # Requires Intel MKL to be installed.
        # Install: ] add LinearSolve Pardiso
        # Note: also try :pardiso_iter for iterative Pardiso variant.
        @isdefined(LinearSolve) ||
            error(":pardiso requires `using LinearSolve` and Intel MKL")
        n = size(A, 1)
        b_dummy = zeros(Float64, n)
        prob  = LinearSolve.LinearProblem(A, b_dummy)
        cache = LinearSolve.init(prob, LinearSolve.MKLPardisoFactorize())
        LinearSolve.solve!(cache)
        return LinearSolveFactor(cache, n)


    elseif solver == :superlu_dist
        # SuperLU_DIST: distributed sparse LU with static pivoting.
        # Better memory efficiency than MUMPS for unsymmetric matrices.
        # No separate symbolic analysis phase — avoids MUMPS analysis bottleneck.
        # Install: ] add SuperLUDIST
        # Requires system SuperLU_DIST libraries:
        #   SuperLUDIST.set_libraries!(libsuperlu_dist_Int32="/path/to/lib")
        # Usage: using SuperLUDIST
        @isdefined(SuperLUDIST) ||
            error(":superlu_dist requires `using SuperLUDIST` and system SuperLU_DIST")
        return SuperLUDIST.splu(A)

    else
        error("Unknown solver: $solver. Options:\n" *
              "  Local:   :klu, :splu, :lu, :rcmsplu, :rcmilu, :rcmilu0, :ilu, :ilu0, :spilu\n" *
              "  Parallel:\n" *
              "    :paru         — ParU (OpenMP) — add LinearSolve ParU_jll\n" *
              "    :pardiso      — MKL Pardiso (OpenMP) — add LinearSolve Pardiso\n" *
              "    :superlu_dist — SuperLU_DIST (MPI+OpenMP) — add SuperLUDIST\n" *
              "    via :mumps_dist preconditioner for MUMPS distributed factorization")
    end
end

# Note: MUMPS.jl (JuliaSmoothOptimizers) defines ldiv! for Mumps objects automatically.
# LinearSolveFactor.ldiv! handles ParU/Pardiso via LinearSolve.jl cache.

# ─────────────────────────────────────────────────────────────────────────────
#  build_asm_preconditioner
#
#  Local Additive Schwarz preconditioner.
#  Each rank factorizes its own completed local matrix.
#  Best for: conforming parallel (~27 iters), non-conforming serial (~6 iters).
#  Non-conforming parallel: ~400 iterations (limited by cross-rank coupling).
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
#  analyze_operator_for_preconditioning
#
#  Data-driven diagnostic answering: "would a Jacobi/Gauss-Seidel/ILU-family
#  iterative preconditioner even have a chance of working on this operator?"
#  This replaces guesswork with the textbook criterion: for the (weighted)
#  Jacobi splitting A = D - R, the stationary iteration
#      x^{k+1} = D^{-1}b + (I - D^{-1}A) x^k
#  converges iff ρ(I - D^{-1}A) < 1, where ρ is the spectral radius (largest
#  eigenvalue magnitude). If ρ ≥ 1, NO amount of tuning a Jacobi, Gauss-Seidel,
#  SOR, or (to first order) ILU-family smoother will make it converge — this
#  is a hard mathematical fact about the operator, not a parameter choice.
#  ρ is estimated cheaply via power iteration (O(n_power_iter) matvecs, no
#  factorization), reusing the same global matvec used by GMRES itself.
#
#  Also reports:
#    - diagonal range (min/max/#zero/#negative) — same check done many times
#      earlier in this codebase's development by hand, now automated
#    - row-wise diagonal dominance ratio r_i = |a_ii| / Σ_{j≠i}|a_ij|,
#      computed from LOCAL row entries only (a conservative-in-name-only,
#      cheap proxy — at partition-boundary rows this OVERESTIMATES r_i since
#      it misses ghost-side off-diagonal contributions, so treat "looks
#      dominant" near boundaries with suspicion; "looks non-dominant" is
#      trustworthy everywhere since adding more terms can only lower r_i)
#
#  This is meant to be called ONCE (e.g. at problem setup), not every
#  timestep — it costs ~n_power_iter matvecs (~O(nnz) each) plus a handful
#  of global Allreduces, which is cheap relative to a single MUMPS
#  factorization but not free.
# ─────────────────────────────────────────────────────────────────────────────

function analyze_operator_for_preconditioning(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        npoin_g      :: Int         = length(ip2gip),
        g_ip2gip     :: Vector{Int} = Int[],
        g_gip2ip     :: Dict{Int,Int} = Dict{Int,Int}(),
        comm         :: MPI.Comm    = MPI.COMM_WORLD,
        n_power_iter :: Int         = 30,
        seed         :: Int         = 12345)

    rank  = MPI.Comm_rank(comm)
    npoin = length(ip2gip)

    # ── Diagonal extraction (global, Allreduce-accumulated) ───────────────
    rows_AL = rowvals(A_local); vals_AL = nonzeros(A_local)
    d_global = zeros(Float64, gnpoin)
    n_zero_diag = 0
    n_neg_diag  = 0
    diag_min = Inf
    diag_max = -Inf
    for j = 1:npoin
        for idx in nzrange(A_local, j)
            if rows_AL[idx] == j
                dv = vals_AL[idx]
                d_global[ip2gip[j]] += dv
                break
            end
        end
    end
    MPI.Allreduce!(d_global, +, comm)
    for j = 1:npoin
        dv = d_global[ip2gip[j]]
        dv == 0.0 && (n_zero_diag += 1)
        dv <  0.0 && (n_neg_diag  += 1)
        diag_min = min(diag_min, dv)
        diag_max = max(diag_max, dv)
    end
    n_zero_diag = MPI.Allreduce(n_zero_diag, +, comm)
    n_neg_diag  = MPI.Allreduce(n_neg_diag,  +, comm)
    diag_min    = MPI.Allreduce(diag_min, MPI.MIN, comm)
    diag_max    = MPI.Allreduce(diag_max, MPI.MAX, comm)

    # ── Row-wise diagonal dominance ratio (local, conservative caveat above) ─
    ratios = Float64[]
    for j = 1:npoin
        diagv = 0.0
        offsum = 0.0
        for idx in nzrange(A_local, j)
            row = rows_AL[idx]
            row > npoin && continue
            if row == j
                diagv = vals_AL[idx]
            else
                offsum += abs(vals_AL[idx])
            end
        end
        if offsum > 0
            push!(ratios, abs(diagv) / offsum)
        end
    end
    n_rows_local        = length(ratios)
    n_not_dominant_local = count(<(1.0), ratios)
    n_rows_total         = MPI.Allreduce(n_rows_local, +, comm)
    n_not_dominant_total = MPI.Allreduce(n_not_dominant_local, +, comm)
    frac_not_dominant = n_rows_total > 0 ? n_not_dominant_total / n_rows_total : NaN
    median_ratio_local = isempty(ratios) ? NaN : sort(ratios)[max(1,length(ratios)÷2)]

    # ── Spectral radius of Jacobi iteration matrix M_J = I - D^{-1}A ──────
    # via power iteration on the GLOBAL matvec (same operator GMRES sees).
    A_parallel = create_parallel_linear_operator(
        A_local, ip2gip, gip2owner, npoin, gnpoin,
        npoin_g, g_ip2gip, g_gip2ip, comm)

    d_safe = [abs(d_global[i]) < 1e-30 ? 1.0 : d_global[i] for i = 1:gnpoin]

    rng = MersenneTwister(seed)   # same seed on every rank → same x on every rank
    x = randn(rng, gnpoin)
    xn = sqrt(MPI.Allreduce(sum(abs2, x), +, comm))
    x ./= xn

    rho_est = NaN
    for it = 1:n_power_iter
        Ax = A_parallel * x
        y  = x .- Ax ./ d_safe          # M_J * x
        yn = sqrt(MPI.Allreduce(sum(abs2, y), +, comm))
        rho_est = yn                     # ||x||=1, so ||M_J x|| ≈ ρ after convergence
        x = y ./ yn
    end

    # ── Report ──────────────────────────────────────────────────────────
    if rank == 0
        @info "═══ Operator diagnostic for preconditioning ═══"
        @info "  Diagonal: min=$diag_min max=$diag_max zero_count=$n_zero_diag neg_count=$n_neg_diag (of $gnpoin)"
        @info "  Diagonal dominance: $(round(100*frac_not_dominant,digits=1))% of rows are NOT diagonally dominant " *
              "(|a_ii| < Σ|a_ij|); rank-0 local median ratio ≈ $(round(median_ratio_local,sigdigits=3)) " *
              "(ratio<1 means row is NOT dominant; local estimate, see caveat in source)"
        @info "  Jacobi iteration spectral radius estimate ρ(I - D⁻¹A) ≈ $(round(rho_est,sigdigits=4)) " *
              "(after $n_power_iter power iterations)"
        if rho_est < 1.0
            @info "  → ρ < 1: Jacobi/Gauss-Seidel/SOR/ILU-family stationary preconditioners are " *
                  "THEORETICALLY VIABLE for this operator. Worth pursuing ILU/SOR tuning."
        else
            @info "  → ρ ≥ 1: Jacobi/Gauss-Seidel/SOR/ILU-family stationary preconditioners CANNOT " *
                  "converge on this operator, by construction, regardless of tuning. This is a hard " *
                  "mathematical property of the operator (streaming/advection dominates absorption on " *
                  "the diagonal), not a parameter-search problem. Viable directions: (a) exact/near-exact " *
                  "block solves that exploit the angle-decoupled structure (what :angular_direct already " *
                  "does), or (b) a genuine physics-based transport sweep / DSA scheme (see earlier " *
                  "discussion — non-trivial to implement for CG-SEM, but the textbook remedy for this " *
                  "operator class)."
        end
        flush(stdout)
    end

    return (diag_min=diag_min, diag_max=diag_max, n_zero_diag=n_zero_diag, n_neg_diag=n_neg_diag,
            frac_not_dominant=frac_not_dominant, rho_jacobi=rho_est)
end


function build_asm_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        solver    :: Symbol      = :klu,
        npoin_g   :: Int         = 0,
        g_ip2gip  :: Vector{Int} = Int[],
        ilu_tau   :: Float64     = 0.01,
        comm      :: MPI.Comm    = MPI.COMM_WORLD,
        analyze   :: Bool        = false,   # run analyze_operator_for_preconditioning on entry
        g_gip2ip  :: Dict{Int,Int} = Dict{Int,Int}())

    rank  = MPI.Comm_rank(comm)
    npoin = length(ip2gip)

    if analyze
        analyze_operator_for_preconditioning(
            A_local, ip2gip, gip2owner, gnpoin;
            npoin_g=(npoin_g > 0 ? npoin_g : npoin), g_ip2gip=g_ip2gip,
            g_gip2ip=g_gip2ip, comm=comm)
    end

    has_ghosts = (npoin_g > npoin) && !isempty(g_ip2gip)
    ghost_gips = has_ghosts ? g_ip2gip : Int[]

    A_complete = complete_preconditioner_matrix(
        A_local, ip2gip, gip2owner, gnpoin, npoin_g,
        has_ghosts ? g_ip2gip : Int[], comm)

    row_has_nz = falses(npoin)
    rows_C = rowvals(A_complete); vals_C = nonzeros(A_complete)
    for j = 1:npoin
        for idx in nzrange(A_complete, j)
            r = rows_C[idx]
            if r <= npoin && vals_C[idx] != 0.0
                row_has_nz[r] = true
            end
        end
    end

    solve_local   = findall(row_has_nz)
    zerorow_local = findall(.!row_has_nz)
    solve_l2g     = ip2gip[solve_local]
    zerorow_l2g   = ip2gip[zerorow_local]
    zero_l2g      = ghost_gips

    rank == 0 && @info "ASM DOF classification: $(length(solve_l2g)) solve, " *
                       "$(length(zerorow_l2g)) zero-row, $(length(zero_l2g)) true ghosts"

    A_reduced = A_complete[solve_local, solve_local]

    rank == 0 && @info "ASM A_reduced: $(size(A_reduced,1))×$(size(A_reduced,2)), " *
                       "nnz=$(nnz(A_reduced))"

    t_factor = @elapsed local_factor = _factorize_matrix(A_reduced, solver, ilu_tau)
    rank == 0 && @info "ASM factorization ($solver): $(round(t_factor,digits=3))s"

    x_sub = zeros(Float64, length(solve_local))
    y_sub = zeros(Float64, length(solve_local))

    prec = ASMPreconditioner(local_factor, solve_l2g, zerorow_l2g, zero_l2g,
                             x_sub, y_sub, gnpoin)

    return LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> _apply_asm!(y, x, prec, comm))
end

# ─────────────────────────────────────────────────────────────────────────────
#  _apply_asm!
# ─────────────────────────────────────────────────────────────────────────────

function _apply_asm!(y    :: AbstractVector{Float64},
                     x    :: AbstractVector{Float64},
                     prec :: ASMPreconditioner,
                     comm :: MPI.Comm)
    fill!(y, 0.0)
    @inbounds for (i, gip) in enumerate(prec.solve_l2g)
        prec.x_sub[i] = x[gip]
    end
    LinearAlgebra.ldiv!(prec.y_sub, prec.local_factor, prec.x_sub)
    @inbounds for (i, gip) in enumerate(prec.solve_l2g)
        y[gip] = prec.y_sub[i]
    end
    MPI.Allreduce!(y, +, comm)
    return y
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_global_factorization_preconditioner
#
#  Assembles the complete global matrix on rank 0 and factorizes it once.
#  Guarantees the same preconditioner quality as serial.
#
#  Solver options for the global system:
#    :klu      — KLU (recommended: fastest for structured physics operators)
#    :rcmsplu  — AMD + exact LU (best quality, highest setup cost)
#    :rcmilu   — AMD + ILU(τ) (fast setup, ~10-15 iterations)
#    :rcmilu0  — AMD + ILU(0) (fastest setup, more iterations)
#    :splu     — UMFPACK without reordering
#
#  Cost per apply: two actual_gnpoin-size collectives + one ldiv! on rank 0.
#  Setup: Allgatherv of all local entries + factorization on rank 0.
#
#  Scalability note: rank 0 holds and factorizes the complete global matrix.
#  For very large problems consider distributed solvers (MUMPS, SuperLU_DIST).
# ─────────────────────────────────────────────────────────────────────────────

function build_global_factorization_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        solver    :: Symbol      = :klu,
        ilu_tau   :: Float64     = 0.01,
        npoin_g   :: Int         = length(ip2gip),
        g_ip2gip  :: Vector{Int} = Int[],
        comm      :: MPI.Comm    = MPI.COMM_WORLD)

    rank   = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    npoin  = length(ip2gip)

    # Actual global DOF count from ip2gip (avoids phantom GID gaps in gnpoin)
    local_max_gid = maximum(ip2gip)
    actual_gnpoin = MPI.Allreduce(local_max_gid, MPI.MAX, comm)

    has_ghosts_g = npoin_g > npoin
    ip2gip_g_loc = Vector{Int}(undef, npoin_g)
    ip2gip_g_loc[1:npoin] .= ip2gip
    if has_ghosts_g; ip2gip_g_loc[npoin+1:npoin_g] .= g_ip2gip; end

    # Gather all entries from A_local
    local_I = Int[]; local_J = Int[]; local_V = Float64[]
    rows_sp = rowvals(A_local); vals_sp = nonzeros(A_local)
    for col_ip = 1:npoin_g
        col_gid = ip2gip_g_loc[col_ip]
        col_gid > actual_gnpoin && continue
        for idx in nzrange(A_local, col_ip)
            row_ip  = rows_sp[idx]
            row_gid = ip2gip_g_loc[row_ip]
            row_gid > actual_gnpoin && continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end

    # Allgatherv to all ranks
    n_local   = Int32(length(local_V))
    counts    = MPI.Allgather([n_local], comm)
    send_ij   = Vector{Int}(undef, 2*n_local)
    for k = 1:n_local
        send_ij[2k-1] = local_I[k]; send_ij[2k] = local_J[k]
    end
    all_ij = MPI.Gatherv(send_ij, Int32.(counts .* 2), 0, comm)
    all_v  = MPI.Gatherv(local_V, counts,              0, comm)
    send_ij = Int[]

    # Rank 0 assembles and factorizes
    global_factor = if rank == 0
        n_total = length(all_v)
        GI = [all_ij[2k-1] for k=1:n_total]
        GJ = [all_ij[2k]   for k=1:n_total]
        all_ij = Int[]
        A_global = sparse(GI, GJ, all_v, actual_gnpoin, actual_gnpoin)
        GI = Int[]; GJ = Int[]; all_v = Float64[]
        @info "Global matrix: $(actual_gnpoin)×$(actual_gnpoin), nnz=$(nnz(A_global))"
        t = @elapsed fac = _factorize_matrix(A_global, solver, ilu_tau)
        @info "Global factorization ($solver): $(round(t,digits=3))s"
        fac
    else
        nothing
    end

    x_global = zeros(Float64, actual_gnpoin)
    y_global = zeros(Float64, gnpoin)

    return LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> begin
            fill!(x_global, 0.0)
            @inbounds for ip = 1:npoin
                gip2owner[ip] == rank || continue
                gid = ip2gip[ip]
                gid <= actual_gnpoin || continue
                x_global[gid] = x[gid]
            end
            MPI.Allreduce!(x_global, MPI.SUM, comm)

            fill!(y_global, 0.0)
            y_sub = view(y_global, 1:actual_gnpoin)
            # Check if global_factor is a MUMPS object (MPI-collective ldiv!)
            # If so all ranks must call ldiv! together.
            # For non-MPI solvers (UMFPACK, KLU, ParU etc.) only rank 0 solves.
            is_mumps = @isdefined(MUMPS) && global_factor isa MUMPS.Mumps
            if is_mumps
                # All ranks participate in MUMPS solve
                LinearAlgebra.ldiv!(y_sub, global_factor, x_global)
            else
                if rank == 0
                    LinearAlgebra.ldiv!(y_sub, global_factor, x_global)
                end
                MPI.Bcast!(y_global, 0, comm)
            end
            copyto!(y, y_global)
            return y
        end)
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_mumps_preconditioner
#
#  Fully distributed MUMPS preconditioner — each MPI rank contributes its
#  local matrix directly. MUMPS handles distributed assembly, symbolic
#  analysis, and numeric factorization internally across all ranks.
#  No gather-to-rank-0 step — true distributed-memory parallel factorization.
#
#  Matrix gathered to rank 0 via Allgatherv; MUMPS factorizes across all ranks.
#
#  Requires: using MUMPS (] add MUMPS)
#  MUMPS must be initialised within the same MPI communicator.
#  NOTE: Known MPI communicator conflicts on macOS with OpenMPI.
#        Recommended for cluster use with system MPI only.
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
#  MUMPSHandle — reusable MUMPS factorization with cached symbolic analysis
#
#  The MUMPS analysis phase (symbolic factorization + METIS ordering) is
#  expensive and doesn't scale well with MPI ranks. At high rank counts it
#  dominates total solver time. Since the sparsity pattern only changes when
#  the mesh is adapted, the analysis can be cached and reused across timesteps.
#
#  Usage:
#    # First call (full analysis + factorization):
#    handle = build_mumps_preconditioner(A_local, ...; full_analysis=true)
#
#    # Subsequent calls with same sparsity (numeric refactorization only):
#    mumps_refactorize!(handle, A_local_new, ip2gip, gip2owner, gnpoin;
#                       npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)
#
#    # After mesh adaptation (full analysis needed again):
#    mumps_refactorize!(handle, A_local_new, ...; full_analysis=true)
# ─────────────────────────────────────────────────────────────────────────────

mutable struct MUMPSHandle
    mumps_obj       :: Any              # MUMPS.Mumps object
    precond_op      :: LinearOperator   # LinearOperator wrapping the solve
    actual_gnpoin   :: Int
    gnpoin          :: Int
    x_global        :: Vector{Float64}
    y_global        :: Vector{Float64}
    analysis_done   :: Bool
end

"""
    mumps_finalize!(handle::MUMPSHandle)

Explicitly free MUMPS C-side buffers. Call before discarding a handle or
before building a new one. MUMPS allocates large internal arrays invisible
to Julia's GC — without this, repeated solver calls accumulate memory.

Example:
    # Each timestep where sparsity is unchanged:
    mumps_refactorize!(handle, A_new, ip2gip, gip2owner, gnpoin; ...)

    # When handle is no longer needed:
    mumps_finalize!(handle)
"""
function mumps_finalize!(handle::MUMPSHandle)
    try
        @isdefined(MUMPS) && MUMPS.finalize(handle.mumps_obj)
    catch e
        @warn "MUMPS finalize: $e"
    end
    return nothing
end

"""
    mumps_refactorize!(handle, A_local, ip2gip, gip2owner, gnpoin; kwargs...)

Update a MUMPSHandle with new matrix values, reusing the symbolic analysis
if `full_analysis=false` (default). Use `full_analysis=true` after mesh
adaptation when the sparsity pattern has changed.

Typical usage: call once per timestep with `full_analysis=false` to skip
the expensive METIS ordering and symbolic factorization.
"""
function mumps_refactorize!(
        handle    :: MUMPSHandle,
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        npoin_g   :: Int         = length(ip2gip),
        g_ip2gip  :: Vector{Int} = Int[],
        full_analysis :: Bool    = false,
        comm      :: MPI.Comm    = MPI.COMM_WORLD)

    rank  = MPI.Comm_rank(comm)
    npoin = length(ip2gip)

    has_ghosts_g = npoin_g > npoin
    ip2gip_g_loc = Vector{Int}(undef, npoin_g)
    ip2gip_g_loc[1:npoin] .= ip2gip
    if has_ghosts_g; ip2gip_g_loc[npoin+1:npoin_g] .= g_ip2gip; end

    local_max_gid = maximum(ip2gip)
    actual_gnpoin = MPI.Allreduce(local_max_gid, MPI.MAX, comm)

    # Collect all entries
    local_I = Int[]; local_J = Int[]; local_V = Float64[]
    rows_sp = rowvals(A_local); vals_sp = nonzeros(A_local)
    for col_ip = 1:npoin_g
        col_gid = ip2gip_g_loc[col_ip]
        col_gid > actual_gnpoin && continue
        for idx in nzrange(A_local, col_ip)
            row_ip  = rows_sp[idx]
            row_gid = ip2gip_g_loc[row_ip]
            row_gid > actual_gnpoin && continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end

    n_local_r = Int32(length(local_V))
    counts_r  = MPI.Allgather([n_local_r], comm)
    send_ij_r = Vector{Int}(undef, 2*n_local_r)
    for k = 1:n_local_r
        send_ij_r[2k-1] = local_I[k]; send_ij_r[2k] = local_J[k]
    end
    all_ij_r = MPI.Gatherv(send_ij_r, Int32.(counts_r .* 2), 0, comm)
    all_v_r  = MPI.Gatherv(local_V, counts_r, 0, comm)
    local_I = Int[]; local_J = Int[]; local_V = Float64[]; send_ij_r = Int[]

    obj = handle.mumps_obj
    if rank == 0
        n_tot = length(all_v_r)
        GI_r  = [all_ij_r[2k-1] for k=1:n_tot]
        GJ_r  = [all_ij_r[2k]   for k=1:n_tot]
        all_ij_r = Int[]
        A_new = sparse(GI_r, GJ_r, all_v_r, actual_gnpoin, actual_gnpoin)
        GI_r = Int[]; GJ_r = Int[]; all_v_r = Float64[]
        MUMPS.associate_matrix!(obj, A_new)
    end

    if full_analysis || !handle.analysis_done
        t_analysis = @elapsed begin
            MUMPS.set_job!(obj, 1)
            MUMPS.invoke_mumps!(obj)
            if obj.infog[1] == 0
                mem_est   = obj.infog[16]
                mem_limit = ceil(Int32, mem_est * 1.5)
                MUMPS.set_icntl!(obj, 23, mem_limit)
            end
        end
        rank == 0 && @info "MUMPS analysis: $(round(t_analysis,digits=3))s"
        handle.analysis_done = true
    end

    MUMPS.set_job!(obj, 2)
    t_factor = @elapsed MUMPS.invoke_mumps!(obj)
    infog1 = obj.infog[1]
    if infog1 == -10
        mem_limit2 = ceil(Int32, obj.infog[16] * 3.0)
        rank == 0 && @warn "MUMPS out of memory, retrying with $(mem_limit2)MB/proc"
        MUMPS.set_icntl!(obj, 23, mem_limit2)
        MUMPS.set_job!(obj, 2)
        MUMPS.invoke_mumps!(obj)
        infog1 = obj.infog[1]
    end
    infog1 == 0 || @warn "MUMPS refactorize failed INFO(1)=$infog1"
    rank == 0 && @info "MUMPS $(full_analysis ? "analysis+" : "")factorization: $(round(t_factor,digits=3))s"
    return handle
end


function build_mumps_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        npoin_g   :: Int         = length(ip2gip),
        g_ip2gip  :: Vector{Int} = Int[],
        comm      :: MPI.Comm    = MPI.COMM_WORLD)

    @isdefined(MUMPS) ||
        error(":mumps_dist requires `using MUMPS` (] add MUMPS)")

    rank   = MPI.Comm_rank(comm)
    npoin  = length(ip2gip)

    has_ghosts_g = npoin_g > npoin

    # Actual global DOF count
    local_max_gid = maximum(ip2gip)
    actual_gnpoin = MPI.Allreduce(local_max_gid, MPI.MAX, comm)

    # Send ALL entries from A_local (including ghost rows) in global GID indexing.
    # Ghost rows carry parent-parent coupling from child elements — essential
    # for non-conforming case. Same approach as build_global_factorization_preconditioner.
    ip2gip_g_loc = Vector{Int}(undef, npoin_g)
    ip2gip_g_loc[1:npoin] .= ip2gip
    if has_ghosts_g; ip2gip_g_loc[npoin+1:npoin_g] .= g_ip2gip; end

    local_I = Int[]; local_J = Int[]; local_V = Float64[]
    rows_sp = rowvals(A_local); vals_sp = nonzeros(A_local)
    for col_ip = 1:npoin_g
        col_gid = ip2gip_g_loc[col_ip]
        col_gid > actual_gnpoin && continue
        for idx in nzrange(A_local, col_ip)
            row_ip  = rows_sp[idx]
            row_gid = ip2gip_g_loc[row_ip]
            row_gid > actual_gnpoin && continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end


    # Gather all entries to rank 0 ONLY using Gatherv
    n_local_m = Int32(length(local_V))
    counts_m  = MPI.Allgather([n_local_m], comm)
    send_ij_m = Vector{Int}(undef, 2*n_local_m)
    for k = 1:n_local_m
        send_ij_m[2k-1] = local_I[k]; send_ij_m[2k] = local_J[k]
    end
    all_ij_m = MPI.Gatherv(send_ij_m, Int32.(counts_m .* 2), 0, comm)
    all_v_m  = MPI.Gatherv(local_V, counts_m, 0, comm)
    local_I = Int[]; local_J = Int[]; local_V = Float64[]; send_ij_m = Int[]

    # Initialise MUMPS — all ranks participate
    icntl_init = MUMPS.default_icntl[:]
    mumps_obj = MUMPS.Mumps{Float64}(MUMPS.mumps_unsymmetric,
                                      icntl_init,
                                      MUMPS.default_cntl64)
    MUMPS.set_icntl!(mumps_obj, 1, -1; displaylevel=0)
    MUMPS.set_icntl!(mumps_obj, 2, -1; displaylevel=0)
    MUMPS.set_icntl!(mumps_obj, 3, -1; displaylevel=0)
    MUMPS.set_icntl!(mumps_obj, 4,  1; displaylevel=0)
    MUMPS.set_icntl!(mumps_obj, 27, -2; displaylevel=0)  # single precision
    MUMPS.set_icntl!(mumps_obj, 35,  2; displaylevel=0)  # BLR: 2=compress factors+contrib
    MUMPS.set_icntl!(mumps_obj, 36,  1; displaylevel=0)  # BLR variant
    MUMPS.set_icntl!(mumps_obj, 37,  1; displaylevel=0)  # BLR during analysis
    MUMPS.set_icntl!(mumps_obj, 38, 200; displaylevel=0) # BLR accuracy: higher=more compression
    MUMPS.set_icntl!(mumps_obj, 14, 20; displaylevel=0)  # workspace: reduced from 50→20
    MUMPS.set_icntl!(mumps_obj, 48,  1; displaylevel=0)  # multithreaded factorization
    try MUMPS.set_cntl!(mumps_obj, 7, 1e-2; displaylevel=0) catch end  # BLR drop: 1e-4→1e-2

    if rank == 0
        n_total_m = length(all_v_m)
        GI_m = [all_ij_m[2k-1] for k=1:n_total_m]
        GJ_m = [all_ij_m[2k]   for k=1:n_total_m]
        all_ij_m = Int[]
        A_global_m = sparse(GI_m, GJ_m, all_v_m, actual_gnpoin, actual_gnpoin)
        GI_m = Int[]; GJ_m = Int[]; all_v_m = Float64[]
        @info "MUMPS: global matrix $(actual_gnpoin)×$(actual_gnpoin), nnz=$(nnz(A_global_m))"
        MUMPS.associate_matrix!(mumps_obj, A_global_m)
    end

    # Analysis phase (JOB=1) — sequential METIS ordering + symbolic factorization
    # Parallel analysis (ParMETIS/PT-Scotch) incompatible with centralized input mode
    # The analysis cost is amortized via handle caching (full_analysis=false on reuse)
    t_analysis = @elapsed begin
        MUMPS.set_job!(mumps_obj, 1)
        MUMPS.invoke_mumps!(mumps_obj)
        if mumps_obj.infog[1] == 0
            mem_est   = mumps_obj.infog[16]
            mem_limit = ceil(Int32, mem_est * 1.2)
            rank == 0 && @info "MUMPS analysis: ~$(mem_est)MB/proc → limit $(mem_limit)MB/proc"
            MUMPS.set_icntl!(mumps_obj, 23, mem_limit)
        end
    end
    if rank == 0
        @info "MUMPS analysis: $(round(t_analysis,digits=3))s INFOG(1)=$(mumps_obj.infog[1]) INFOG(16)=$(mumps_obj.infog[16])MB est"
        flush(stdout)
    end

    # Numeric factorization (JOB=2)
    MUMPS.set_job!(mumps_obj, 2)
    t_factor = @elapsed MUMPS.invoke_mumps!(mumps_obj)
    infog1 = mumps_obj.infog[1]
    if infog1 == -10
        mem_limit2 = ceil(Int32, mumps_obj.infog[16] * 3.0)
        rank == 0 && @warn "MUMPS out of memory, retrying with $(mem_limit2)MB/proc"
        MUMPS.set_icntl!(mumps_obj, 23, mem_limit2)
        MUMPS.set_job!(mumps_obj, 2)
        t_factor = @elapsed MUMPS.invoke_mumps!(mumps_obj)
        infog1 = mumps_obj.infog[1]
    end
    if infog1 != 0
        @warn "MUMPS factorization failed INFO(1)=$infog1 INFO(2)=$(mumps_obj.infog[2])"
        # INFO(1)=-3: integer workarray too small → increase ICNTL(14)
        # INFO(1)=-8: integer workarray too small → increase ICNTL(14)
        # INFO(1)=-9: real workarray too small → increase ICNTL(14)
        # INFO(1)=-10: numerically singular matrix
    end
    rank == 0 && @info "MUMPS numeric factorization: $(round(t_factor,digits=3))s INFO(1)=$infog1"


    # Work vectors
    rhs_buf  = zeros(Float64, actual_gnpoin)
    y_global = zeros(Float64, gnpoin)

    op = LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> begin
            fill!(rhs_buf, 0.0)
            @inbounds for ip = 1:npoin
                gip2owner[ip] == rank || continue
                gid = ip2gip[ip]
                gid <= actual_gnpoin || continue
                rhs_buf[gid] = x[gid]
            end
            MPI.Allreduce!(rhs_buf, MPI.SUM, comm)

            # MUMPS solve — collective, all ranks participate
            # Uses cached factorization from mumps_obj (analysis reused
            # across timesteps via mumps_refactorize! when mesh unchanged)
            MUMPS.associate_rhs!(mumps_obj, rhs_buf)
            MUMPS.solve!(mumps_obj)

            fill!(y_global, 0.0)
            if rank == 0
                sol = MUMPS.get_solution(mumps_obj)
                copyto!(view(y_global, 1:actual_gnpoin), sol)
            end
            MPI.Bcast!(y_global, 0, comm)
            copyto!(y, y_global)
            return y
        end)

    handle = MUMPSHandle(mumps_obj, op, actual_gnpoin, gnpoin,
                         rhs_buf, y_global, true)
    # Register Julia finalizer so GC can clean up if mumps_finalize! isn't called
    finalizer(h -> begin
        try; @isdefined(MUMPS) && MUMPS.finalize(h.mumps_obj); catch; end
    end, handle)
    return handle
end


# ─────────────────────────────────────────────────────────────────────────────
#  Jacobi preconditioner
# ─────────────────────────────────────────────────────────────────────────────

function build_jacobi_preconditioner(A_local, ip2gip, gnpoin, comm=MPI.COMM_WORLD)
    npoin = length(ip2gip)
    d_global = zeros(Float64, gnpoin)
    for ip = 1:npoin; d_global[ip2gip[ip]] += A_local[ip,ip]; end
    MPI.Allreduce!(d_global, +, comm)
    d_min = minimum(abs.(d_global[d_global .!= 0.0]))
    for i = 1:gnpoin; abs(d_global[i]) < 1e-14 && (d_global[i] = d_min); end
    d_inv = 1.0 ./ d_global
    return LinearOperator(Float64, gnpoin, gnpoin, true, true,
                          (y, x) -> y .= d_inv .* x)
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_inner_gmres_preconditioner
# ─────────────────────────────────────────────────────────────────────────────

function build_inner_gmres_preconditioner(A_parallel::LinearOperator, gnpoin::Int;
        inner_tol::Float64=1e-2, inner_itmax::Int=10)
    return LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> begin
            y_inner, _ = Krylov.gmres(A_parallel, x;
                atol=inner_tol, rtol=inner_tol, itmax=inner_itmax, verbose=0)
            copyto!(y, y_inner); return y
        end)
end



# ─────────────────────────────────────────────────────────────────────────────
#  SuperLUHandle — reusable SuperLU_DIST factorization
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SuperLUHandle
    F               :: Any
    A_dist          :: Any
    grid            :: Any
    actual_gnpoin   :: Int
    gnpoin          :: Int
    first_row       :: Int
    local_nrows     :: Int
    ip2gip          :: Vector{Int}
    gip2owner       :: Vector{Int}
    x_global        :: Vector{Float64}
    y_global        :: Vector{Float64}
    precond_op      :: LinearOperator
end

function superlu_finalize!(handle::SuperLUHandle)
    handle.F = nothing; handle.A_dist = nothing
    return nothing
end

function build_superlu_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        npoin_g   :: Int         = length(ip2gip),
        g_ip2gip  :: Vector{Int} = Int[],
        comm      :: MPI.Comm    = MPI.COMM_WORLD)

    @isdefined(SuperLUDIST) ||
        error("build_superlu_preconditioner requires `using SuperLUDIST, SparseBase, CIndices`")
    @isdefined(CIndices) ||
        error("build_superlu_preconditioner requires `using CIndices` (dependency of SuperLUDIST)")

    # Load SuperLU_DIST library explicitly — bypasses preference system which
    # may not be set on compute nodes. Try both lib and lib64 paths.
    # RTLD_GLOBAL ensures MPI symbols resolve from system MPI, not Julia's.
    let superlu_lib = get(ENV, "SUPERLU_DIST_LIB",
                          "/home/tissaoui/local/lib64/libsuperlu_dist.so")
        try
            Libdl.dlopen(superlu_lib, Libdl.RTLD_GLOBAL)
            rank == 0 && @info "SuperLU_DIST: loaded $superlu_lib"
        catch e
            error("SuperLU_DIST: failed to load $superlu_lib: $e\n" *
                  "Set SUPERLU_DIST_LIB env var to the correct path.")
        end
    end

    rank   = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    npoin  = length(ip2gip)
    has_ghosts_g = npoin_g > npoin

    local_max_gid = maximum(ip2gip)
    actual_gnpoin = MPI.Allreduce(local_max_gid, MPI.MAX, comm)

    # Build ip2gip_g
    ip2gip_g_loc = Vector{Int}(undef, npoin_g)
    ip2gip_g_loc[1:npoin] .= ip2gip
    if has_ghosts_g; ip2gip_g_loc[npoin+1:npoin_g] .= g_ip2gip; end

    # Gather all entries
    local_I = Int[]; local_J = Int[]; local_V = Float64[]
    rows_sp = rowvals(A_local); vals_sp = nonzeros(A_local)
    for col_ip = 1:npoin_g
        col_gid = ip2gip_g_loc[col_ip]
        col_gid > actual_gnpoin && continue
        for idx in nzrange(A_local, col_ip)
            row_ip  = rows_sp[idx]
            row_gid = ip2gip_g_loc[row_ip]
            row_gid > actual_gnpoin && continue
            push!(local_I, row_gid)
            push!(local_J, col_gid)
            push!(local_V, vals_sp[idx])
        end
    end

    n_loc = Int32(length(local_V))
    cnts  = MPI.Allgather([n_loc], comm)
    sij   = Vector{Int}(undef, 2*n_loc)
    for k = 1:n_loc; sij[2k-1] = local_I[k]; sij[2k] = local_J[k]; end
    aij   = MPI.Allgatherv(sij, Int32.(cnts .* 2), comm)
    av    = MPI.Allgatherv(local_V, cnts, comm)

    # Build process grid — try square, fall back to 1xN
    nprow, npcol = 1, nprocs
    for p = isqrt(nprocs):-1:1
        if nprocs % p == 0; nprow = p; npcol = nprocs ÷ p; break; end
    end
    rank == 0 && @info "SuperLU_DIST: process grid $(nprow)x$(npcol)"
    grid = SuperLUDIST.Grid{Int64}(nprow, npcol, comm)

    # Scatter to DistributedSuperMatrix
    A_dist = SuperLUDIST.DistributedSuperMatrix{Float64,Int64}(grid)
    csr_g  = nothing
    chunks = nothing
    if rank == 0
        n_tot = length(av)
        GI = [aij[2k-1] for k=1:n_tot]
        GJ = [aij[2k]   for k=1:n_tot]
        A_global = sparse(GI, GJ, av, actual_gnpoin, actual_gnpoin)
        @info "SuperLU_DIST: global matrix $(actual_gnpoin)x$(actual_gnpoin), nnz=$(nnz(A_global))"
        # Build CSR store — SuperLUDIST uses CIndex internally for indices
        # CSC of A^T gives colptr=row ptrs, rowval=col indices in CSR format
        A_csr_sp = SparseMatrixCSC(A_global')
        csr_g = SparseBase.SinglyCompressedStore(
            SparseBase.RowMajor(),
            CIndices.CIndex{Int64}.(A_csr_sp.colptr),
            CIndices.CIndex{Int64}.(A_csr_sp.rowval),
            A_csr_sp.nzval,
            (size(A_global,1), size(A_global,2)))
        chunks = SuperLUDIST.Communication.distribute_evenly(actual_gnpoin, nprocs)
    end

    t_f = @elapsed begin
        SuperLUDIST.Communication.scatterstore!(A_dist, csr_g, chunks; root=0)
        # Disable ILU (ILU_DropRule=0) — default options trigger ILU mode
        # which is not implemented in SuperLU_DIST. Use direct LU only.
        opts_slu = SuperLUDIST.Options()
        opts_slu.ILU_DropRule = Int32(0)
        opts_slu.PrintStat    = SuperLUBase.Common.NO
        F = lu!(A_dist; options=opts_slu)
    end
    rank == 0 && @info "SuperLU_DIST factorization: $(round(t_f,digits=3))s"

    # Set OMP threads to avoid oversubscription
    omp_t = parse(Int, get(ENV, "OMP_NUM_THREADS", "1"))
    SuperLUDIST.superlu_set_num_threads(Int64, omp_t)

    first_row   = A_dist.first_row
    local_nrows = SuperLUDIST.Communication.localsize(A_dist, 1)

    x_global = zeros(Float64, actual_gnpoin)
    y_global = zeros(Float64, gnpoin)

    op = LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> begin
            # Gather owned DOF values
            fill!(x_global, 0.0)
            @inbounds for ip = 1:npoin
                gip2owner[ip] == rank || continue
                gid = ip2gip[ip]
                gid <= actual_gnpoin || continue
                x_global[gid] = x[gid]
            end
            MPI.Allreduce!(x_global, MPI.SUM, comm)

            # Extract local rows, solve, scatter back
            b_local = reshape(x_global[first_row:first_row+local_nrows-1], :, 1)
            x_local = F \ b_local

            fill!(y_global, 0.0)
            y_global[first_row:first_row+local_nrows-1] .= vec(x_local)
            MPI.Allreduce!(y_global, MPI.SUM, comm)
            copyto!(y, y_global)
            return y
        end)

    handle = SuperLUHandle(F, A_dist, grid, actual_gnpoin, gnpoin,
                           first_row, local_nrows, ip2gip, gip2owner,
                           x_global, y_global, op)
    finalizer(handle) do h
        try
            # SuperLU_DIST manages its own memory — nulling references
            # allows GC to collect Julia-side objects; C-side cleanup
            # is handled by SuperLUDIST.jl's own finalizers on F and A_dist
            h.F = nothing
            h.A_dist = nothing
        catch
        end
    end
    return handle
end


# ─────────────────────────────────────────────────────────────────────────────
#  build_spatial_amg_preconditioner
#
#  Physics-based two-level preconditioner for 5D radiation transport:
#    1. Project 5D operator onto 3D spatial subspace via P (5D→3D)
#    2. Build BoomerAMG on projected 3D spatial operator A_spa = P*A*P^T
#    3. Apply: e_5d = P^T * (A_spa \ (P * r_5d))
#
#  Works for both conforming and non-conforming angular meshes.
#
#  Conforming case:  ip_5d = (ip_3d-1)*extra_npoin + ip_ext
#                   → ip_3d = (ip_5d-1) ÷ extra_npoin + 1
#
#  Non-conforming:  ip_3d extracted from connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
#                   and mesh.connijk[iel,i,j,k]
#
#  The projected A_spa captures the spatial elliptic structure (diffusion-like
#  after angular averaging) that AMG handles well, while GMRES handles the
#  angular streaming.
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_spatial_projection(npoin_5d, npoin_3d, extra_npoin; connijk_spa=nothing, mesh=nothing)

Build the sparse projection matrix P of size (npoin_3d × npoin_5d) that maps
5D DOFs to 3D spatial DOFs. P[ip_3d, ip_5d] = 1 if ip_5d belongs to ip_3d.

Conforming case (connijk_spa===nothing): uses formula ip_3d = (ip_5d-1)÷extra_npoin+1
Non-conforming case: iterates over mesh elements using connijk and connijk_spa.
"""
function build_spatial_projection(
        npoin_5d    :: Int,
        npoin_3d    :: Int,
        extra_npoin :: Int;
        connijk_spa  = nothing,
        mesh         = nothing,
        extra_mesh   = nothing,
        extra_meshes = nothing)

    I_P = Int[]
    J_P = Int[]
    V_P = Float64[]

    if connijk_spa === nothing
        # ── Conforming case ───────────────────────────────────────────────
        # A_local is already DSS'd so node sharing is correctly encoded.
        # P just needs the angular quadrature weight for each angular DOF.
        sizehint!(I_P, npoin_5d)
        sizehint!(J_P, npoin_5d)
        sizehint!(V_P, npoin_5d)

        if extra_mesh !== nothing
            w_ext = zeros(Float64, extra_npoin)
            for e_ext = 1:extra_mesh.extra_nelem
                nop = extra_mesh.extra_nop[e_ext]
                for iϕ = 1:nop+1, iθ = 1:nop+1
                    ip_ext = extra_mesh.extra_connijk[e_ext, iθ, iϕ]
                    Je     = extra_mesh.extra_metrics.Je[e_ext, iθ, iϕ]
                    w_ext[ip_ext] += Je * extra_mesh.ωθ[iθ] * extra_mesh.ωθ[iϕ]
                end
            end
            for ip_5d = 1:npoin_5d
                ip_3d  = (ip_5d - 1) ÷ extra_npoin + 1
                ip_ext = (ip_5d - 1) % extra_npoin + 1
                push!(I_P, ip_3d)
                push!(J_P, ip_5d)
                push!(V_P, w_ext[ip_ext])
            end
        else
            # No quadrature info — uniform weights
            for ip_5d = 1:npoin_5d
                ip_3d = (ip_5d - 1) ÷ extra_npoin + 1
                push!(I_P, ip_3d)
                push!(J_P, ip_5d)
                push!(V_P, 1.0)
            end
        end

    else
        # ── Non-conforming case ───────────────────────────────────────────
        # A_local is DSS'd — no node_div needed in P.
        w_acc = Dict{Tuple{Int,Int}, Float64}()
        nelem = mesh.nelem
        ngl   = mesh.ngl
        for iel = 1:nelem
            em      = (extra_meshes !== nothing) ? extra_meshes[iel] : extra_mesh
            n_e_ext = size(connijk_spa[iel], 4)
            nop_i   = size(connijk_spa[iel], 5)
            nop_j   = size(connijk_spa[iel], 6)
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                ip_3d = mesh.connijk[iel, i, j, k]
                for e_ext = 1:n_e_ext
                    for jθ = 1:nop_j, iθ = 1:nop_i
                        ip_5d = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        Je    = em.extra_metrics.Je[e_ext, iθ, jθ]
                        w     = Je * em.ωθ[iθ] * em.ωθ[jθ]
                        key   = (ip_3d, ip_5d)
                        w_acc[key] = get(w_acc, key, 0.0) + w
                    end
                end
            end
        end
        for ((ip_3d, ip_5d), w) in w_acc
            push!(I_P, ip_3d)
            push!(J_P, ip_5d)
            push!(V_P, w)
        end
    end

    return sparse(I_P, J_P, V_P, npoin_3d, npoin_5d)
end


"""
    build_spatial_amg_preconditioner(A_local, ip2gip_5d, ip2gip_3d, gip2owner_3d,
                                      gnpoin_3d, gnpoin_5d, extra_npoin; kwargs...)

Build the spatial AMG preconditioner for the 5D radiation transport operator.
Returns a LinearOperator that applies P^T * (A_spa_inv * (P * r)).
"""
function build_spatial_amg_preconditioner(
        A_local      :: AbstractSparseMatrix{Float64},
        ip2gip_5d    :: AbstractVector{Int},   # 5D global IPs
        ip2gip_3d   :: AbstractVector{Int},   # 3D spatial global IPs (mesh.ip2gip)
        gip2owner_3d :: AbstractVector{Int},
        gnpoin_3d   :: Int,
        gnpoin_5d   :: Int,
        extra_npoin :: Int;                   # angular DOFs per spatial node (conforming)
        connijk_spa  = nothing,               # non-conforming angular connectivity
        mesh         = nothing,               # mesh object (needed for non-conforming)
        extra_mesh   = nothing,               # conforming: single extra_mesh object
        extra_meshes = nothing,               # non-conforming: per-element extra_mesh[iel]
        npoin_g_5d  :: Int         = length(ip2gip_5d),
        g_ip2gip_5d :: Vector{Int} = Int[],
        npoin_g_3d  :: Int         = length(ip2gip_3d),
        g_ip2gip_3d :: Vector{Int} = Int[],
        gip2owner_5d :: AbstractVector{Int}   = Int[],
        comm        :: MPI.Comm    = MPI.COMM_WORLD,
        # BoomerAMG kwargs
        CoarsenType     :: Int     = 6,
        RelaxType       :: Int     = 6,
        MaxLevels       :: Int     = 25,
        NumSweeps       :: Int     = 1,
        StrongThreshold :: Float64 = 0.25,
        PrintLevel      :: Int     = 0,
        local_weight    :: Float64 = 1.0)   # weight for additive Jacobi smoother (0=disable)

    rank  = MPI.Comm_rank(comm)
    npoin_5d = length(ip2gip_5d)
    npoin_3d = length(ip2gip_3d)

    # ── Step 1: Extract diagonal for Jacobi smoother ─────────────────────
    # Simple diagonal (Jacobi) scaling as additive correction.
    # D^{-1} is applied additively: y += local_weight * D^{-1} * x
    # This handles the badly-scaled DOFs directly without amplifying A.
    rank == 0 && @info "SpatialAMG: extracting diagonal..."
    flush(stdout)
    rows_AL = rowvals(A_local); vals_AL = nonzeros(A_local)
    d = ones(Float64, npoin_5d)
    for j = 1:npoin_5d
        for idx in nzrange(A_local, j)
            if rows_AL[idx] == j
                d[j] = abs(vals_AL[idx]) < 1e-30 ? 1.0 : vals_AL[idx]
                break
            end
        end
    end

    # ── Step 2: Build projection P (for restriction/prolongation) ────────
    rank == 0 && @info "SpatialAMG: building projection P ($npoin_3d×$npoin_5d)..."
    flush(stdout)
    P = build_spatial_projection(npoin_5d, npoin_3d, extra_npoin;
                                  connijk_spa=connijk_spa, mesh=mesh,
                                  extra_mesh=extra_mesh, extra_meshes=extra_meshes)

    # ── Step 3: Extract per-angular-direction spatial submatrices ─────────
    # For conforming case: direction ip_ext has spatial submatrix
    #   A_spa[ip_ext][ip_3d, jp_3d] = A_local[(ip_3d-1)*nang+ip_ext, (jp_3d-1)*nang+ip_ext]
    # Single pass over A_local: route each entry to its angular direction.
    rank == 0 && @info "SpatialAMG: extracting $extra_npoin per-direction spatial submatrices..."
    flush(stdout)

    # COO buffers for each direction
    I_dirs = [Int[]     for _ = 1:extra_npoin]
    J_dirs = [Int[]     for _ = 1:extra_npoin]
    V_dirs = [Float64[] for _ = 1:extra_npoin]

    rows_AL = rowvals(A_local); vals_AL = nonzeros(A_local)
    for j = 1:npoin_5d
        jp_ext = (j - 1) % extra_npoin + 1   # angular direction of column j
        jp_3d  = (j - 1) ÷ extra_npoin + 1   # spatial node of column j
        for idx in nzrange(A_local, j)
            row = rows_AL[idx]
            row > npoin_5d && continue        # skip ghost rows
            ip_ext = (row - 1) % extra_npoin + 1
            ip_ext == jp_ext || continue      # only same-direction entries
            ip_3d = (row - 1) ÷ extra_npoin + 1
            push!(I_dirs[ip_ext], ip_3d)
            push!(J_dirs[ip_ext], jp_3d)
            push!(V_dirs[ip_ext], vals_AL[idx])
        end
    end

    A_dirs = [sparse(I_dirs[k], J_dirs[k], V_dirs[k], npoin_3d, npoin_3d)
              for k = 1:extra_npoin]
    I_dirs = nothing; J_dirs = nothing; V_dirs = nothing

    diag_vals_all = [extrema(diag(A_dirs[k])) for k = 1:extra_npoin]
    rank == 0 && @info "SpatialAMG: per-direction diag range: " *
        "min=$(minimum(x[1] for x in diag_vals_all)) " *
        "max=$(maximum(x[2] for x in diag_vals_all))"
    flush(stdout)

    # ── Verify block-diagonal structure ──────────────────────────────────
    # Check that A_local has no entries coupling different angular directions.
    # If truly block-diagonal in angle, off-angular-block entries should be zero.
    if rank == 0
        n_cross = 0
        rows_AL2 = rowvals(A_local); vals_AL2 = nonzeros(A_local)
        for j = 1:npoin_5d
            jp_ext = (j - 1) % extra_npoin + 1
            for idx in nzrange(A_local, j)
                row = rows_AL2[idx]
                row > npoin_5d && continue
                ip_ext = (row - 1) % extra_npoin + 1
                if ip_ext != jp_ext && abs(vals_AL2[idx]) > 1e-15
                    n_cross += 1
                end
            end
        end
        @info "Block-diagonal check: $n_cross cross-angular entries (should be 0 for exact decoupling)"
        flush(stdout)
    end

    # ── Step 4: Factorize each per-direction spatial matrix ───────────────
    # Use sparse LU (KLU) — handles indefinite matrices unlike AMG.
    # Each matrix is npoin_3d × npoin_3d ≈ 68921 — small enough for direct solve.
    rank == 0 && @info "SpatialAMG: factorizing $extra_npoin per-direction matrices (KLU)..."
    flush(stdout)
    t_factor = @elapsed begin
        dir_factors = Vector{Any}(undef, extra_npoin)
        for k = 1:extra_npoin
            dir_factors[k] = lu(A_dirs[k])
            if rank == 0 && k % 10 == 0
                @info "SpatialAMG: factorized $k/$extra_npoin"
                flush(stdout)
            end
        end
    end
    A_dirs = nothing
    rank == 0 && @info "SpatialAMG: all factorizations ready ($(round(t_factor,digits=2))s)"
    flush(stdout)

    # ── Step 5: Build LinearOperator for apply ────────────────────────────
    # Distribute angular directions across MPI ranks.
    # Each rank owns ceil(extra_npoin/nprocs) directions and solves them locally.
    # Communication: 2 Allreduces per apply (broadcast residual, gather solution).
    nprocs = MPI.Comm_size(comm)

    # Assign directions to this rank: rank r owns directions r+1, r+1+nprocs, ...
    my_dirs = collect(rank+1:nprocs:extra_npoin)  # 0-based rank, 1-based directions

    r_g5d = zeros(Float64, gnpoin_5d)   # global 5D residual (broadcast)
    e_g5d = zeros(Float64, gnpoin_5d)   # global 5D solution (accumulate)
    r_loc = zeros(Float64, npoin_3d)    # local 3D residual for one direction
    e_loc = zeros(Float64, npoin_3d)    # local 3D solution for one direction

    op = LinearOperator(Float64, gnpoin_5d, gnpoin_5d, false, false,
        (y, x) -> begin
            # Step 1: broadcast full 5D residual to all ranks
            copyto!(r_g5d, x)
            MPI.Allreduce!(r_g5d, +, comm)

            # Step 2: each rank solves its assigned directions
            fill!(e_g5d, 0.0)
            for ip_ext in my_dirs
                # Extract local 3D residual for direction ip_ext
                @inbounds for ip_3d = 1:npoin_3d
                    ip_5d = (ip_3d - 1) * extra_npoin + ip_ext
                    ip_5d > npoin_5d && (r_loc[ip_3d] = 0.0; continue)
                    r_loc[ip_3d] = r_g5d[ip2gip_5d[ip_5d]]
                end

                # Direct solve for this direction (local, no communication)
                ldiv!(e_loc, dir_factors[ip_ext], r_loc)

                # Accumulate solution into global 5D vector
                @inbounds for ip_3d = 1:npoin_3d
                    ip_5d = (ip_3d - 1) * extra_npoin + ip_ext
                    ip_5d > npoin_5d && continue
                    e_g5d[ip2gip_5d[ip_5d]] = e_loc[ip_3d]
                end
            end

            # Step 3: gather solutions from all ranks
            MPI.Allreduce!(e_g5d, +, comm)
            copyto!(y, e_g5d)
            y
        end)

    return op
end


# ─────────────────────────────────────────────────────────────────────────────
#  solve_angular_direct
#
#  Direct block-diagonal solver for radiation transport with inexact integration.
#  With inexact integration, A is exactly block-diagonal in angle — each angular
#  direction ip_ext has an independent 3D spatial system. This function solves
#  all extra_npoin systems directly using KLU, bypassing GMRES entirely.
#
#  Parallelism: angular directions distributed across MPI ranks.
#  Each rank factorizes and solves its assigned directions independently.
#  Communication: 2 Allreduces (gather RHS, scatter solution).
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
#  AngularDirectCache
#
#  Caches symbolic factorizations per sparsity pattern across calls.
#  On first call: analysis + numeric factorization for each pattern group.
#  On subsequent calls: numeric refactorization only (reuse symbolic).
#  Reduces setup time from O(n_pats × analysis + n_dirs × refactor)
#  to O(n_dirs × refactor) on reuse.
# ─────────────────────────────────────────────────────────────────────────────

mutable struct AngularDirectCache
    # Persistent MUMPS instances keyed by SPARSITY PATTERN hash — not by
    # direction index. Empirically there are only ~3-6 distinct sparsity
    # patterns across all extra_npoin directions (structure appears to
    # depend on the streaming-derivative basis-function structure, not on
    # κ/σ values), so this bounds the number of simultaneously-alive MUMPS
    # instances to a handful, unlike a per-direction cache (up to
    # extra_npoin instances), which was observed to crash Intel MPI's
    # shared-memory collective module once enough concurrent MUMPS
    # instances were held alive at once.
    #
    # Analysis (JOB=1) is done once per unique pattern hash and never
    # repeated for that hash again, in this call or any later one — the
    # sparsity STRUCTURE is assumed not to depend on κ/σ or on which
    # specific direction produced it, only on which streaming-derivative
    # terms are structurally present. Refactorization (JOB=2, new numeric
    # values) is still done every single time a pattern is encountered,
    # since VALUES differ direction-to-direction and call-to-call even when
    # the STRUCTURE repeats.
    #
    # The caller is responsible for calling free_angular_direct_cache!
    # whenever the underlying MESH changes (which could in principle change
    # which sparsity patterns occur) — there is no automatic invalidation,
    # deliberately, matching the caller-controlled design used throughout
    # this file (every attempt at implicit same-matrix/same-mesh detection
    # during development caused a correctness or deadlock bug).
    mumps_by_pattern :: Dict{UInt64, MUMPS.Mumps{Float64}}
end

AngularDirectCache() = AngularDirectCache(Dict{UInt64, MUMPS.Mumps{Float64}}())

function free_angular_direct_cache!(cache::AngularDirectCache)
    for (_, m) in cache.mumps_by_pattern
        try MUMPS.finalize!(m) catch end
    end
    empty!(cache.mumps_by_pattern)
end

function solve_angular_direct(
        A_local     :: AbstractSparseMatrix{Float64},
        b           :: AbstractVector{Float64},
        ip2gip_5d   :: AbstractVector{Int},
        ip2gip_3d   :: AbstractVector{Int},
        gnpoin_5d   :: Int,
        gnpoin_3d   :: Int,
        extra_npoin :: Int;
        comm        :: MPI.Comm = MPI.COMM_WORLD,
        cache       :: Union{AngularDirectCache, Nothing} = nothing,
        gather_batch_size :: Int = 20)

    rank   = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    npoin_5d = length(ip2gip_5d)
    npoin_3d = length(ip2gip_3d)
    # Global 5D IDs use dense_gip = rank of spatial GID in sorted unique list,
    # not the original sparse GID. Replicate the same mapping here.
    local_gips     = collect(ip2gip_3d)
    n_loc          = Int32(length(local_gips))
    counts_g       = MPI.Allgather([n_loc], comm)
    all_gips_raw   = MPI.Allgatherv(local_gips, counts_g, comm)
    all_gips_sorted = sort(unique(all_gips_raw))
    gid_to_dense   = Dict{Int,Int}(gid => i for (i,gid) in enumerate(all_gips_sorted))
    gnpoin_3d_dense = length(all_gips_sorted)

    # True global 5D dense size = gnpoin_3d_dense × extra_npoin
    # (overrides gnpoin_5d argument which may be incorrect)
    gnpoin_5d = gnpoin_3d_dense * extra_npoin

    rank == 0 && @info "AngularDirect: gnpoin_3d=$gnpoin_3d dense=$gnpoin_3d_dense extra_npoin=$extra_npoin gnpoin_5d=$gnpoin_5d"
    flush(stdout)

    # ── Step 2: Extract local COO per OWNED direction only ──────────────
    # No Allgatherv — each rank only extracts its owned directions.
    # For MUMPS solve: all ranks contribute their local entries to rank 0
    # direction-by-direction using Reduce, not bulk Allgatherv.
    rank == 0 && @info "AngularDirect: extracting $extra_npoin per-direction matrices..."
    flush(stdout)

    nprocs      = MPI.Comm_size(comm)
    my_dirs     = collect(rank+1:nprocs:extra_npoin)
    my_dirs_set = Set(my_dirs)

    rows_AL = rowvals(A_local); vals_AL = nonzeros(A_local)

    # Extract local COO for ALL directions (needed for Reduce)
    # Store as sparse matrices directly using local entries only
    # These are INCOMPLETE (missing other ranks' entries)
    A_local_dirs = Dict{Int, SparseMatrixCSC{Float64,Int}}()
    I_loc = [Int[]     for _ = 1:extra_npoin]
    J_loc = [Int[]     for _ = 1:extra_npoin]
    V_loc = [Float64[] for _ = 1:extra_npoin]
    for j = 1:npoin_5d
        jp_ext   = (j - 1) % extra_npoin + 1
        jp_3d    = (j - 1) ÷ extra_npoin + 1
        jp_dense = gid_to_dense[ip2gip_3d[jp_3d]]
        for idx in nzrange(A_local, j)
            row = rows_AL[idx]
            row > npoin_5d && continue
            ip_ext = (row - 1) % extra_npoin + 1
            ip_ext == jp_ext || continue
            ip_3d    = (row - 1) ÷ extra_npoin + 1
            ip_dense = gid_to_dense[ip2gip_3d[ip_3d]]
            push!(I_loc[ip_ext], ip_dense)
            push!(J_loc[ip_ext], jp_dense)
            push!(V_loc[ip_ext], vals_AL[idx])
        end
    end

    # Build incomplete local sparse matrices — one per direction
    # Complete matrices assembled via Allreduce on dense value arrays
    for k = 1:extra_npoin
        A_local_dirs[k] = sparse(I_loc[k], J_loc[k], V_loc[k],
                                  gnpoin_3d_dense, gnpoin_3d_dense)
    end
    I_loc = nothing; J_loc = nothing; V_loc = nothing

    if rank == 0
        d1 = diag(A_local_dirs[my_dirs[1]])
        @info "AngularDirect: local dir $(my_dirs[1]): nnz=$(nnz(A_local_dirs[my_dirs[1]])) zero_diag=$(sum(d1.==0))"
        flush(stdout)
    end


    t_total = @elapsed begin

    # ── Step 5: Gather full RHS ───────────────────────────────────────────
    b_global = zeros(Float64, gnpoin_5d)
    for ip = 1:npoin_5d
        b_global[ip2gip_5d[ip]] += b[ip]
    end
    MPI.Allreduce!(b_global, +, comm)

    # ── Step 6: Solve all directions using MUMPS ───────────────────────────
    # All ranks participate in every MUMPS solve — no deadlock.
    # NOTE: pattern-based grouping/reuse was removed. It was computed from
    # each rank's LOCAL partial matrix (A_local_dirs[k]), which differs
    # across ranks (different ranks own different spatial nodes), causing
    # different ranks to group/iterate directions in different orders —
    # a collective-call desync (MUMPS/MPI hang) since invoke_mumps! requires
    # identical JOB sequencing on every rank. Iterating k=1:extra_npoin in a
    # single fixed order, identical on every rank, with a full analysis+
    # factorization per direction, guarantees no rank divergence. Slower,
    # but correct. Pattern reuse can be reintroduced later, but must be
    # decided from data that is identical on every rank (e.g. computed once
    # from the gathered/complete matrix on rank 0 and broadcast).
    my_dirs_set = Set(my_dirs)

    x_global = zeros(Float64, gnpoin_5d)
    b_3d     = zeros(Float64, gnpoin_3d_dense)
    x_3d     = zeros(Float64, gnpoin_3d_dense)

    # ── Helper: gather complete matrices for a BATCH of directions at once ──
    # Packs [k, row, col, val] quadruples for every direction in the batch
    # into ONE send buffer per rank, then ONE Allgather(counts)+Gatherv(data)
    # round trip for the whole batch — instead of one round trip PER
    # direction. Each MPI collective carries fixed per-call overhead
    # regardless of payload size, so batching cuts that fixed cost by
    # roughly the batch size, at the cost of a larger (but still bounded,
    # not all-extra_npoin-at-once) transient buffer on rank 0.
    function gather_dir_matrix_batch(krange)
        send_b = Float64[]
        for k in krange
            A_k_loc = A_local_dirs[k]
            rows_k  = rowvals(A_k_loc); vals_k = nonzeros(A_k_loc)
            for col = 1:gnpoin_3d_dense
                for idx in nzrange(A_k_loc, col)
                    push!(send_b, Float64(k), Float64(rows_k[idx]), Float64(col), vals_k[idx])
                end
            end
        end
        nb     = Int32(length(send_b))
        nc_b   = MPI.Allgather([nb], comm)
        recv_b = MPI.Gatherv(send_b, nc_b, 0, comm)

        A_batch = Dict{Int, SparseMatrixCSC{Float64,Int}}()
        if rank == 0
            I_b = Dict{Int, Vector{Int}}(k => Int[] for k in krange)
            J_b = Dict{Int, Vector{Int}}(k => Int[] for k in krange)
            V_b = Dict{Int, Vector{Float64}}(k => Float64[] for k in krange)
            ne_b = length(recv_b) ÷ 4
            for m = 1:ne_b
                base = (m-1)*4
                kk = Int(recv_b[base+1])
                ig = Int(recv_b[base+2])
                jg = Int(recv_b[base+3])
                v  = recv_b[base+4]
                push!(I_b[kk], ig); push!(J_b[kk], jg); push!(V_b[kk], v)
            end
            for k in krange
                A_batch[k] = sparse(I_b[k], J_b[k], V_b[k], gnpoin_3d_dense, gnpoin_3d_dense)
            end
        else
            for k in krange
                A_batch[k] = spzeros(Float64, gnpoin_3d_dense, gnpoin_3d_dense)
            end
        end
        return A_batch
    end

    t_setup = @elapsed begin
    _t_start = time()

    t_gather = 0.0; t_job1 = 0.0; t_job2 = 0.0; t_job3 = 0.0
    t_rhs    = 0.0; t_bcast = 0.0
    n_new_patterns    = 0
    n_reused_patterns = 0

    # Pattern-hash-keyed dict of MUMPS instances, bounded to the number of
    # DISTINCT sparsity patterns (empirically ~3-6), not the number of
    # directions (up to 600) — this is the same lookup whether or not a
    # persistent cache was supplied by the caller:
    #   - cache === nothing: a fresh, call-local dict, discarded (all
    #     instances finalized) at the end of this call — matches the
    #     within-call-only reuse already confirmed stable, just written once
    #     instead of duplicated into a separate code path.
    #   - cache !== nothing: the caller's persistent dict — instances survive
    #     past the end of this call for reuse on a later call (e.g. the next
    #     scattering outer-iteration, or LW after SW if patterns happen to
    #     coincide — though SW/LW should normally use SEPARATE cache objects
    #     since their κ/σ, and so likely their patterns, generally differ).
    pattern_dict = cache !== nothing ? cache.mumps_by_pattern :
                                        Dict{UInt64, MUMPS.Mumps{Float64}}()

    for batch_start in 1:gather_batch_size:extra_npoin
        batch_end   = min(batch_start + gather_batch_size - 1, extra_npoin)
        batch_range = batch_start:batch_end

        t0 = time()
        A_batch = gather_dir_matrix_batch(batch_range)
        t_gather += time() - t0

    for k in batch_range
        # Fixed, globally identical iteration order — every rank does the
        # exact same sequence of collective calls in the exact same order.
        A_k = A_batch[k]

        # ── Pattern hash decided on rank 0 (from the just-gathered, complete
        # matrix), then broadcast — the ONLY place this decision is made,
        # never from per-rank partial local data, which is what caused the
        # earlier hang.
        pat_hash_buf = zeros(UInt64, 1)
        if rank == 0
            pat_hash_buf[1] = hash(A_k.colptr) ⊻ hash(A_k.rowval)
        end
        MPI.Bcast!(pat_hash_buf, 0, comm)
        pat_hash = pat_hash_buf[1]

        is_new_pattern = !haskey(pattern_dict, pat_hash)
        if is_new_pattern
            n_new_patterns += 1
            mumps_k = MUMPS.Mumps{Float64}(MUMPS.mumps_unsymmetric,
                                            MUMPS.default_icntl[:],
                                            MUMPS.default_cntl64[:])
            MUMPS.set_icntl!(mumps_k, 1, -1; displaylevel=0)
            MUMPS.set_icntl!(mumps_k, 2, -1; displaylevel=0)
            MUMPS.set_icntl!(mumps_k, 3, -1; displaylevel=0)
            MUMPS.set_icntl!(mumps_k, 4,  1; displaylevel=0)
            pattern_dict[pat_hash] = mumps_k
        else
            n_reused_patterns += 1
            mumps_k = pattern_dict[pat_hash]
        end

        if rank == 0
            MUMPS.associate_matrix!(mumps_k, A_k)
        end
        if is_new_pattern
            t0 = time()
            MUMPS.set_job!(mumps_k, 1)
            MUMPS.invoke_mumps!(mumps_k)
            t_job1 += time() - t0
            if rank == 0
                MUMPS.associate_matrix!(mumps_k, A_k)
            end
        end
        t0 = time()
        MUMPS.set_job!(mumps_k, 2)
        MUMPS.invoke_mumps!(mumps_k)
        t_job2 += time() - t0

        t0 = time()
        fill!(b_3d, 0.0)
        if k in my_dirs_set
            for gid_5d = k:extra_npoin:gnpoin_5d
                dense_gip       = (gid_5d - k) ÷ extra_npoin + 1
                b_3d[dense_gip] = b_global[gid_5d]
            end
        end
        MPI.Allreduce!(b_3d, +, comm)
        t_rhs += time() - t0

        if rank == 0
            MUMPS.associate_rhs!(mumps_k, b_3d)
        end
        t0 = time()
        MUMPS.set_job!(mumps_k, 3)
        MUMPS.invoke_mumps!(mumps_k)
        t_job3 += time() - t0

        t0 = time()
        if rank == 0
            x_3d .= vec(MUMPS.get_solution(mumps_k))
        end
        MPI.Bcast!(x_3d, 0, comm)
        t_bcast += time() - t0

        if k in my_dirs_set
            for gid_5d = k:extra_npoin:gnpoin_5d
                dense_gip        = (gid_5d - k) ÷ extra_npoin + 1
                x_global[gid_5d] = x_3d[dense_gip]
            end
        end

        if rank == 0 && (k == 1 || k % 10 == 0 || k == extra_npoin)
            @info "AngularDirect: direction $k/$extra_npoin done " *
                  "(new_patterns=$n_new_patterns, reused_patterns=$n_reused_patterns), " *
                  "elapsed=$(round(time() - _t_start, digits=1))s"
            flush(stdout)
        end
    end
    end # for batch_start

    if cache === nothing
        for (_, m) in pattern_dict
            try MUMPS.finalize!(m) catch end
        end
    end

    if rank == 0
        @info "AngularDirect: timing breakdown — gather=$(round(t_gather,digits=2))s " *
              "JOB1(analysis)=$(round(t_job1,digits=2))s JOB2(factorize)=$(round(t_job2,digits=2))s " *
              "JOB3(solve)=$(round(t_job3,digits=2))s rhs_assemble=$(round(t_rhs,digits=2))s " *
              "bcast=$(round(t_bcast,digits=2))s"
        @info "AngularDirect: $(length(pattern_dict)) distinct sparsity patterns; " *
              "this call did $n_new_patterns new JOB=1 analyses and $n_reused_patterns " *
              "reused-pattern JOB=2-only refactorizations" *
              (cache !== nothing ? " (persistent cache — instances kept alive for next call)" :
                                    " (call-local — instances finalized at end of this call)")
        flush(stdout)
    end
    end # t_setup
    rank == 0 && @info "AngularDirect: factorize+solve done in $(round(t_setup,digits=2))s"
    flush(stdout)

    A_local_dirs = nothing

    # Combine solutions across ranks
    MPI.Allreduce!(x_global, +, comm)

    end # t_total
    rank == 0 && @info "AngularDirect: total solve time $(round(t_total,digits=2))s"
    flush(stdout)

    rank == 0 && @info "AngularDirect: solve complete"
    flush(stdout)

    # Return global vector indexed by dense 5D global IDs
    return x_global
end

# ─────────────────────────────────────────────────────────────────────────────
#  build_boomeramg_preconditioner
#
#  BoomerAMG algebraic multigrid via direct ccall to libHYPRE.
#  Uses hypre_boomeramg.jl module — no HYPRE.jl package needed.
#  Fully distributed — each rank provides owned rows directly to HYPRE.
#
#  Requires:
#    1. Build HYPRE from source: ./configure --prefix=$HOME/local
#         --with-MPI --enable-shared --with-openmp --disable-fortran
#    2. include("hypre_boomeramg.jl") before calling
#    3. HypreAMG.set_library!(get(ENV, "HYPRE_LIB", "/home/tissaoui/local/lib/libHYPRE.so"))
#    4. HypreAMG.init()
#    5. Add to job script: export HYPRE_LIB=$HOME/local/lib/libHYPRE.so
#                          export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
#
#  Tuning kwargs: CoarsenType, RelaxType, StrongThreshold, NumSweeps
# ─────────────────────────────────────────────────────────────────────────────

function build_boomeramg_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        npoin_g   :: Int         = length(ip2gip),
        g_ip2gip  :: Vector{Int} = Int[],
        comm      :: MPI.Comm    = MPI.COMM_WORLD,
        CoarsenType     :: Int     = 6,
        RelaxType       :: Int     = 6,
        MaxLevels       :: Int     = 25,
        NumSweeps       :: Int     = 1,
        StrongThreshold :: Float64 = 0.25,
        Restriction     :: Int     = 0,
        InterpType      :: Int     = -1,
        NonGalerkinTol  :: Float64 = 0.0,
        PrintLevel      :: Int     = 0)

    @isdefined(HypreAMG) ||
        error(":boomeramg requires:\n" *
              "  include(\"hypre_boomeramg.jl\")\n" *
              "  HypreAMG.set_library!(\"/home/tissaoui/local/lib/libHYPRE.so\")\n" *
              "  HypreAMG.init()")

    rank = MPI.Comm_rank(comm)
    t = @elapsed handle = HypreAMG.setup(
        A_local, ip2gip, gip2owner, gnpoin;
        npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm,
        CoarsenType=CoarsenType, RelaxType=RelaxType,
        MaxLevels=MaxLevels, NumSweeps=NumSweeps,
        StrongThreshold=StrongThreshold, PrintLevel=PrintLevel,
        Restriction=Restriction, InterpType=InterpType,
        NonGalerkinTol=NonGalerkinTol)

    rank == 0 && @info "BoomerAMG setup: $(round(t,digits=3))s"

    op = LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> HypreAMG.apply!(y, handle, x))

    # Register finalizer
    finalizer(_ -> HypreAMG.destroy!(handle), op)
    return op
end


# ─────────────────────────────────────────────────────────────────────────────
#  solve_parallel_gmres_asm
#
#  Preconditioner options:
#
#  LOCAL (ASM — each rank factorizes its own completed submatrix):
#    :klu       — KLU local factorization
#    :rcmsplu   — AMD + exact LU local
#    :rcmilu    — AMD + ILU(τ) local  [best for conforming parallel]
#    :rcmilu0   — AMD + ILU(0) local  [fastest setup, conforming parallel]
#    :splu/:lu  — other exact local solvers
#    :ilu/:ilu0/:spilu — other incomplete local solvers
#
#  GLOBAL (rank 0 holds and factorizes complete global matrix):
#    :global_klu    — KLU on global matrix [fastest setup, ~10 iters]
#    :global_lu     — AMD + exact LU on global matrix [best quality]
#    :global_ilu    — AMD + ILU(τ) on global matrix [fast setup, ~15 iters]
#    :global_ilu0   — AMD + ILU(0) on global matrix [fastest setup]
#
#  OTHER:
#    :jacobi        — diagonal scaling
#    :inner_gmres   — inner GMRES (diagnostic use)
#    :none          — no preconditioning
# ─────────────────────────────────────────────────────────────────────────────

function solve_parallel_gmres_asm(
        ip2gip, gip2owner, A_local, b, gnpoin, npoin, x_prev;
        precond     :: Symbol  = :rcmilu,
        restart     :: Int     = 50,
        tol         :: Float64 = 1e-6,
        itmax       :: Int     = 10000,
        npoin_g     :: Int     = 0,
        g_ip2gip    :: Vector{Int} = Int[],
        g_gip2ip    = Int[],
        asm_solver  :: Symbol  = :rcmilu,
        asm_ilu_tau :: Float64 = 0.01,
        inner_tol   :: Float64 = 1e-2,
        inner_itmax :: Int     = 10,
        A_pre_rp    :: Union{AbstractSparseMatrix{Float64}, Nothing} = nothing,
        b_pre_rp    :: Union{AbstractVector{Float64}, Nothing} = nothing,
        R_mat       = nothing,
        P_mat       = nothing,
        mumps_handle_ref   :: Union{Ref{Union{MUMPSHandle,Nothing}},  Nothing} = nothing,
        superlu_handle_ref :: Union{Ref{Union{SuperLUHandle,Nothing}}, Nothing} = nothing,
        prev_nnz_ref       :: Union{Ref{Int}, Nothing}                          = nothing,
        connijk_spa   = nothing,
        extra_nelem   = nothing,
        extra_nops    = nothing,
        extra_connijk = nothing,
        nelem         :: Int   = 0,
        ngl           :: Int   = 0,
        ladaptive     :: Bool  = false,
        npoin_ang     = nothing,
        npoin_space   :: Int   = 0,
        comm          :: MPI.Comm = MPI.COMM_WORLD,
        # BoomerAMG tuning (used when precond=:boomeramg)
        amg_coarsen      :: Int     = 6,
        amg_relax        :: Int     = 6,
        amg_threshold    :: Float64 = 0.25,
        amg_sweeps       :: Int     = 1,
        amg_levels       :: Int     = 25,
        amg_restriction  :: Int     = 0,    # 0=classical, 1=AIR-1, 2=AIR-2
        amg_interp       :: Int     = -1,   # -1=default, 4=pAIR-1, 5=pAIR-2
        amg_nongalerkin  :: Float64 = 0.0,  # >0 enables non-Galerkin coarsening
        amg_print        :: Int     = 0,
        # Spatial AMG kwargs (used when precond=:spatial_amg)
        spa_ip2gip       :: AbstractVector{Int}   = Int[],  # 3D spatial ip2gip (mesh.ip2gip)
        spa_gip2owner    :: AbstractVector{Int}   = Int[],  # 3D spatial gip2owner
        spa_gnpoin       :: Int                   = 0,      # 3D global DOF count
        spa_extra_npoin  :: Int                   = 0,      # angular DOFs per spatial node
        spa_connijk_spa  = nothing,                          # non-conforming angular connijk
        spa_mesh         = nothing,                          # mesh object
        spa_extra_mesh   = nothing,                          # conforming: extra_mesh object
        spa_extra_meshes = nothing,                          # non-conforming: per-element
        spa_npoin_g      :: Int                   = 0,
        spa_g_ip2gip     :: Vector{Int}           = Int[],
        spa_local_weight :: Float64               = 1.0,  # Jacobi smoother weight (0=disable)
        spa_angular_cache :: Union{AngularDirectCache, Nothing} = nothing)

    rank = MPI.Comm_rank(comm)

    # ── Angular direct solve (bypasses GMRES for block-diagonal systems) ──
    if precond == :angular_direct
        isempty(spa_ip2gip)   && error(":angular_direct requires spa_ip2gip (mesh.ip2gip)")
        spa_extra_npoin == 0  && error(":angular_direct requires spa_extra_npoin")
        x_sol = solve_angular_direct(
            A_local, b, ip2gip, spa_ip2gip,
            0,
            spa_gnpoin, spa_extra_npoin;
            comm=comm, cache=spa_angular_cache)
        # x_sol[gid] = solution at dense global 5D ID (1:gnpoin_3d_dense*extra_npoin)
        # ip2gip[ip] gives the dense global 5D ID for local IP ip
        gnpoin_5d_dense = spa_gnpoin  # placeholder — actual size from x_sol
        x_local = zeros(Float64, max(npoin_g, npoin))
        for ip = 1:npoin
            gid = ip2gip[ip]
            1 <= gid <= length(x_sol) && (x_local[ip] = x_sol[gid])
        end
        return x_local
    end

    A_parallel = create_parallel_linear_operator(
        A_local, ip2gip, gip2owner, npoin, gnpoin,
        npoin_g, g_ip2gip, g_gip2ip, comm)

    b_global = zeros(Float64, gnpoin)
    for ip = 1:npoin; b_global[ip2gip[ip]] += b[ip]; end
    MPI.Allreduce!(b_global, +, comm)

    x0 = zeros(Float64, gnpoin)
    if !isempty(x_prev)
        for ip = 1:npoin; x0[ip2gip[ip]] = x_prev[ip]; end
        MPI.Allreduce!(x0, +, comm)
        r0_rel = norm(b_global - A_parallel * x0) / max(norm(b_global), 1e-30)
        rank == 0 && @info "Warm start: residual=$(round(r0_rel,sigdigits=3))"
        if r0_rel < tol
            rank == 0 && @info "Warm start: already converged."
            x_local = zeros(Float64, max(npoin_g, npoin))
            for ip = 1:npoin; x_local[ip] = x0[ip2gip[ip]]; end
            return x_local
        end
    end

    # ── Route preconditioner symbol ───────────────────────────────────────────
    # Local ASM solvers — route through build_asm_preconditioner
    local_solvers = (:klu, :splu, :lu, :ilu, :ilu0, :spilu,
                     :rcmsplu, :rcmilu, :rcmilu0)
    # Global solvers — route through build_global_factorization_preconditioner
    global_solver_map = Dict(
        :global_klu      => :klu,
        :global_lu       => :rcmsplu,
        :global_ilu      => :rcmilu,
        :global_ilu0     => :rcmilu0,
        :global_paru        => :paru,         # ParU (OpenMP)
        :global_pardiso     => :pardiso,      # MKL Pardiso (OpenMP)
        :global_superlu     => :superlu_dist, # SuperLU_DIST (MPI+OpenMP)
    )

    effective_precond = precond
    effective_solver  = asm_solver

    if precond in local_solvers
        effective_precond = :asm
        effective_solver  = precond
    elseif haskey(global_solver_map, precond)
        effective_precond = :global
        effective_solver  = global_solver_map[precond]
    elseif precond == :mumps_dist
        effective_precond = :mumps_dist
    elseif precond == :boomeramg
        effective_precond = :boomeramg
    elseif precond == :spatial_amg
        effective_precond = :spatial_amg
    elseif precond == :diag_jacobi
        effective_precond = :diag_jacobi
    elseif precond == :superlu_dist
        effective_precond = :superlu_dist
    end

    # ── Build preconditioner ──────────────────────────────────────────────────
    N = if effective_precond == :asm
        build_asm_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
            solver=effective_solver, npoin_g=npoin_g, g_ip2gip=g_ip2gip,
            ilu_tau=asm_ilu_tau, comm=comm)

    elseif effective_precond == :global
        build_global_factorization_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
            solver=effective_solver, ilu_tau=asm_ilu_tau,
            npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)

    elseif effective_precond == :boomeramg
        # BoomerAMG algebraic multigrid — fully distributed, no gather step.
        # Requires: using HYPRE; HYPRE.Init() before first call.
        # Tuning: pass CoarsenType, RelaxType, StrongThreshold as kwargs.
        # For advection-dominated (radiation transport): try RelaxType=16
        # For strongly non-symmetric: try CoarsenType=8, StrongThreshold=0.5
        build_boomeramg_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
            npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm,
            CoarsenType=amg_coarsen, RelaxType=amg_relax,
            StrongThreshold=amg_threshold, NumSweeps=amg_sweeps,
            MaxLevels=amg_levels, Restriction=amg_restriction,
            InterpType=amg_interp,
            NonGalerkinTol=amg_nongalerkin, PrintLevel=amg_print)

    elseif effective_precond == :spatial_amg
        isempty(spa_ip2gip)    && error(":spatial_amg requires spa_ip2gip (mesh.ip2gip, 3D spatial)")
        spa_gnpoin == 0        && error(":spatial_amg requires spa_gnpoin (3D global DOF count)")
        spa_extra_npoin == 0   && spa_connijk_spa === nothing &&
            error(":spatial_amg requires spa_extra_npoin (conforming) or spa_connijk_spa (non-conforming)")
        _npoin_g_3d = spa_npoin_g > 0 ? spa_npoin_g : length(spa_ip2gip)
        build_spatial_amg_preconditioner(
            A_local, ip2gip, spa_ip2gip, spa_gip2owner, spa_gnpoin, gnpoin,
            spa_extra_npoin;
            connijk_spa=spa_connijk_spa, mesh=spa_mesh,
            extra_mesh=spa_extra_mesh, extra_meshes=spa_extra_meshes,
            npoin_g_5d=npoin_g, g_ip2gip_5d=g_ip2gip,
            npoin_g_3d=_npoin_g_3d, g_ip2gip_3d=spa_g_ip2gip,
            comm=comm,
            CoarsenType=amg_coarsen, RelaxType=amg_relax,
            StrongThreshold=amg_threshold, NumSweeps=amg_sweeps,
            MaxLevels=amg_levels, PrintLevel=amg_print,
            local_weight=spa_local_weight)

    elseif effective_precond == :diag_jacobi
        # Jacobi preconditioner with strengthened diagonal.
        # For each row: d_eff = max(|diag|, threshold * max|off-diag|)
        # This handles badly conditioned rows where streaming >> absorption.
        rows_AL = rowvals(A_local); vals_AL = nonzeros(A_local)
        d_eff = zeros(Float64, gnpoin)
        for j = 1:npoin
            gid_j = ip2gip[j]
            diag_j = 0.0; max_off = 0.0
            for idx in nzrange(A_local, j)
                ip_row = rows_AL[idx]
                ip_row > npoin && continue
                if ip_row == j
                    diag_j = vals_AL[idx]
                else
                    max_off = max(max_off, abs(vals_AL[idx]))
                end
            end
            d_eff[gid_j] = abs(diag_j) < spa_local_weight * max_off ?
                           sign(diag_j) * spa_local_weight * max_off : diag_j
            abs(d_eff[gid_j]) < 1e-30 && (d_eff[gid_j] = 1.0)
        end
        MPI.Allreduce!(d_eff, +, comm)
        rank == 0 && @info "Jacobi: d_eff range [$(minimum(abs.(d_eff[d_eff.!=0]))), $(maximum(abs.(d_eff)))]"
        flush(stdout)
        LinearOperator(Float64, gnpoin, gnpoin, false, false,
            (y, x) -> begin
                @inbounds for i = 1:gnpoin; y[i] = x[i] / d_eff[i]; end
                y
            end)

    elseif effective_precond == :superlu_dist
        # SuperLU_DIST with optional handle caching for timestep reuse.
        # No separate analysis phase — cost same regardless of mesh change.
        _cur_nnz_slu = nnz(A_local)
        if superlu_handle_ref === nothing || superlu_handle_ref[] === nothing
            _sh = build_superlu_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
                npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)
            superlu_handle_ref !== nothing && (superlu_handle_ref[] = _sh)
            # If no ref, GC will clean up via the finalizer in build_superlu_preconditioner
            _sh.precond_op
        else
            _info = (prev_nnz_ref !== nothing && prev_nnz_ref[] != _cur_nnz_slu) ?
                "mesh changed" : "mesh unchanged"
            rank == 0 && @info "SuperLU_DIST: refactorizing ($_info)"
            superlu_refactorize!(superlu_handle_ref[], A_local, ip2gip, gip2owner, gnpoin;
                npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)
            superlu_handle_ref[].precond_op
        end

    elseif effective_precond == :mumps_dist
        # MUMPS with optional handle caching for timestep reuse.
        # Pass mumps_handle_ref and prev_nnz_ref (declared outside timestepping loop):
        #   - First call / no ref: full analysis + factorization
        #   - Mesh unchanged: numeric refactorization only (skips METIS)
        #   - Mesh changed: full analysis + factorization
        _cur_nnz_m = nnz(A_local)
        _mesh_changed_m = prev_nnz_ref !== nothing && prev_nnz_ref[] != _cur_nnz_m
        if mumps_handle_ref === nothing || mumps_handle_ref[] === nothing
            _mh = build_mumps_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
                npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)
            if mumps_handle_ref !== nothing
                mumps_handle_ref[] = _mh
            else
                # No ref — register finalizer so GC cleans up C-side memory
                finalizer(h -> try; @isdefined(MUMPS) && MUMPS.finalize(h.mumps_obj); catch; end, _mh)
            end
            _mh.precond_op
        elseif _mesh_changed_m
            rank == 0 && @info "MUMPS: mesh changed (nnz $(prev_nnz_ref[])→$_cur_nnz_m), full analysis"
            mumps_refactorize!(mumps_handle_ref[], A_local, ip2gip, gip2owner, gnpoin;
                npoin_g=npoin_g, g_ip2gip=g_ip2gip, full_analysis=true, comm=comm)
            mumps_handle_ref[].precond_op
        else
            rank == 0 && @info "MUMPS: mesh unchanged, numeric refactorization only"
            mumps_refactorize!(mumps_handle_ref[], A_local, ip2gip, gip2owner, gnpoin;
                npoin_g=npoin_g, g_ip2gip=g_ip2gip, full_analysis=false, comm=comm)
            mumps_handle_ref[].precond_op
        end

    elseif effective_precond == :inner_gmres
        build_inner_gmres_preconditioner(A_parallel, gnpoin;
            inner_tol=inner_tol, inner_itmax=inner_itmax)

    elseif precond == :jacobi
        build_jacobi_preconditioner(A_local, ip2gip, gnpoin, comm)

    else
        # :none or unknown — identity preconditioner
        LinearOperator(Float64, gnpoin, gnpoin, true, true, (y,x)->copyto!(y,x))
    end

    rank == 0 && @info "GMRES: precond=$precond (solver=$effective_solver) " *
                       "restart=$restart tol=$tol DOFs=$gnpoin ranks=$(MPI.Comm_size(comm))"

    # Global preconditioners are nonlinear (contain MPI collectives) → fgmres
    use_fgmres = effective_precond in (:global, :inner_gmres, :mumps_dist, :superlu_dist, :boomeramg, :spatial_amg, :diag_jacobi)
    krylov_solve = use_fgmres ? Krylov.fgmres : Krylov.gmres

    t_solve = @elapsed begin
        x, stats = if isempty(x_prev)
            krylov_solve(A_parallel, b_global; N=N, memory=restart, restart=true,
                         atol=tol, rtol=tol, itmax=itmax, verbose=(rank==0) ? 1 : 0)
        else
            krylov_solve(A_parallel, b_global, x0; N=N, memory=restart, restart=true,
                         atol=tol, rtol=tol, itmax=itmax, verbose=(rank==0) ? 1 : 0)
        end
    end

    if rank == 0
        final_res = isempty(stats.residuals) ? 0.0 : stats.residuals[end]
        @info "GMRES done: $(stats.niter) iters, $(round(t_solve,digits=2))s, " *
              "res=$(round(final_res,sigdigits=4)), converged=$(stats.solved)"
        !stats.solved && @warn "GMRES did not converge."
    end

    prev_nnz_ref !== nothing && (prev_nnz_ref[] = nnz(A_local))

    x_local = zeros(Float64, max(npoin_g, npoin))
    for ip = 1:npoin; x_local[ip] = x[ip2gip[ip]]; end
    return x_local
end

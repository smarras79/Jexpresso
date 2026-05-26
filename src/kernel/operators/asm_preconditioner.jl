using SparseArrays
using LinearAlgebra
using KLU
using ILUZero
using IncompleteLU
using Krylov
using LinearOperators
using AMD

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
#  Used for :paru, :pardiso, :strumpack solvers.
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

    elseif solver == :strumpack
        # STRUMPACK: sparse direct solver with optional low-rank compression.
        # Supports both shared (OpenMP) and distributed (MPI) parallelism.
        # Requires STRUMPACK_jll to be imported BEFORE using LinearSolve
        # so the extension is loaded:
        #   import STRUMPACK_jll; using LinearSolve
        # Install: ] add LinearSolve STRUMPACK_jll
        @isdefined(LinearSolve) ||
            error(":strumpack requires `import STRUMPACK_jll; using LinearSolve`")
        isdefined(LinearSolve, :STRUMPACKFactorization) ||
            error(":strumpack: STRUMPACKFactorization not found in LinearSolve.
" *
                  "  1. Ensure STRUMPACK_jll is installed: ] add STRUMPACK_jll
" *
                  "  2. Import BEFORE using LinearSolve: import STRUMPACK_jll; using LinearSolve
" *
                  "  3. Verify libstrumpack is discoverable (check with Libdl.find_library(["libstrumpack"]))")
        n = size(A, 1)
        b_dummy = zeros(Float64, n)
        prob  = LinearSolve.LinearProblem(A, b_dummy)
        # STRUMPACKFactorization targets single-node multithreaded sparse LU
        # Thread count controlled by OMP_NUM_THREADS environment variable
        cache = LinearSolve.init(prob, LinearSolve.STRUMPACKFactorization())
        LinearSolve.solve!(cache)
        return LinearSolveFactor(cache, n)

    elseif solver == :mumps
        # MUMPS: multifrontal parallel sparse direct solver (JuliaSmoothOptimizers).
        # Install: ] add MUMPS;  Usage: using MUMPS
        # NOTE: Known MPI communicator conflicts on macOS with OpenMPI.
        #       Recommended for cluster use with system MPI only.
        @isdefined(MUMPS) ||
            error(":mumps requires `using MUMPS` (] add MUMPS)")
        # High-level API: mumps_factorize takes SparseMatrixCSC directly
        mumps_obj = MUMPS.Mumps{Float64}(MUMPS.mumps_unsymmetric,
                                          MUMPS.default_icntl,
                                          MUMPS.default_cntl64)
        # Suppress most MUMPS output
        MUMPS.set_icntl!(mumps_obj, 1, -1)
        MUMPS.set_icntl!(mumps_obj, 2, -1)
        MUMPS.set_icntl!(mumps_obj, 3, -1)
        MUMPS.set_icntl!(mumps_obj, 4, 1)

        # ── Dynamic memory allocation ────────────────────────────────────────
        # Run analysis phase only (JOB=1) to get memory estimates,
        # then set ICNTL(23) (memory per process in MB) from the estimate
        # with a safety margin. This scales automatically with problem size.
        MUMPS.associate_matrix!(mumps_obj, A)
        MUMPS.set_job!(mumps_obj, 1)   # analysis only
        MUMPS.invoke_mumps!(mumps_obj)
        analysis_err = MUMPS.get_infog(mumps_obj, 1)
        if analysis_err == 0
            # INFOG(16) = estimated max memory per process (MB)
            # INFOG(17) = estimated total memory (MB)
            mem_per_proc = MUMPS.get_infog(mumps_obj, 16)
            mem_total    = MUMPS.get_infog(mumps_obj, 17)
            # Add 50% safety margin and set as memory limit per process
            mem_limit = ceil(Int, mem_per_proc * 1.5)
            MUMPS.set_icntl!(mumps_obj, 23, mem_limit)
            @info "MUMPS analysis: ~$(mem_per_proc)MB/proc, ~$(mem_total)MB total → allocating $(mem_limit)MB/proc"
        else
            # Analysis failed — fall back to generous fixed relaxation
            @warn "MUMPS analysis phase failed (INFO(1)=$analysis_err), using ICNTL(14)=200"
            MUMPS.set_icntl!(mumps_obj, 14, 200)
        end

        # ── Factorization (JOB=2) ────────────────────────────────────────────
        MUMPS.set_job!(mumps_obj, 2)
        MUMPS.invoke_mumps!(mumps_obj)
        infog1 = MUMPS.get_infog(mumps_obj, 1)
        if infog1 == -10
            # Still out of memory — retry with doubled limit
            mem_limit2 = ceil(Int, MUMPS.get_infog(mumps_obj, 16) * 3.0)
            @warn "MUMPS out of memory (INFO(1)=-10), retrying with $(mem_limit2)MB/proc"
            MUMPS.set_icntl!(mumps_obj, 23, mem_limit2)
            MUMPS.set_job!(mumps_obj, 2)
            MUMPS.invoke_mumps!(mumps_obj)
            infog1 = MUMPS.get_infog(mumps_obj, 1)
        end
        infog1 == 0 || error("MUMPS factorization failed INFO(1)=$infog1 INFO(2)=$(MUMPS.get_infog(mumps_obj,2))")
        return mumps_obj

    else
        error("Unknown solver: $solver. Options:
" *
              "  Local:   :klu, :splu, :lu, :rcmsplu, :rcmilu, :rcmilu0, :ilu, :ilu0, :spilu
" *
              "  Parallel (require extra packages):
" *
              "    :paru      — ParU (OpenMP, SuiteSparse) — add LinearSolve ParU_jll
" *
              "    :pardiso   — MKL Pardiso (OpenMP, Intel MKL) — add LinearSolve
" *
              "    :strumpack — STRUMPACK (OpenMP/MPI) — add LinearSolve STRUMPACK_jll
" *
              "    :mumps     — MUMPS (MPI distributed) — add MUMPS_jll")
    end
end

# Note: MUMPS.jl (JuliaSmoothOptimizers) defines ldiv! for Mumps objects automatically.
# LinearSolveFactor.ldiv! handles ParU/Pardiso/STRUMPACK via LinearSolve.jl cache.

# ─────────────────────────────────────────────────────────────────────────────
#  build_asm_preconditioner
#
#  Local Additive Schwarz preconditioner.
#  Each rank factorizes its own completed local matrix.
#  Best for: conforming parallel (~27 iters), non-conforming serial (~6 iters).
#  Non-conforming parallel: ~400 iterations (limited by cross-rank coupling).
# ─────────────────────────────────────────────────────────────────────────────

function build_asm_preconditioner(
        A_local   :: AbstractSparseMatrix{Float64},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        gnpoin    :: Int;
        solver    :: Symbol      = :klu,
        npoin_g   :: Int         = 0,
        g_ip2gip  :: Vector{Int} = Int[],
        ilu_tau   :: Float64     = 0.01,
        comm      :: MPI.Comm    = MPI.COMM_WORLD)

    rank  = MPI.Comm_rank(comm)
    npoin = length(ip2gip)

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
    all_ij = MPI.Allgatherv(send_ij, Int32.(counts .* 2), comm)
    all_v  = MPI.Allgatherv(local_V, counts, comm)

    # Rank 0 assembles and factorizes
    global_factor = if rank == 0
        n_total = length(all_v)
        GI = [all_ij[2k-1] for k=1:n_total]
        GJ = [all_ij[2k]   for k=1:n_total]
        A_global = sparse(GI, GJ, all_v, actual_gnpoin, actual_gnpoin)
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
#  Contrast with :global_mumps which gathers to rank 0 first (cheaper for
#  small rank counts but loses distributed scalability).
#
#  Requires: using MUMPS (] add MUMPS)
#  MUMPS must be initialised within the same MPI communicator.
#  NOTE: Known MPI communicator conflicts on macOS with OpenMPI.
#        Recommended for cluster use with system MPI only.
# ─────────────────────────────────────────────────────────────────────────────

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

    # Build global GID-indexed local matrix contribution.
    # Each rank provides its owned rows in global (GID) indexing.
    # MUMPS assembles the global system from all ranks' contributions.
    has_ghosts_g = npoin_g > npoin
    ip2gip_g_loc = Vector{Int}(undef, npoin_g)
    ip2gip_g_loc[1:npoin] .= ip2gip
    if has_ghosts_g; ip2gip_g_loc[npoin+1:npoin_g] .= g_ip2gip; end

    # Actual global DOF count
    local_max_gid = maximum(ip2gip)
    actual_gnpoin = MPI.Allreduce(local_max_gid, MPI.MAX, comm)

    # Build GID-indexed local matrix (owned rows only, all cols)
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

    # Assemble local sparse matrix in global GID indexing
    A_local_gid = sparse(local_I, local_J, local_V,
                         actual_gnpoin, actual_gnpoin)

    # Initialise MUMPS — all ranks participate
    mumps_handle = MUMPS.Mumps{Float64}(
        MUMPS.mumps_unsymmetric,
        MUMPS.default_icntl,
        MUMPS.default_cntl64)

    # Suppress MUMPS output
    MUMPS.set_icntl!(mumps_handle, 1, -1)
    MUMPS.set_icntl!(mumps_handle, 2, -1)
    MUMPS.set_icntl!(mumps_handle, 3, -1)
    MUMPS.set_icntl!(mumps_handle, 4, 1)

    # Each rank provides its local contribution; MUMPS assembles globally
    MUMPS.associate_matrix!(mumps_handle, A_local_gid)

    # ── Dynamic memory: analysis phase first ─────────────────────────────────
    MUMPS.set_job!(mumps_handle, 1)
    MUMPS.invoke_mumps!(mumps_handle)
    analysis_err = MUMPS.get_infog(mumps_handle, 1)
    if analysis_err == 0
        mem_per_proc = MUMPS.get_infog(mumps_handle, 16)
        mem_total    = MUMPS.get_infog(mumps_handle, 17)
        mem_limit    = ceil(Int, mem_per_proc * 1.5)
        MUMPS.set_icntl!(mumps_handle, 23, mem_limit)
        rank == 0 && @info "MUMPS dist analysis: ~$(mem_per_proc)MB/proc, ~$(mem_total)MB total → $(mem_limit)MB/proc"
    else
        rank == 0 && @warn "MUMPS dist analysis failed (INFO(1)=$analysis_err), using ICNTL(14)=200"
        MUMPS.set_icntl!(mumps_handle, 14, 200)
    end

    # ── Factorization ─────────────────────────────────────────────────────────
    MUMPS.set_job!(mumps_handle, 2)
    t = @elapsed MUMPS.invoke_mumps!(mumps_handle)
    infog1 = MUMPS.get_infog(mumps_handle, 1)
    if infog1 == -10
        mem_limit2 = ceil(Int, MUMPS.get_infog(mumps_handle, 16) * 3.0)
        rank == 0 && @warn "MUMPS out of memory, retrying with $(mem_limit2)MB/proc"
        MUMPS.set_icntl!(mumps_handle, 23, mem_limit2)
        MUMPS.set_job!(mumps_handle, 2)
        t = @elapsed MUMPS.invoke_mumps!(mumps_handle)
        infog1 = MUMPS.get_infog(mumps_handle, 1)
    end
    if infog1 != 0
        rank == 0 && @warn "MUMPS distributed factorization failed INFO(1)=$infog1 — preconditioner will be identity"
    else
        rank == 0 && @info "MUMPS distributed factorization: $(round(t,digits=3))s"
    end

    # Work vector for RHS/solution
    rhs_buf = zeros(Float64, actual_gnpoin)
    y_global = zeros(Float64, gnpoin)

    return LinearOperator(Float64, gnpoin, gnpoin, false, false,
        (y, x) -> begin
            # Gather owned DOF values into complete global RHS on all ranks
            fill!(rhs_buf, 0.0)
            @inbounds for ip = 1:npoin
                gip2owner[ip] == rank || continue
                gid = ip2gip[ip]
                gid <= actual_gnpoin || continue
                rhs_buf[gid] = x[gid]
            end
            MPI.Allreduce!(rhs_buf, MPI.SUM, comm)

            # MUMPS solve — all ranks participate
            MUMPS.associate_rhs!(mumps_handle, rhs_buf)
            MUMPS.solve!(mumps_handle)

            # Retrieve solution on rank 0, broadcast to all
            fill!(y_global, 0.0)
            if rank == 0
                sol = MUMPS.get_solution(mumps_handle)
                copyto!(view(y_global, 1:actual_gnpoin), sol)
            end
            MPI.Bcast!(y_global, 0, comm)
            copyto!(y, y_global)
            return y
        end)
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
        connijk_spa   = nothing,
        extra_nelem   = nothing,
        extra_nops    = nothing,
        extra_connijk = nothing,
        nelem         :: Int   = 0,
        ngl           :: Int   = 0,
        ladaptive     :: Bool  = false,
        npoin_ang     = nothing,
        npoin_space   :: Int   = 0,
        comm          :: MPI.Comm = MPI.COMM_WORLD)

    rank = MPI.Comm_rank(comm)

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
        :global_paru     => :paru,      # ParU parallel (OpenMP)
        :global_pardiso  => :pardiso,   # MKL Pardiso parallel (OpenMP)
        :global_strumpack => :strumpack, # STRUMPACK (OpenMP/MPI)
        :global_mumps    => :mumps,     # MUMPS distributed (MPI)
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

    elseif effective_precond == :mumps_dist
        # Fully distributed MUMPS — no gather, true distributed factorization
        build_mumps_preconditioner(A_local, ip2gip, gip2owner, gnpoin;
            npoin_g=npoin_g, g_ip2gip=g_ip2gip, comm=comm)

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
    use_fgmres = effective_precond in (:global, :inner_gmres, :mumps_dist)
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

    x_local = zeros(Float64, max(npoin_g, npoin))
    for ip = 1:npoin; x_local[ip] = x[ip2gip[ip]]; end
    return x_local
end
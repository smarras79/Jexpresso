
 
# ─────────────────────────────────────────────────────────────────────────────
#  MatvecCacheSparse
#
#  Sparse point-to-point matvec assembler for nprocs > 4.
#  Sends only nonzero contributions using GLOBAL indices directly —
#  eliminates the proc_ip2gip table (nprocs × npoin per rank) entirely.
#
#  Memory per rank:
#    Old: O(nprocs × npoin)  for proc_ip2gip
#    New: O(nprocs + nnz_local) — only counts and nonzero buffers
#
#  The nnzeros, request arrays, and send/recv buffers are pre-allocated
#  and reused across matvec calls to avoid allocation in the hot path.
# ─────────────────────────────────────────────────────────────────────────────
 
mutable struct MatvecCacheSparse
    # Number of nonzero send contributions from this rank (updated each call)
    nnzero_local    :: Int
    # Number of nonzero contributions received from each rank (updated each call)
    nnzeros         :: Vector{Int}
    # Send buffers: global indices and values of nonzero contributions
    send_gip        :: Vector{Int}       # global indices (pre-allocated to npoin_g)
    send_val        :: Vector{Float64}   # values (pre-allocated to npoin_g)
    # Receive buffers: resized each call based on nnzeros exchange
    recv_gip        :: Vector{Vector{Int}}
    recv_val        :: Vector{Vector{Float64}}
    # Pre-allocated MPI request arrays (2*(nprocs-1) requests per round)
    req_counts      :: MPI.MultiRequest
    req_indices     :: MPI.MultiRequest
    req_values      :: MPI.MultiRequest
    nprocs          :: Int
    rank            :: Int
end
 
function setup_MatvecSparse_assembler(ip2gip_g::AbstractVector{Int},
                                      gip2owner::AbstractVector{Int},
                                      comm::MPI.Comm = MPI.COMM_WORLD)
    rank    = MPI.Comm_rank(comm)
    nprocs  = MPI.Comm_size(comm)
    npoin_g = length(ip2gip_g)
 
    nreq = 2 * (nprocs - 1)
 
    return MatvecCacheSparse(
        0,                                              # nnzero_local
        zeros(Int, nprocs),                            # nnzeros
        zeros(Int,     npoin_g),                       # send_gip (max size)
        zeros(Float64, npoin_g),                       # send_val (max size)
        [Int[]     for _ in 1:nprocs],                 # recv_gip
        [Float64[] for _ in 1:nprocs],                 # recv_val
        MPI.MultiRequest(nreq),                        # req_counts
        MPI.MultiRequest(nreq),                        # req_indices
        MPI.MultiRequest(nreq),                        # req_values
        nprocs, rank)
end
 
function assemble_mpi_matvec_sparse!(y_local  :: AbstractVector{Float64},
                                     y_global :: AbstractVector{Float64},
                                     ip2gip_g :: AbstractVector{Int},
                                     cache    :: MatvecCacheSparse,
                                     comm     :: MPI.Comm = MPI.COMM_WORLD)
    nprocs = cache.nprocs
    rank   = cache.rank
    npoin_g = length(ip2gip_g)
 
    # ── Pack nonzero contributions as (global_index, value) pairs ────────────
    # Send global indices directly — no proc_ip2gip lookup needed on receiver.
    nnz = 0
    @inbounds for ip = 1:npoin_g
        val = y_local[ip]
        val == zero(Float64) && continue
        nnz += 1
        cache.send_gip[nnz] = ip2gip_g[ip]
        cache.send_val[nnz] = val
    end
    cache.nnzero_local = nnz
 
    # ── Round 1: exchange nnzero counts ──────────────────────────────────────
    req_idx = 1
    @inbounds for i in 0:nprocs-1
        i == rank && continue
        MPI.Isend(nnz, comm, cache.req_counts[req_idx]; dest=i, tag=10)
        req_idx += 1
        MPI.Irecv!(@view(cache.nnzeros[i+1:i+1]), comm,
                   cache.req_counts[req_idx]; source=i, tag=10)
        req_idx += 1
    end
    cache.nnzeros[rank+1] = nnz
    MPI.Waitall(cache.req_counts)
 
    # ── Resize receive buffers to match incoming counts ───────────────────────
    # Only reallocates when size changes — stable after warmup iterations.
    @inbounds for i in 0:nprocs-1
        i == rank && continue
        n = cache.nnzeros[i+1]
        length(cache.recv_gip[i+1]) != n && resize!(cache.recv_gip[i+1], n)
        length(cache.recv_val[i+1]) != n && resize!(cache.recv_val[i+1], n)
    end
 
    # ── Round 2: exchange global indices ─────────────────────────────────────
    req_idx = 1
    @inbounds for i in 0:nprocs-1
        i == rank && continue
        send_view = @view(cache.send_gip[1:nnz])
        MPI.Isend(send_view, comm, cache.req_indices[req_idx]; dest=i, tag=11)
        req_idx += 1
        n = cache.nnzeros[i+1]
        recv_view = n > 0 ? cache.recv_gip[i+1] : @view(cache.recv_gip[i+1][1:0])
        MPI.Irecv!(recv_view, comm, cache.req_indices[req_idx]; source=i, tag=11)
        req_idx += 1
    end
    MPI.Waitall(cache.req_indices)
 
    # ── Round 3: exchange values ──────────────────────────────────────────────
    req_idx = 1
    @inbounds for i in 0:nprocs-1
        i == rank && continue
        send_view = @view(cache.send_val[1:nnz])
        MPI.Isend(send_view, comm, cache.req_values[req_idx]; dest=i, tag=12)
        req_idx += 1
        n = cache.nnzeros[i+1]
        recv_view = n > 0 ? cache.recv_val[i+1] : @view(cache.recv_val[i+1][1:0])
        MPI.Irecv!(recv_view, comm, cache.req_values[req_idx]; source=i, tag=12)
        req_idx += 1
    end
    MPI.Waitall(cache.req_values)
 
    # ── Accumulate: self contributions + received ─────────────────────────────
    fill!(y_global, zero(Float64))
 
    # Self contributions
    @inbounds for j = 1:nnz
        y_global[cache.send_gip[j]] += cache.send_val[j]
    end
 
    # Remote contributions — global index used directly, no lookup needed
    @inbounds for i in 0:nprocs-1
        i == rank && continue
        n = cache.nnzeros[i+1]
        gips = cache.recv_gip[i+1]
        vals = cache.recv_val[i+1]
        for j = 1:n
            y_global[gips[j]] += vals[j]
        end
    end
end
 
# ─────────────────────────────────────────────────────────────────────────────
#  create_parallel_linear_operator
#
#  Builds a parallel matrix-vector product operator for distributed GMRES.
#  Uses Allreduce for nprocs <= 4 (simpler, lower latency at small scale).
#  Uses sparse point-to-point for nprocs > 4 (avoids O(gnpoin) Allreduce).
#
#  Memory improvements over original:
#  - No proc_ip2gip table (was nprocs × npoin per rank)
#  - Send/recv buffers sized to max local nnz, not full npoin
#  - Global indices sent directly, eliminating indirect lookup on receive
# ─────────────────────────────────────────────────────────────────────────────
 
function create_parallel_linear_operator(
        A_local   :: AbstractMatrix{T},
        ip2gip    :: AbstractVector{Int},
        gip2owner :: AbstractVector{Int},
        npoin     :: Int,
        gnpoin    :: Int,
        npoin_g   :: Int,
        g_ip2gip  :: AbstractVector{Int},
        g_gip2ip,
        comm      :: MPI.Comm = MPI.COMM_WORLD) where T
 
    nprocs = MPI.Comm_size(comm)
    rank   = MPI.Comm_rank(comm)
 
    # ── Pre-allocate fixed-size work arrays ───────────────────────────────────
    y_local  = zeros(Float64, gnpoin)
    y_temp   = zeros(Float64, npoin_g)
    x_local  = zeros(Float64, npoin_g)
 
    # Build extended ip2gip covering ghost nodes
    ip2gip_g = zeros(Int, npoin_g)
    ip2gip_g[1:npoin] .= ip2gip
    if npoin_g > npoin
        @inbounds for ip = npoin+1:npoin_g
            ip2gip_g[ip] = g_ip2gip[ip - npoin]
        end
    end
 
    # Sparse assembler cache (only needed for nprocs > 4)
    pM = nprocs > 4 ? setup_MatvecSparse_assembler(ip2gip_g, gip2owner, comm) :
                      nothing
 
    # ── Matvec ────────────────────────────────────────────────────────────────
    function matvec!(y, x)
        # Gather x into local extended vector
        @inbounds for ip = 1:npoin_g
            x_local[ip] = x[ip2gip_g[ip]]
        end
 
        # Local sparse matvec
        mul!(y_temp, A_local, x_local)
 
        if nprocs <= 4
            # Small rank count: Allreduce is simpler and lower latency
            fill!(y_local, zero(T))
            @inbounds for ip = 1:npoin_g
                y_local[ip2gip_g[ip]] += y_temp[ip]
            end
            MPI.Allreduce!(y_local, +, comm)
            copyto!(y, y_local)
        else
            # Large rank count: sparse point-to-point assembly
            assemble_mpi_matvec_sparse!(y_temp, y_local, ip2gip_g, pM, comm)
            copyto!(y, y_local)
        end
        return y
    end
 
    function rmatvec!(y, x)
        @inbounds for ip = 1:npoin_g
            x_local[ip] = x[ip2gip_g[ip]]
        end
        mul!(y_temp, A_local', x_local)
        if nprocs <= 4
            fill!(y_local, zero(T))
            @inbounds for ip = 1:npoin_g
                y_local[ip2gip_g[ip]] += y_temp[ip]
            end
            MPI.Allreduce!(y_local, +, comm)
            copyto!(y, y_local)
        else
            assemble_mpi_matvec_sparse!(y_temp, y_local, ip2gip_g, pM, comm)
            copyto!(y, y_local)
        end
        return y
    end
 
    return LinearOperator{T}(gnpoin, gnpoin, false, false, matvec!, rmatvec!, rmatvec!)
end


function build_block_jacobi_preconditioner_adaptive(
        A_local       :: AbstractMatrix{Float64},
        ip2gip        :: AbstractVector{Int},
        gnpoin        :: Int,
        connijk_spa,
        extra_nelem,
        extra_nops,
        extra_connijk,
        nelem         :: Int,
        ngl           :: Int,
        ladaptive     :: Bool,
        npoin         :: Int,
        npoin_ang     :: Int,        # used only for non-adaptive case
        comm          :: MPI.Comm = MPI.COMM_WORLD)

    rank = MPI.Comm_rank(comm)

    # ── Step 1: build per-spatial-node list of (local_ipg, global_ipg) ───────
    # Key insight: ip2gip maps local combined index -> global combined index
    # for both adaptive and non-adaptive cases. We collect all local ip_g
    # indices that belong to each spatial node ip, then translate to global.

    # local_ipg_sets[ip] = sorted vector of local combined DOF indices
    # global_ipg_sets[ip] = corresponding global combined DOF indices
    local_ipg_sets  = [Int[] for _ = 1:npoin]
    global_ipg_sets = [Int[] for _ = 1:npoin]

    if ladaptive
        for iel = 1:nelem
            n_ang_elem = extra_nelem[iel]
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                # Get spatial node for this (iel,i,j,k)
                # connijk_spa[iel][i,j,k,e_ext,iθ,iϕ] gives local ip_g
                # We need the spatial ip — use the spatial connijk
                ip = mesh_local.connijk[iel, i, j, k]
                for e_ext = 1:n_ang_elem
                    nop = extra_nops[iel][e_ext]
                    for iθ = 1:nop+1, iϕ = 1:nop+1
                        ip_g_local  = connijk_spa[iel][i, j, k, e_ext, iθ, iϕ]
                        ip_g_global = ip2gip[ip_g_local]
                        # Avoid duplicates from shared element nodes
                        if ip_g_local ∉ local_ipg_sets[ip]
                            push!(local_ipg_sets[ip],  ip_g_local)
                            push!(global_ipg_sets[ip], ip_g_global)
                        end
                    end
                end
            end
        end
    else
        # Non-adaptive: ip_g_local = (ip-1)*npoin_ang + ip_ext
        for ip = 1:npoin
            for ip_ext = 1:npoin_ang
                ip_g_local  = (ip - 1) * npoin_ang + ip_ext
                ip_g_global = ip2gip[ip_g_local]
                push!(local_ipg_sets[ip],  ip_g_local)
                push!(global_ipg_sets[ip], ip_g_global)
            end
        end
    end

    # Sort by local index for consistent block extraction
    for ip = 1:npoin
        if !isempty(local_ipg_sets[ip])
            sort_idx = sortperm(local_ipg_sets[ip])
            local_ipg_sets[ip]  = local_ipg_sets[ip][sort_idx]
            global_ipg_sets[ip] = global_ipg_sets[ip][sort_idx]
        end
    end

    # ── Step 2: extract and invert each spatial node's angular block ──────────
    block_inv       = Vector{Union{Matrix{Float64},Nothing}}(nothing, npoin)
    block_sizes     = Int[]

    t_build = @elapsed begin
        for ip = 1:npoin
            local_ipg = local_ipg_sets[ip]
            isempty(local_ipg) && continue

            n_ang = length(local_ipg)
            push!(block_sizes, n_ang)

            # Extract dense block from sparse A_local using local indices
            block = zeros(Float64, n_ang, n_ang)
            for (bi, row) in enumerate(local_ipg)
                for (bj, col) in enumerate(local_ipg)
                    v = A_local[row, col]
                    v != 0.0 && (block[bi, bj] = v)
                end
            end

            try
                block_inv[ip] = inv(block)
            catch e
                @warn "Rank $rank: block inversion failed at ip=$ip " *
                      "(size $n_ang×$n_ang): $e. Using diagonal fallback."
                block_inv[ip] = diagm(1.0 ./ max.(abs.(diag(block)), 1e-14))
            end
        end
    end

    if rank == 0
        @info "Adaptive Block Jacobi preconditioner:"
        @info "  Local blocks:     $(count(!isnothing, block_inv))"
        @info "  Block size range: $(extrema(block_sizes))"
        @info "  Block size mean:  $(round(mean(block_sizes), digits=1))"
        @info "  Build time:       $(round(t_build, digits=3))s"
    end

    # ── Step 3: apply function using global indices ───────────────────────────
    # The Krylov solver operates on the global vector x of length gnpoin.
    # For each spatial node ip on this rank, apply the block inverse to
    # the subvector of x at the global combined indices for that node.
    # Nodes not owned by this rank are left as zero — the Allreduce in
    # the Krylov solver accumulates contributions across ranks.

    function apply_block_jacobi!(y, x)
        fill!(y, 0.0)
        for ip = 1:npoin
            isnothing(block_inv[ip]) && continue
            global_ipg = global_ipg_sets[ip]
            isempty(global_ipg) && continue
            @inbounds y[global_ipg] .= block_inv[ip] * x[global_ipg]
        end
        return y
    end

    return LinearOperator(Float64, gnpoin, gnpoin, false, false,
                          apply_block_jacobi!)
end

function solve_parallel_lsqr(ip2gip, gip2owner, A_local, b, gnpoin, npoin; npoin_g=0, g_ip2gip=Int[], g_gip2ip=Int[], tol::Float64 = 1e-7)
   
     # Create parallel linear operator (better than AbstractMatrix)
    A_parallel = create_parallel_linear_operator(A_local, ip2gip, gip2owner, npoin, gnpoin, npoin_g, g_ip2gip, g_gip2ip)

    # Gather the full RHS vector with proper Allgatherv
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    nprocs = MPI.Comm_size(MPI.COMM_WORLD)

    #build contributions to global RHS
    b_global = zeros(Float64,gnpoin)
    for ip=1:npoin
        gip = ip2gip[ip]
        b_global[gip] = b[ip]
    end
    #Use allreduce to build global rhs vector
    MPI.Allreduce!(b_global, +, MPI.COMM_WORLD)
    # Solve using LSQR
    maxiter = gnpoin
    x, _ = Krylov.lsqr(A_parallel, b_global;
                   atol = tol,
                   rtol = tol,
                   itmax = maxiter,
                   verbose = (rank == 0) ? 1 : 0)  # Only print on rank 0

    #return to local indexing
    x_local = zeros(Float64,npoin_g)
    for ip =1:npoin
        gip = ip2gip[ip]
        x_local[ip] = x[gip]
    end

    return x_local
end

# ── GMRES solver ──────────────────────────────────────────────────────────────
"""
    solve_parallel_gmres(ip2gip, gip2owner, A_local, b, gnpoin, npoin;
                         precond=:jacobi, restart=50, tol=1e-6,
                         npoin_g, g_ip2gip, g_gip2ip,
                         npoin_ang=nothing)

Solve the RT system using restarted GMRES with optional preconditioning.

# Preconditioner options
- `:none`        — no preconditioning (identity)
- `:jacobi`      — diagonal scaling. Recommended for LW (diagonally dominant).
- `:block_jacobi`— angular block inverse. Recommended for SW diffuse
                   (near-conservative scattering). Requires `npoin_ang`.

# Restart
`restart=50` is a good default. Increase to 100-200 if convergence stalls.
Memory cost scales as O(restart × gnpoin).

# Tolerance
- LW:  tol=1e-5 gives heating rate accuracy ~0.01 K/day
- SW:  tol=1e-4 is usually sufficient (direct beam dominates)
"""
function solve_parallel_gmres(ip2gip, gip2owner, A_local, b, gnpoin, npoin, x_prev;
                               precond   :: Symbol  = :jacobi,
                               restart   :: Int     = 50,
                               tol       :: Float64 = 1e-6,
                               itmax     :: Int     = 10000,
                               npoin_g   :: Int     = 0,
                               g_ip2gip  :: Vector{Int} = Int[],
                               g_gip2ip  = Int[],
                               # Angular mesh connectivity for block Jacobi
                               connijk_spa   = nothing,
                               extra_nelem   = nothing,
                               extra_nops    = nothing,
                               extra_connijk = nothing,
                               nelem         :: Int     = 0,
                               ngl           :: Int     = 0,
                               ladaptive     :: Bool    = false,
                               npoin_ang     = nothing,
                               npoin_space   :: Int     = 0,
                               comm      :: MPI.Comm = MPI.COMM_WORLD)

    rank = MPI.Comm_rank(comm)

    A_parallel = create_parallel_linear_operator(
        A_local, ip2gip, gip2owner, npoin, gnpoin,
        npoin_g, g_ip2gip, g_gip2ip, comm)

    b_global = zeros(Float64, gnpoin)
    for ip = 1:npoin
        b_global[ip2gip[ip]] += b[ip]
    end
    MPI.Allreduce!(b_global, +, comm)

    # Assemble global initial guess from previous timestep solution
    x0 = zeros(Float64, gnpoin)
    if !isempty(x_prev)
        for ip = 1:npoin
            x0[ip2gip[ip]] = x_prev[ip]
        end
        MPI.Allreduce!(x0, +, comm)

        # Check if warm start is useful — compute initial residual
        # If ||b - A x0|| / ||b|| is already below tol, skip solve entirely
        r0     = b_global - A_parallel * x0
        r0_rel = norm(r0) / max(norm(b_global), 1e-30)
        if rank == 0
            @info "Warm start: initial residual = $(round(r0_rel, sigdigits=3))"
        end
        if r0_rel < tol
            if rank == 0
                @info "Warm start: initial guess already converged, skipping solve."
            end
            x_local = zeros(Float64, max(npoin_g, npoin))
            for ip = 1:npoin
                x_local[ip] = x0[ip2gip[ip]]
            end
            return x_local
        end
    end
    
    # Build preconditioner
    N = if precond == :jacobi
        build_jacobi_preconditioner(A_local, ip2gip, gnpoin, comm)
    elseif precond == :block_jacobi
        isnothing(npoin_ang) && error("block_jacobi requires npoin_ang")
        n_spatial_local = npoin ÷ npoin_ang
        build_block_jacobi_preconditioner(
            A_local, ip2gip, gnpoin, n_spatial_local, npoin_ang, comm)
    elseif precond == :ilu
        if rank == 0
            @info "Building ILU preconditioner on each rank (local block)"
        end
        ilu_factor = IncompleteLU.ilu(A_local, τ = 0.01)
        # Apply locally: scatter global x to local, apply ILU, scatter back
        LinearOperator(Float64, gnpoin, gnpoin, true, true,
            (y, x) -> begin
                x_local = [x[ip2gip[ip]] for ip = 1:npoin]
                y_local = similar(x_local)
                LinearAlgebra.ldiv!(y_local, ilu_factor, x_local)
                fill!(y, 0.0)
                for ip = 1:npoin
                    y[ip2gip[ip]] += y_local[ip]
                end
                MPI.Allreduce!(y, +, comm)
            end)
    elseif precond == :lu
        if rank == 0
            @info "Building full LU preconditioner on each rank (local block)"
        end
        lu_factor = lu(A_local)
        # Apply locally: scatter global x to local, apply LU, scatter back
        LinearOperator(Float64, gnpoin, gnpoin, true, true,
            (y, x) -> begin
                x_local = [x[ip2gip[ip]] for ip = 1:npoin]
                y_local = similar(x_local)
                LinearAlgebra.ldiv!(y_local, lu_factor, x_local)
                fill!(y, 0.0)
                for ip = 1:npoin
                    y[ip2gip[ip]] += y_local[ip]
                end
                MPI.Allreduce!(y, +, comm)
            end)
    elseif precond == :klu
        # KLU multifrontal sparse direct solver as preconditioner
        # KLU uses a fill-reducing ordering and is typically faster than UMFPACK
        # for circuit-simulation and similar sparse structures
        if rank == 0
            @info "Building KLU multifrontal preconditioner on each rank (local block)"
        end
        klu_factor = klu(A_local)
        LinearOperator(Float64, gnpoin, gnpoin, true, true,
            (y, x) -> begin
                x_local = [x[ip2gip[ip]] for ip = 1:npoin]
                y_local = similar(x_local)
                LinearAlgebra.ldiv!(y_local, klu_factor, x_local)
                fill!(y, 0.0)
                for ip = 1:npoin
                    y[ip2gip[ip]] += y_local[ip]
                end
                MPI.Allreduce!(y, +, comm)
            end)
    else
        LinearOperator(Float64, gnpoin, gnpoin, true, true,
                       (y, x) -> copyto!(y, x))
    end

    t_solve = @elapsed begin
        if isempty(x_prev)
            x, stats = Krylov.gmres(A_parallel, b_global;
                                     N=N,
                                     memory  = restart,
                                     restart = true,
                                     atol    = tol,
                                     rtol    = tol,
                                     itmax   = itmax,
                                     verbose = (rank == 0) ? 1 : 0)
        else
            x, stats = Krylov.gmres(A_parallel, b_global, x0;
                                     N=N,
                                     memory  = restart,
                                     restart = true,
                                     atol    = tol,
                                     rtol    = tol,
                                     itmax   = itmax,
                                     verbose = (rank == 0) ? 1 : 0)
        end
    end

    if rank == 0
        @info "GMRES: $(stats.niter) iterations, " *
              "$(round(t_solve, digits=2))s, " *
              "residual=$(stats.residuals), " *
              "converged=$(stats.solved)"
        if !stats.solved
            @warn "GMRES did not converge in $itmax iterations. " *
                  "Consider increasing memory or loosening tol."
        end
    end

    x_local = zeros(Float64, max(npoin_g, npoin))
    for ip = 1:npoin
        x_local[ip] = x[ip2gip[ip]]
    end
    return x_local
end


# ── BiCGSTAB solver ───────────────────────────────────────────────────────────
"""
    solve_parallel_bicgstab(ip2gip, gip2owner, A_local, b, gnpoin, npoin;
                             precond=:jacobi, tol=1e-6, ...)

Solve using BiCGSTAB with optional preconditioning.

BiCGSTAB uses 2 matvecs per iteration (vs GMRES restart/iter amortized).
It has lower memory cost than GMRES (O(n) vs O(restart × n)) but can
stall on highly non-symmetric systems.

Recommended when:
- Memory is constrained and restart GMRES is expensive
- The matrix is only mildly non-symmetric
- As a cross-check against GMRES

Not recommended for the SW diffuse problem with near-conservative
scattering — GMRES + Block Jacobi will typically outperform it there.
"""
function solve_parallel_bicgstab(ip2gip, gip2owner, A_local, b, gnpoin, npoin;
                                  precond   :: Symbol  = :jacobi,
                                  tol       :: Float64 = 1e-6,
                                  itmax     :: Int     = 500,
                                  npoin_g   :: Int     = 0,
                                  g_ip2gip  :: Vector{Int} = Int[],
                                  g_gip2ip  = Int[],
                                  npoin_ang :: Union{Int,Nothing} = nothing,
                                  comm      :: MPI.Comm = MPI.COMM_WORLD)

    rank   = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    A_parallel = create_parallel_linear_operator(
        A_local, ip2gip, gip2owner, npoin, gnpoin,
        npoin_g, g_ip2gip, g_gip2ip, comm)

    b_global = zeros(Float64, gnpoin)
    for ip = 1:npoin
        b_global[ip2gip[ip]] = b[ip]
    end
    MPI.Allreduce!(b_global, +, comm)

    N = if precond == :jacobi
        build_jacobi_preconditioner(A_local, ip2gip, gnpoin, comm)
    elseif precond == :block_jacobi
        isnothing(npoin_ang) && error("block_jacobi requires npoin_ang")
        n_spatial_local = npoin ÷ npoin_ang
        build_block_jacobi_preconditioner(
            A_local, ip2gip, gnpoin, n_spatial_local, npoin_ang, comm)
    else
        LinearOperator(Float64, gnpoin, gnpoin, true, true,
                       (y, x) -> copyto!(y, x))
    end

    if rank == 0
        @info "Starting parallel BiCGSTAB solve..."
        @info "  Global DOFs:    $gnpoin"
        @info "  Preconditioner: $precond"
        @info "  Tolerance:      $tol"
        @info "  Max iterations: $itmax"
        @info "  Processes:      $nprocs"
    end

    t_solve = @elapsed begin
        x, stats = Krylov.bicgstab(A_parallel, b_global;
                                    N      = N,
                                    atol   = tol,
                                    rtol   = tol,
                                    itmax  = itmax,
                                    verbose = (rank == 0) ? 1 : 0)
    end

    if rank == 0
        @info "BiCGSTAB converged: $(stats.solved) in $(stats.niter) iterations, " *
              "$(round(t_solve, digits=3))s, " *
              "residual=$(round(stats.residuals[end], sigdigits=3))"
        if !stats.solved
            @warn "BiCGSTAB did not converge. Switch to GMRES with larger restart."
        end
    end

    x_local = zeros(Float64, npoin_g > 0 ? npoin_g : npoin)
    for ip = 1:npoin
        x_local[ip] = x[ip2gip[ip]]
    end

    return x_local
end

function diagnose_rt_matrix(A_local, ip2gip, gnpoin, label, comm=MPI.COMM_WORLD)
    rank = MPI.Comm_rank(comm)

    n = size(A_local, 1)
    d = abs.(diag(A_local))
    
    # Row sum of off-diagonal absolute values
    offdiag_sum = vec(sum(abs.(A_local), dims=2)) .- d
    
    # Diagonal dominance ratio per row
    dd = d ./ max.(offdiag_sum, 1e-30)
    
    # Fraction of rows that are diagonally dominant
    frac_dd = mean(dd .> 1.0)
    
    # Distribution of dd ratio
    pct = [1, 10, 25, 50, 75, 90, 99]
    dd_sorted = sort(dd)
    n_pct = length(dd_sorted)
    dd_pct = [dd_sorted[clamp(round(Int, p/100 * n_pct), 1, n_pct)] for p in pct]

    if rank == 0
        @info "$label matrix diagonal dominance:"
        @info "  Rows:                    $n"
        @info "  Diag dominant fraction:  $(round(frac_dd*100, digits=1))%"
        @info "  DD ratio percentiles:"
        for (p, v) in zip(pct, dd_pct)
            @info "    $(lpad(p,3))th: $(round(v, sigdigits=3))"
        end
        @info "  Diagonal sign: $(all(diag(A_local) .> 0) ? "all positive ✓" : "mixed ✗")"
        @info "  Diag range:    $(round.(extrema(d), sigdigits=4))"
        @info "  Offdiag sum range: $(round.(extrema(offdiag_sum), sigdigits=4))"
    end
end

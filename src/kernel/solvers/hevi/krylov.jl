#=============================================================================
 krylov.jl -- a right-preconditioned restarted GMRES over MPI-distributed
 spectral-element fields.

 WHY THIS AND NOT Krylov.jl / IterativeSolvers.jl
 ------------------------------------------------
 Both are in the manifest and both want an `AbstractVector` with a working
 `dot`. The thing being solved for here is not one:

   * it is stored as an `npoin x nimp` matrix, not a vector;
   * it is REDUNDANT. A continuous-Galerkin field keeps a private copy of
     every node that sits on a rank interface, so `dot(x, y)` summed
     rank-locally and reduced counts those nodes once per holder. The result
     is a bilinear form, but it is not the inner product of the space the
     operator is self-consistent in, and Arnoldi built on it produces a basis
     that is not orthogonal in any norm the residual is measured in.

 Wrapping the field in a custom `AbstractVector` to satisfy those packages
 means implementing the redundancy-aware `dot` anyway, plus `similar`,
 broadcasting and a `mul!`, and then debugging through two layers of
 abstraction when a rank disagrees. The Arnoldi loop below is ~80 lines and
 says exactly what it does.

 THE INNER PRODUCT
 -----------------
 Weight each local node by the reciprocal of the number of ranks holding it:

     <x, y> = Allreduce( Σ_ip Σ_q  w[ip] · x[ip,q] · y[ip,q] )      w = 1/mult

 and the sum over holders telescopes back to one contribution per GLOBAL node.
 The multiplicity is measured, not derived: `assemble_mpi!` leaves every holder
 of a node with the sum over all holders, so assembling a field of ones returns
 the multiplicity in the field's own indexing -- including local duplicates
 from periodicity, which a `gip2owner .== rank` mask would get wrong.

 Everything the iteration touches stays single-valued across ranks. The
 operator ends in DSS + `assemble_mpi!` + the inverse mass matrix, and the
 preconditioner ends in `column_scatter!`; both hand back a field on which
 every holder of a node agrees. So the weighted form really is an inner
 product on the Krylov subspace, not merely on paper.

 CLASSICAL GRAM-SCHMIDT, TWICE, AND WHY NOT MODIFIED
 ---------------------------------------------------
 Modified Gram-Schmidt is the textbook choice and costs `j` separate
 `Allreduce`s at Arnoldi step `j` -- each one a full latency round trip, and
 the reason a Krylov method that looks fine on one rank falls over at scale.
 Classical Gram-Schmidt batches all `j` projections into ONE reduction, and
 run twice ("CGS2") it is as stable as MGS to working precision. Three
 reductions per step, independent of `j`.
=============================================================================#

"""
    DistributedInner(w, comm)

The multiplicity-weighted inner product of a redundantly-stored CG field.
`w[ip]` is `1/(number of ranks holding node ip)`.
"""
struct DistributedInner
    w::Vector{Float64}
    comm::MPI.Comm
    parallel::Bool
end

"""
    build_distributed_inner(dss_cache, npoin, comm) -> DistributedInner

Measure the node multiplicities and build the inner product from them.

`dss_cache === nothing` means there is no assembler (a serial 1D case, which
this scheme refuses anyway) and every node is held once.
"""
function build_distributed_inner(dss_cache, npoin::Int, comm::MPI.Comm)
    w = ones(Float64, npoin)
    if dss_cache !== nothing
        m = ones(Float64, npoin, 1)
        assemble_mpi!(m, dss_cache)
        @inbounds for ip = 1:npoin
            # A multiplicity below one cannot happen; if the assembler ever
            # produces one, dividing by it would silently blow the norm up
            # instead of failing, so clamp and let the setup self-check speak.
            w[ip] = 1.0 / max(m[ip, 1], 1.0)
        end
    end
    return DistributedInner(w, comm, MPI.Comm_size(comm) > 1)
end

"""
    ddot(di, X, Y) -> Float64

`<X, Y>` over the global mesh.
"""
function ddot(di::DistributedInner, X::AbstractMatrix, Y::AbstractMatrix)
    w = di.w
    s = 0.0
    @inbounds for q = 1:size(X, 2), ip = 1:size(X, 1)
        s += w[ip] * X[ip, q] * Y[ip, q]
    end
    return di.parallel ? MPI.Allreduce(s, MPI.SUM, di.comm) : s
end

"""
    ddots!(out, di, X, V, j)

`out[i] = <X, V[i]>` for `i = 1:j`, in ONE reduction rather than `j` of them.
This is the whole reason the Arnoldi loop below uses classical Gram-Schmidt.
"""
function ddots!(out::Vector{Float64}, di::DistributedInner,
                X::AbstractMatrix, V::Vector{Matrix{Float64}}, j::Int)
    w = di.w
    n, m = size(X)
    fill!(out, 0.0)
    @inbounds for i = 1:j
        Vi = V[i]
        s = 0.0
        for q = 1:m, ip = 1:n
            s += w[ip] * X[ip, q] * Vi[ip, q]
        end
        out[i] = s
    end
    # The whole buffer is reduced, not `view(out, 1:j)`: the tail is zeroed
    # above, the extra bytes are a handful of doubles, and it keeps this off
    # the question of which SubArrays MPI.jl will accept as a buffer.
    di.parallel && MPI.Allreduce!(out, MPI.SUM, di.comm)
    return out
end

dnorm(di::DistributedInner, X::AbstractMatrix) = sqrt(max(ddot(di, X, X), 0.0))

"""
    GMRESWorkspace

Everything the stage solve allocates, once, at setup.

`V` holds `m+1` Arnoldi basis vectors, each `npoin x nimp`. The preconditioned
basis is deliberately NOT stored: recovering the update at the end of a cycle
costs one extra preconditioner application instead of `m` more full fields,
which on a production mesh is the difference between tens and hundreds of
megabytes per rank.
"""
mutable struct GMRESWorkspace
    m::Int                            # restart length
    maxiter::Int                      # total Arnoldi steps across restarts
    rtol::Float64
    atol::Float64
    V::Vector{Matrix{Float64}}
    W::Matrix{Float64}                # matvec / orthogonalisation scratch
    Z::Matrix{Float64}                # preconditioner scratch
    C::Matrix{Float64}                # correction accumulator
    H::Matrix{Float64}                # (m+1) x m Hessenberg, Givens-reduced
    cs::Vector{Float64}
    sn::Vector{Float64}
    g::Vector{Float64}
    y::Vector{Float64}
    hs::Vector{Float64}               # batched projection coefficients
    inner::DistributedInner
    # running totals, for the profile line
    nsolve::Int
    niter::Int
    nunconverged::Int
    last_iters::Int
    last_relres::Float64
    worst_relres::Float64
end

function GMRESWorkspace(npoin::Int, nimp::Int, inner::DistributedInner;
                        m::Int = 20, maxiter::Int = 200,
                        rtol::Float64 = 1.0e-8, atol::Float64 = 1.0e-30)
    m >= 1 || error("GMRES: the restart length must be at least 1; got $m.")
    maxiter >= m ||
        error("GMRES: :imex_maxiter ($maxiter) is below the restart length ($m), so ",
              "the iteration would be cut off mid-cycle and could never reach a ",
              "restart. Raise :imex_maxiter or lower :imex_restart.")
    z() = zeros(Float64, npoin, nimp)
    return GMRESWorkspace(m, maxiter, rtol, atol,
                          [z() for _ = 1:m+1], z(), z(), z(),
                          zeros(Float64, m + 1, m),
                          zeros(Float64, m), zeros(Float64, m), zeros(Float64, m + 1),
                          zeros(Float64, m), zeros(Float64, m),
                          inner, 0, 0, 0, 0, 0.0, 0.0)
end

"""
    gmres_solve!(X, B, ws, matvec!, precon!) -> (iters, relres, converged)

Solve `A X = B` for `X`, with `matvec!(out, in)` applying `A` and
`precon!(Z)` applying `M^-1` in place.

Right-preconditioned: the iteration runs on `A M^-1 y = B` and reports the
residual of the ORIGINAL system, so the tolerance means what it says. Left
preconditioning would minimise `‖M^-1 (B - A X)‖` instead, and with a
preconditioner as strong as a column solve that norm is a poor proxy for the
one the time integrator cares about.

`X` is overwritten and used as the initial guess only in the sense that it
starts at zero: the stage solve's right-hand side is a deviation field that is
already small, and a zero guess makes the reported relative residual
unambiguous.
"""
function gmres_solve!(X::AbstractMatrix, B::AbstractMatrix, ws::GMRESWorkspace,
                      matvec!, precon!)

    di = ws.inner
    m  = ws.m
    V, W, Z, C = ws.V, ws.W, ws.Z, ws.C
    H, cs, sn, g, y, hs = ws.H, ws.cs, ws.sn, ws.g, ws.y, ws.hs

    fill!(X, 0.0)
    bnorm = dnorm(di, B)
    ws.nsolve += 1

    # An exactly-zero right-hand side is not a degenerate case here, it is the
    # NORMAL one at t = 0 on a state at rest: the operator acts on u - qe and
    # the stage vector is qe. X = 0 is the exact answer and 0/0 is not.
    if !(bnorm > 0.0)
        ws.last_iters = 0; ws.last_relres = 0.0
        return (0, 0.0, true)
    end

    tol   = max(ws.rtol * bnorm, ws.atol)
    total = 0
    resid = bnorm

    while true
        # r = B - A X. On the first cycle X is zero and the matvec is skipped;
        # on a restart it is not, and recomputing the residual explicitly
        # (rather than carrying the Arnoldi estimate across the restart) is
        # what keeps the reported number the true residual.
        if total == 0
            copyto!(V[1], B)
        else
            matvec!(W, X)
            @inbounds @. V[1] = B - W
        end
        β = dnorm(di, V[1])
        if β <= tol || total >= ws.maxiter
            resid = β
            break
        end

        @inbounds @. V[1] = V[1] * (1.0 / β)
        fill!(g, 0.0); g[1] = β
        fill!(H, 0.0)
        j = 0
        breakdown = false

        for jj = 1:m
            j = jj
            copyto!(Z, V[jj])
            precon!(Z)                 # z = M^-1 v_j
            matvec!(W, Z)              # w = A z

            # --- CGS2: two batched projection passes, three reductions ------
            ddots!(hs, di, W, V, jj)
            @inbounds for i = 1:jj
                H[i, jj] = hs[i]
                axpy_field!(W, -hs[i], V[i])
            end
            ddots!(hs, di, W, V, jj)
            @inbounds for i = 1:jj
                H[i, jj] += hs[i]
                axpy_field!(W, -hs[i], V[i])
            end

            hnext = dnorm(di, W)
            H[jj+1, jj] = hnext

            # --- Givens: apply the previous rotations, then the new one -----
            @inbounds for i = 1:jj-1
                t = cs[i] * H[i, jj] + sn[i] * H[i+1, jj]
                H[i+1, jj] = -sn[i] * H[i, jj] + cs[i] * H[i+1, jj]
                H[i, jj]   = t
            end
            ρ = hypot(H[jj, jj], H[jj+1, jj])
            if ρ == 0.0
                cs[jj] = 1.0; sn[jj] = 0.0
            else
                cs[jj] = H[jj, jj] / ρ; sn[jj] = H[jj+1, jj] / ρ
            end
            H[jj, jj]   = ρ
            H[jj+1, jj] = 0.0
            g[jj+1] = -sn[jj] * g[jj]
            g[jj]   =  cs[jj] * g[jj]

            total += 1
            resid = abs(g[jj+1])

            # A zero subdiagonal is a LUCKY breakdown: the Krylov space is
            # invariant and the current iterate is exact. Normalising by it
            # would produce NaNs instead.
            if hnext <= 1.0e-14 * max(β, 1.0)
                breakdown = true
                break
            end
            @inbounds @. V[jj+1] = W * (1.0 / hnext)

            (resid <= tol || total >= ws.maxiter) && break
        end

        # --- back-substitute and form the update -----------------------------
        @inbounds for i = j:-1:1
            t = g[i]
            for k = i+1:j
                t -= H[i, k] * y[k]
            end
            y[i] = H[i, i] == 0.0 ? 0.0 : t / H[i, i]
        end
        fill!(C, 0.0)
        @inbounds for i = 1:j
            y[i] != 0.0 && axpy_field!(C, y[i], V[i])
        end
        precon!(C)                     # right preconditioning: X += M^-1 (V y)
        @inbounds @. X = X + C

        (resid <= tol || total >= ws.maxiter || breakdown) && break
    end

    ws.niter += total
    ws.last_iters = total
    rel = resid / bnorm
    ws.last_relres = rel
    ws.worst_relres = max(ws.worst_relres, rel)
    converged = resid <= tol
    converged || (ws.nunconverged += 1)
    return (total, rel, converged)
end

"""
    axpy_field!(Y, a, X)

`Y .+= a .* X`, written out because `@.` over two `Matrix{Float64}`s in this
loop is on the stage-solve hot path and this makes the bounds check and the
broadcast machinery go away.
"""
@inline function axpy_field!(Y::AbstractMatrix, a::Float64, X::AbstractMatrix)
    @inbounds for q = 1:size(Y, 2), ip = 1:size(Y, 1)
        Y[ip, q] += a * X[ip, q]
    end
    return Y
end

#=============================================================================
 tools/poisson3d_benchmark/poisson3d.jl

 Solver comparison for the fully periodic Poisson problem in d = 2 or 3
 dimensions,
     -∇²u = f   on [0,2π]^d,  u periodic in every direction,
 at sizes where the asymptotic cost of the solvers shows (millions of
 unknowns in 3D). It is the d-dimensional generalisation of the 2D benchmark
 in tools/periodic_poisson_benchmark, with the same exact solution family,
 the same seven solvers and the same timing protocol.

 WHY A SEPARATE ASSEMBLY. Jexpresso's linear-solve path (Laplacian assembly,
 periodic reduction, static condensation) is 2D only, and stores every
 element matrix densely. On a Cartesian mesh of affine elements with GLL
 collocation (Jexpresso's default "inexact" quadrature, Q = N) the SEM
 stiffness matrix is EXACTLY a sum of Kronecker products of assembled 1D
 matrices,
     3D:  K = Mz⊗My⊗Kx + Mz⊗Ky⊗Mx + Kz⊗My⊗Mx,     M = Mz⊗My⊗Mx (lumped),
 because ψ_i(ξ_k) = δ_ik makes the quadrature separate by direction. This is
 the same discrete operator Jexpresso assembles (and the one of the CEED
 bake-off problem BP5), built in seconds for millions of unknowns.
 verify_2d.jl checks that in 2D it reproduces Jexpresso's periodic SEM
 solution (same errors to round-off).

 The 1D ingredients come from Jexpresso itself: LGL nodes and weights
 (basis_structs_ξ_ω!, Kopriva's algorithm), the Lagrange basis derivatives
 (LagrangeInterpolatingPolynomials_classic); the AMG and skeleton solves are
 Jexpresso's (jx_amg_setup / jx_amg_solve, el_skeleton_solve); the spectral
 solvers are Jexpresso's (_collocation_axis with FourierDerivativeMatrix,
 FFTPoissonSolver).

 SOLVERS (all return the zero-mean solution; the singular SEM systems pin
 their first unknown and shift to zero M-weighted mean afterwards):
   :sem         sparse Cholesky (CHOLMOD) of the full SEM system
   :sem_amg     smoothed-aggregation AMG + CG, full system
   :sem_jacobi  Jacobi (diagonal) preconditioned CG, full system (= BP5)
   :sc_direct   element-level static condensation, skeleton by Cholesky
   :sc_amg      element-level static condensation, skeleton by AMG + CG
   :ps          pseudo-spectral Fourier collocation (matrix diagonalisation,
                O(N_g^(d+1)))
   :fft         FFT (FFTW, planned with FFTW.ESTIMATE)
=============================================================================#
module P3D

using Jexpresso, SparseArrays, LinearAlgebra, Printf
const JX = Jexpresso
const KA = Jexpresso.KernelAbstractions
# METIS (nested-dissection ordering for the sparse Cholesky factorisations):
# installed with Jexpresso as a dependency of Gridap's partitioning, loaded by
# package id since it is not a direct dependency of the project.
const Metis = Base.require(Base.PkgId(Base.UUID("2679e427-3c69-5b7f-982b-ece356f1e94b"), "Metis"))
using Base.Threads: @threads, nthreads

export run_config, SOLVERS

const SOLVERS = (:sem, :sem_amg, :sem_jacobi, :sc_direct, :sc_amg, :ps, :fft)
const SEM_SOLVERS = (:sem, :sem_amg, :sem_jacobi, :sc_direct, :sc_amg)
const LBOX = 2π

# ---------------------------------------------------------------------------
# Exact solution: product of periodic Poisson kernels (d-dimensional version
# of problems/Elliptic/poisson_periodic_sem/user_source.jl)
#   u = A ( Π_k p(x_k) - (c²-1)^(-d/2) ),  p(s) = 1/(c - cos s),
#   c = (r + 1/r)/2; A scales the peak to u(0,…,0) = 1.
# Fourier coefficients decay like r^(|k_1|+…+|k_d|): not band-limited. r sets
# how hard the problem is: the peak has width ~√(2(c-1)) in each direction.
# The 2D benchmark uses r = 0.8; in 3D the product of three such peaks is not
# resolved at feasible sizes, so the 3D default is r = 0.5 (c = 1.25).
# ---------------------------------------------------------------------------
const R_DEFAULT = Dict(2 => 0.8, 3 => 0.5)
_c(r) = (r + 1 / r) / 2
_p(s, c)   = 1 / (c - cos(s))
_pss(s, c) = (q = c - cos(s); -cos(s) / q^2 + 2 * sin(s)^2 / q^3)
amp(d, c) = 1 / ((c - 1)^(-d) - (c^2 - 1)^(-d / 2))
umean(d, c) = (c^2 - 1)^(-d / 2)

# u and f on the tensor grid of the 1D coordinates x (the same in every
# direction), as d-dimensional arrays with x_1 fastest
function exact_on_grid(x::Vector{Float64}, d::Int, r::Float64)
    c = _c(r); A = amp(d, c); p = _p.(x, c); q = _pss.(x, c); n = length(x)
    if d == 2
        U = A .* (p .* p' .- umean(2, c))
        F = -A .* (q .* p' .+ p .* q')
    else
        pz = reshape(p, 1, 1, n); qz = reshape(q, 1, 1, n)
        U = A .* (p .* p' .* pz .- umean(3, c))
        F = -A .* (q .* p' .* pz .+ p .* q' .* pz .+ p .* p' .* qz)
    end
    return U, F
end

# ---------------------------------------------------------------------------
# SEM operators
# ---------------------------------------------------------------------------
"""
    sem1d(N, ne) -> (; K1, m1, x, Ke, me)

Periodic 1D SEM operators on `ne` elements of order `N` over [0, 2π):
assembled stiffness K1 (exactly symmetric), lumped mass m1, node coordinates
x, and the element stiffness Ke / element mass diagonal me.
"""
function sem1d(N::Int, ne::Int)
    JX.MPI.Initialized() || JX.MPI.Init()           # Jexpresso's LGL routine asks for the MPI rank
    lgl = JX.basis_structs_ξ_ω!(JX.LGL(), N, KA.CPU())
    ξ = collect(Float64, lgl.ξ); ω = collect(Float64, lgl.ω)
    _, dψ = JX.LagrangeInterpolatingPolynomials_classic(ξ, ξ, Float64, KA.CPU())
    dψ = Matrix{Float64}(dψ)                       # dψ[i,k] = ψ_i'(ξ_k)
    h = LBOX / ne; J = h / 2
    Ke = [sum(ω[k] * dψ[i, k] * dψ[l, k] for k in 1:N+1) / J for i in 1:N+1, l in 1:N+1]
    Ke = (Ke + Ke') ./ 2                           # exactly symmetric
    me = J .* ω
    n1 = ne * N
    I = Int[]; Jv = Int[]; V = Float64[]
    m1 = zeros(n1); x = zeros(n1)
    for e in 0:ne-1, i in 0:N, l in 0:N
        push!(I, mod(e * N + i, n1) + 1); push!(Jv, mod(e * N + l, n1) + 1); push!(V, Ke[i+1, l+1])
    end
    for e in 0:ne-1, i in 0:N
        m1[mod(e * N + i, n1) + 1] += me[i+1]
        i < N && (x[e * N + i + 1] = e * h + (ξ[i+1] + 1) * J)
    end
    return (; K1 = sparse(I, Jv, V, n1, n1), m1, x, Ke, me, N, ne, n1)
end

"""
    sem_system(N, ne, d) -> (; K, m, b, fmean, x, U, o1)

The d-dimensional periodic SEM system K u = b, b = M f projected to zero sum,
on ne^d elements of order N; U is the exact solution at the nodes.
"""
function sem_system(N::Int, ne::Int, d::Int, r::Float64)
    o1 = sem1d(N, ne)
    Md = spdiagm(o1.m1); K1 = o1.K1
    K = d == 2 ? kron(Md, K1) + kron(K1, Md) :
                 kron(Md, kron(Md, K1)) + kron(Md, kron(K1, Md)) + kron(K1, kron(Md, Md))
    m = d == 2 ? kron(o1.m1, o1.m1) : kron(o1.m1, kron(o1.m1, o1.m1))
    U, F = exact_on_grid(o1.x, d, r)
    b = m .* vec(F)
    fmean = sum(b) / sum(m)
    b .-= fmean .* m
    return (; K, m, b, fmean, x = o1.x, U = vec(U), o1)
end

gauge!(u, m) = (u .-= dot(m, u) / sum(m); u)

function sem_error(u, U, m)
    e = u .- U
    return (linf = maximum(abs, e), l2rel = sqrt(dot(m, e .^ 2) / dot(m, U .^ 2)))
end

# ---------------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------------
macro phase(T, key, ex)
    quote
        local t0 = time_ns()
        local v = $(esc(ex))
        $(esc(T))[$(esc(key))] = get($(esc(T)), $(esc(key)), 0.0) + (time_ns() - t0) / 1e9
        v
    end
end

# ---------------------------------------------------------------------------
# Full-system SEM solvers
# ---------------------------------------------------------------------------
# fill-reducing ordering for CHOLMOD: METIS nested dissection (optimal
# O(n^(4/3)) fill in 3D; CHOLMOD's own default here is AMD, whose fill grows
# markedly faster in 3D), or CHOLMOD's default (:amd)
fill_ordering(A, ordering) = ordering === :metis ? Int.(first(Metis.permutation(A))) : nothing
_chol(A, p) = p === nothing ? cholesky(Symmetric(A)) : cholesky(Symmetric(A); perm = p)

function solve_sem!(T, info, sys, solver; rtol = 1e-12, itmax = 100_000, amg_method = :sa, ordering = :metis)
    Kp = @phase T :setup sys.K[2:end, 2:end]
    bp = sys.b[2:end]
    u = zeros(length(sys.b))
    if solver === :sem
        F = @phase T :setup _chol(Kp, fill_ordering(Kp, ordering))
        info[:factor_nnz] = nnz(F)
        u[2:end] = @phase T :solve F \ bp
    elseif solver === :sem_amg
        S = @phase T :setup JX.jx_amg_setup(Kp; method = amg_method)
        u[2:end] = @phase T :solve JX.jx_amg_solve(S, bp; rtol = rtol, itmax = itmax)
        info[:iters] = JX.JX_AMG_STATS[].iters
        info[:amg_levels] = JX.JX_AMG_STATS[].levels
    elseif solver === :sem_jacobi
        Dinv = @phase T :setup Diagonal(1 ./ diag(Kp))
        x, st = @phase T :solve JX.Krylov.cg(Kp, bp; M = Dinv, rtol = rtol, atol = 0.0, itmax = itmax)
        u[2:end] = x
        info[:iters] = st.niter
        st.solved || @warn "Jacobi-CG did not converge in $itmax iterations"
    end
    return gauge!(u, sys.m)
end

# ---------------------------------------------------------------------------
# Element-level static condensation (d-dimensional)
#
# For every element e, with its (N-1)^d interior nodes o and its boundary
# nodes b (the element's share of the skeleton), the element stiffness is
# split into A_oo, A_ob, A_bo, A_bb and condensed:
#     S^e = A_bb - A_bo A_oo⁻¹ A_ob,     t^e = A_oo⁻¹ f_o,
# the skeleton system  B u_s = b_s - Σ_e A_bo t^e  with B = Σ_e S^e assembled,
# then the interiors are recovered element by element:
#     u_o = A_oo⁻¹ (f_o - A_ob u_b).
# This is the static condensation of elementLearning_Axb! (same equations),
# done at element level: an interior node belongs to one element, so the
# interior rows of the global matrix are that element's rows. A_oo⁻¹ is
# applied through a dense Cholesky factorisation of A_oo, computed per element
# (as on a general mesh, even though on this uniform mesh all element matrices
# are equal), once for the condensation and again for the recovery, as
# elementLearning_Axb! does. The skeleton system goes to Jexpresso's
# el_skeleton_solve (symmetrisation, then CHOLMOD Cholesky or AMG + CG), with
# one unknown pinned (singular: no Dirichlet boundary).
# ---------------------------------------------------------------------------
function element_matrix(o1, d)
    Ke = o1.Ke; Me = Diagonal(o1.me)
    return d == 2 ? kron(Me, Ke) + kron(Ke, Me) :
                    kron(Me, kron(Me, Ke)) + kron(Me, kron(Ke, Me)) + kron(Ke, kron(Me, Me))
end

function solve_sc!(T, info, sys, solver; rtol = 1e-12, amg_method = "sa", d::Int, ordering = :metis)
    o1 = sys.o1; N = o1.N; ne = o1.ne; n1 = o1.n1
    nl = N + 1
    Ael = @phase T :setup element_matrix(o1, d)
    # local nodes (a_1,…,a_d) ∈ 0:N, a_1 fastest; interior = all a_k in 1:N-1
    locs = d == 2 ? [(a, b) for b in 0:N for a in 0:N] : [(a, b, c) for c in 0:N for b in 0:N for a in 0:N]
    isint = [all(1 .<= t .<= N - 1) for t in locs]
    io = findall(isint); ib = findall(.!isint)
    Aoo = Symmetric(Matrix(Ael[io, io])); Aob = Matrix(Ael[io, ib]); Abo = Matrix(Ael[ib, io]); Abb = Matrix(Ael[ib, ib])
    no, nb = length(io), length(ib)
    # skeleton numbering: nodes with a global 1D index ≡ 0 (mod N) in some direction
    n = length(sys.b)
    spos = zeros(Int32, n); ns = 0
    on1 = [(i - 1) % N == 0 for i in 1:n1]         # 1D node on an element boundary
    @phase T :setup begin
        I = 0
        if d == 2
            for j in 1:n1, i in 1:n1
                I += 1
                (on1[i] || on1[j]) && (ns += 1; spos[I] = ns)
            end
        else
            for k in 1:n1, j in 1:n1, i in 1:n1
                I += 1
                (on1[i] || on1[j] || on1[k]) && (ns += 1; spos[I] = ns)
            end
        end
    end
    ne_tot = ne^d
    lin(t) = d == 2 ? (t[1] + n1 * t[2] + 1) : (t[1] + n1 * (t[2] + n1 * t[3]) + 1)
    function gnodes(e)       # global linear indices of the element's local nodes
        ec = d == 2 ? (e % ne, e ÷ ne) : (e % ne, (e ÷ ne) % ne, e ÷ (ne * ne))
        return [lin(ntuple(k -> mod(ec[k] * N + t[k], n1), d)) for t in locs]
    end
    b = sys.b
    # --- condensation, in chunks of elements (bounded triplet memory) ---
    B = spzeros(ns, ns); fs = b[findall(!iszero, spos)]
    fsk = zeros(ns)                                 # Σ_e A_bo t^e
    chunk = max(nthreads() * 64, 1 + 40_000_000 ÷ (nb * nb))
    @phase T :sc_condense begin
        for c0 in 0:chunk:ne_tot-1
            c1 = min(c0 + chunk, ne_tot) - 1; m = c1 - c0 + 1
            Ii = Vector{Int32}(undef, m * nb * nb); Jj = similar(Ii); Vv = Vector{Float64}(undef, m * nb * nb)
            dfs = [zeros(ns) for _ in 1:nthreads()]
            @threads :static for e in c0:c1
                g = gnodes(e); gb = spos[g[ib]]; fo = b[g[io]]
                Ch = cholesky(Aoo)
                Tm = Ch \ Aob                        # A_oo⁻¹ A_ob
                te = Ch \ fo                         # A_oo⁻¹ f_o
                Se = Abb - Abo * Tm
                off = (e - c0) * nb * nb
                @inbounds for jj in 1:nb, ii in 1:nb
                    k = off + ii + (jj - 1) * nb
                    Ii[k] = gb[ii]; Jj[k] = gb[jj]; Vv[k] = Se[ii, jj]
                end
                df = dfs[threadid_static()]
                v = Abo * te
                @inbounds for ii in 1:nb
                    df[gb[ii]] += v[ii]
                end
            end
            B += sparse(Ii, Jj, Vv, ns, ns)
            for df in dfs; fsk .+= df; end
        end
    end
    info[:skeleton] = ns
    info[:skeleton_nnz] = nnz(B)
    get(info, :keep_skeleton, false) && (info[:B] = B)
    rhs = fs .- fsk
    # ordering of the pinned skeleton matrix (METIS wants an exactly symmetric
    # pattern; B is symmetric up to round-off, so symmetrise first)
    perm = solver === :sc_direct ? (@phase T :setup (Bp = B[2:end, 2:end]; fill_ordering((Bp + Bp') ./ 2, ordering))) : nothing
    empty!(JX.JX_TIMINGS)
    us = JX.el_skeleton_solve(B, rhs; solver = solver === :sc_amg ? :amg : :direct,
                              amg_method = amg_method, amg_rtol = rtol, singular = true, perm = perm)
    T[:setup] = get(T, :setup, 0.0) + get(JX.JX_TIMINGS, :sc_factor, 0.0)
    T[:solve] = get(T, :solve, 0.0) + get(JX.JX_TIMINGS, :sc_skeleton, 0.0)
    solver === :sc_amg && (info[:iters] = JX.JX_AMG_STATS[].iters)
    # --- interior recovery ---
    u = zeros(n)
    sk = findall(!iszero, spos)
    u[sk] = us
    @phase T :sc_recover begin
        @threads :static for e in 0:ne_tot-1
            g = gnodes(e)
            Ch = cholesky(Aoo)
            u[g[io]] = Ch \ (b[g[io]] .- Aob * u[g[ib]])
        end
    end
    T[:setup] += T[:sc_condense]
    T[:solve] += T[:sc_recover]
    return gauge!(u, sys.m)
end

# thread slot for :static scheduling (threadid is stable within a :static loop)
threadid_static() = Threads.threadid()

# ---------------------------------------------------------------------------
# Spectral solvers on the uniform N_g^d grid, N_g = ne*N
# ---------------------------------------------------------------------------
fourier_grid(Ng) = [(i - 1) * LBOX / Ng for i in 1:Ng]

function fourier_error(u, U)
    u = u .+ (sum(U) - sum(u)) / length(u)
    e = u .- U
    return (linf = maximum(abs, e), l2rel = norm(e) / norm(U))
end

# apply the matrix A along dimension k of the d-dimensional array X (n^d)
function mode_mul!(Y, A, X, k, d, n)
    if k == 1
        mul!(reshape(Y, n, :), A, reshape(X, n, :))
    elseif k == d
        mul!(reshape(Y, :, n), reshape(X, :, n), transpose(A))
    else                                  # k == 2, d == 3
        for c in 1:n
            mul!(view(Y, :, :, c), view(X, :, :, c), transpose(A))
        end
    end
    return Y
end

function solve_ps!(T, info, Ng, d, r)
    x = fourier_grid(Ng)
    U, F = @phase T :rhs exact_on_grid(x, d, r)
    Q, invλ = @phase T :setup begin
        Q, λ, isc = JX._collocation_axis(Ng, LBOX, JX.FourierDerivativeMatrix)
        inv = d == 2 ? [(isc[i] && isc[j]) ? 0.0 : 1 / (λ[i] + λ[j]) for i in 1:Ng, j in 1:Ng] :
                       [(isc[i] && isc[j] && isc[k]) ? 0.0 : 1 / (λ[i] + λ[j] + λ[k]) for i in 1:Ng, j in 1:Ng, k in 1:Ng]
        Q, inv
    end
    u = @phase T :solve begin
        W1 = copy(F); W2 = similar(F)
        for k in 1:d                          # Û = F ×_k Qᵀ
            mode_mul!(W2, transpose(Q), W1, k, d, Ng); W1, W2 = W2, W1
        end
        W1 .*= invλ
        for k in 1:d                          # U = Û ×_k Q
            mode_mul!(W2, Q, W1, k, d, Ng); W1, W2 = W2, W1
        end
        W1
    end
    return vec(u), vec(U)
end

function solve_fft!(T, info, Ng, d, r)
    x = fourier_grid(Ng)
    U, F = @phase T :rhs exact_on_grid(x, d, r)
    S = @phase T :setup JX.FFTPoissonSolver(ntuple(_ -> Ng, d), ntuple(_ -> LBOX, d); flags = JX.FFTW.ESTIMATE)
    u = similar(F)
    @phase T :solve JX.fft_poisson_solve!(u, S, F)
    return vec(u), vec(U)
end

# ---------------------------------------------------------------------------
# One configuration
# ---------------------------------------------------------------------------
"""
    run_config(solver, d, ne, N; r = R_DEFAULT[d], rtol = 1e-12, ordering = :metis) -> row

Solve the d-dimensional periodic problem once with `solver` on ne^d elements
of order N (the spectral solvers on the (ne N)^d grid) and return the row:
errors, phase timings, iterations, sizes. `r` is the Fourier decay rate of
the exact solution (the difficulty of the problem); `ordering` the
fill-reducing ordering of the sparse Cholesky factorisations (:metis nested
dissection, or :amd for CHOLMOD's default).
"""
function run_config(solver::Symbol, d::Int, ne::Int, N::Int; r::Float64 = R_DEFAULT[d],
                    rtol = 1e-12, amg_method = "sa", ordering::Symbol = :metis)
    solver in SOLVERS || error("unknown solver $solver")
    T = Dict{Symbol, Float64}(); info = Dict{Symbol, Any}()
    Ng = ne * N; n = Ng^d
    if solver in SEM_SOLVERS
        sys = @phase T :assembly sem_system(N, ne, d, r)
        info[:nnz] = nnz(sys.K)
        u = solver in (:sc_direct, :sc_amg) ?
            solve_sc!(T, info, sys, solver; rtol = rtol, amg_method = amg_method, d = d, ordering = ordering) :
            solve_sem!(T, info, sys, solver; rtol = rtol, amg_method = Symbol(amg_method), ordering = ordering)
        err = sem_error(u, sys.U, sys.m)
    else
        u, U = solver === :ps ? solve_ps!(T, info, Ng, d, r) : solve_fft!(T, info, Ng, d, r)
        err = fourier_error(u, U)
    end
    g(k) = get(T, k, 0.0)
    return (solver = solver, d = d, r = r, ne = ne, nop = N, Ng = Ng, n = n,
            solved = get(info, :skeleton, n - (solver in SEM_SOLVERS ? 1 : 0)),
            linf = err.linf, l2rel = err.l2rel,
            assembly = g(:assembly), rhs = g(:rhs), setup = g(:setup), solve = g(:solve),
            total = g(:assembly) + g(:rhs) + g(:setup) + g(:solve),
            iters = get(info, :iters, 0), nnz = get(info, :nnz, 0),
            factor_nnz = get(info, :factor_nnz, 0), skeleton_nnz = get(info, :skeleton_nnz, 0),
            ordering = ordering)
end

end # module

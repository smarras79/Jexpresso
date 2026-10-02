#=============================================================================
 pmultigrid.jl — p-multigrid preconditioners for the conjugate-gradient
 solve of the 3D periodic SEM system (affine hexahedra), in two flavours:

   :linsolve_pmg => "amg"   p-levels N → ⌊N/2⌋ → … → 1, then ALGEBRAIC
                            multigrid (smoothed aggregation) on the p = 1 level
   :linsolve_pmg => "gmg"   the same p-levels, then GEOMETRIC h-multigrid on the
                            p = 1 grid (element vertices): trilinear transfers,
                            Galerkin coarse operators, coarsening by 2 down to a
                            small grid solved by sparse Cholesky

 The levels are spectral-element discretisations built by Jexpresso's SEM:
   • order N: the matrix K of the solve itself (sem_setup → DSS_laplace_sparse_3D
     → periodic reduction);
   • order p < N: the same operator rediscretised on the same elements with the
     LGL basis of order p (basis_structs_ξ_ω!, build_Interpolation_basis!), from
     the element geometry of the metrics: on an affine (box) element of the
     collocation SEM,  A_e = c_x ω̂⊗ω̂⊗K̂ + c_y ω̂⊗K̂⊗ω̂ + c_z K̂⊗ω̂⊗ω̂  (i fastest),
     c_x = J ξ_x², c_y = J η_y², c_z = J ζ_z²,  K̂, ω̂ the order-p GLL stiffness
     and weights. The order-N form is checked against the diagonal of K.
 Transfers between orders are the tensor-product Lagrange interpolation
 ℓ_a^{(p)}(ξ_k^{(q)}), applied element by element (three 1D passes): the
 prolongation writes each node from one owning element, the restriction (its
 transpose) accumulates element by element, elements coloured so that the
 threads never write the same node.

 Smoother on every level: Chebyshev–Jacobi of degree `degree` on the upper
 part [lower·λ, 1.1·λ] of the spectrum of D⁻¹A, λ estimated by power
 iteration: a polynomial in D⁻¹A, the same before and after the coarse
 correction, so the V-cycle is a symmetric preconditioner. Its work is
 sparse matvecs (threaded, JXSymCSC) and vector updates (threaded): no
 factorisation, no Gauss–Seidel sweep.

 The system is the singular periodic one (constants in the null space, a
 compatible right-hand side): CG runs on K itself; the preconditioner
 projects out the mean; the coarsest solves pin one node (a compatible
 singular system is then solved exactly by the pinned one).
=============================================================================#

# One multigrid level: the operator and its smoother's data and work vectors
struct _MGLevel
    A    :: SparseMatrixCSC{Float64, Int}
    op   :: JXSymCSC{Int}
    dinv :: Vector{Float64}
    λ    :: Float64                       # estimate of λ_max(D⁻¹A)
    x    :: Vector{Float64}
    b    :: Vector{Float64}
    r    :: Vector{Float64}
    d    :: Vector{Float64}
    w    :: Vector{Float64}
end
function _MGLevel(A::SparseMatrixCSC{Float64, Int})
    n = size(A, 1)
    dinv = Vector{Float64}(undef, n)
    @inbounds for j = 1:n
        q = _sc_findrow(A, j, j)
        (q > 0 && A.nzval[q] > 0) || error(" # p-multigrid: non-positive diagonal at row $j.")
        dinv[j] = 1 / A.nzval[q]
    end
    L = _MGLevel(A, JXSymCSC(A), dinv, 0.0, zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))
    return _MGLevel(A, L.op, dinv, _mg_lambda_max(L), L.x, L.b, L.r, L.d, L.w)
end

# Transfer between p-levels: element-wise tensor-product interpolation
struct _PTransfer
    Q     :: Matrix{Float64}            # (pf+1) × (pc+1): ℓ_a^{(pc)}(ξ_k^{(pf)})
    Qt    :: Matrix{Float64}
    connf :: Matrix{Int}                # (nelem, (pf+1)³), i fastest
    connc :: Matrix{Int}                # (nelem, (pc+1)³)
    own   :: Matrix{Bool}               # ((pf+1)³, nelem): element owns that fine node
    cols  :: _SCColours                 # no two elements of a colour share a coarse node
    bufs  :: Vector{NTuple{4, Vector{Float64}}}   # per thread: fine, coarse, two scratch
end

# Transfer between h-levels (p = 1 grid): sparse trilinear interpolation
struct _HTransfer
    P  :: SparseMatrixCSC{Float64, Int}           # n_fine × n_coarse
    Pt :: SparseMatrixCSC{Float64, Int}           # its transpose (column dots: threaded)
end

# Coarsest solves (pinned node 1)
struct _AMGCoarse{P}
    P  :: P
    rk :: Vector{Float64}
    zk :: Vector{Float64}
end
struct _CholCoarse{F}
    F  :: F
    rk :: Vector{Float64}
end

"""
    JXPMG

p-multigrid V-cycle preconditioner (jx_pmg_setup); apply with
`ldiv!(z, M, r)`.
"""
struct JXPMG{C}
    levels  :: Vector{_MGLevel}
    orders  :: Vector{Int}               # p of the p-levels (levels 1:length(orders))
    ptrans  :: Vector{_PTransfer}        # ptrans[ℓ]: between p-levels ℓ and ℓ+1
    htrans  :: Vector{_HTransfer}        # then between the h-levels
    coarse  :: C
    degree  :: Int
    lower   :: Float64
    kind    :: Symbol
end

# ── threaded vector loops ────────────────────────────────────────────────────
@inline _mg_nt(n) = n < 20_000 ? 1 : Threads.nthreads()

# λ_max(D⁻¹A) by power iteration (deterministic start, constants removed)
function _mg_lambda_max(L::_MGLevel; its::Int = 20)
    n = length(L.x);  v = L.d;  w = L.w
    @inbounds for i = 1:n;  v[i] = sin(0.7 * i) + 0.3 * cos(1.3 * i);  end
    v .-= sum(v) / n
    λ = 0.0
    for _ = 1:its
        mul!(w, L.op, v)
        w .*= L.dinv
        λ = sqrt(dot(w, w) / dot(v, v))
        v .= w ./ norm(w)
    end
    return λ
end

# Chebyshev–Jacobi smoothing of A x = b (x = 0 on entry if `zero`)
function _mg_cheb!(L::_MGLevel, deg::Int, lower::Float64, zero::Bool)
    x = L.x;  b = L.b;  r = L.r;  d = L.d;  w = L.w;  dinv = L.dinv
    n = length(x);  nt = _mg_nt(n)
    hi = 1.1 * L.λ;  lo = lower * L.λ
    θ = (hi + lo) / 2;  δ = (hi - lo) / 2;  σ = θ / δ;  ρ = 1 / σ
    zero || mul!(w, L.op, x)
    _sc_chunks(n, nt) do _, a, z
        @inbounds for i = a:z
            ri = dinv[i] * (zero ? b[i] : b[i] - w[i])
            r[i] = ri;  d[i] = ri / θ
            x[i] = zero ? d[i] : x[i] + d[i]
        end
        return nothing
    end
    for k = 2:deg
        mul!(w, L.op, d)
        ρn = 1 / (2σ - ρ);  c1 = ρn * ρ;  c2 = 2ρn / δ
        _sc_chunks(n, nt) do _, a, z
            @inbounds for i = a:z
                r[i] -= dinv[i] * w[i]
                d[i]  = c1 * d[i] + c2 * r[i]
                x[i] += d[i]
            end
            return nothing
        end
        ρ = ρn
    end
    return x
end

# r = b − A x
function _mg_residual!(L::_MGLevel)
    mul!(L.w, L.op, L.x)
    n = length(L.x)
    _sc_chunks(n, _mg_nt(n)) do _, a, z
        @inbounds for i = a:z;  L.r[i] = L.b[i] - L.w[i];  end
        return nothing
    end
    return L.r
end

# Y (mo³) = (Q⊗Q⊗Q) X (mi³), Q mo × mi, i fastest; T1 (mo·mi²), T2 (mo²·mi) scratch
@inline function _tp3!(Y, Q, X, T1, T2, mo::Int, mi::Int)
    @inbounds for jk = 0:mi*mi-1                                     # direction i
        oi = jk * mi;  oo = jk * mo
        for a = 1:mo
            s = 0.0
            for c = 1:mi;  s += Q[a, c] * X[oi+c];  end
            T1[oo+a] = s
        end
    end
    @inbounds for k = 0:mi-1, a = 1:mo                               # direction j
        oa = (k * mo + (a - 1)) * mo
        for i = 1:mo;  T2[oa+i] = 0.0;  end
        for c = 1:mi
            q = Q[a, c];  oc = (k * mi + (c - 1)) * mo
            for i = 1:mo;  T2[oa+i] += q * T1[oc+i];  end
        end
    end
    m2 = mo * mo
    @inbounds for a = 1:mo                                           # direction k
        oa = (a - 1) * m2
        for ij = 1:m2;  Y[oa+ij] = 0.0;  end
        for c = 1:mi
            q = Q[a, c];  oc = (c - 1) * m2
            for ij = 1:m2;  Y[oa+ij] += q * T2[oc+ij];  end
        end
    end
    return Y
end

# x_fine += P x_coarse (p-levels)
function _mg_prolong_add!(xf, xc, T::_PTransfer)
    nelem = size(T.connf, 1);  mf = size(T.Q, 1);  mc = size(T.Q, 2)
    nt = min(Threads.nthreads(), nelem)
    _sc_chunks(nelem, nt) do t, lo, hi
        F, C, T1, T2 = T.bufs[t]
        @inbounds for e = lo:hi
            for a = 1:mc^3;  C[a] = xc[T.connc[e, a]];  end
            _tp3!(F, T.Q, C, T1, T2, mf, mc)
            for a = 1:mf^3
                T.own[a, e] && (xf[T.connf[e, a]] += F[a])
            end
        end
        return nothing
    end
    return xf
end

# b_coarse = Pᵀ r_fine (p-levels)
function _mg_restrict!(bc, rf, T::_PTransfer)
    mf = size(T.Q, 1);  mc = size(T.Q, 2)
    fill!(bc, 0.0)
    for (c1, c2) in T.cols
        ne = c2 - c1 + 1
        _sc_chunks(ne, min(Threads.nthreads(), ne)) do t, lo, hi
            F, C, T1, T2 = T.bufs[t]
            @inbounds for q = c1+lo-1:c1+hi-1
                e = Int(T.cols.elems[q])
                for a = 1:mf^3;  F[a] = T.own[a, e] ? rf[T.connf[e, a]] : 0.0;  end
                _tp3!(C, T.Qt, F, T1, T2, mc, mf)
                for a = 1:mc^3;  bc[T.connc[e, a]] += C[a];  end
            end
            return nothing
        end
    end
    return bc
end

# y = M x for a CSC matrix via the column dots of its transpose (threaded)
@inline _mg_csc_apply!(y, Mt::SparseMatrixCSC, x) = mul!(y, JXSymCSCView(Mt), x)

# Mᵀ-column dot products: y[j] = Σ_p Mt[rowval[p], j] x[rowval[p]] = (M x)_j
struct JXSymCSCView{Ti}
    A :: SparseMatrixCSC{Float64, Ti}
end
function LinearAlgebra.mul!(y::AbstractVector, S::JXSymCSCView, x::AbstractVector)
    A = S.A;  n = size(A, 2)
    _sc_chunks(n, _mg_nt(n)) do _, lo, hi
        @inbounds for j = lo:hi
            s = 0.0
            for p = A.colptr[j]:A.colptr[j+1]-1;  s += A.nzval[p] * x[A.rowval[p]];  end
            y[j] = s
        end
        return nothing
    end
    return y
end

_coarse_solve!(x, b, C::_AMGCoarse) = begin
    n = length(b)
    @inbounds for i = 2:n;  C.rk[i-1] = b[i];  end
    fill!(C.zk, 0.0)
    LinearAlgebra.ldiv!(C.zk, C.P, C.rk)
    x[1] = 0.0
    @inbounds for i = 2:n;  x[i] = C.zk[i-1];  end
    x
end
_coarse_solve!(x, b, C::_CholCoarse) = begin
    n = length(b)
    @inbounds for i = 2:n;  C.rk[i-1] = b[i];  end
    z = C.F \ C.rk
    x[1] = 0.0
    @inbounds for i = 2:n;  x[i] = z[i-1];  end
    x
end

# V-cycle on level ℓ: L.x ≈ A⁻¹ L.b
function _mg_vcycle!(M::JXPMG, ℓ::Int)
    L = M.levels[ℓ]
    if ℓ == length(M.levels)
        _coarse_solve!(L.x, L.b, M.coarse)
        return nothing
    end
    _mg_cheb!(L, M.degree, M.lower, true)
    _mg_residual!(L)
    Lc = M.levels[ℓ+1];  np = length(M.orders)
    if ℓ < np
        _mg_restrict!(Lc.b, L.r, M.ptrans[ℓ])
    else
        _mg_csc_apply!(Lc.b, M.htrans[ℓ-np+1].P, L.r)         # Pᵀ r: column dots of P
    end
    _mg_vcycle!(M, ℓ + 1)
    if ℓ < np
        _mg_prolong_add!(L.x, Lc.x, M.ptrans[ℓ])
    else
        _mg_csc_apply!(L.r, M.htrans[ℓ-np+1].Pt, Lc.x)       # P x: column dots of Pᵀ
        L.x .+= L.r
    end
    _mg_cheb!(L, M.degree, M.lower, false)
    return nothing
end

function LinearAlgebra.ldiv!(z::AbstractVector, M::JXPMG, r::AbstractVector)
    L = M.levels[1];  n = length(r)
    μ = sum(r) / n
    @inbounds for i = 1:n;  L.b[i] = r[i] - μ;  end
    _mg_vcycle!(M, 1)
    μ = sum(L.x) / n
    @inbounds for i = 1:n;  z[i] = L.x[i] - μ;  end
    return z
end

# ── setup ───────────────────────────────────────────────────────────────────
# 1D GLL data of order p from Jexpresso's basis routines: nodes, weights, stiffness
function _pmg_gll(p::Int)
    lgl = basis_structs_ξ_ω!(LGL(), p, CPU())
    ξ = Vector{Float64}(lgl.ξ);  ω = Vector{Float64}(lgl.ω)
    B = build_Interpolation_basis!(LagrangeBasis(), ξ, ξ, Float64, CPU())
    return ξ, ω, sc_gll_stiffness(Matrix{Float64}(B.dψ), ω)
end

# interpolation from order pc to order pf: Q[k, a] = ℓ_a^{(pc)}(ξ_k^{(pf)})
function _pmg_interp(pc::Int, pf::Int)
    ξc, _, _ = _pmg_gll(pc);  ξf, _, _ = _pmg_gll(pf)
    B = build_Interpolation_basis!(LagrangeBasis(), ξc, ξf, Float64, CPU())
    return Matrix{Float64}(transpose(B.ψ))
end

# the order-p SEM operator Σ_e A_e on the periodic numbering conn (i fastest)
function _pmg_assemble(conn::Matrix{Int}, c::Matrix{Float64}, K̂, ω, p::Int, n::Int)
    nelem = size(conn, 1);  m = p + 1;  per = m^3 * 3m
    I = Vector{Int}(undef, nelem * per);  J = similar(I);  V = Vector{Float64}(undef, nelem * per)
    lin(i, j, k) = i + (j - 1) * m + (k - 1) * m * m
    Threads.@threads :static for e = 1:nelem
        q = (e - 1) * per
        cx = c[1, e];  cy = c[2, e];  cz = c[3, e]
        @inbounds for k = 1:m, j = 1:m, i = 1:m
            r = conn[e, lin(i, j, k)]
            for a = 1:m
                q += 1;  I[q] = r;  J[q] = conn[e, lin(a, j, k)];  V[q] = cx * K̂[i, a] * ω[j] * ω[k]
                q += 1;  I[q] = r;  J[q] = conn[e, lin(i, a, k)];  V[q] = cy * ω[i] * K̂[j, a] * ω[k]
                q += 1;  I[q] = r;  J[q] = conn[e, lin(i, j, a)];  V[q] = cz * ω[i] * ω[j] * K̂[k, a]
            end
        end
    end
    A = sparse(I, J, V, n, n)
    return (A + sparse(transpose(A))) ./ 2          # exactly symmetric for JXSymCSC
end

# greedy element colouring: no two elements of a colour share a node of conn
function _mg_colour(conn::Matrix{Int}, n::Int)
    nelem, npl = size(conn)
    eptr = zeros(Int, n + 1)
    @inbounds for e = 1:nelem, a = 1:npl;  eptr[conn[e, a]+1] += 1;  end
    eptr[1] = 1
    for i = 1:n;  eptr[i+1] += eptr[i];  end
    elist = Vector{Int32}(undef, eptr[n+1] - 1);  nxt = eptr[1:n]
    @inbounds for e = 1:nelem, a = 1:npl
        g = conn[e, a];  elist[nxt[g]] = e;  nxt[g] += 1
    end
    pos = Int32.(1:n)                                         # every node "active"
    act = trues(npl, nelem)
    return _sc_colour(conn, pos, act, eptr, elist, npl, nelem)
end

function _pmg_transfer(connf, connc, pf, pc, nf)
    nelem = size(connf, 1);  mf = pf + 1;  mc = pc + 1
    own  = Matrix{Bool}(undef, mf^3, nelem)
    seen = falses(nf)
    @inbounds for e = 1:nelem, a = 1:mf^3
        g = connf[e, a];  own[a, e] = !seen[g];  seen[g] = true
    end
    Q = _pmg_interp(pc, pf)
    mx = max(mf, mc)
    bufs = [(zeros(mf^3), zeros(mc^3), zeros(mx^3), zeros(mx^3)) for _ = 1:Threads.nthreads()]
    return _PTransfer(Q, Matrix(transpose(Q)), connf, connc, own,
                      _mg_colour(connc, maximum(connc)), bufs)
end

"""
    jx_pmg_setup(sem, K, cls; coarse = :amg, degree = 3, lower = 0.25, amg_method = "sa") -> JXPMG

p-multigrid preconditioner for the 3D periodic SEM system K (class numbering
`cls` of the nodes of `sem.mesh`). Requires uniform, axis-aligned box
elements on a fully periodic box (the 3D benchmark deck); errors otherwise.
"""
function jx_pmg_setup(sem, K::SparseMatrixCSC, cls::Vector{Int}; coarse::Symbol = :amg,
                      degree::Int = 3, lower::Float64 = 0.25, amg_method = "sa")
    coarse in (:amg, :gmg) || error(" # :linsolve_pmg => \"$coarse\"; expected \"amg\" or \"gmg\".")
    mesh = sem.mesh;  met = sem.metrics
    Int(mesh.nsd) == 3 || error(" # p-multigrid: 3D only.")
    ngl = Int(mesh.ngl);  N = ngl - 1;  nelem = Int(mesh.nelem)
    X = mesh.coords;  cj = mesh.connijk
    t0 = time_ns()

    # ── element geometry: uniform axis-aligned boxes, integer position ────────
    # local axis d of element e runs along physical axis ax[d, e], in the
    # direction sg[d, e] (Jexpresso's elements need not be oriented alike)
    lo = (minimum(view(X, 1, :)), minimum(view(X, 2, :)), minimum(view(X, 3, :)))
    h  = zeros(3);  E = zeros(Int, 3, nelem);  c = zeros(3, nelem)
    ax = zeros(Int, 3, nelem);  sg = zeros(Int, 3, nelem)
    for e = 1:nelem
        g0 = cj[e, 1, 1, 1]
        corner = (cj[e, ngl, 1, 1], cj[e, 1, ngl, 1], cj[e, 1, 1, ngl])
        mn = [X[k, g0] for k = 1:3]
        for d = 1:3
            v = X[:, corner[d]] .- X[:, g0]
            k = argmax(abs.(v));  hd = abs(v[k])
            all(abs(v[q]) <= 1e-10 * hd for q = 1:3 if q != k) ||
                error(" # p-multigrid: element $e is not an axis-aligned box (local axis $d: $v).")
            ax[d, e] = k;  sg[d, e] = v[k] > 0 ? 1 : -1
            h[k] == 0 && (h[k] = hd)
            abs(hd - h[k]) <= 1e-10 * h[k] || error(" # p-multigrid: elements of different sizes.")
            mn[k] = min(mn[k], X[k, corner[d]])
        end
        sort(ax[:, e]) == [1, 2, 3] || error(" # p-multigrid: element $e has two local axes along one direction.")
        for k = 1:3;  E[k, e] = round(Int, (mn[k] - lo[k]) / h[k]);  end
        J = met.Je[e, 1, 1, 1]
        c[1, e] = J * (met.dξdx[e, 1, 1, 1]^2 + met.dξdy[e, 1, 1, 1]^2 + met.dξdz[e, 1, 1, 1]^2)
        c[2, e] = J * (met.dηdx[e, 1, 1, 1]^2 + met.dηdy[e, 1, 1, 1]^2 + met.dηdz[e, 1, 1, 1]^2)
        c[3, e] = J * (met.dζdx[e, 1, 1, 1]^2 + met.dζdy[e, 1, 1, 1]^2 + met.dζdz[e, 1, 1, 1]^2)
    end
    ne = [maximum(E[d, :]) + 1 for d = 1:3]
    prod(ne) == nelem || error(" # p-multigrid: the elements do not tile a box.")
    prod(ne .* N) == size(K, 1) ||
        error(" # p-multigrid: the system is not periodic in x, y and z ($(size(K,1)) unknowns).")

    # ── fine level: Jexpresso's K, its numbering; the tensor form checked on diag(K)
    lin(i, j, k, m) = i + (j - 1) * m + (k - 1) * m * m
    connN = Matrix{Int}(undef, nelem, ngl^3)
    for e = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        connN[e, lin(i, j, k, ngl)] = cls[cj[e, i, j, k]]
    end
    K̂N = sc_gll_stiffness(Matrix{Float64}(sem.basis.dψ), Vector{Float64}(sem.ω))
    ωN = Vector{Float64}(sem.ω)
    dg = zeros(size(K, 1))
    for e = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        dg[connN[e, lin(i, j, k, ngl)]] += c[1, e] * K̂N[i, i] * ωN[j] * ωN[k] +
            c[2, e] * ωN[i] * K̂N[j, j] * ωN[k] + c[3, e] * ωN[i] * ωN[j] * K̂N[k, k]
    end
    dK = Vector{Float64}(diag(K))
    maximum(abs, dg .- dK) <= 1e-9 * maximum(abs, dK) ||
        error(" # p-multigrid: K is not the tensor-product SEM operator of these elements " *
              "(max diagonal mismatch $(maximum(abs, dg .- dK))).")
    Kf = K isa SparseMatrixCSC{Float64, Int} ? K : SparseMatrixCSC{Float64, Int}(K)
    JX_TIMINGS[:pmg_geometry] = (time_ns() - t0) / 1e9;  t0 = time_ns()

    # ── p-levels ──────────────────────────────────────────────────────────────
    orders = [N];  while orders[end] > 1;  push!(orders, orders[end] ÷ 2);  end
    levels = [_MGLevel(Kf)];  conns = [connN];  ns = [size(K, 1)]
    for p in orders[2:end]
        m = p + 1
        conn = Matrix{Int}(undef, nelem, m^3)
        nd = ne .* p
        g = zeros(Int, 3)
        for e = 1:nelem, k = 0:p, j = 0:p, i = 0:p
            for (d, l) in ((1, i), (2, j), (3, k))           # local index -> physical grid index
                q = ax[d, e]
                g[q] = mod(E[q, e] * p + (sg[d, e] > 0 ? l : p - l), nd[q])
            end
            conn[e, lin(i + 1, j + 1, k + 1, m)] = 1 + g[1] + nd[1] * (g[2] + nd[2] * g[3])
        end
        _, ω, K̂ = _pmg_gll(p)
        push!(levels, _MGLevel(_pmg_assemble(conn, c, K̂, ω, p, prod(nd))))
        push!(conns, conn);  push!(ns, prod(nd))
    end
    JX_TIMINGS[:pmg_levels] = (time_ns() - t0) / 1e9;  t0 = time_ns()
    ptrans = [_pmg_transfer(conns[l], conns[l+1], orders[l], orders[l+1], ns[l])
              for l = 1:length(orders)-1]

    JX_TIMINGS[:pmg_transfers] = (time_ns() - t0) / 1e9;  t0 = time_ns()
    # ── coarse: AMG on p = 1, or geometric h-levels on the p = 1 grid ─────────
    htrans = _HTransfer[]
    if coarse === :amg
        A1 = levels[end].A
        S  = jx_amg_setup(A1[2:end, 2:end]; method = amg_method)
        C  = _AMGCoarse(S.P, zeros(size(A1, 1) - 1), zeros(size(A1, 1) - 1))
    else
        dims = copy(ne)
        while all(iseven, dims) && prod(dims) > 512
            P  = _hmg_interp(dims)
            A  = levels[end].A
            Ac = sparse(transpose(P)) * (A * P)
            Ac = (Ac + sparse(transpose(Ac))) ./ 2
            push!(htrans, _HTransfer(P, sparse(transpose(P))))
            push!(levels, _MGLevel(Ac))
            dims = dims .÷ 2
        end
        Ac = levels[end].A
        F  = cholesky(Symmetric(Ac[2:end, 2:end]))
        C  = _CholCoarse(F, zeros(size(Ac, 1) - 1))
    end
    JX_TIMINGS[:pmg_coarse] = (time_ns() - t0) / 1e9
    return JXPMG(levels, orders, ptrans, htrans, C, degree, lower, coarse)
end

# trilinear interpolation on the periodic p = 1 grid: dims → dims/2
function _hmg_interp(dims::Vector{Int})
    nf = prod(dims);  dc = dims .÷ 2;  nc = prod(dc)
    I = Int[];  J = Int[];  V = Float64[]
    sizehint!(I, 8nf);  sizehint!(J, 8nf);  sizehint!(V, 8nf)
    par(i, n) = iseven(i) ? ((i ÷ 2, 1.0),) : (((i - 1) ÷ 2, 0.5), (mod((i + 1) ÷ 2, n), 0.5))
    for k = 0:dims[3]-1, j = 0:dims[2]-1, i = 0:dims[1]-1
        row = 1 + i + dims[1] * (j + dims[2] * k)
        for (ci, wi) in par(i, dc[1]), (cj, wj) in par(j, dc[2]), (ck, wk) in par(k, dc[3])
            push!(I, row);  push!(J, 1 + ci + dc[1] * (cj + dc[2] * ck));  push!(V, wi * wj * wk)
        end
    end
    return sparse(I, J, V, nf, nc)
end

"""
    jx_pmg_cg(K, b, M; rtol = 1e-12, itmax = 1000) -> x

CG on the (singular, compatible) periodic system K x = b, preconditioned by
the p-multigrid V-cycle M; the iterations and residual go to JX_AMG_STATS.
"""
function jx_pmg_cg(K::SparseMatrixCSC, b::AbstractVector, M::JXPMG; rtol::Real = 1e-12, itmax::Int = 1000)
    x, st = Krylov.cg(M.levels[1].op, Vector{Float64}(b); M = M, ldiv = true,
                      rtol = Float64(rtol), atol = 0.0, itmax = itmax, history = true)
    r0 = isempty(st.residuals) ? NaN : first(st.residuals)
    rel = isempty(st.residuals) || r0 == 0 ? 0.0 : last(st.residuals) / r0
    JX_AMG_STATS[] = (iters = st.niter, rel_resid = rel, levels = length(M.levels),
                      method = Symbol("pmg_", M.kind))
    st.solved || @warn " # p-multigrid CG did not converge to rtol=$rtol in $itmax iterations " *
                       "(final relative residual $rel)."
    return x
end

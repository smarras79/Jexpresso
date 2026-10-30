#---------------------------------------------------------------------------------
# Split time integration: du/dt = L u + N(u), with a stiff LINEAR part L
# integrated exactly by an exponential integrator (or implicitly by an IMEX
# one) and the nonlinear rest N(u) by Jexpresso's ordinary rhs!.
#
# Why. A higher-order linear term makes an explicit step useless: on a mesh
# whose -∂ₓ² has largest eigenvalue λ, ∂ₓ⁴ has λ², so an explicit Δt has to
# be ~1/λ² (Kuramoto-Sivashinsky on 64 p=4 elements: λ² ≈ 5500, Δt ≈ 2e-4,
# 750 000 steps to t = 150). Taking L out of the explicit part removes that
# limit altogether; the step is then set by the nonlinear dynamics alone
# (Δt = 1/4 there: 600 steps). This is the ETDRK4 approach of Cox &
# Matthews (2002) and Kassam & Trefethen (2005), on spectral elements.
#
# Deck (user_inputs.jl), 1D continuous Galerkin only:
#
#     :lsplit_linear => true,
#     :split_d2      => [c2],       # one coefficient per equation
#     :split_d4      => [c4],
#     :ode_solver    => ETDRK4(),   # or any SplitODEProblem solver
#
# gives, for every equation ieq,
#
#     L = c2[ieq] D2 + c4[ieq] D2²,     D2 = -M⁻¹ K,
#
# where K is the weak Laplacian (∫ψᵢ' ψⱼ' dx, LGL quadrature) and M the
# lumped LGL mass: D2 is exactly the second derivative that the viscous
# operator applies (_expansion_visc!), and D2² is the C⁰ mixed (Ciarlet-
# Raviart) fourth derivative. N(u) is whatever user_flux!/user_source! give,
# so the deck keeps :lvisc => false: the derivatives live in L.
#
# THE φ-FUNCTIONS. Exponential integrators need φₖ(hL). ExponentialUtilities
# computes them for a general dense matrix one column at a time, each column
# an exponential of an (n+k)×(n+k) matrix: seconds for n = 256 with ‖hL‖ ~
# 10³. Here they are exact and cheap, because L is similar to a SYMMETRIC
# matrix: with S = M^{1/2} and B = S⁻¹ K S⁻¹ (symmetric),
#
#     D2 = S⁻¹ (-B) S,    L = S⁻¹ (-c2 B + c4 B²) S,
#
# so every function of hL is S⁻¹ g(B) S for a scalar function g, and all it
# takes is the eigendecomposition of B. A function of L is never formed: it
# is kept as its diagonal g(Λ) in the eigenbasis of B, so sums, products,
# powers and multiplication by L (the algebra OrdinaryDiffEq's
# expRK_operators does on them) act on diagonals, and setting up ETDRK4's
# six operators costs O(n). Applying one is a change of basis, a scaling and
# the change back, and there are two ways to change basis:
#
# - UNIFORM PERIODIC MESH (SEMBloch), the common case. Every element is the
#   same, so with the nodes numbered element by element (ngl - 1 per element,
#   the right end node being the next element's first) B is block-circulant,
#   and the DFT over the elements block-diagonalizes it: one Hermitian
#   (ngl-1)×(ngl-1) eigenproblem per wavenumber (Bloch waves). A product is a
#   real FFT over the elements, an (ngl-1)×(ngl-1) block per wavenumber and
#   the inverse FFT: O(n log nel), the cost of a Fourier pseudospectral code.
#   Kassam & Trefethen mesh (n = 256): 1.2 μs per product. L u itself is
#   applied the same way (1.2 μs; the sparse product takes 2.0 μs).
# - ANY OTHER MESH (SEMEigen): one dense symmetric eigendecomposition
#   B = Q Λ Qᵀ, and Qᵀ, Q per product with the one dense Q that all the
#   operators share (it stays in cache): O(n²), 11 μs at n = 256 (forming the
#   six operators instead, as dense n×n matrices, was memory-bound: 47 μs).
#   L u is applied with the sparse matrix (14 nonzeros per row at p = 4).
#
# build_split_linear_operator checks the block-circulant structure entry by
# entry (to 1e-12 relative) and falls back to the dense eigenvectors when it
# does not hold, e.g. on a stretched mesh.
#
# Scope: periodic 1D (:lperiodic_1d => true), CPU, ContGal, one rank. The
# dense fallback is one n³ eigensolve (~40 ms for n = 256), the right trade
# for 1D problems of up to a few thousand nodes; past that, an IMEX solver on
# the sparse L (KenCarp4, which receives L.A) is the scalable choice.
#---------------------------------------------------------------------------------

"""
    SEMSplitBasis{T}

Eigendata of B = S⁻¹ K S⁻¹ shared by the linear operator of a split problem
and every function of it (see the header of split_linear.jl): [`SEMBloch`](@ref)
on a uniform periodic mesh, [`SEMEigen`](@ref) on any other. `Λ` holds the
eigenvalues of B in the basis' own order, and an operator's diagonal `d` has
`length(Λ)` entries per equation in that same order.
"""
abstract type SEMSplitBasis{T} end

"""
    SEMEigen{T}

B = Q diag(Λ) Qᵀ for one equation's npoin nodes (the same for all
equations), and the diagonal of S = M^{1/2}.
"""
struct SEMEigen{T} <: SEMSplitBasis{T}
    Q::Matrix{T}
    Λ::Vector{T}
    s::Vector{T}
    sinv::Vector{T}
    npoin::Int
    neqs::Int
end

"""
    SEMBloch{T}

B on a uniform periodic mesh, block-diagonalized by the DFT over the
elements. Slot p = (e-1)·nb + i is the i-th of the nb = ngl - 1 nodes that
element e owns (the elements taken in the order of the periodic chain), and
node `perm[p]` of the mesh. With Cⱼ the nb×nb block that couples an element
to the j-th one after it,

    B̂ₖ = Σⱼ Cⱼ e^{2πi jk/nel} = Uₖ diag(Λₖ) Uₖᴴ,    k = 0 … nel÷2,

(Hermitian; wavenumber nel - k is the complex conjugate of k, which the real
transform does not store). `Λ[(k)·nb + m]` is the m-th eigenvalue of B̂ₖ. `X`,
`Xh`, `Yh` and the FFTW plans `pf` (rfft over the elements) and `pb` (brfft
back) are the workspace of every product: these operators are applied one
at a time (one rank, one thread), never concurrently.
"""
struct SEMBloch{T,PF,PB} <: SEMSplitBasis{T}
    perm::Vector{Int}
    sp::Vector{T}               # s[perm]
    sinvp::Vector{T}            # sinv[perm]
    U::Array{Complex{T},3}      # nb × nb × nk
    Λ::Vector{T}                # nb·nk
    nb::Int
    nel::Int
    nk::Int
    npoin::Int
    neqs::Int
    X::Matrix{T}                # nb × nel
    Xh::Matrix{Complex{T}}      # nb × nk
    Yh::Matrix{Complex{T}}      # nb × nk
    pf::PF
    pb::PB
end

"""
    SEMLinearOperator1D{T,B,PT}

The linear part α·L of a split problem (see the header of split_linear.jl).
Products use `S`, L itself as a function of B, on a uniform periodic mesh
([`SEMBloch`](@ref)), and the sparse `A` on any other (`S === nothing`: the
sparse product beats the dense eigenvectors). `exp`,
`ExponentialUtilities.phi` and products with functions of L return
[`SEMSpectralOperator1D`](@ref)s. `α` is the scalar factor of the expression
this object stands for: `dt*L` is the same data with α = dt, which is how
OrdinaryDiffEq's exponential integrators ask for φ(dt L).
"""
struct SEMLinearOperator1D{T,B<:SEMSplitBasis{T},PT} <: AbstractMatrix{T}
    A::SparseMatrixCSC{T,Int}
    E::B
    λ::Vector{T}            # eigenvalues of L, equation by equation, in E's order
    α::T
    S::PT
end

"""
    SEMSpectralOperator1D{T,B,C}

A function of the linear operator, S⁻¹ g(B) S equation by equation, stored as
its diagonal `d` = g(Λ) in the basis `E` (see the header of split_linear.jl).
`c` is what a product needs besides `E`: the two work vectors of the dense
basis, or, for [`SEMBloch`](@ref), the nb×nb blocks Uₖ diag(dₖ) Uₖᴴ of every
wavenumber and equation, formed once as static matrices (so that their
products unroll). `mul!` allocates nothing.
"""
struct SEMSpectralOperator1D{T,B<:SEMSplitBasis{T},C} <: AbstractMatrix{T}
    E::B
    d::Vector{T}
    c::C
end
function SEMSpectralOperator1D(E::SEMSplitBasis{T}, d::Vector{T}) where {T}
    c = _spectral_cache(E, d)
    return SEMSpectralOperator1D{T,typeof(E),typeof(c)}(E, d, c)
end
_spectral_cache(E::SEMEigen{T}, d) where {T} = (w1 = zeros(T, E.npoin), w2 = zeros(T, E.npoin))
function _spectral_cache(E::SEMBloch{T}, d) where {T}
    nb = E.nb; nk = E.nk; nλ = length(E.Λ)
    Gk = zeros(Complex{T}, nb, nb)
    G  = Vector{SMatrix{nb,nb,Complex{T},nb*nb}}(undef, nk*E.neqs)
    for ieq = 1:E.neqs, k = 1:nk
        fill!(Gk, zero(Complex{T}))
        for j = 1:nb, i = 1:nb, m = 1:nb
            Gk[i, j] += E.U[i, m, k]*d[(ieq - 1)*nλ + (k - 1)*nb + m]*conj(E.U[j, m, k])
        end
        G[(ieq - 1)*nk + k] = Gk
    end
    return G
end

#------------------------------ the linear operator ------------------------------
Base.size(L::SEMLinearOperator1D) = size(L.A)
Base.getindex(L::SEMLinearOperator1D, i::Int, j::Int) = L.α*L.A[i, j]

# y = a·L x + b·y: through the Bloch form S when there is one, the sparse A otherwise
_linear_apply!(y, ::Nothing, L, x, a, b) = mul!(y, L.A, x, a, b)
_linear_apply!(y, S::SEMSpectralOperator1D, L, x, a, b) = _spectral_apply!(y, S, x, a, b)
LinearAlgebra.mul!(y::AbstractVector, L::SEMLinearOperator1D, x::AbstractVector) =
    _linear_apply!(y, L.S, L, x, L.α, false)
LinearAlgebra.mul!(y::AbstractVector, L::SEMLinearOperator1D, x::AbstractVector, a::Number, b::Number) =
    _linear_apply!(y, L.S, L, x, a*L.α, b)
Base.:*(L::SEMLinearOperator1D, x::AbstractVector) = L.α*(L.A*x)
Base.:*(L::SEMLinearOperator1D, X::AbstractMatrix) = L.α*(L.A*X)
Base.:*(a::Number, L::SEMLinearOperator1D{T}) where {T} = SEMLinearOperator1D(L.A, L.E, L.λ, T(a*L.α), L.S)
Base.:*(L::SEMLinearOperator1D, a::Number) = a*L
Base.:/(L::SEMLinearOperator1D, a::Number) = inv(a)*L

Base.exp(L::SEMLinearOperator1D) = SEMSpectralOperator1D(L.E, exp.(L.α .* L.λ))

function ExponentialUtilities.phi(L::SEMLinearOperator1D{T}, k::Integer; kwargs...) where {T}
    z = L.α .* L.λ
    return [SEMSpectralOperator1D(L.E, _phi_scalar.(z, j)) for j = 0:k]
end

#--------------------------- functions of the operator ---------------------------
Base.size(P::SEMSpectralOperator1D) = (P.E.npoin*P.E.neqs, P.E.npoin*P.E.neqs)
function Base.getindex(P::SEMSpectralOperator1D{T,<:SEMEigen}, i::Int, j::Int) where {T}
    n = P.E.npoin
    ei, ii = divrem(i - 1, n); ej, jj = divrem(j - 1, n)
    ei == ej || return zero(T)
    Q = P.E.Q; off = ei*n
    acc = zero(T)
    @inbounds for k = 1:n
        acc += Q[ii+1, k]*P.d[off+k]*Q[jj+1, k]
    end
    return P.E.sinv[ii+1]*acc*P.E.s[jj+1]
end
# an entry of a Bloch-represented operator: its j-th column (inspection only)
function Base.getindex(P::SEMSpectralOperator1D{T,<:SEMBloch}, i::Int, j::Int) where {T}
    x = zeros(T, size(P, 2)); x[j] = one(T)
    return (P*x)[i]
end

# y = a·P x + b·y, equation by equation: S x, Qᵀ, diag(d), Q, S⁻¹
function _spectral_apply!(y::AbstractVector, P::SEMSpectralOperator1D{T,<:SEMEigen}, x::AbstractVector, a, b) where {T}
    E = P.E; n = E.npoin; w1 = P.c.w1; w2 = P.c.w2
    for ieq = 1:E.neqs
        off = (ieq - 1)*n
        @inbounds for i = 1:n
            w1[i] = E.s[i]*x[off+i]
        end
        mul!(w2, transpose(E.Q), w1)
        @inbounds for i = 1:n
            w2[i] *= P.d[off+i]
        end
        mul!(w1, E.Q, w2)
        if iszero(b)
            @inbounds for i = 1:n
                y[off+i] = a*E.sinv[i]*w1[i]
            end
        else
            @inbounds for i = 1:n
                y[off+i] = a*E.sinv[i]*w1[i] + b*y[off+i]
            end
        end
    end
    return y
end

# y = a·P x + b·y, equation by equation: S x gathered element by element, rfft
# over the elements, the nb×nb block of every wavenumber, brfft (which is nel
# times the inverse), S⁻¹ and scatter
function _spectral_apply!(y::AbstractVector, P::SEMSpectralOperator1D{T,<:SEMBloch}, x::AbstractVector, a, b) where {T}
    E = P.E; G = P.c
    n = E.npoin; nk = E.nk
    perm = E.perm; X = E.X; Xh = E.Xh; Yh = E.Yh
    c = a/E.nel
    for ieq = 1:E.neqs
        off = (ieq - 1)*n
        @inbounds for p = 1:n
            X[p] = E.sp[p]*x[off+perm[p]]
        end
        mul!(Xh, E.pf, X)
        _bloch_blocks!(Yh, G, Xh, (ieq - 1)*nk)
        mul!(X, E.pb, Yh)
        if iszero(b)
            @inbounds for p = 1:n
                y[off+perm[p]] = c*E.sinvp[p]*X[p]
            end
        else
            @inbounds for p = 1:n
                q = off + perm[p]
                y[q] = c*E.sinvp[p]*X[p] + b*y[q]
            end
        end
    end
    return y
end
# Ŷₖ = Gₖ X̂ₖ, wavenumber by wavenumber, with the static nb×nb block of equation
# block `off` (fully unrolled: 0.39 μs for the 33 4×4 blocks of n = 256,
# against 0.59 μs with loops over a dense array)
function _bloch_blocks!(Yh, G::Vector{SMatrix{N,N,Complex{T},L}}, Xh, off) where {N,T,L}
    @inbounds for k in axes(Xh, 2)
        xk = SVector{N,Complex{T}}(ntuple(j -> @inbounds(Xh[j, k]), Val(N)))
        yk = G[off+k]*xk
        for i = 1:N
            Yh[i, k] = yk[i]
        end
    end
    return Yh
end
LinearAlgebra.mul!(y::AbstractVector, P::SEMSpectralOperator1D, x::AbstractVector) =
    _spectral_apply!(y, P, x, true, false)
LinearAlgebra.mul!(y::AbstractVector, P::SEMSpectralOperator1D, x::AbstractVector, a::Number, b::Number) =
    _spectral_apply!(y, P, x, a, b)
Base.:*(P::SEMSpectralOperator1D, x::AbstractVector) = _spectral_apply!(similar(x), P, x, true, false)

_same_basis(P1, P2) = P1.E === P2.E ||
    error("functions of two different split linear operators cannot be combined")
Base.:+(P1::SEMSpectralOperator1D, P2::SEMSpectralOperator1D) = (_same_basis(P1, P2); SEMSpectralOperator1D(P1.E, P1.d .+ P2.d))
Base.:-(P1::SEMSpectralOperator1D, P2::SEMSpectralOperator1D) = (_same_basis(P1, P2); SEMSpectralOperator1D(P1.E, P1.d .- P2.d))
Base.:-(P::SEMSpectralOperator1D) = SEMSpectralOperator1D(P.E, -P.d)
Base.:*(P1::SEMSpectralOperator1D, P2::SEMSpectralOperator1D) = (_same_basis(P1, P2); SEMSpectralOperator1D(P1.E, P1.d .* P2.d))
Base.:*(a::Number, P::SEMSpectralOperator1D) = SEMSpectralOperator1D(P.E, a .* P.d)
Base.:*(P::SEMSpectralOperator1D, a::Number) = a*P
Base.:/(P::SEMSpectralOperator1D, a::Number) = SEMSpectralOperator1D(P.E, P.d ./ a)
Base.:^(P::SEMSpectralOperator1D, p::Integer) = SEMSpectralOperator1D(P.E, P.d .^ p)
Base.:*(L::SEMLinearOperator1D, P::SEMSpectralOperator1D) = (_same_basis(L, P); SEMSpectralOperator1D(P.E, (L.α .* L.λ) .* P.d))
Base.:*(P::SEMSpectralOperator1D, L::SEMLinearOperator1D) = L*P

# φⱼ(z) = Σₘ zᵐ/(m+j)!  (φ₀ = exp). Taylor for |z| < 1, where the closed form
# (eᶻ - Σ_{m<j} zᵐ/m!)/zʲ cancels; the closed form elsewhere, where it does
# not (z ≤ -1 for the dissipative modes, 0 < z ≤ 1 at most for the unstable
# ones of a well-resolved run).
function _phi_scalar(z::T, j::Integer) where {T<:AbstractFloat}
    j == 0 && return exp(z)
    if abs(z) < one(T)
        s    = zero(T)
        term = one(T)
        for m = 1:j
            term /= m                       # 1/j!
        end
        for m = 0:30
            s    += term
            term *= z/(m + j + 1)
        end
        return s
    else
        p = zero(T)
        t = one(T)
        for m = 0:j-1
            p += t
            t *= z/(m + 1)
        end
        return (exp(z) - p)/z^j
    end
end

_split_coeffs(inputs, key, neqs) = begin
    c = get(inputs, key, 0.0)
    c isa Number ? fill(Float64(c), neqs) : Float64.(collect(c))
end

"""
    _sem_bloch(Bs, s, sinv, connijk, nelem, ngl, npoin, neqs) -> SEMBloch or nothing

The Bloch representation of the symmetric B = S⁻¹ K S⁻¹ (sparse `Bs`), or
`nothing` when B is not block-circulant over the elements (a non-uniform
mesh), in which case the caller falls back to the dense eigenvectors. See
[`SEMBloch`](@ref).
"""
function _sem_bloch(Bs::SparseMatrixCSC{T,Int}, s, sinv, connijk, nelem, ngl, npoin, neqs) where {T}

    nb = ngl - 1
    nb >= 1 || return nothing

    # the elements in the order of the periodic chain: each starts at the node
    # where the previous one ends, and the last one ends where the first starts
    first_of = Dict{Int,Int}()
    for e = 1:nelem
        haskey(first_of, connijk[e, 1, 1, 1]) && return nothing
        first_of[connijk[e, 1, 1, 1]] = e
    end
    order = zeros(Int, nelem); order[1] = 1
    for b = 2:nelem
        order[b] = get(first_of, connijk[order[b-1], ngl, 1, 1], 0)
        order[b] == 0 && return nothing
    end
    connijk[order[nelem], ngl, 1, 1] == connijk[order[1], 1, 1, 1] || return nothing
    perm = vec([connijk[order[b], i, 1, 1] for i = 1:nb, b = 1:nelem])
    (length(perm) == npoin && isperm(perm)) || return nothing
    slot = invperm(perm)

    # Cⱼ from the first block row, then every stored entry of B checked against
    # them, and every nonzero of the Cⱼ found in every block row: B is then
    # block-circulant to 1e-12 relative
    rows = rowvals(Bs); vals = nonzeros(Bs)
    C = zeros(T, nb, nb, nelem)
    for col = 1:npoin, q in nzrange(Bs, col)
        br, ir = divrem(slot[rows[q]] - 1, nb); bc, ic = divrem(slot[col] - 1, nb)
        br == 0 && (C[ir+1, ic+1, bc+1] = vals[q])
    end
    tol  = 1.0e-12*maximum(abs, vals; init = zero(T))
    nnzC = count(c -> abs(c) > tol, C)
    seen = zeros(Int, nelem)
    for col = 1:npoin, q in nzrange(Bs, col)
        br, ir = divrem(slot[rows[q]] - 1, nb); bc, ic = divrem(slot[col] - 1, nb)
        c = C[ir+1, ic+1, mod(bc - br, nelem)+1]
        abs(vals[q] - c) <= tol || return nothing
        abs(c) > tol && (seen[br+1] += 1)
    end
    all(==(nnzC), seen) || return nothing

    # one Hermitian nb×nb eigenproblem per wavenumber
    nk = nelem ÷ 2 + 1
    U  = zeros(Complex{T}, nb, nb, nk)
    Λ  = zeros(T, nb*nk)
    Bk = zeros(Complex{T}, nb, nb)
    for k = 0:nk-1
        fill!(Bk, zero(Complex{T}))
        for j = 0:nelem-1
            Cj = view(C, :, :, j+1)
            any(!iszero, Cj) || continue
            Bk .+= Cj .* cispi(T(2*mod(j*k, nelem))/nelem)
        end
        F = eigen(Hermitian(Bk))
        U[:, :, k+1] .= F.vectors
        Λ[k*nb+1:(k+1)*nb] .= F.values
    end

    X  = zeros(T, nb, nelem)
    Xh = zeros(Complex{T}, nb, nk)
    Yh = zeros(Complex{T}, nb, nk)
    pf = FFTW.plan_rfft(X, 2; flags = FFTW.ESTIMATE)
    pb = FFTW.plan_brfft(Xh, nelem, 2; flags = FFTW.ESTIMATE)
    return SEMBloch{T,typeof(pf),typeof(pb)}(perm, s[perm], sinv[perm], U, Λ, nb, nelem, nk,
                                             npoin, neqs, X, Xh, Yh, pf, pb)
end

"""
    build_split_linear_operator(params, inputs) -> SEMLinearOperator1D

Assemble L = c2 D2 + c4 D2² (block-diagonal over the equations, in Jexpresso's
column-major (npoin, neqs) layout of u) and the eigendecomposition of B: its
Bloch form on a uniform periodic mesh, the dense one otherwise. Built once
per run and kept in params.split_cache: the precompile warm-up, the
integrator warm-up and the real solve all reuse it. See the header of
split_linear.jl.
"""
function build_split_linear_operator(params, inputs)

    params.split_cache[] isa SEMLinearOperator1D && return params.split_cache[]

    params.SD isa NSD_1D ||
        error(":lsplit_linear is implemented for 1D only (this case is $(typeof(params.SD))).")
    params.AD isa ContGal ||
        error(":lsplit_linear needs the continuous Galerkin discretization (this case is $(typeof(params.AD))).")
    inputs[:backend] == CPU() ||
        error(":lsplit_linear runs on the CPU backend only.")
    inputs[:lperiodic_1d] == true ||
        error(":lsplit_linear needs :lperiodic_1d => true: the linear operator carries no boundary conditions.")
    MPI.Comm_size(get_mpi_comm()) == 1 ||
        error(":lsplit_linear runs on one MPI rank (the 1D path is serial).")

    T     = Float64
    m     = params.mesh1d
    ω     = params.ω
    dψ    = m.dψ
    npoin = m.npoin
    nelem = Int(params.mesh.nelem)
    ngl   = Int(params.mesh.ngl)
    neqs  = Int(params.neqs)
    Minv  = params.Minv
    Minv isa AbstractVector ||
        error(":lsplit_linear needs the lumped (diagonal) mass matrix.")

    c2 = _split_coeffs(inputs, :split_d2, neqs)
    c4 = _split_coeffs(inputs, :split_d4, neqs)
    (length(c2) == neqs && length(c4) == neqs) ||
        error(":split_d2 and :split_d4 need one coefficient per equation ($neqs).")

    # K: the weak Laplacian, element by element exactly as _expansion_visc!
    # applies it (LGL quadrature at the nodes, ωJac·(dξ/dx)² weight)
    nnzK = nelem*ngl*ngl*ngl
    Ki = Vector{Int}(undef, nnzK); Kj = Vector{Int}(undef, nnzK); Kv = Vector{T}(undef, nnzK)
    n = 0
    @inbounds for iel = 1:nelem, k = 1:ngl
        wk = ω[k]*m.Je[iel, k]*m.dξdx[iel, k]^2
        for j = 1:ngl, i = 1:ngl
            n += 1
            Ki[n] = m.connijk[iel, i, 1, 1]
            Kj[n] = m.connijk[iel, j, 1, 1]
            Kv[n] = dψ[i, k]*dψ[j, k]*wk
        end
    end
    K  = sparse(Ki, Kj, Kv, npoin, npoin)          # duplicates are summed: the DSS
    D2 = -(Diagonal(Minv)*K)

    # symmetric similarity: B = S⁻¹ K S⁻¹, S = M^{1/2}; its eigenbasis by
    # Bloch waves on a uniform periodic mesh, dense otherwise
    s    = sqrt.(one(T) ./ Minv)
    sinv = sqrt.(Minv)
    Bs   = (K .* sinv) .* transpose(sinv)
    E    = _sem_bloch(Bs, s, sinv, m.connijk, nelem, ngl, npoin, neqs)
    if E === nothing
        F = eigen(Symmetric(Matrix(Bs)))
        E = SEMEigen{T}(Matrix(F.vectors), F.values, s, sinv, npoin, neqs)
    end

    nλ = length(E.Λ)
    λ  = zeros(T, nλ*neqs)
    blocks = SparseMatrixCSC{T,Int}[]
    for ieq = 1:neqs
        λ[(ieq - 1)*nλ + 1 : ieq*nλ] .= -c2[ieq] .* E.Λ .+ c4[ieq] .* E.Λ.^2
        push!(blocks, c2[ieq]*D2 + c4[ieq]*(D2*D2))
    end
    A = blockdiag(blocks...)
    dropzeros!(A)

    S = E isa SEMBloch ? SEMSpectralOperator1D(E, copy(λ)) : nothing
    L = SEMLinearOperator1D{T,typeof(E),typeof(S)}(A, E, λ, one(T), S)
    params.split_cache[] = L
    return L
end

"""
    jexpresso_ode_problem(inputs, params, u, tspan)

The ODE problem of the time loop: `ODEProblem(rhs!, …)`, or, with
:lsplit_linear => true, `SplitODEProblem(L, rhs!, …)` with the SEM linear
operator of build_split_linear_operator. Exponential integrators receive the
operator itself (exact φ-functions); every other split solver (IMEX) the
sparse matrix, which is what its linear solves want. FullSpecialize either
way: see time_loop!.
"""
function jexpresso_ode_problem(inputs, params, u, tspan)
    if get(inputs, :lsplit_linear, false) == true
        L  = build_split_linear_operator(params, inputs)
        Lf = inputs[:ode_solver] isa OrdinaryDiffEq.OrdinaryDiffEqCore.OrdinaryDiffEqExponentialAlgorithm ? L : L.A
        # ODEFunction around each part, as SplitFunction{iip}(f1, f2) does, but
        # FullSpecialize on rhs! too: the default AutoSpecialize wraps it in a
        # FunctionWrapper that erases the type of params (see time_loop!)
        f = SplitFunction{true, FullSpecialize}(ODEFunction(MatrixOperator(Lf)),
                                                ODEFunction{true, FullSpecialize}(rhs!))
        return SplitODEProblem(f, u, tspan, params)
    end
    return ODEProblem{true, FullSpecialize}(rhs!, u, tspan, params)
end

"""
    check_split_time_grid(inputs, tstops)

Exponential integrators build their φ operators once, for one Δt; a step
shortened to land on an output time would silently use the wrong ones. So
every output time must be a whole number of steps from :tinit.
"""
function check_split_time_grid(inputs, tstops)
    get(inputs, :lsplit_linear, false) == true || return nothing
    inputs[:ode_solver] isa OrdinaryDiffEq.OrdinaryDiffEqCore.OrdinaryDiffEqExponentialAlgorithm || return nothing
    Δt = Float64(inputs[:Δt]); t0 = Float64(inputs[:tinit])
    for t in vcat(collect(Float64, tstops), Float64(inputs[:tend]))
        n = (t - t0)/Δt
        abs(n - round(n)) <= 1.0e-9*max(1.0, abs(n)) ||
            error(":lsplit_linear with $(nameof(typeof(inputs[:ode_solver]))): output time $t is not a whole number of steps Δt = $Δt from :tinit = $t0.")
    end
    return nothing
end

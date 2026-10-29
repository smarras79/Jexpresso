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
# so the deck keeps :lvisc => false: the derivatives live in L. L u itself is
# applied with the sparse matrix (14 nonzeros per row at p = 4).
#
# THE φ-FUNCTIONS. Exponential integrators need φₖ(hL). ExponentialUtilities
# computes them for a general dense matrix one column at a time, each column
# an exponential of an (n+k)×(n+k) matrix: seconds for n = 256 with ‖hL‖ ~
# 10³. Here they are exact and cheap, because L is similar to a SYMMETRIC
# matrix: with S = M^{1/2} and B = S⁻¹ K S⁻¹ (symmetric),
#
#     D2 = S⁻¹ (-B) S,    L = S⁻¹ Q diag(λ) Qᵀ S,    λ = -c2 Λ + c4 Λ²,
#
# from one symmetric eigendecomposition B = Q Λ Qᵀ, and then every function
# of hL is S⁻¹ Q diag(f(hλ)) Qᵀ S. Those operators are never formed: they are
# kept as their diagonal f(hλ) and applied as Qᵀ, a scaling and Q, all with
# the ONE matrix Q. Formed, ETDRK4's six coefficient matrices are six distinct
# dense n×n matrices, and their products were memory-bound (measured: 47 μs
# per product at n = 256, 45 % of the run); through the shared Q, which stays
# in cache, each product is two ~9 μs ones and setting the operators up is
# O(n) instead of O(n³). Sums, products, powers and multiplication by L, which
# is the algebra OrdinaryDiffEq's expRK_operators does on them, act on the
# diagonals.
#
# Scope: periodic 1D (:lperiodic_1d => true), CPU, ContGal, one rank. Q is
# dense (one eigensolve, n³, ~40 ms for n = 256), the right trade for 1D
# problems of up to a few thousand nodes; past that, an IMEX solver on the
# sparse L (KenCarp4, which receives L.A) is the scalable choice.
#---------------------------------------------------------------------------------

"""
    SEMEigen{T}

Eigendata shared by the linear operator of a split problem and every function
of it: B = S⁻¹ K S⁻¹ = Q diag(Λ) Qᵀ for one equation's npoin nodes (the same
for all equations), and the diagonal of S = M^{1/2}.
"""
struct SEMEigen{T}
    Q::Matrix{T}
    Λ::Vector{T}
    s::Vector{T}
    sinv::Vector{T}
    npoin::Int
    neqs::Int
end

"""
    SEMLinearOperator1D{T}

The linear part α·L of a split problem (see the header of split_linear.jl).
Products use the sparse `A`. `exp`, `ExponentialUtilities.phi` and products
with functions of L return [`SEMSpectralOperator1D`](@ref)s. `α` is the scalar
factor of the expression this object stands for: `dt*L` is the same data with
α = dt, which is how OrdinaryDiffEq's exponential integrators ask for φ(dt L).
"""
struct SEMLinearOperator1D{T} <: AbstractMatrix{T}
    A::SparseMatrixCSC{T,Int}
    E::SEMEigen{T}
    λ::Vector{T}            # eigenvalues of L, equation by equation
    α::T
end

"""
    SEMSpectralOperator1D{T}

A function of the linear operator, S⁻¹ Q diag(d) Qᵀ S equation by equation,
stored as its diagonal `d` and applied through the shared `Q` (see the header
of split_linear.jl). `w1`, `w2` are its own workspace, so `mul!` allocates
nothing.
"""
struct SEMSpectralOperator1D{T} <: AbstractMatrix{T}
    E::SEMEigen{T}
    d::Vector{T}
    w1::Vector{T}
    w2::Vector{T}
end
SEMSpectralOperator1D(E::SEMEigen{T}, d::Vector{T}) where {T} =
    SEMSpectralOperator1D{T}(E, d, zeros(T, E.npoin), zeros(T, E.npoin))

#------------------------------ the linear operator ------------------------------
Base.size(L::SEMLinearOperator1D) = size(L.A)
Base.getindex(L::SEMLinearOperator1D, i::Int, j::Int) = L.α*L.A[i, j]

LinearAlgebra.mul!(y::AbstractVector, L::SEMLinearOperator1D, x::AbstractVector) =
    mul!(y, L.A, x, L.α, false)
LinearAlgebra.mul!(y::AbstractVector, L::SEMLinearOperator1D, x::AbstractVector, a::Number, b::Number) =
    mul!(y, L.A, x, a*L.α, b)
Base.:*(L::SEMLinearOperator1D, x::AbstractVector) = L.α*(L.A*x)
Base.:*(L::SEMLinearOperator1D, X::AbstractMatrix) = L.α*(L.A*X)
Base.:*(a::Number, L::SEMLinearOperator1D) = SEMLinearOperator1D(L.A, L.E, L.λ, a*L.α)
Base.:*(L::SEMLinearOperator1D, a::Number) = a*L
Base.:/(L::SEMLinearOperator1D, a::Number) = inv(a)*L

Base.exp(L::SEMLinearOperator1D) = SEMSpectralOperator1D(L.E, exp.(L.α .* L.λ))

function ExponentialUtilities.phi(L::SEMLinearOperator1D{T}, k::Integer; kwargs...) where {T}
    z = L.α .* L.λ
    return [SEMSpectralOperator1D(L.E, _phi_scalar.(z, j)) for j = 0:k]
end

#--------------------------- functions of the operator ---------------------------
Base.size(P::SEMSpectralOperator1D) = (P.E.npoin*P.E.neqs, P.E.npoin*P.E.neqs)
function Base.getindex(P::SEMSpectralOperator1D{T}, i::Int, j::Int) where {T}
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

# y = a·P x + b·y, equation by equation: S x, Qᵀ, diag(d), Q, S⁻¹
function _spectral_apply!(y::AbstractVector, P::SEMSpectralOperator1D, x::AbstractVector, a, b)
    E = P.E; n = E.npoin; w1 = P.w1; w2 = P.w2
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
    build_split_linear_operator(params, inputs) -> SEMLinearOperator1D

Assemble L = c2 D2 + c4 D2² (block-diagonal over the equations, in Jexpresso's
column-major (npoin, neqs) layout of u) and its eigendecomposition. Built once
per run and kept in params.split_cache: the precompile warm-up, the integrator
warm-up and the real solve all reuse it. See the header of split_linear.jl.
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

    # symmetric similarity: B = S⁻¹ K S⁻¹ = Q Λ Qᵀ, S = M^{1/2}
    s    = sqrt.(one(T) ./ Minv)
    sinv = sqrt.(Minv)
    F    = eigen(Symmetric((Matrix(K) .* sinv) .* transpose(sinv)))
    E    = SEMEigen{T}(Matrix(F.vectors), F.values, s, sinv, npoin, neqs)

    λ = zeros(T, npoin*neqs)
    blocks = SparseMatrixCSC{T,Int}[]
    for ieq = 1:neqs
        λ[(ieq - 1)*npoin + 1 : ieq*npoin] .= -c2[ieq] .* E.Λ .+ c4[ieq] .* E.Λ.^2
        push!(blocks, c2[ieq]*D2 + c4[ieq]*(D2*D2))
    end
    A = blockdiag(blocks...)
    dropzeros!(A)

    L = SEMLinearOperator1D{T}(A, E, λ, one(T))
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

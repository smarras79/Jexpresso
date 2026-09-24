#=============================================================================
 fourier_collocation.jl — pseudo-spectral (Fourier collocation) solver for the
 periodic Poisson equation

     -∇²u = f      on [x0, x0+Lx) × [y0, y0+Ly),   u periodic,

 on the uniform collocation grid x_i = x0 + i Lx/Nx, i = 0..Nx-1 (and the same
 in y), with Kopriva's Fourier collocation machinery ("Implementing Spectral
 Methods for Partial Differential Equations", ch. 4; the routines are in
 src/kernel/infrastructure/Kopriva_functions.jl):

   • the 1-D Fourier derivative matrix D is FourierDerivativeMatrix(N)
     (Algorithm 18, nodes 2πj/N), scaled by 2π/L for a period L;
   • second derivatives are built from D applied twice, D·D, as Kopriva's
     Fourier collocation diffusion operator does (Algorithm 41), with ONE
     correction on the Nyquist mode (below) that makes the result the exact
     Fourier second-derivative collocation matrix D⁽²⁾;
   • the collocation Laplacian on the tensor grid is  D⁽²⁾ₓ ⊗ I + I ⊗ D⁽²⁾ᵧ.

 Unlike the FFT solver (fft_poisson_core.jl), which works on Fourier
 coefficients, everything here happens in PHYSICAL space with dense derivative
 matrices: this is the pseudo-spectral method, and its cost reflects that.

 DIRECT SOLVE by matrix diagonalisation. D⁽²⁾ₓ is real symmetric (D is
 antisymmetric), so -D⁽²⁾ₓ = Qₓ Λₓ Qₓᵀ with orthogonal Qₓ, and the same in y.
 The collocation system  -(D⁽²⁾ₓ U + U D⁽²⁾ᵧ) = F  (U[i,j] = u(x_i, y_j)) becomes

       Û = Qₓᵀ F Qᵧ ,   Û_ij ← Û_ij / (λˣ_i + λʸ_j) ,   U = Qₓ Û Qᵧᵀ ,

 with λ ≥ 0: four dense N×N matrix products per solve, O(N³), no
 N²×N² matrix ever formed. This is the exact solution of the collocation
 system, not an approximation to it.

 THE NYQUIST MODE AND THE NULL SPACE. For even N, D·D annihilates the
 constant AND the Nyquist mode cos(N x/2): D maps the Nyquist mode to a
 multiple of sin(N x/2), which vanishes at every node. Taken literally, D·D
 would then divide a mode that is Nyquist in x and k_y in y by k_y² instead of
 (N/2)² + k_y², amplifying the Nyquist-line content of f by up to (N/2)²/k_y².
 Measured on a smooth, non-band-limited f at 40×36 (vs the FFT's 9.5e-10):
     D·D taken literally ................ L∞ error 4.3e-6
     Nyquist lines dropped from u ....... L∞ error 3.1e-8   (u's own Nyquist
                                          content is lost)
     Nyquist given its true eigenvalue .. agrees with the FFT to round-off
 so the last is what this solver does:
   • the double eigenvalue 0 of each 1-D D·D is given its EXACT eigenvectors —
     the constant 1/√N and the Nyquist vector (-1)^j/√N — instead of the
     arbitrary basis an eigensolver returns for a degenerate eigenvalue;
   • the Nyquist vector gets its true eigenvalue (πN/L)² of -d²/dx²
     (cos(N x/2) at the nodes), which turns D·D into D⁽²⁾;
   • the constant×constant mode is the only null mode left; it is set to zero,
     which fixes the ZERO-MEAN solution (the gauge of the FFT and the periodic
     SEM solve) and drops the mean of f, which has no periodic solution
     (reported in `S.fmean[]`).
 The pseudo-spectral and the FFT solutions are then the same discrete
 solution (the Fourier interpolant of f, differentiated exactly); what differs
 is how it is computed: dense physical-space matrices here, O(N³) per solve,
 against the FFT's O(N² log N).

 Only 2-D grids with EVEN Nx, Ny (Algorithm 18). This file depends on
 LinearAlgebra and on FourierDerivativeMatrix only, so it is unit-tested
 standalone by test/poisson_periodic/test_pseudospectral_poisson.jl.
=============================================================================#

"""
    FourierCollocationPoissonSolver((Nx, Ny), (Lx, Ly); derivative_matrix = FourierDerivativeMatrix)

Eigen-decompositions of the 1-D collocation operators -D⁽²⁾ (from Kopriva's
Fourier derivative matrix D, squared, with the Nyquist correction described in
the file header) and the work arrays for the periodic Poisson
problem on an Nx × Ny grid of periods (Lx, Ly). Build once (O(N³)), then call
[`fourier_collocation_poisson_solve!`](@ref) as often as needed.
"""
struct FourierCollocationPoissonSolver
    dims  :: NTuple{2, Int}
    Ls    :: NTuple{2, Float64}
    Qx    :: Matrix{Float64}
    Qy    :: Matrix{Float64}
    invλ  :: Matrix{Float64}          # 1/(λˣ_i + λʸ_j); 0 on the mean mode
    W1    :: Matrix{Float64}          # work arrays
    W2    :: Matrix{Float64}
    fmean :: Base.RefValue{Float64}   # mean of the last f (removed by the solve)
end

# -d²/dx² on one periodic axis: eigenpairs of -D⁽²⁾ built from Kopriva's D·D,
# with the constant and the Nyquist eigenvectors set exactly and the Nyquist
# eigenvalue corrected to (πN/L)². `isconst` flags the constant eigenvector.
function _collocation_axis(N::Int, L::Float64, derivative_matrix)
    N >= 2 && iseven(N) ||
        throw(ArgumentError("Fourier collocation needs an even number of points per axis, got $N"))
    L > 0 || throw(ArgumentError("period $L must be positive"))
    D  = derivative_matrix(N) .* (2π / L)       # Kopriva Alg. 18, scaled to period L
    A  = -(D * D)                                 # -D·D, symmetric positive semi-definite
    E  = eigen(Symmetric((A + A') ./ 2))
    Q, λ = E.vectors, E.values
    # the non-zero eigenvalues are ≥ (2π/L)² (mode k = 1); the null ones are round-off
    null = findall(abs.(λ) .< 1e-6 * (2π / L)^2)
    length(null) == 2 ||
        error("Fourier collocation: expected a 2-dimensional null space of D·D, found $(length(null))")
    Q[:, null[1]] .= 1 / sqrt(N)                            # constant
    Q[:, null[2]] .= [(-1)^j / sqrt(N) for j in 0:N-1]      # Nyquist: cos(N x/2) at the nodes
    λ[null[1]] = 0.0
    λ[null[2]] = (π * N / L)^2                              # its true -d²/dx² eigenvalue
    isconst = falses(N)
    isconst[null[1]] = true
    return Q, λ, isconst
end

function FourierCollocationPoissonSolver(dims::NTuple{2, Integer}, Ls::NTuple{2, Real};
                                         derivative_matrix = FourierDerivativeMatrix)
    Nx, Ny = Int(dims[1]), Int(dims[2])
    Lx, Ly = Float64(Ls[1]), Float64(Ls[2])
    Qx, λx, cx = _collocation_axis(Nx, Lx, derivative_matrix)
    Qy, λy, cy = _collocation_axis(Ny, Ly, derivative_matrix)
    invλ = Matrix{Float64}(undef, Nx, Ny)
    @inbounds for j = 1:Ny, i = 1:Nx
        invλ[i, j] = (cx[i] && cy[j]) ? 0.0 : 1.0 / (λx[i] + λy[j])
    end
    return FourierCollocationPoissonSolver((Nx, Ny), (Lx, Ly), Qx, Qy, invλ,
                                           zeros(Nx, Ny), zeros(Nx, Ny), Ref(0.0))
end

"""
    fourier_collocation_poisson_solve!(u, S, f) -> u

Overwrite the Nx × Ny array `u` with the zero-mean pseudo-spectral solution of
-∇²u = f (see the file header). The mean removed from `f` is stored in
`S.fmean[]`. Allocation-free; `u` and `f` may alias.
"""
function fourier_collocation_poisson_solve!(u::AbstractMatrix{Float64},
                                            S::FourierCollocationPoissonSolver,
                                            f::AbstractMatrix{Float64})
    size(f) == S.dims || throw(DimensionMismatch("f has size $(size(f)), solver expects $(S.dims)"))
    size(u) == S.dims || throw(DimensionMismatch("u has size $(size(u)), solver expects $(S.dims)"))
    S.fmean[] = sum(f) / length(f)
    mul!(S.W1, transpose(S.Qx), f)          # Qₓᵀ F
    mul!(S.W2, S.W1, S.Qy)                  #      … Qᵧ
    S.W2 .*= S.invλ                         # divide by the eigenvalues
    mul!(S.W1, S.Qx, S.W2)                  # Qₓ Û
    mul!(u, S.W1, transpose(S.Qy))          #      … Qᵧᵀ
    return u
end

"""
    fourier_collocation_poisson_solve(f, Ls) -> u

One-shot wrapper: build the solver for `size(f)` and `Ls`, solve once.
"""
function fourier_collocation_poisson_solve(f::AbstractMatrix{<:Real}, Ls::NTuple{2, Real})
    S = FourierCollocationPoissonSolver(size(f), Ls)
    return fourier_collocation_poisson_solve!(Matrix{Float64}(undef, size(f)), S, Matrix{Float64}(f))
end

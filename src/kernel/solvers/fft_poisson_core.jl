#=============================================================================
 fft_poisson_core.jl — FFTW-based solver for the periodic Poisson equation

     -∇²u = f      on the box  Π_d [x0_d, x0_d + L_d),   u periodic,

 in 1, 2 or 3 dimensions, on the uniform grid  x_d,i = x0_d + (i-1) L_d / n_d,
 i = 1..n_d (the periodic image i = n_d+1 ≡ 1 is NOT stored).

 The Fourier modes exp(i k·x) diagonalise every translation-invariant operator
 on that grid, so the solve is three passes and no iteration:

     f̂ = rfft(f)             real-to-complex forward transform  (FFTW)
     û = f̂ / λ(k)            pointwise, λ = eigenvalue of the discrete -∇²
     u = brfft(û) / Π n_d    complex-to-real backward transform  (FFTW)

 THE NULL SPACE. On a periodic domain the constants are in the kernel of -∇²
 (λ(0) = 0), so a solution exists only if ∫f = 0 and is then unique only up to
 a constant. The solver
   • fixes the constant by returning the ZERO-MEAN solution (û(0) = 0), and
   • solves the PROJECTED problem  -∇²u = f - mean(f)  when f has a non-zero
     mean. The removed mean is kept in `S.fmean[]` so the caller can decide
     whether it was round-off or a genuinely incompatible right-hand side.
 Laplace's equation (f ≡ 0) therefore returns u ≡ 0, the only zero-mean
 periodic harmonic function.

 TWO DISCRETISATIONS of -∇², selected by `laplacian`:

   :spectral   λ = Σ_d k_d²                          (Fourier–Galerkin)
               Exact for every mode the grid resolves: a band-limited f gives
               u to round-off, a smooth periodic f converges exponentially.

   :fd2        λ = Σ_d (4/h_d²) sin²(k_d h_d / 2)     (second-order centred
               finite differences, the 3/5/7-point stencil). The FFT then
               returns, to round-off, the solution of the SPARSE linear system
               A u = f - mean(f) with A the periodic FD Laplacian. This is the
               mode an iterative / multigrid solver of the same system is to be
               compared against: the two must agree to solver tolerance, not
               merely to truncation error.

 The wavenumbers are the signed ones, k = 2π m / L with m ∈ (-n/2, n/2]. For
 even n the Nyquist mode m = n/2 is its own alias; λ depends on k² only, so
 the sign ambiguity is harmless.

 Any grid size works (FFTW handles all n; powers of 2, 3, 5, 7 are fastest).

 USAGE
     S = FFTPoissonSolver((nx, ny), (Lx, Ly); laplacian = :spectral)
     u = similar(f)
     fft_poisson_solve!(u, S, f)          # allocation-free, plan reused
     S.fmean[]                            # mean removed from f
 or, one-shot,
     u = fft_poisson_solve(f, (Lx, Ly); laplacian = :fd2)

 This file depends on FFTW and LinearAlgebra only (no Jexpresso types), so it
 is unit-tested standalone by test/poisson_periodic/test_fft_poisson.jl. The
 Jexpresso driver that feeds it from a case deck is fft_laplace.jl.
=============================================================================#

const FFT_POISSON_LAPLACIANS = (:spectral, :fd2)

"""
    FFTPoissonSolver(dims, Ls; laplacian = :spectral, flags = FFTW.ESTIMATE)

Precomputed FFTW plans, work array and inverse eigenvalues for the periodic
Poisson problem `-∇²u = f` on an `dims[1] × … × dims[ND]` grid of periods `Ls`.
Build it once, then call [`fft_poisson_solve!`](@ref) as many times as needed.

`flags` is passed to FFTW's planner: `FFTW.ESTIMATE` plans instantly;
`FFTW.MEASURE` spends time planning to run faster when the solver is reused.
"""
struct FFTPoissonSolver{ND, TF, TB}
    dims      :: NTuple{ND, Int}
    Ls        :: NTuple{ND, Float64}
    laplacian :: Symbol
    fhat      :: Array{ComplexF64, ND}   # half-spectrum work array (rfft layout)
    invλ      :: Array{Float64, ND}      # 1/(λ Π n_d); 0 on the null mode
    fwd       :: TF                      # rfft  plan
    bwd       :: TB                      # brfft plan (unnormalised inverse)
    fmean     :: Base.RefValue{Float64}  # mean of the last f (removed by the solve)
end

function FFTPoissonSolver(dims::NTuple{ND, Integer}, Ls::NTuple{ND, Real};
                          laplacian::Symbol = :spectral,
                          flags::UInt32 = FFTW.ESTIMATE) where {ND}
    laplacian in FFT_POISSON_LAPLACIANS ||
        throw(ArgumentError("laplacian = :$laplacian; expected one of $(FFT_POISSON_LAPLACIANS)"))
    all(>(0), dims) || throw(ArgumentError("grid size $dims must be positive"))
    all(>(0), Ls)   || throw(ArgumentError("periods $Ls must be positive"))

    n  = ntuple(d -> Int(dims[d]), ND)
    L  = ntuple(d -> Float64(Ls[d]), ND)
    nh = ntuple(d -> d == 1 ? n[1] ÷ 2 + 1 : n[d], ND)      # rfft halves dim 1

    # Plan on scratch arrays: FFTW.MEASURE overwrites the arrays it plans on.
    fhat = zeros(ComplexF64, nh)
    fwd  = FFTW.plan_rfft(zeros(Float64, n); flags = flags)
    bwd  = FFTW.plan_brfft(fhat, n[1]; flags = flags)

    # Per-axis eigenvalue factors on the rfft index layout.
    λ1d = ntuple(ND) do d
        m = d == 1 ? (0:nh[1]-1) : [j <= n[d] ÷ 2 ? j : j - n[d] for j in 0:n[d]-1]
        k = (2π / L[d]) .* m
        if laplacian === :spectral
            k .^ 2
        else
            h = L[d] / n[d]
            (4 / h^2) .* sin.(k .* (h / 2)) .^ 2
        end
    end

    scale = 1.0 / prod(n)             # brfft is unnormalised
    invλ  = Array{Float64}(undef, nh)
    @inbounds for I in CartesianIndices(invλ)
        λ = 0.0
        for d = 1:ND
            λ += λ1d[d][I[d]]
        end
        # λ vanishes only on the constant mode (every other k has k_d ≠ 0 for
        # some d, and 0 < |k_d h_d/2| ≤ π/2 keeps the fd2 factor positive).
        invλ[I] = I == first(CartesianIndices(invλ)) ? 0.0 : scale / λ
    end

    return FFTPoissonSolver{ND, typeof(fwd), typeof(bwd)}(n, L, laplacian, fhat, invλ,
                                                          fwd, bwd, Ref(0.0))
end

"""
    fft_poisson_solve!(u, S::FFTPoissonSolver, f) -> u

Overwrite `u` with the zero-mean periodic solution of `-∇²u = f - mean(f)`.
`u` and `f` are `Float64` arrays of size `S.dims` (they may alias). The mean
that was removed from `f` is stored in `S.fmean[]`. Allocation-free.
"""
function fft_poisson_solve!(u::Array{Float64, ND}, S::FFTPoissonSolver{ND},
                            f::Array{Float64, ND}) where {ND}
    size(f) == S.dims || throw(DimensionMismatch("f has size $(size(f)), solver expects $(S.dims)"))
    size(u) == S.dims || throw(DimensionMismatch("u has size $(size(u)), solver expects $(S.dims)"))
    mul!(S.fhat, S.fwd, f)
    S.fmean[] = real(S.fhat[1]) / length(f)
    S.fhat .*= S.invλ
    mul!(u, S.bwd, S.fhat)            # brfft may destroy fhat; it is rebuilt every call
    return u
end

"""
    fft_poisson_solve(f, Ls; laplacian = :spectral) -> u

One-shot convenience wrapper: builds an [`FFTPoissonSolver`](@ref) for
`size(f)` and `Ls`, solves once and returns the zero-mean solution.
"""
function fft_poisson_solve(f::AbstractArray{<:Real, ND}, Ls::NTuple{ND, Real};
                           laplacian::Symbol = :spectral) where {ND}
    S = FFTPoissonSolver(size(f), Ls; laplacian = laplacian)
    return fft_poisson_solve!(Array{Float64, ND}(undef, size(f)), S, Array{Float64, ND}(f))
end

"""
    periodic_grid_lines(dims, Ls, x0s) -> NTuple of coordinate vectors

The uniform periodic grid the solver works on: `x0_d + (i-1) L_d / n_d`,
`i = 1..n_d` (the right-hand periodic image is not included).
"""
periodic_grid_lines(dims::NTuple{ND, Integer}, Ls::NTuple{ND, Real},
                    x0s::NTuple{ND, Real} = ntuple(_ -> 0.0, ND)) where {ND} =
    ntuple(d -> [Float64(x0s[d]) + (i - 1) * Float64(Ls[d]) / dims[d] for i in 1:dims[d]], ND)

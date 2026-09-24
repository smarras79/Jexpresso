#=============================================================================
 test/poisson_periodic/test_pseudospectral_poisson.jl — the pseudo-spectral
 (Fourier collocation) periodic Poisson solver
 (src/kernel/solvers/fourier_collocation.jl), built on Kopriva's
 FourierDerivativeMatrix (src/kernel/infrastructure/Kopriva_functions.jl).

     julia --project=test/poisson_periodic -e 'using Pkg; Pkg.instantiate()'
     julia --project=test/poisson_periodic test/poisson_periodic/test_pseudospectral_poisson.jl

 NO `using Jexpresso`: Kopriva_functions.jl is included directly; it needs
 only KernelAbstractions (for its one GPU kernel) and the TInt/TFloat aliases
 that Jexpresso normally defines.

 WHAT IS CHECKED
   1. Kopriva's Fourier derivative matrix itself: exact on trigonometric
      polynomials below the Nyquist mode, antisymmetric, kills the Nyquist
      mode (which is why the solver corrects D·D on that mode).
   2. The solver-agnostic benchmark (benchmark.jl), on the grids a 2-D,
      even-size solver supports: Laplace equation, linearity, incompatible
      RHS, band-limited exactness, exponential convergence.
   3. Agreement with the FFTW solver (same Fourier interpolation space, so
      the two must agree to round-off on resolved data).
   4. Allocation-free reuse, the reported mean, argument checks.
=============================================================================#
using Test, LinearAlgebra
using KernelAbstractions
import FFTW

const TInt   = Int64          # aliases Jexpresso defines before including Kopriva_functions.jl
const TFloat = Float64
include(joinpath(@__DIR__, "..", "..", "src", "kernel", "infrastructure", "Kopriva_functions.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "fourier_collocation.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "fft_poisson_core.jl"))
include(joinpath(@__DIR__, "benchmark.jl"))
using .PeriodicPoissonBenchmark

ps_solve(f, Ls)  = fourier_collocation_poisson_solve(f, Ls)
fft_solve(f, Ls) = fft_poisson_solve(f, Ls)
ps_supports(dims) = length(dims) == 2 && all(iseven, dims)

@testset verbose = true "Pseudo-spectral (Fourier collocation) periodic Poisson solver" begin

    @testset "Kopriva FourierDerivativeMatrix (Algorithm 18)" begin
        for N in (8, 16, 30)
            D = FourierDerivativeMatrix(N)
            x = [2π * j / N for j in 0:N-1]
            @test maximum(abs, D .+ D') < 1e-13                    # antisymmetric
            for k in 1:(N ÷ 2 - 1)                                  # exact below Nyquist
                @test maximum(abs, D * sin.(k .* x) .- k .* cos.(k .* x)) < 1e-11 * k
            end
            @test maximum(abs, D * cos.((N ÷ 2) .* x)) < 1e-12      # Nyquist is annihilated
        end
        @test_throws ArgumentError FourierDerivativeMatrix(9)
    end

    verify_periodic_poisson_solver("pseudo-spectral", ps_solve; supports = ps_supports)

    @testset "comparison harness: pseudo-spectral vs FFTW" begin
        compare_periodic_poisson_solvers("pseudo-spectral", ps_solve, "FFTW", fft_solve;
            cases = (("modes2d", (64, 64)), ("aniso2d", (96, 48)),
                     ("smooth2d", (64, 64)), ("smooth2d", (40, 36))),
            rtol = 1e-10, reps = 2)
    end

    @testset "plan reuse is allocation-free and repeatable" begin
        p = problem("aniso2d")
        f, uex = sample_problem(p, (64, 48))
        S = FourierCollocationPoissonSolver(size(f), p.Ls)
        u = similar(f)
        fourier_collocation_poisson_solve!(u, S, f)
        u1 = copy(u)
        fourier_collocation_poisson_solve!(u, S, f)
        @test u == u1
        @test (@allocated fourier_collocation_poisson_solve!(u, S, f)) == 0
        @test error_norms(u, uex).linf < 1e-12
    end

    @testset "reported RHS mean, zero-mean answer" begin
        p = problem("modes2d")
        f, _ = sample_problem(p, (32, 32))
        S = FourierCollocationPoissonSolver(size(f), p.Ls)
        u = similar(f)
        fourier_collocation_poisson_solve!(u, S, f);         @test abs(S.fmean[]) < 1e-13
        fourier_collocation_poisson_solve!(u, S, f .+ 0.75); @test S.fmean[] ≈ 0.75
        @test abs(sum(u)) / length(u) < 1e-13
    end

    @testset "argument checks" begin
        @test_throws ArgumentError FourierCollocationPoissonSolver((9, 8), (1.0, 1.0))
        @test_throws ArgumentError FourierCollocationPoissonSolver((8, 8), (1.0, -1.0))
        S = FourierCollocationPoissonSolver((8, 8), (1.0, 1.0))
        @test_throws DimensionMismatch fourier_collocation_poisson_solve!(zeros(8, 8), S, zeros(8, 10))
    end
end

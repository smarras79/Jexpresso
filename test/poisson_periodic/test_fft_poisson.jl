#=============================================================================
 test/poisson_periodic/test_fft_poisson.jl — the FFTW periodic Poisson solver
 (src/kernel/solvers/fft_poisson_core.jl).

     julia --project=test/poisson_periodic -e 'using Pkg; Pkg.instantiate()'
     julia --project=test/poisson_periodic test/poisson_periodic/test_fft_poisson.jl

 NO `using Jexpresso`: the solver core depends on FFTW and LinearAlgebra only,
 so it is included directly and tested in seconds. The small environment in
 this directory is all it needs.

 WHAT IS CHECKED
   1. The solver-agnostic benchmark (benchmark.jl): the Laplace equation,
      linearity, the null space / incompatible RHS, band-limited exactness to
      round-off and exponential convergence for a non-band-limited solution.
   2. The comparison harness: the FFT vs an independent dense-DFT spectral
      solve (reference_spectral_solve) — same mathematics, no FFT library, so
      they must agree to round-off.
   3. The nodal scoring path (nodal_error_norms) that the SEM-mesh solvers —
      direct SEM, AlgebraicMultigrid on the SEM system, static condensation —
      will be scored with: on the FFT grid it must reduce to error_norms.
   4. FFT-specific behaviour: plan reuse without allocation, the reported RHS
      mean, the 1D case against a closed form, argument checks.
=============================================================================#
using Test, LinearAlgebra
import FFTW

include(joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "fft_poisson_core.jl"))
include(joinpath(@__DIR__, "benchmark.jl"))
using .PeriodicPoissonBenchmark

# The benchmark's calling convention: solve(f, Ls) -> u
fft_solve(f, Ls) = fft_poisson_solve(f, Ls)

@testset verbose = true "FFTW periodic Poisson solver" begin

    verify_periodic_poisson_solver("FFTW", fft_solve)

    @testset "comparison harness: FFTW vs dense-DFT reference" begin
        compare_periodic_poisson_solvers("FFTW", fft_solve,
                                         "dense DFT", reference_spectral_solve;
                                         rtol = 1e-11, reps = 2)
    end

    @testset "nodal error norms agree with the grid ones (SEM comparison path)" begin
        # The SEM solvers of step 3 are scored through nodal_error_norms; on
        # the uniform grid with uniform weights it must reduce to error_norms.
        for (pname, dims) in (("smooth2d", (40, 36)), ("modes3d", (10, 12, 8)))
            p = problem(pname)
            f, uex = sample_problem(p, dims)
            u  = fft_solve(f, p.Ls)
            eg = error_norms(u, uex)
            en = nodal_error_norms(p, grid_points(p, dims), vec(u))
            @test en.linf  ≈ eg.linf  rtol = 1e-10
            @test en.l2rel ≈ eg.l2rel rtol = 1e-10
            # the gauge is weight-consistent: a constant shift of u changes
            # nothing (up to the round-off of adding and removing the constant)
            @test abs(nodal_error_norms(p, grid_points(p, dims), vec(u) .+ 7.0;
                                        w = fill(0.3, length(u))).linf - eg.linf) < 1e-12
        end
    end

    @testset "plan reuse is allocation-free and repeatable" begin
        p = problem("aniso2d")
        f, uex = sample_problem(p, (64, 48))
        S = FFTPoissonSolver(size(f), p.Ls; flags = FFTW.MEASURE)
        u = similar(f)
        fft_poisson_solve!(u, S, f)
        u1 = copy(u)
        fft_poisson_solve!(u, S, f)
        @test u == u1
        @test (@allocated fft_poisson_solve!(u, S, f)) == 0
        @test error_norms(u, uex).linf < 1e-12
    end

    @testset "reported RHS mean" begin
        p = problem("modes2d")
        f, _ = sample_problem(p, (32, 32))
        S = FFTPoissonSolver(size(f), p.Ls)
        u = similar(f)
        fft_poisson_solve!(u, S, f);          @test abs(S.fmean[]) < 1e-13
        fft_poisson_solve!(u, S, f .+ 0.75);  @test S.fmean[] ≈ 0.75
        @test abs(sum(u)) / length(u) < 1e-14  # zero-mean answer either way
    end

    @testset "1D: u = sin(3x) + cos(5x) on [0,2π)" begin
        n  = 32
        x  = periodic_grid_lines((n,), (2π,))[1]
        u  = fft_poisson_solve(9sin.(3x) .+ 25cos.(5x), (2π,))
        @test maximum(abs, u .- (sin.(3x) .+ cos.(5x))) < 1e-13
    end

    @testset "argument checks" begin
        @test_throws ArgumentError FFTPoissonSolver((8, 0), (1.0, 1.0))
        @test_throws ArgumentError FFTPoissonSolver((8, 8), (1.0, -1.0))
        S = FFTPoissonSolver((8, 8), (1.0, 1.0))
        @test_throws DimensionMismatch fft_poisson_solve!(zeros(8, 8), S, zeros(8, 9))
    end
end

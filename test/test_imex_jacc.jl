# Standalone correctness test for the IMEX JACC linear-algebra kernels
# (src/kernel/solvers/imex_jacc.jl). It is self-contained: the file only
# depends on JACC, SparseArrays and LinearAlgebra, so it can be checked without
# building the whole Jexpresso module.
#
# Run from the repo root with the project active:
#
#     julia --project test/test_imex_jacc.jl
#
# On a machine with a JACC GPU backend configured (CUDA/AMDGPU), the same test
# exercises the device path because `JACC.Array` then returns a device array.

using Test
using SparseArrays
using LinearAlgebra

const PROJECT_ROOT = dirname(@__DIR__)
include(joinpath(PROJECT_ROOT, "src", "kernel", "solvers", "imex_jacc.jl"))

@testset "IMEX JACC kernels" begin

    @testset "jacc_spmv! matches sparse mul!" begin
        for n in (16, 257)
            A = sprand(n, n, 0.2) + 2.0 * I        # well-conditioned, nonsymmetric
            x = rand(n)
            Aj = JaccSparseCSR(A)
            xj = JACC.Array(copy(x))
            yj = JACC.Array(zeros(n))
            jacc_spmv!(yj, Aj, xj)
            @test Array(yj) ≈ A * x rtol = 1e-12
        end
    end

    @testset "jacc_bicgstab! solves a Helmholtz-like operator" begin
        # Mimic the IMEX stage operator A = I - λ L with L = -μ (Minv .* K):
        # K symmetric positive (1D Laplacian), Minv a positive diagonal scaling,
        # so A is nonsymmetric but has a real, positive spectrum.
        n  = 400
        K  = spdiagm(-1 => fill(-1.0, n-1), 0 => fill(2.0, n), 1 => fill(-1.0, n-1))
        Minv = 0.5 .+ rand(n)                       # positive diagonal
        μ  = 1.0e-2
        λ  = 1.0e-3
        L  = -μ * (Minv .* K)                       # row scaling -> nonsymmetric
        A  = sparse(1.0I, n, n) - λ * L

        xtrue = rand(n)
        b     = A * xtrue

        Aj = JaccSparseCSR(A)
        bj = JACC.Array(copy(b))
        xj = JACC.Array(zeros(n))
        work = ntuple(_ -> JACC.Array(zeros(n)), 6)

        conv, iters, resnorm = jacc_bicgstab!(xj, Aj, bj, work;
                                              rtol = 1e-10, atol = 1e-14, itmax = 1000)
        @test conv
        @test Array(xj) ≈ xtrue rtol = 1e-7
        @test resnorm ≤ 1e-8 * norm(b)
    end

    @testset "jacc_bicgstab! handles a zero RHS" begin
        n  = 32
        A  = sprand(n, n, 0.3) + 2.0 * I
        Aj = JaccSparseCSR(A)
        bj = JACC.Array(zeros(n))
        xj = JACC.Array(ones(n))
        work = ntuple(_ -> JACC.Array(zeros(n)), 6)
        conv, iters, _ = jacc_bicgstab!(xj, Aj, bj, work)
        @test conv
        @test iters == 0
        @test all(Array(xj) .== 0.0)
    end
end

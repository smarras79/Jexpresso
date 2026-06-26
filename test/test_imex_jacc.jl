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

    @testset "jacc_spmv! is allocation-free on the CPU path" begin
        # Regression guard: the CPU Vector dispatch must NOT go through
        # JACC.parallel_for, whose threaded/serial CPU backend allocates per call
        # (152 KB/call at 1 thread) and doubled this case's allocation.
        n  = 1000
        A  = sprand(n, n, 0.05) + 2.0 * I
        Aj = JaccSparseCSR(A)
        x  = JACC.Array(rand(n))
        y  = JACC.Array(zeros(n))
        jacc_spmv!(y, Aj, x)                       # warm up / compile
        if y isa Vector                            # CPU JACC backend only
            @test (@allocated jacc_spmv!(y, Aj, x)) == 0
        end
    end

    @testset "jacc_bicgstab! reproduces the sparse-LU solve (CPU path parity)" begin
        # The constant-operator CPU path (_imex_rk_run_const!) solves each stage
        # with a cached sparse LU (lu/ldiv!), i.e. exactly `A \ b`. The JACC path
        # (_imex_rk_run_const_jacc!) replaces that with jacc_bicgstab!. This
        # checks the two agree to ~machine precision on the real stage operator
        # A = I - λL, so swapping the solver does not change the stage solution.
        n  = 2000
        K  = spdiagm(-1 => fill(-1.0, n-1), 0 => fill(2.0, n), 1 => fill(-1.0, n-1))
        Minv = 0.5 .+ rand(n)
        μ, λ = 1.0e-2, 1.0e-3
        L  = -μ * (Minv .* K)
        A  = sparse(1.0I, n, n) - λ * L
        b  = A * rand(n)

        x_direct = A \ b                       # what lu/ldiv! computes on the CPU path

        Aj = JaccSparseCSR(A)
        xj = JACC.Array(zeros(n))
        work = ntuple(_ -> JACC.Array(zeros(n)), 6)
        conv, _, _ = jacc_bicgstab!(xj, Aj, JACC.Array(copy(b)), work;
                                    rtol = 1e-10, atol = 1e-14, itmax = 1000)
        @test conv
        @test norm(Array(xj) .- x_direct) ≤ 1e-9 * norm(x_direct)
    end

    @testset "offload round-trip: host rhs -> device solve -> host (parity)" begin
        # Exercises exactly what _imex_rk_run_const_jacc_offload! does per stage:
        # assemble the operator + rhs on the host, upload to device arrays
        # (JACC.Array), solve with jacc_bicgstab!, copy the solution back. On the
        # JACC CPU backend the "device" arrays are host Vectors, so this both
        # validates the mechanics and matches the device path bit-for-bit.
        n  = 1500
        K  = spdiagm(-1 => fill(-1.0, n-1), 0 => fill(2.0, n), 1 => fill(-1.0, n-1))
        Minv = 0.5 .+ rand(n)
        μ, λ = 1.0e-2, 1.0e-3
        A  = sparse(1.0I, n, n) - λ * (-μ * (Minv .* K))

        rhs_h = A * rand(n)                       # host rhs
        x_ref = A \ rhs_h                         # host LU reference

        A_d   = JaccSparseCSR(A)                  # operator uploaded once
        b_d   = JACC.Array(zeros(n))
        x_d   = JACC.Array(zeros(n))
        work  = ntuple(_ -> JACC.Array(zeros(n)), 6)
        x_h   = zeros(n)

        copyto!(b_d, rhs_h)                        # host -> device
        conv, _, _ = jacc_bicgstab!(x_d, A_d, b_d, work; rtol = 1e-10, atol = 1e-14, itmax = 1000)
        copyto!(x_h, x_d)                          # device -> host
        @test conv
        @test norm(x_h .- x_ref) ≤ 1e-9 * norm(x_ref)
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

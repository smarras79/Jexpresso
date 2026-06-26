# GPU equivalent of test/test_imex_jacc.jl.
#
# Same correctness checks for the IMEX JACC linear-algebra kernels
# (src/kernel/solvers/imex_jacc.jl), but exercised on an NVIDIA GPU by
# switching JACC to its CUDA backend before the kernels run. When no CUDA
# device / package is available the whole test set is skipped (marked broken)
# so it is safe to include from CI on a CPU-only machine.
#
# Run on a GPU box with the project active:
#
#     julia --project test/test_imex_jacc_gpu.jl
#
# For AMD GPUs, swap `using CUDA` for `using AMDGPU` and
# `JACC.set_backend("cuda")` for `JACC.set_backend("amdgpu")`.

using Test
using SparseArrays
using LinearAlgebra

const PROJECT_ROOT_GPU = dirname(@__DIR__)

# Bring up a CUDA-backed JACC, or report why we cannot and skip.
const _GPU_READY = Ref(false)
try
    @eval using CUDA
    if CUDA.functional()
        @eval using JACC
        JACC.set_backend("cuda")
        _GPU_READY[] = true
    else
        @info "test_imex_jacc_gpu: CUDA is installed but not functional; skipping GPU tests."
    end
catch err
    @info "test_imex_jacc_gpu: no usable CUDA/JACC GPU backend; skipping GPU tests." exception = err
end

if _GPU_READY[]
    include(joinpath(PROJECT_ROOT_GPU, "src", "kernel", "solvers", "imex_jacc.jl"))

    @testset "IMEX JACC kernels (CUDA)" begin

        @testset "jacc_spmv! matches sparse mul! on GPU" begin
            for n in (16, 257)
                A = sprand(n, n, 0.2) + 2.0 * I
                x = rand(n)
                Aj = JaccSparseCSR(A)                 # nzval/colind/rowptr -> CuArray
                xj = JACC.array(copy(x))              # -> CuArray
                yj = JACC.array(zeros(n))
                jacc_spmv!(yj, Aj, xj)
                @test Array(yj) ≈ A * x rtol = 1e-10
            end
        end

        @testset "jacc_bicgstab! solves I - λL on GPU" begin
            n  = 400
            K  = spdiagm(-1 => fill(-1.0, n-1), 0 => fill(2.0, n), 1 => fill(-1.0, n-1))
            Minv = 0.5 .+ rand(n)
            μ  = 1.0e-2
            λ  = 1.0e-3
            L  = -μ * (Minv .* K)
            A  = sparse(1.0I, n, n) - λ * L

            xtrue = rand(n)
            b     = A * xtrue

            Aj = JaccSparseCSR(A)
            bj = JACC.array(copy(b))
            xj = JACC.array(zeros(n))
            work = ntuple(_ -> JACC.array(zeros(n)), 6)

            conv, iters, resnorm = jacc_bicgstab!(xj, Aj, bj, work;
                                                  rtol = 1e-10, atol = 1e-14, itmax = 1000)
            @test conv
            @test Array(xj) ≈ xtrue rtol = 1e-6
            @test resnorm ≤ 1e-7 * norm(b)
        end
    end
else
    @testset "IMEX JACC kernels (CUDA) [skipped - no GPU]" begin
        @test_broken false
    end
end

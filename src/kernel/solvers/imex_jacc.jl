#----------------------------------------------------------------------------
# imex_jacc.jl
#
# Hardware-agnostic (CPU / NVIDIA / AMD) linear-algebra kernels used by the
# IMEX time integrator, written with JACC.jl so that the *implicit* per-stage
# solve runs on whatever device JACC is configured for.
#
# Why this file exists
# --------------------
# The constant-operator IMEX-RK fast path (`_imex_rk_run_const!` in
# IMEXTimeIntegrators.jl) solves, once per stage,
#
#     (I - λ L) x = b
#
# on the CPU with a cached sparse LU factorization (`lu` / `ldiv!`). Sparse LU
# has no portable GPU equivalent, so it is the single piece that pins that path
# to the host. This file replaces it with:
#
#   * `JaccSparseCSR`     - a CSR sparse operator whose three arrays live on the
#                           JACC device (so the same struct works on CPU and GPU).
#   * `jacc_spmv!`        - a portable sparse mat-vec  y = A x  (one JACC thread
#                           per row), used both for the residual and inside the
#                           iterative solve.
#   * `jacc_bicgstab!`    - a matrix-free BiCGSTAB that only needs `jacc_spmv!`
#                           plus dot / norm / broadcast, all of which are already
#                           defined for the device array type (Array, CuArray,
#                           ROCArray, ...). BiCGSTAB (not CG) because the operator
#                           I - λ L = I + λμ (M⁻¹K) is *not* symmetric (M⁻¹ is a
#                           diagonal scaling of the symmetric stiffness K), even
#                           though its spectrum is real and positive.
#
# The surrounding vector algebra (`axpy!`, broadcasting, `norm`) and the
# explicit RHS (`rhs!` inside `S_fun!`) are intentionally left on the native
# KernelAbstractions backend: they already dispatch on the backend that `u`
# lives on, so on a GPU run they execute on the GPU unchanged. JACC is used
# *only* for the sparse operator and the iterative solve - the parts that the
# direct-LU path could not take to the device.
#
# Selection: a case opts in with `:limex_jacc => true` (and, for an actual GPU
# run, `:backend => CUDABackend()` plus a JACC backend configured for CUDA).
#----------------------------------------------------------------------------

using JACC
using SparseArrays
using LinearAlgebra

#----------------------------------------------------------------------------
# JACC backend selection
#
# Map the Jexpresso KernelAbstractions backend onto a JACC backend. JACC's
# active backend is a global/preference setting; we only *attempt* to switch it
# at run time and swallow any error so that a stock CPU build (the default,
# always-available JACC backend) keeps working without CUDA.jl/AMDGPU.jl
# installed. For a real GPU run the user is expected to have JACC configured
# for the matching device (LocalPreferences.toml or `JACC.set_backend`).
#----------------------------------------------------------------------------
function imex_jacc_init_backend!(inputs)
    be = get(inputs, :backend, nothing)
    name = nothing
    if be !== nothing
        bestr = string(typeof(be))
        if occursin("CUDA", bestr)
            name = "cuda"
        elseif occursin("ROC", bestr) || occursin("AMDGPU", bestr)
            name = "amdgpu"
        end
    end
    if name !== nothing
        try
            JACC.set_backend(name)
        catch err
            @warn "imex_jacc: could not switch JACC backend to $name; using JACC default. ($err)"
        end
    end
    return nothing
end

#----------------------------------------------------------------------------
# CSR sparse operator with device-resident arrays.
#
# Julia's SparseMatrixCSC is column-major (CSC). `transpose(A)` materialized
# with `sparse` gives the CSC of Aᵀ, whose `colptr/rowval/nzval` are exactly the
# CSR `rowptr/colind/nzval` of A - the natural layout for a row-parallel SpMV.
#----------------------------------------------------------------------------
struct JaccSparseCSR{VI, VF}
    n::Int          # number of rows (= columns; operator is square)
    rowptr::VI      # length n+1
    colind::VI      # length nnz
    nzval::VF       # length nnz
end

"""
    JaccSparseCSR(A::SparseMatrixCSC; index_type, float_type)

Build a device-resident CSR copy of `A`. `index_type`/`float_type` let a GPU
run use 32-bit indices/precision; they default to the host matrix's own types.
"""
function JaccSparseCSR(A::SparseMatrixCSC;
                       index_type::Type = Int,
                       float_type::Type = eltype(A))
    n = size(A, 1)
    @assert size(A, 2) == n "JaccSparseCSR expects a square operator"
    At     = sparse(transpose(A))                 # CSC(Aᵀ) == CSR(A)
    rowptr = JACC.Array(convert(Vector{index_type}, At.colptr))
    colind = JACC.Array(convert(Vector{index_type}, At.rowval))
    nzval  = JACC.Array(convert(Vector{float_type}, At.nzval))
    return JaccSparseCSR{typeof(rowptr), typeof(nzval)}(n, rowptr, colind, nzval)
end

# One JACC thread per row: y[i] = Σ_k nzval[k] * x[colind[k]]. Kept a top-level
# function (not a closure) so the GPU backends can compile it.
function _imex_jacc_spmv_row!(i, y, rowptr, colind, nzval, x)
    acc = zero(eltype(y))
    @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
        acc += nzval[k] * x[colind[k]]
    end
    @inbounds y[i] = acc
    return nothing
end

"""
    jacc_spmv!(y, A::JaccSparseCSR, x)

Portable sparse mat-vec `y .= A * x`. `x`, `y` must be device arrays matching the
active JACC backend (the same type as the integrator's `similar(u)` buffers).
"""
function jacc_spmv!(y, A::JaccSparseCSR, x)
    JACC.parallel_for(A.n, _imex_jacc_spmv_row!, y, A.rowptr, A.colind, A.nzval, x)
    return y
end

#----------------------------------------------------------------------------
# Matrix-free BiCGSTAB on the device.
#
# Solves A x = b for a general (non-symmetric) operator with real, positive
# spectrum - exactly the I - λL Helmholtz operator of the IMEX implicit stage.
# Only `jacc_spmv!` plus dot/norm/broadcast are used, so it runs wherever the
# device array type supports those (CPU Array, CuArray, ROCArray, ...).
#
# `x` is overwritten with the solution (started from x = 0). All scratch vectors
# are passed in (`work` = (r, rhat, p, v, s, t)) so the time loop stays
# allocation-light. Returns (converged::Bool, iters::Int, resnorm).
#----------------------------------------------------------------------------
function jacc_bicgstab!(x, A::JaccSparseCSR, b, work;
                        rtol = 1.0e-8, atol = 1.0e-12, itmax = 500)
    r, rhat, p, v, s, t = work
    T = eltype(x)

    fill!(x, zero(T))
    copyto!(r, b)            # r = b - A*0 = b
    copyto!(rhat, r)
    fill!(p, zero(T))
    fill!(v, zero(T))

    bnorm = norm(b)
    if bnorm == zero(T)
        fill!(x, zero(T))
        return (true, 0, zero(real(T)))
    end
    tol = max(atol, rtol * bnorm)

    ρ_old = one(T); α = one(T); ω = one(T)
    resnorm = norm(r)
    iters = 0
    converged = resnorm <= tol

    while !converged && iters < itmax
        iters += 1
        ρ_new = dot(rhat, r)
        if ρ_new == zero(T)
            break                       # breakdown; restart would be needed
        end
        β = (ρ_new / ρ_old) * (α / ω)
        # p = r + β (p - ω v)
        @. p = r + β * (p - ω * v)
        jacc_spmv!(v, A, p)
        α = ρ_new / dot(rhat, v)
        # s = r - α v
        @. s = r - α * v
        resnorm = norm(s)
        if resnorm <= tol
            @. x = x + α * p
            converged = true
            break
        end
        jacc_spmv!(t, A, s)
        tt = dot(t, t)
        ω = (tt == zero(T)) ? zero(T) : dot(t, s) / tt
        # x = x + α p + ω s ;  r = s - ω t
        @. x = x + α * p + ω * s
        @. r = s - ω * t
        resnorm = norm(r)
        converged = resnorm <= tol
        if ω == zero(T)
            break
        end
        ρ_old = ρ_new
    end

    return (converged, iters, resnorm)
end

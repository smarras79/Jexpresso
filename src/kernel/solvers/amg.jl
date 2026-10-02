#=============================================================================
 amg.jl — algebraic multigrid (AlgebraicMultigrid.jl) as a preconditioner of
 the conjugate-gradient method (Krylov.jl) for the SPD systems of the SEM
 linear solves:

   • the full SEM system (standard_linsolve! with :linsolve_amg => true)
   • the statically condensed skeleton system of element learning
     (elementLearning_Axb! / elementLearning_infer! with
      :EL_skeleton_solver => "amg")

 Two AMG flavours, :amg_method =>
   "sa"  smoothed aggregation (default) — robust for the high-order SEM
         stiffness, whose off-diagonal entries are not all ≤ 0
   "rs"  classical Ruge–Stüben
 CG stops at the relative residual :amg_rtol (default 1e-12), or after
 :amg_itmax iterations (default 1000; full-system solves). The hierarchy
 is built once (setup) and then applied as a V-cycle preconditioner by every
 CG iteration (solve). The iteration count and final residual of the last
 solve are kept in JX_AMG_STATS.

 A singular periodic system (constants in the null space) is passed here with
 one unknown pinned, so AMG always sees an SPD matrix.
=============================================================================#

const JX_AMG_STATS = Ref((iters = 0, rel_resid = NaN, levels = 0, method = :none))

_amg_method(s) = (m = Symbol(lowercase(string(s)));
                  m in (:sa, :rs) || error(" # AMG: :amg_method => \"$s\"; expected \"sa\" or \"rs\".");
                  m)

"""
    jx_amg_setup(A; method = :sa) -> S

Build the AMG hierarchy of the SPD sparse matrix `A` and its V-cycle
preconditioner (the one-time setup of the AMG solve).
"""
function jx_amg_setup(A::SparseMatrixCSC; method = :sa)
    m   = _amg_method(method)
    # AlgebraicMultigrid wants Int indices (no copy when A already has them)
    Ai  = A isa SparseMatrixCSC{Float64, Int} ? A : SparseMatrixCSC{Float64, Int}(A)
    ml  = m === :sa ? AlgebraicMultigrid.smoothed_aggregation(Ai) :
                      AlgebraicMultigrid.ruge_stuben(Ai)
    # CG's matrix-vector product: threaded column dot products when A is
    # exactly symmetric (A x = Aᵀ x), else SparseArrays' serial product
    op  = issymmetric(Ai) ? JXSymCSC(Ai) : Ai
    return (A = Ai, op = op, ml = ml, P = AlgebraicMultigrid.aspreconditioner(ml), method = m)
end

# An exactly symmetric CSC matrix as a CG operator: (A x)_j = Σ_p A[rowval[p], j] x[rowval[p]],
# one column per output entry, threaded without write conflicts
struct JXSymCSC{Ti}
    A :: SparseMatrixCSC{Float64, Ti}
end
Base.size(S::JXSymCSC) = size(S.A)
Base.size(S::JXSymCSC, d) = size(S.A, d)
Base.eltype(::JXSymCSC) = Float64
function LinearAlgebra.mul!(y::AbstractVector, S::JXSymCSC, x::AbstractVector)
    A = S.A;  n = size(A, 2);  nt = Threads.nthreads()
    Threads.@threads :static for c = 1:nt
        lo = div((c - 1) * n, nt) + 1;  hi = div(c * n, nt)
        @inbounds for j = lo:hi
            s = 0.0
            for p = A.colptr[j]:A.colptr[j+1]-1;  s += A.nzval[p] * x[A.rowval[p]];  end
            y[j] = s
        end
    end
    return y
end

"""
    jx_amg_solve(S, b; rtol = 1e-12, itmax = 1000) -> x

AMG-preconditioned CG solve of `S.A x = b` (S from `jx_amg_setup`).
"""
function jx_amg_solve(S, b::AbstractVector; rtol::Real = 1e-12, itmax::Int = 1000)
    x, st = Krylov.cg(S.op, Vector{Float64}(b); M = S.P, ldiv = true,
                      rtol = Float64(rtol), atol = 0.0, itmax = itmax, history = true)
    r0 = isempty(st.residuals) ? NaN : first(st.residuals)
    rel = isempty(st.residuals) || r0 == 0 ? 0.0 : last(st.residuals) / r0
    JX_AMG_STATS[] = (iters = st.niter, rel_resid = rel,
                      levels = length(S.ml.levels) + 1, method = S.method)
    st.solved || @warn " # AMG-CG did not converge to rtol=$rtol in $itmax iterations " *
                       "(final relative residual $rel)."
    return x
end

# AMG options from a case deck.
jx_amg_options(inputs) = (method = get(inputs, :amg_method, "sa"),
                          rtol   = Float64(get(inputs, :amg_rtol, 1e-12)),
                          itmax  = Int(get(inputs, :amg_itmax, 1000)))

# ── Fill-reducing ordering of the sparse Cholesky factorisations ──────────────
# :sparse_ordering => "cholmod" (default: CHOLMOD's own, AMD) | "metis" (METIS
# nested dissection: near-optimal fill on 3D meshes, n^(4/3)).
function jx_sparse_ordering(inputs)
    o = Symbol(lowercase(string(get(inputs, :sparse_ordering, "cholmod"))))
    o in (:cholmod, :metis) ||
        error(" # :sparse_ordering => \"$o\"; expected \"cholmod\" or \"metis\".")
    return o
end

# METIS is installed with Jexpresso (a dependency of Gridap's partitioning)
# but is not a direct dependency: it is loaded by package id on first use.
const _JX_METIS_ID = Base.PkgId(Base.UUID("2679e427-3c69-5b7f-982b-ece356f1e94b"), "Metis")

"""
    jx_metis_perm(A) -> Vector{Int}

METIS nested-dissection fill-reducing permutation of the exactly symmetric
sparse matrix `A`, for `cholesky(Symmetric(A); perm = p)`.
"""
function jx_metis_perm(A::SparseMatrixCSC)
    Metis = Base.require(_JX_METIS_ID)
    p, _ = Base.invokelatest(Metis.permutation, A)
    return Vector{Int}(p)
end

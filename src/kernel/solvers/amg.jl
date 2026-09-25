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
    Ai  = SparseMatrixCSC{Float64, Int}(A)          # AlgebraicMultigrid wants Int indices
    ml  = m === :sa ? AlgebraicMultigrid.smoothed_aggregation(Ai) :
                      AlgebraicMultigrid.ruge_stuben(Ai)
    return (A = Ai, ml = ml, P = AlgebraicMultigrid.aspreconditioner(ml), method = m)
end

"""
    jx_amg_solve(S, b; rtol = 1e-12, itmax = 1000) -> x

AMG-preconditioned CG solve of `S.A x = b` (S from `jx_amg_setup`).
"""
function jx_amg_solve(S, b::AbstractVector; rtol::Real = 1e-12, itmax::Int = 1000)
    x, st = Krylov.cg(S.A, Vector{Float64}(b); M = S.P, ldiv = true,
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

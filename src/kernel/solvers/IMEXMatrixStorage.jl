#-------------------------------------------------------------------------------
# Persistent storage for the IMEX implicit operator.
#
# The matrix-free IMEX path applies the implicit operator `L u` through a
# user-supplied callback (`L_fun!`). That path avoids global matrix storage,
# but it also blocks any preconditioner or mixed-precision strategy that needs
# direct access to the assembled sparse operator (ILU, Ruge–Stuben AMG,
# smoothed-aggregation AMG with a coarsened copy in reduced precision, etc.).
#
# This file gives the IMEX time integrators an opt-in alternative: an
# assembled-matrix path that stores the implicit operator
#
#     L_curr = I − λ L
#
# once (per time step, or per stage, depending on the integrator), keeps it
# around across nonlinear iterations, and carries along mixed-precision copies
# plus a cached preconditioner. The user selects between the two paths through
# the new `:matrix_storage` input flag:
#
#     :assembled   — assemble and store the sparse operator (this file)
#     :matrix_free — apply the operator through `L_fun!` (no global matrix)
#
# The three precision slots in `StoredIMEXMatrix` are:
#
#     L_full   — master copy (nonlinear-residual precision)
#     L_solver — copy handed to the Krylov solver
#     L_prec   — copy used to build the preconditioner
#
# When the three precisions coincide (the default) all three fields alias the
# same matrix, so the default path incurs no extra storage.
#-------------------------------------------------------------------------------

mutable struct StoredIMEXMatrix{Tf, Ts, Tp}
    L_full   :: SparseMatrixCSC{Tf, Int}
    L_solver :: SparseMatrixCSC{Ts, Int}
    L_prec   :: SparseMatrixCSC{Tp, Int}
    prec     :: Any
    prec_sp  :: Dict
end

"""
    StoredIMEXMatrix(L; solver_precision, prec_sp)

Wrap a freshly assembled sparse `L` in the storage struct, building the
mixed-precision copies and the preconditioner once.

`prec_sp[:precision]` sets the preconditioner-copy precision.
`prec_sp[:prec_type]` is one of `"AMG"`, `"ILU"`, or `"Jacobi"`.
"""
function StoredIMEXMatrix(L::SparseMatrixCSC;
                          solver_precision::Type = eltype(L),
                          prec_sp::Dict = Dict(:maxiter   => 1,
                                               :abstol    => 1e-8,
                                               :precision => eltype(L),
                                               :prec_type => "AMG"))
    Tf = eltype(L)
    Ts = solver_precision
    Tp = get(prec_sp, :precision, Tf)

    L_solver = Ts === Tf ? L : _imex_cast_sparse(L, Ts)
    L_prec   = Tp === Tf ? L : _imex_cast_sparse(L, Tp)
    prec     = _imex_build_preconditioner(L_prec, prec_sp)

    return StoredIMEXMatrix{Tf, Ts, Tp}(L, L_solver, L_prec, prec, prec_sp)
end

"""
    update_stored_matrix!(S, L)

Replace the stored operator with a freshly assembled `L`, refreshing the
mixed-precision copies and rebuilding the preconditioner. Call this whenever
the caller re-assembles the implicit operator (e.g. `upd_L = true`, or at the
start of a new RK stage with a different `a_tilde_ii`).
"""
function update_stored_matrix!(S::StoredIMEXMatrix{Tf, Ts, Tp},
                               L::SparseMatrixCSC) where {Tf, Ts, Tp}
    S.L_full   = eltype(L) === Tf ? L : _imex_cast_sparse(L, Tf)
    S.L_solver = eltype(L) === Ts ? L : _imex_cast_sparse(L, Ts)
    S.L_prec   = eltype(L) === Tp ? L : _imex_cast_sparse(L, Tp)
    S.prec     = _imex_build_preconditioner(S.L_prec, S.prec_sp)
    return S
end

"""
    rebuild_stored_preconditioner!(S)

Rebuild the preconditioner from the current `L_prec` without touching the
matrix copies. Useful when tuning preconditioner options between iterations.
"""
function rebuild_stored_preconditioner!(S::StoredIMEXMatrix)
    S.prec = _imex_build_preconditioner(S.L_prec, S.prec_sp)
    return S
end

_imex_cast_sparse(L::SparseMatrixCSC, ::Type{T}) where {T} =
    SparseMatrixCSC{T, Int}(L)

function _imex_build_preconditioner(L_prec::SparseMatrixCSC, prec_sp::Dict)
    prec_type = lowercase(string(get(prec_sp, :prec_type, "AMG")))
    if prec_type == "amg"
        return MyPrecClass.MyPrec(L_prec, RugeStubenAMG(), prec_sp)
    elseif prec_type == "ilu"
        droptol = get(prec_sp, :droptol, get(prec_sp, :ilu_tol, 0.1))
        return MyPrecClass.MyPrec(L_prec, IncompleteLU.ilu(L_prec; τ = droptol), prec_sp)
    elseif prec_type == "jacobi"
        return MyPrecClass.MyPrec(L_prec, Diagonal(diag(L_prec)), prec_sp)
    else
        error("StoredIMEXMatrix: unknown :prec_type = \"$prec_type\" " *
              "(expected \"AMG\", \"ILU\", or \"Jacobi\").")
    end
end

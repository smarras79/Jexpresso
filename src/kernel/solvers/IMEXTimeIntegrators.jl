#----------------------------------------------------------------------------
# IMEXTimeIntegrators.jl
#
# Jexpresso's native IMEX (implicit-explicit) time integrator.
#
# This is the "proper Jexpresso" home for the IMEX time advancement that used
# to live on a parallel, copy-pasted entry chain (IMEXJexpresso.jl / IMEXrun.jl
# / IMEXdrivers.jl). Instead of duplicating the driver, the integrator is
# selected exactly like every other time scheme: through `inputs[:ode_solver]`.
# When `inputs[:ode_solver] isa IMEX`, `time_loop!` (TimeIntegrators.jl)
# dispatches straight to `imex_time_loop!` below.
#
# Splitting:
#   du/dt = S(u) + L u
#     - S(u)  : the (stiff-free) EXPLICIT part, evaluated through the standard
#               `rhs!` (advection / sources). Supplied as `inputs[:S_fun]`.
#     - L u   : the IMPLICIT part, a linear operator (e.g. the viscous
#               -μ M⁻¹ K Laplacian). Supplied as `inputs[:L_fun]` (apply) and
#               `inputs[:build_L]` (assemble the sparse operator).
#
# Schemes (`inputs[:method]`):
#   "RK"        : IMEX additive Runge-Kutta. Butcher tableaux in `inputs[:coeff]`
#                 as (A_RK,b_RK,c_RK) for the explicit part and
#                 (A_RK_tilde,b_RK_tilde,c_RK_tilde) for the implicit (DIRK) part.
#   "multistep" : IMEX multistep, coefficients (alpha,beta,xi) in `inputs[:coeff]`.
#
# `inputs[:delta]`:
#   0 -> pure explicit (no implicit solve), 1 -> IMEX.
#
# Linear solve per stage/step:
#   Each implicit stage forms  L_curr = I - λ L_temp  and solves
#   L_curr x = nonl_res inside a fixed-point nonlinear loop. The solve uses
#   `inputs[:lsolve]` when provided (e.g. a direct sparse `L_curr \ b`), and
#   otherwise falls back to a direct sparse factorization. Krylov / AMG
#   preconditioning is intentionally left out of this lean integration; cases
#   that need it can plug a custom `lsolve` closure.
#----------------------------------------------------------------------------

# Pull in the portable JACC linear-algebra kernels used by the GPU/JACC path
# (`:limex_jacc => true`). Included here, next to the only code that uses them,
# so the integrator and its solver always load together regardless of the
# top-level include list in Jexpresso.jl. `@__DIR__` is this file's directory
# (src/kernel/solvers), where imex_jacc.jl also lives.
include(joinpath(@__DIR__, "imex_jacc.jl"))

# Precision study hooks for the JACC offload path. When
# `IMEX_JACC_SOLVE_PRECISION[]` is a concrete type (Float16/Float32/Float64) it
# overrides `inputs[:imex_jacc_solve_precision]` for the device solve; the
# offload runner stashes that run's diagnostics in `IMEX_JACC_LAST_DIAG` so
# `run_imex_precision_study` can collect them across runs.
const IMEX_JACC_SOLVE_PRECISION = Ref{Union{DataType,Nothing}}(nothing)
const IMEX_JACC_LAST_DIAG       = Ref{Any}(nothing)

function imex_time_loop!(inputs, params, u)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    println_rank(" #   IMEX integrator ........................... "; msg_rank = rank)

    #--------------------------------------------------------------------------
    # Problem sizes
    #--------------------------------------------------------------------------
    mesh   = params.mesh
    npoint = mesh.npoin
    neqs   = params.neqs
    unkwn  = neqs * npoint

    # The user S_fun!/bcs_fun! closures carry the legacy (params, sem) /
    # (..., sem, qp) signatures. Everything they need now lives in `params`,
    # so the mesh/SEM bundle is passed through `params` and `sem`/`qp` are
    # only kept for signature compatibility.
    sem = nothing
    qp  = params.qp

    #--------------------------------------------------------------------------
    # Time control
    #--------------------------------------------------------------------------
    Δt  = inputs[:Δt]
    t_n = inputs[:tinit]

    #--------------------------------------------------------------------------
    # Precision and solver knobs
    #--------------------------------------------------------------------------
    lsolve           = inputs[:lsolve]
    solver_precision = inputs[:solver_precision]
    nl_precision     = inputs[:nl_precision]

    nl_atol     = inputs[:nl_atol]
    nl_rtol     = inputs[:nl_rtol]
    max_nl_iter = inputs[:max_nl_iter]

    #--------------------------------------------------------------------------
    # IMEX configuration
    #--------------------------------------------------------------------------
    delta  = inputs[:delta]          # 0 -> explicit, 1 -> IMEX
    method = inputs[:method]         # "RK" or "multistep"
    k      = inputs[:k]              # # of RK stages, or # of multistep steps
    if k < 1
        error("IMEX: k should be a positive integer (got $k)")
    end
    coeff   = inputs[:coeff]

    S_fun!  = inputs[:S_fun]
    bcs_fun! = inputs[:bcs_fun]
    # `upd_L` is always defined (false in the pure-explicit case) so the hot
    # loops can test it without a `delta == 1` guard.
    upd_L = (delta == 1) ? inputs[:upd_L] : false
    if delta == 1
        L_fun!  = inputs[:L_fun]
        build_L = inputs[:build_L]
    end

    backend = inputs[:backend]

    #--------------------------------------------------------------------------
    # Linear solve for  L_curr x = b
    #--------------------------------------------------------------------------
    function imex_linsolve(L_curr, b)
        if lsolve === nothing
            return L_curr \ b           # direct sparse factorization fallback
        else
            return lsolve(L_curr, b)    # user-provided solve (e.g. direct / Krylov)
        end
    end

    #--------------------------------------------------------------------------
    # Assemble the per-stage implicit operator  L_curr = I - λ L_temp
    #--------------------------------------------------------------------------
    Is_cache = sparse(one(solver_precision) * I, Int64(unkwn), Int64(unkwn))

    function L_update(uloc, tloc, λ)
        if delta == 0
            return solver_precision.(Is_cache)
        else
            # build_L must already fold in the inverse mass matrix scaling.
            L_temp = build_L(uloc, tloc, params)
            return solver_precision.(Is_cache - λ * L_temp)
        end
    end

    # When the implicit operator is constant in time (`upd_L == false`) the
    # per-stage system matrix `I - λ L` never changes, so it is assembled and
    # LU-factorized ONCE and the factorization is reused for every step. This
    # is the decisive memory/time win: without it, a fresh sparse LU (with full
    # 2-D fill-in) would be built and thrown away on every stage of every step
    # (~k·nsteps factorizations), which is what made this case allocate tens of
    # GiB. The factorization is passed to the solve in place of the raw matrix;
    # `F \ b` reuses it, so no re-factorization happens in the hot loop.
    operator_is_constant = (delta == 1 && !upd_L)

    # Avoid an extra residual copy when the residual and solver precisions match
    # (the common case, e.g. both Float64).
    same_precision = (solver_precision === nl_precision)
    solver_rhs(res) = same_precision ? res : solver_precision.(res)

    #--------------------------------------------------------------------------
    # RHS builders (in-place into pre-allocated buffers)
    #--------------------------------------------------------------------------
    # Runge-Kutta: rhs = u_prev + Σ_{j<i} Δt[ A_ex_ij S_j + (A_im_ij - A_ex_ij) L_j ]
    function construct_rhs_rk!(rhs, s_j, l_j, u_prev, U_stages, t_n, i,
                               A_RK, A_RK_tilde, c_RK)
        copyto!(rhs, u_prev)
        for j = 1:i-1
            fill!(s_j, zero(eltype(s_j)))
            time = t_n + c_RK[j] * Δt
            S_fun!(s_j, U_stages[j], time, params, sem)
            if delta == 1
                fill!(l_j, zero(eltype(l_j)))
                L_fun!(l_j, U_stages[j], time, params)
                α = Δt * A_RK[i, j]
                β = Δt * (A_RK_tilde[i, j] - A_RK[i, j])
                @. rhs += α * s_j + β * l_j
            else
                α = Δt * A_RK[i, j]
                @. rhs += α * s_j
            end
        end
        return rhs
    end

    # Multistep: rhs = Σ_n α_n u_{n} + λ Σ_n β_n ( S_n - L_n )
    function construct_rhs_multistep!(rhs, s_j, l_j, u_prev, t_n,
                                      alpha, beta, lambda)
        fill!(rhs, zero(eltype(rhs)))
        for n = 1:k
            @. rhs += alpha[n] * u_prev[n]
            fill!(s_j, zero(eltype(s_j)))
            time = t_n - (n - 1) * Δt
            S_fun!(s_j, u_prev[n], time, params, sem)
            if delta == 1
                fill!(l_j, zero(eltype(l_j)))
                L_fun!(l_j, u_prev[n], time, params)
                γ = lambda * beta[n]
                @. rhs += γ * (s_j - l_j)
            else
                γ = lambda * beta[n]
                @. rhs += γ * s_j
            end
        end
        return rhs
    end

    #--------------------------------------------------------------------------
    # Diagnostics output schedule
    #--------------------------------------------------------------------------
    dosetimes = collect(TFloat.(inputs[:diagnostics_at_times]))
    sort!(dosetimes)
    next_out_idx = 1
    iout = 0
    n_step = 0

    function maybe_write_output!()
        # Write at the first step whose time has reached the next requested
        # diagnostic time (robust to Δt round-off accumulation).
        wrote = false
        while next_out_idx <= length(dosetimes) && t_n + 1.0e-10 >= dosetimes[next_out_idx]
            iout += 1
            u2uaux!(params.uaux, u, neqs, npoint)
            println_rank(@sprintf(" #   IMEX: t = %.6f   step = %d", t_n, n_step); msg_rank = rank)
            write_output(params.SD, u, params.uaux, t_n, iout,
                         params.mesh, params.mp,
                         params.connijk_original, params.poin_in_bdy_face_original,
                         params.x_original, params.y_original, params.z_original,
                         inputs[:output_dir], inputs,
                         params.qp.qvars, params.qp.qoutvars,
                         inputs[:outformat];
                         nvar = params.qp.neqs, qexact = params.qp.qe)
            next_out_idx += 1
            wrote = true
        end
        return wrote
    end

    #--------------------------------------------------------------------------
    # Pre-allocated work buffers (allocation-free hot loops)
    #--------------------------------------------------------------------------
    rhs_buf      = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
    s_j_buf      = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
    l_j_buf      = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
    nonl_res_buf = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
    Lu_buf       = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))

    #----------------------------------------------------------------------------
    # Runge-Kutta IMEX
    #----------------------------------------------------------------------------
    if method == "RK"
        for key in (:A_RK, :b_RK, :c_RK, :A_RK_tilde, :b_RK_tilde, :c_RK_tilde)
            haskey(coeff, key) || error("IMEX/RK: missing coefficient $(key) in :coeff")
        end
        A_RK       = coeff[:A_RK];       b_RK       = coeff[:b_RK];       c_RK       = coeff[:c_RK]
        A_RK_tilde = coeff[:A_RK_tilde]; b_RK_tilde = coeff[:b_RK_tilde]; c_RK_tilde = coeff[:c_RK_tilde]

        # Fast path for the common IMEX case: a time-independent implicit
        # operator (upd_L == false) and matched residual/solver precision. The
        # whole hot loop runs inside a type-stable function barrier so that the
        # stage vectors handed to `rhs!` are concretely typed - otherwise the
        # Any-typed locals here would make every `rhs!` evaluation type-unstable
        # and box, allocating GiB over a run.
        if operator_is_constant && same_precision
            # GPU/JACC opt-in (`:limex_jacc => true`): the per-stage implicit
            # solve and the residual sparse mat-vec run on the JACC device via a
            # portable BiCGSTAB instead of the host sparse LU. Everything else
            # (vector algebra, explicit `rhs!`) stays on the native backend that
            # `u` lives on, so a GPU run only needs `:backend => CUDABackend()`
            # and a JACC backend configured for the matching device.
            if get(inputs, :limex_jacc, false)
                if get(inputs, :limex_jacc_offload, false)
                    # HYBRID offload: the Jexpresso pipeline (mesh / rhs! / vector
                    # algebra) stays on the host (`:backend => CPU()`), and ONLY
                    # the per-stage implicit solve (I - λL) x = b is run on the GPU
                    # via JACC — the operator is uploaded once and each stage's
                    # small rhs/solution vector is copied host↔device per step.
                    _imex_rk_run_const_jacc_offload!(u, build_L, L_update, S_fun!, bcs_fun!,
                                                     params, qp, sem, inputs,
                                                     A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                                                     Δt, k, nl_atol, nl_rtol, max_nl_iter,
                                                     t_n, inputs[:tend], dosetimes, neqs, npoint, rank)
                    println_rank(" #   IMEX integrator (JACC offload) ............ DONE"; msg_rank = rank)
                    return u
                end
                _imex_rk_run_const_jacc!(u, build_L, L_update, S_fun!, bcs_fun!,
                                         params, qp, sem, inputs,
                                         A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                                         Δt, k, nl_atol, nl_rtol, max_nl_iter,
                                         t_n, inputs[:tend], dosetimes, neqs, npoint, rank)
                println_rank(" #   IMEX integrator (JACC) .................... DONE"; msg_rank = rank)
                return u
            end
            # The constant-operator barrier prints a per-region allocation
            # breakdown at the end when `:lalloc_summary => true` (runtime flag).
            _imex_rk_run_const!(u, L_update, S_fun!, L_fun!, bcs_fun!,
                                params, qp, sem, inputs,
                                A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                                Δt, k, nl_atol, nl_rtol, max_nl_iter,
                                t_n, inputs[:tend], dosetimes, neqs, npoint, rank)
            println_rank(" #   IMEX integrator ........................... DONE"; msg_rank = rank)
            return u
        end

        U_stages = Dict{Int,Any}()
        for j = 1:k
            U_stages[j] = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
        end

        # Constant-operator fast path: assemble and factorize each stage's
        # `I - λ_i L` once and reuse the factorization for every time step.
        Lmat_stage = Vector{Any}(undef, k)   # sparse operator   (for residual mul!)
        Lfac_stage = Vector{Any}(undef, k)   # its LU factorization (for the solve)
        if operator_is_constant
            for i = 1:k
                λi            = Δt * A_RK_tilde[i, i]
                Lmat_stage[i] = L_update(u, t_n, λi)
                Lfac_stage[i] = lu(Lmat_stage[i])
            end
        end

        while (abs(t_n - inputs[:tend]) > 1.0e-14 && t_n < inputs[:tend])
            for i = 1:k
                time_tilde = t_n + c_RK_tilde[i] * Δt
                λ          = Δt * A_RK_tilde[i, i]

                # Implicit operator for this stage. `L_mat` is the matrix (used
                # for the residual); `L_solve` is what is handed to the linear
                # solve - a cached factorization when the operator is constant,
                # otherwise the freshly-built matrix.
                if operator_is_constant
                    L_mat   = Lmat_stage[i]
                    L_solve = Lfac_stage[i]
                else
                    L_mat   = L_update(U_stages[i], time_tilde, λ)
                    L_solve = L_mat
                end

                # Stage RHS
                rhs = construct_rhs_rk!(rhs_buf, s_j_buf, l_j_buf,
                                        u, U_stages, t_n, i,
                                        A_RK, A_RK_tilde, c_RK)
                if bcs_fun! !== nothing
                    bcs_fun!(rhs, L_mat, time_tilde, params, sem, qp)
                end

                # Nonlinear fixed-point loop: solve L_curr U_i = rhs
                nonl_res = nonl_res_buf
                copyto!(nonl_res, rhs)
                nl_norm_0 = norm(nonl_res)
                nl_norm_k = nl_norm_0
                nl_iter   = 1
                while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                    x = imex_linsolve(L_solve, solver_rhs(nonl_res))
                    U_stages[i] .+= x

                    if upd_L
                        L_mat   = L_update(U_stages[i], time_tilde, λ)
                        L_solve = L_mat
                        bcs_fun! !== nothing && bcs_fun!(rhs, L_mat, time_tilde, params, sem, qp)
                    end

                    mul!(Lu_buf, L_mat, U_stages[i])
                    @. nonl_res = rhs - Lu_buf
                    nl_norm_k = norm(nonl_res)
                    nl_iter += 1
                end
            end

            # Solution update: u += Σ_i Δt b_RK[i] S(U_i)
            for i = 1:k
                fill!(s_j_buf, zero(eltype(s_j_buf)))
                time = t_n + c_RK[i] * Δt
                S_fun!(s_j_buf, U_stages[i], time, params, sem)
                axpy!(Δt * b_RK[i], s_j_buf, u)
            end

            t_n   += Δt
            n_step += 1

            maybe_write_output!()

            for j = 1:k
                fill!(U_stages[j], zero(eltype(U_stages[j])))
            end
        end

    #----------------------------------------------------------------------------
    # Multistep IMEX
    #----------------------------------------------------------------------------
    elseif method == "multistep"
        for key in (:alpha, :beta, :xi)
            haskey(coeff, key) || error("IMEX/multistep: missing coefficient $(key) in :coeff")
        end
        alpha  = coeff[:alpha]
        beta   = coeff[:beta]
        xi     = coeff[:xi]
        lambda = xi * Δt

        # History buffers
        u_prev = Dict{Int,Any}()
        for n = 1:k
            u_prev[n] = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
        end
        u_prev[1] .= u

        # Warm-up the k-step history with explicit Euler micro-steps
        if k > 1
            Δt_expl = inputs[:Δt_expl]
            warm_rhs = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
            for n = 2:k
                copyto!(warm_rhs, u_prev[1])
                fill!(s_j_buf, zero(eltype(s_j_buf)))
                S_fun!(s_j_buf, u_prev[1], t_n, params, sem)
                axpy!(Δt_expl, s_j_buf, warm_rhs)
                for istep = n:-1:2
                    u_prev[istep] .= u_prev[istep-1]
                end
                u_prev[1] .= warm_rhs
                u .= warm_rhs
                t_n += Δt_expl
            end
        end

        u_next = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
        # Constant operator (upd_L == false): assemble + factorize once, reuse.
        L_mat   = L_update(u_next, t_n + Δt, lambda)
        L_solve = operator_is_constant ? lu(L_mat) : L_mat

        while (abs(t_n - inputs[:tend]) > 1.0e-14 && t_n < inputs[:tend])
            rhs = construct_rhs_multistep!(rhs_buf, s_j_buf, l_j_buf,
                                           u_prev, t_n, alpha, beta, lambda)
            if bcs_fun! !== nothing
                bcs_fun!(rhs, L_mat, t_n + Δt, params, sem, qp)
            end

            nonl_res = nonl_res_buf
            copyto!(nonl_res, rhs)
            nl_norm_0 = norm(nonl_res)
            nl_norm_k = nl_norm_0
            nl_iter   = 1
            while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                x = imex_linsolve(L_solve, solver_rhs(nonl_res))
                u_next .+= x

                if upd_L
                    L_mat   = L_update(u_next, t_n + Δt, lambda)
                    L_solve = L_mat
                    bcs_fun! !== nothing && bcs_fun!(rhs, L_mat, t_n + Δt, params, sem, qp)
                end

                mul!(Lu_buf, L_mat, u_next)
                @. nonl_res = rhs - Lu_buf
                nl_norm_k = norm(nonl_res)
                nl_iter += 1
            end

            # Shift the history and commit u_next
            for n = k:-1:2
                u_prev[n] .= u_prev[n-1]
            end
            u_prev[1] .= u_next
            u .= u_next

            t_n   += Δt
            n_step += 1

            maybe_write_output!()

            fill!(u_next, zero(eltype(u_next)))
        end
    else
        error("IMEX: unknown :method = $(repr(method)). Expected \"RK\" or \"multistep\".")
    end

    # Always emit a final snapshot if the schedule did not already cover t_end.
    if iout == 0 || next_out_idx <= length(dosetimes)
        iout += 1
        u2uaux!(params.uaux, u, neqs, npoint)
        write_output(params.SD, u, params.uaux, t_n, iout,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar = params.qp.neqs, qexact = params.qp.qe)
    end

    println_rank(" #   IMEX integrator ........................... DONE"; msg_rank = rank)

    return u
end


# Thin wrapper around write_output for the IMEX diagnostic snapshots.
function _imex_write(params, u, t, iout, inputs)
    write_output(params.SD, u, params.uaux, t, iout,
                 params.mesh, params.mp,
                 params.connijk_original, params.poin_in_bdy_face_original,
                 params.x_original, params.y_original, params.z_original,
                 inputs[:output_dir], inputs,
                 params.qp.qvars, params.qp.qoutvars,
                 inputs[:outformat];
                 nvar = params.qp.neqs, qexact = params.qp.qe)
    return nothing
end


#----------------------------------------------------------------------------
# Type-stable hot loop for the constant-operator IMEX additive Runge-Kutta
# path (delta == 1, upd_L == false, matched residual/solver precision).
#
# Why a separate function: `imex_time_loop!` reads the time step, buffers and
# user closures out of `inputs`/runtime-typed allocations, so inside it the
# stage vectors and `S_fun!`/`L_fun!` are inferred `Any`. Passing an Any-typed
# stage vector into `rhs!` makes the entire RHS evaluation type-unstable and
# box on every array op - which is what made this small case allocate GiB. By
# taking `u` and the closures as positional arguments, Julia specializes this
# function on their concrete runtime types; the stage vectors and work buffers
# are then allocated here with `similar(u)`, so they are concretely typed and
# `rhs!` stays allocation-light.
#
# The constant implicit operator `I - λᵢ L` is assembled and LU-factorized
# once per stage and reused for every step; the per-stage solve is an in-place
# `ldiv!` with the cached factorization (equivalent to the case's direct
# `:lsolve = L \ b`). A custom iterative `:lsolve` is only honoured on the
# general path in `imex_time_loop!`.
#----------------------------------------------------------------------------
function _imex_rk_run_const!(u, L_update, S_fun!, L_fun!, bcs_fun!,
                             params, qp, sem, inputs,
                             A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                             Δt, k, nl_atol, nl_rtol, max_nl_iter,
                             tinit, tend, dosetimes, neqs, npoint, rank)

    # Work buffers and stage vectors. `u` is a concrete positional argument
    # here, so `similar(u)` yields concretely-typed arrays.
    rhs_buf = fill!(similar(u), zero(eltype(u)))
    res_buf = fill!(similar(u), zero(eltype(u)))
    Lu_buf  = fill!(similar(u), zero(eltype(u)))
    x_buf   = fill!(similar(u), zero(eltype(u)))
    U_stages = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]

    # Cached explicit/implicit stage evaluations S(U_i) and L(U_i). Each is
    # computed ONCE, right after stage i is solved, and reused both by the
    # later stages' RHS assembly and by the final solution update. This drops
    # the number of (expensive) `rhs!` evaluations from ~k(k+1)/2 + k per step
    # to exactly k per step.
    S_store = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]
    L_store = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]

    # Constant per-stage operators I - λ_i L, with λ_i = Δt * Ã[i,i]. Stages
    # that share the same λ_i (e.g. ARS(2,3,2) stages 2 and 3) are the SAME
    # matrix, so assemble + LU-factorize once per DISTINCT λ_i and share the
    # factorization. This avoids redundant (and memory-heavy) factorizations.
    λs = [Δt * A_RK_tilde[i, i] for i in 1:k]
    L0   = L_update(u, tinit, λs[1])
    F0   = lu(L0)
    Lmat = Vector{typeof(L0)}(undef, k)
    Lfac = Vector{typeof(F0)}(undef, k)
    Lmat[1] = L0
    Lfac[1] = F0
    for i in 2:k
        reuse = 0
        for j in 1:i-1
            if λs[j] == λs[i]
                reuse = j
                break
            end
        end
        if reuse > 0
            Lmat[i] = Lmat[reuse]
            Lfac[i] = Lfac[reuse]
        else
            Lmat[i] = L_update(u, tinit, λs[i])
            Lfac[i] = lu(Lmat[i])
        end
    end

    t_n          = tinit
    n_step       = 0
    next_out_idx = 1
    iout         = 0

    # Per-region allocation accounting. `@allocated` works at runtime (unlike
    # @timeit_debug, which is gated at module-load time by JEXPRESSO_ALLOC_SUMMARY)
    # so this populates even when only `:lalloc_summary => true` is set in
    # user_inputs.jl. The accounting itself is cheap (a couple of gc_bytes reads
    # per call) and the breakdown is printed once at the end.
    lmeasure = alloc_summary_enabled(inputs)
    a_asm = 0; a_solve = 0; a_mul = 0; a_Sfun = 0; a_Lfun = 0; a_out = 0

    while (abs(t_n - tend) > 1.0e-14 && t_n < tend)
        for i = 1:k
            time_tilde = t_n + c_RK_tilde[i] * Δt
            L_mat = Lmat[i]
            L_fac = Lfac[i]

            # Stage RHS from cached stage evaluations:
            #   u + Σ_{j<i} Δt[ A_ex S_j + (A_im - A_ex) L_j ]
            a_asm += @allocated begin
                copyto!(rhs_buf, u)
                for j = 1:i-1
                    aex = Δt * A_RK[i, j]
                    bim = Δt * (A_RK_tilde[i, j] - A_RK[i, j])
                    @. rhs_buf += aex * S_store[j] + bim * L_store[j]
                end
            end
            if bcs_fun! !== nothing
                bcs_fun!(rhs_buf, L_mat, time_tilde, params, sem, qp)
            end

            # Solve (I - λL) U_i = rhs with the cached factorization. The
            # operator is linear, so the fixed-point loop converges in a single
            # solve; the residual check is kept for robustness.
            copyto!(res_buf, rhs_buf)
            nl_norm_0 = norm(res_buf)
            nl_norm_k = nl_norm_0
            nl_iter   = 1
            while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                a_solve += @allocated ldiv!(x_buf, L_fac, res_buf)
                U_stages[i] .+= x_buf
                a_mul   += @allocated mul!(Lu_buf, L_mat, U_stages[i])
                @. res_buf = rhs_buf - Lu_buf
                nl_norm_k = norm(res_buf)
                nl_iter  += 1
            end

            # Cache S(U_i) and L(U_i) for this stage (reused above and below).
            ti = t_n + c_RK[i] * Δt
            a_Sfun += @allocated S_fun!(S_store[i], U_stages[i], ti, params, sem)
            a_Lfun += @allocated L_fun!(L_store[i], U_stages[i], ti, params)
        end

        # Solution update: u += Σ_i Δt b_RK[i] S(U_i)  (cached)
        for i = 1:k
            axpy!(Δt * b_RK[i], S_store[i], u)
        end

        t_n    += Δt
        n_step += 1

        while next_out_idx <= length(dosetimes) && t_n + 1.0e-10 >= dosetimes[next_out_idx]
            iout += 1
            u2uaux!(params.uaux, u, neqs, npoint)
            println_rank(@sprintf(" #   IMEX: t = %.6f   step = %d", t_n, n_step); msg_rank = rank)
            a_out += @allocated _imex_write(params, u, t_n, iout, inputs)
            next_out_idx += 1
        end

        for j = 1:k
            fill!(U_stages[j], zero(eltype(U_stages[j])))
        end
    end

    # Final snapshot if the schedule did not already cover t_end.
    if iout == 0 || next_out_idx <= length(dosetimes)
        iout += 1
        u2uaux!(params.uaux, u, neqs, npoint)
        a_out += @allocated _imex_write(params, u, t_n, iout, inputs)
    end

    if lmeasure && rank == 0
        gib(b) = b / (1024.0^3)
        total = a_asm + a_solve + a_mul + a_Sfun + a_Lfun + a_out
        println()
        println(" # ===== IMEX per-region allocation (constant-operator RK path) =====")
        @printf("   stage RHS assembly        : %8.3f GiB\n", gib(a_asm))
        @printf("   linear solve (ldiv!)      : %8.3f GiB\n", gib(a_solve))
        @printf("   residual mul!             : %8.3f GiB\n", gib(a_mul))
        @printf("   S(U) [rhs!]               : %8.3f GiB\n", gib(a_Sfun))
        @printf("   L(U) [sparse mul!]        : %8.3f GiB\n", gib(a_Lfun))
        @printf("   output (write_output/VTK) : %8.3f GiB\n", gib(a_out))
        println("   ----------------------------------------------")
        @printf("   total measured            : %8.3f GiB\n", gib(total))
        println()
    end

    return nothing
end


#----------------------------------------------------------------------------
# JACC / GPU constant-operator IMEX additive Runge-Kutta path.
#
# A device-portable twin of `_imex_rk_run_const!`. It is selected by
# `:limex_jacc => true` and differs in exactly the piece the direct-LU path
# could not move off the host: the per-stage implicit solve `(I - λL) x = b`
# and the residual sparse mat-vec are done on the JACC device with
# `jacc_bicgstab!` / `jacc_spmv!` (imex_jacc.jl) instead of `lu` / `ldiv!`.
#
# Everything else is shared with the host path's logic and stays on the native
# KernelAbstractions backend of `u`:
#   * stage RHS assembly, the solution update and the residual combination are
#     plain broadcast / `axpy!` on `similar(u)` buffers (run on GPU when `u` is
#     a device array);
#   * S(U_i) is the explicit `rhs!` through `S_fun!`, which already dispatches on
#     that backend;
#   * L(U_i) is a device SpMV with the raw L operator.
#
# So on a CPU run (the default, and what is exercised by the test suite) this
# uses the JACC CPU backend and is bit-for-bit a BiCGSTAB-instead-of-LU solve;
# on `:backend => CUDABackend()` with a CUDA-configured JACC it runs on the GPU
# with no further code change.
#----------------------------------------------------------------------------
function _imex_rk_run_const_jacc!(u, build_L, L_update, S_fun!, bcs_fun!,
                                  params, qp, sem, inputs,
                                  A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                                  Δt, k, nl_atol, nl_rtol, max_nl_iter,
                                  tinit, tend, dosetimes, neqs, npoint, rank)

    # Configure JACC for the requested device (no-op / default on CPU).
    imex_jacc_init_backend!(inputs)

    index_type = (TInt === Int32) ? Int32 : Int
    float_type = eltype(u)

    # Work buffers and stage vectors live on the native backend of `u`.
    rhs_buf = fill!(similar(u), zero(eltype(u)))
    res_buf = fill!(similar(u), zero(eltype(u)))
    Lu_buf  = fill!(similar(u), zero(eltype(u)))
    x_buf   = fill!(similar(u), zero(eltype(u)))
    U_stages = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]
    S_store  = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]
    L_store  = [fill!(similar(u), zero(eltype(u))) for _ in 1:k]

    # BiCGSTAB scratch (r, rhat, p, v, s, t).
    bicg_work = ntuple(_ -> fill!(similar(u), zero(eltype(u))), 6)
    itmax = get(inputs, :imex_jacc_itmax, max(500, 10 * Int(ceil(sqrt(npoint)))))

    # Raw implicit operator L (for L(U_i) = L * U_i), assembled once on the host
    # and copied to the device as CSR.
    L_raw  = build_L(u, tinit, params)
    L_jacc = JaccSparseCSR(L_raw; index_type = index_type, float_type = float_type)

    # Per-stage system operators I - λ_i L. Stages that share the same λ_i share
    # the (host sparse + device CSR) operator, mirroring the LU-dedup of the host
    # path. The host sparse copy is kept only so a `bcs_fun!` that mutates the
    # matrix still has something to act on; the solve/SpMV use the CSR.
    λs   = [Δt * A_RK_tilde[i, i] for i in 1:k]
    A0_h = L_update(u, tinit, λs[1])
    Ahost = Vector{typeof(A0_h)}(undef, k)
    Ajacc = Vector{typeof(L_jacc)}(undef, k)
    Ahost[1] = A0_h
    Ajacc[1] = JaccSparseCSR(A0_h; index_type = index_type, float_type = float_type)
    for i in 2:k
        reuse = 0
        for j in 1:i-1
            if λs[j] == λs[i]
                reuse = j
                break
            end
        end
        if reuse > 0
            Ahost[i] = Ahost[reuse]
            Ajacc[i] = Ajacc[reuse]
        else
            Ahost[i] = L_update(u, tinit, λs[i])
            Ajacc[i] = JaccSparseCSR(Ahost[i]; index_type = index_type, float_type = float_type)
        end
    end

    t_n          = tinit
    n_step       = 0
    next_out_idx = 1
    iout         = 0

    while (abs(t_n - tend) > 1.0e-14 && t_n < tend)
        for i = 1:k
            time_tilde = t_n + c_RK_tilde[i] * Δt

            # Stage RHS:  u + Σ_{j<i} Δt[ A_ex S_j + (A_im - A_ex) L_j ]
            copyto!(rhs_buf, u)
            for j = 1:i-1
                aex = Δt * A_RK[i, j]
                bim = Δt * (A_RK_tilde[i, j] - A_RK[i, j])
                @. rhs_buf += aex * S_store[j] + bim * L_store[j]
            end
            if bcs_fun! !== nothing
                bcs_fun!(rhs_buf, Ahost[i], time_tilde, params, sem, qp)
            end

            # Solve (I - λL) U_i = rhs on the device. The operator is linear, so
            # the fixed-point loop converges in one BiCGSTAB solve; the residual
            # check below is kept for robustness / parity with the host path.
            copyto!(res_buf, rhs_buf)
            nl_norm_0 = norm(res_buf)
            nl_norm_k = nl_norm_0
            nl_iter   = 1
            while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                conv, bi, _ = jacc_bicgstab!(x_buf, Ajacc[i], res_buf, bicg_work;
                                             rtol = nl_rtol, atol = nl_atol, itmax = itmax)
                if !conv && rank == 0
                    @warn "IMEX/JACC BiCGSTAB did not converge (stage $i, step $n_step): $bi iters"
                end
                U_stages[i] .+= x_buf
                jacc_spmv!(Lu_buf, Ajacc[i], U_stages[i])
                @. res_buf = rhs_buf - Lu_buf
                nl_norm_k = norm(res_buf)
                nl_iter  += 1
            end

            # Cache S(U_i) [explicit rhs!] and L(U_i) [device SpMV].
            ti = t_n + c_RK[i] * Δt
            S_fun!(S_store[i], U_stages[i], ti, params, sem)
            jacc_spmv!(L_store[i], L_jacc, U_stages[i])
        end

        # Solution update: u += Σ_i Δt b_RK[i] S(U_i)
        for i = 1:k
            axpy!(Δt * b_RK[i], S_store[i], u)
        end

        t_n    += Δt
        n_step += 1

        while next_out_idx <= length(dosetimes) && t_n + 1.0e-10 >= dosetimes[next_out_idx]
            iout += 1
            u2uaux!(params.uaux, u, neqs, npoint)
            println_rank(@sprintf(" #   IMEX(JACC): t = %.6f   step = %d", t_n, n_step); msg_rank = rank)
            _imex_write(params, u, t_n, iout, inputs)
            next_out_idx += 1
        end

        for j = 1:k
            fill!(U_stages[j], zero(eltype(U_stages[j])))
        end
    end

    # Final snapshot if the schedule did not already cover t_end.
    if iout == 0 || next_out_idx <= length(dosetimes)
        iout += 1
        u2uaux!(params.uaux, u, neqs, npoint)
        _imex_write(params, u, t_n, iout, inputs)
    end

    return nothing
end


#----------------------------------------------------------------------------
# HYBRID GPU offload of the constant-operator IMEX additive Runge-Kutta path.
#
# Selected by `:limex_jacc => true` AND `:limex_jacc_offload => true`, with
# `:backend => CPU()` (the default). The entire Jexpresso pipeline — mesh,
# sem_setup, the explicit `rhs!` (S_fun!), the stage-RHS assembly and the
# solution update — runs on the HOST exactly as in the proven CPU path. The
# ONLY thing moved onto the GPU is the per-stage implicit linear solve
# (I - λL) x = b, performed with `jacc_bicgstab!` on JACC device arrays:
#
#   * the distinct per-stage operators I - λ_i L are assembled once on the host
#     and uploaded to the device as CSR (`JaccSparseCSR`, whose arrays are
#     `JACC.array`s — CuArrays when JACC's backend is CUDA);
#   * each stage copies its small rhs vector host→device, solves on the device,
#     and copies the solution device→host;
#   * S(U_i) is the host `rhs!`; L(U_i) is a host sparse mat-vec with the cached
#     operator (both feed the next stages' RHS).
#
# This sidesteps the (incomplete) GPU port of the rest of Jexpresso: nothing but
# the solver touches the device, so there are no GPU scalar-indexing / MPI-type
# issues. Requires CUDA loaded and JACC on its CUDA backend — run with
# `Jexpresso.run_case(...; backend = :cuda)` (which loads CUDA for JACC while the
# case keeps `:backend => CPU()`). If JACC is on its CPU backend the device
# arrays are plain host arrays and this simply runs the solve on the CPU.
#----------------------------------------------------------------------------
function _imex_rk_run_const_jacc_offload!(u, build_L, L_update, S_fun!, bcs_fun!,
                                          params, qp, sem, inputs,
                                          A_RK, b_RK, c_RK, A_RK_tilde, c_RK_tilde,
                                          Δt, k, nl_atol, nl_rtol, max_nl_iter,
                                          tinit, tend, dosetimes, neqs, npoint, rank)

    T = eltype(u)                 # host precision (the Jexpresso pipeline runs in T)
    n = length(u)
    index_type = (TInt === Int32) ? Int32 : Int

    # Precision of the device SOLVE only. A global override (set by
    # run_imex_precision_study) wins over the input; otherwise the input, else T.
    S = IMEX_JACC_SOLVE_PRECISION[]
    if S === nothing
        S = get(inputs, :imex_jacc_solve_precision, T)
    end
    S === nothing && (S = T)

    # Host buffers / stage vectors (u is a host array on the CPU backend, in T).
    rhs_h    = fill!(similar(u), zero(T))
    x_h      = fill!(similar(u), zero(T))
    rhs_s    = Vector{S}(undef, n)   # host staging in solve precision (T -> S)
    x_s      = Vector{S}(undef, n)   # host staging in solve precision (S -> T)
    U_stages = [fill!(similar(u), zero(T)) for _ in 1:k]
    S_store  = [fill!(similar(u), zero(T)) for _ in 1:k]
    L_store  = [fill!(similar(u), zero(T)) for _ in 1:k]

    # Raw host operator L (for L(U_i) = L * U_i on the host).
    L_raw = build_L(u, tinit, params)

    # Per-stage operators I - λ_i L assembled on the host (double) and uploaded to
    # the device as CSR in the SOLVE precision S (deduplicated by λ_i). The host
    # sparse copy is kept only so a `bcs_fun!` that mutates the matrix has
    # something to act on; the solve uses the device CSR.
    λs    = [Δt * A_RK_tilde[i, i] for i in 1:k]
    A0_h  = L_update(u, tinit, λs[1])
    A0_d  = JaccSparseCSR(A0_h; index_type = index_type, float_type = S)
    Ahost = Vector{typeof(A0_h)}(undef, k)
    Adev  = Vector{typeof(A0_d)}(undef, k)
    Ahost[1] = A0_h
    Adev[1]  = A0_d
    for i in 2:k
        reuse = 0
        for j in 1:i-1
            if λs[j] == λs[i]
                reuse = j
                break
            end
        end
        if reuse > 0
            Ahost[i] = Ahost[reuse]
            Adev[i]  = Adev[reuse]
        else
            Ahost[i] = L_update(u, tinit, λs[i])
            Adev[i]  = JaccSparseCSR(Ahost[i]; index_type = index_type, float_type = S)
        end
    end

    # Device vectors for the solve, in precision S (CuArrays when JACC backend == CUDA).
    b_d    = JACC.array(zeros(S, n))
    x_d    = JACC.array(zeros(S, n))
    work_d = ntuple(_ -> JACC.array(zeros(S, n)), 6)
    itmax  = get(inputs, :imex_jacc_itmax, max(500, 10 * Int(ceil(sqrt(npoint)))))

    # Where is the GPU load? The solve arrays are a CuArray (GPU) when JACC is on
    # its CUDA backend, or a plain Vector (CPU) otherwise — printing their type is
    # direct proof of where (I - λL)x = b actually runs.
    on_gpu = !(x_d isa Array)
    n_ops  = length(unique(λs))
    if rank == 0
        println()
        println(" #   ┌─ IMEX/JACC — implicit-solve GPU offload ───────────────────────")
        println(" #   │   CPU (host) : mesh, sem_setup, explicit rhs! S(u), stage-RHS")
        println(" #   │                assembly, L(u) mat-vec, solution update, output")
        println(" #   │   GPU (dev)  : per-stage implicit solve (I - λL)x = b  [BiCGSTAB]")
        println(" #   │   solve array type : ", typeof(x_d))
        println(" #   │   solve precision  : ", S, "   (host pipeline in ", T, ")")
        println(" #   │   running on GPU?  : ", on_gpu ? "YES" : "NO  (JACC CPU backend / no GPU)")
        println(" #   │   device operators : ", n_ops, " distinct (I - λL) CSR uploaded once")
        println(" #   └────────────────────────────────────────────────────────────────")
        if !on_gpu
            @warn """
            IMEX/JACC offload is running on the CPU (JACC array type = $(typeof(x_d))).
            JACC picks its backend from a COMPILE-TIME preference, so it cannot be
            switched mid-session — `enable_cuda!()` / `backend = :cuda` only writes the
            preference for the NEXT session. To run the solve on the GPU, do this ONCE
            and then RESTART Julia:
                using CUDA, JACC
                JACC.set_backend("cuda")        # writes LocalPreferences.toml
                # ... exit and restart Julia ...
            then `using CUDA, JACC` again and re-run. Verify with `Jexpresso.jacc_status()`
            that the array type is a CuArray before the run.
            """
        end
        println()
    end

    # Diagnostics accumulators for this run.
    nsolve    = 0
    tot_iters = 0
    nonconv   = 0
    res_sum   = 0.0
    res_max   = 0.0
    solve_ns  = UInt64(0)

    t_n          = tinit
    n_step       = 0
    next_out_idx = 1
    iout         = 0

    while (abs(t_n - tend) > 1.0e-14 && t_n < tend)
        for i = 1:k
            time_tilde = t_n + c_RK_tilde[i] * Δt

            # `trace` prints the CPU↔GPU hand-off for the FIRST step only, so the
            # user can see the alternation without flooding the screen over every
            # step. The actual device residency is the same on every step.
            trace = (rank == 0 && n_step == 0)

            # Stage RHS on the host: u + Σ_{j<i} Δt[ A_ex S_j + (A_im - A_ex) L_j ]
            copyto!(rhs_h, u)
            for j = 1:i-1
                aex = Δt * A_RK[i, j]
                bim = Δt * (A_RK_tilde[i, j] - A_RK[i, j])
                @. rhs_h += aex * S_store[j] + bim * L_store[j]
            end
            if bcs_fun! !== nothing
                bcs_fun!(rhs_h, Ahost[i], time_tilde, params, sem, qp)
            end
            trace && println(@sprintf("   [CPU] step 0 stage %d: assembled stage RHS on host", i))

            # Offload the linear solve (I - λL) U_i = rhs to the device. The
            # operator is linear, so a single BiCGSTAB solve (iterating to
            # tolerance internally) is exact up to that tolerance.
            trace && print(@sprintf("   [GPU] step 0 stage %d: solving (I-λL)x=b on device [%s] ... ", i, string(S)))
            rhs_s .= rhs_h                            # T -> S (host cast)
            copyto!(b_d, rhs_s)                       # host -> device
            t0 = time_ns()
            conv, bi, rn = jacc_bicgstab!(x_d, Adev[i], b_d, work_d;
                                          rtol = nl_rtol, atol = nl_atol, itmax = itmax)
            solve_ns += time_ns() - t0
            copyto!(x_s, x_d)                         # device -> host (S)
            x_h .= x_s                                # S -> T
            nsolve    += 1
            tot_iters += bi
            conv || (nonconv += 1)
            rnf = Float64(rn); res_sum += rnf; res_max = max(res_max, rnf)
            trace && println(@sprintf("done (%d iters, ‖res‖=%.2e), x copied back to host", bi, rnf))
            if !conv && rank == 0 && n_step == 0
                @warn "IMEX/JACC offload BiCGSTAB did not reach tol in $(string(S)) (stage $i, step 0): $bi iters, ‖res‖=$(rnf). Reduced precision may not converge; see end-of-run diagnostics."
            end
            U_stages[i] .+= x_h

            # Cache S(U_i) [host rhs!] and L(U_i) [host sparse mat-vec].
            ti = t_n + c_RK[i] * Δt
            S_fun!(S_store[i], U_stages[i], ti, params, sem)
            mul!(L_store[i], L_raw, U_stages[i])
            trace && println(@sprintf("   [CPU] step 0 stage %d: S(U)=rhs! and L(U) mat-vec on host", i))
        end

        # Solution update: u += Σ_i Δt b_RK[i] S(U_i)
        for i = 1:k
            axpy!(Δt * b_RK[i], S_store[i], u)
        end

        t_n    += Δt
        n_step += 1

        while next_out_idx <= length(dosetimes) && t_n + 1.0e-10 >= dosetimes[next_out_idx]
            iout += 1
            u2uaux!(params.uaux, u, neqs, npoint)
            println_rank(@sprintf(" #   IMEX(JACC offload): t = %.6f   step = %d", t_n, n_step); msg_rank = rank)
            _imex_write(params, u, t_n, iout, inputs)
            next_out_idx += 1
        end

        for j = 1:k
            fill!(U_stages[j], zero(eltype(U_stages[j])))
        end
    end

    # Final snapshot if the schedule did not already cover t_end.
    if iout == 0 || next_out_idx <= length(dosetimes)
        iout += 1
        u2uaux!(params.uaux, u, neqs, npoint)
        _imex_write(params, u, t_n, iout, inputs)
    end

    # Stash and print this run's diagnostics (collected by run_imex_precision_study).
    unorm = sqrt(Float64(sum(abs2, u)))
    diag = (precision     = S,
            on_gpu        = on_gpu,
            nsteps        = n_step,
            nsolve        = nsolve,
            total_iters   = tot_iters,
            avg_iters     = nsolve > 0 ? tot_iters / nsolve : 0.0,
            mean_resnorm  = nsolve > 0 ? res_sum / nsolve : 0.0,
            max_resnorm   = res_max,
            nonconverged  = nonconv,
            solve_seconds = solve_ns / 1.0e9,
            final_unorm   = unorm)
    IMEX_JACC_LAST_DIAG[] = diag

    if rank == 0
        println()
        println(" #   ┌─ IMEX/JACC offload diagnostics ────────────────────────────────")
        @printf(" #   │   solve precision     : %s   (host pipeline %s, device %s)\n",
                string(S), string(T), on_gpu ? "GPU" : "CPU")
        @printf(" #   │   implicit solves     : %d   over %d steps\n", nsolve, n_step)
        @printf(" #   │   BiCGSTAB iterations : %d total, %.2f avg/solve\n", tot_iters, diag.avg_iters)
        @printf(" #   │   residual ‖b-Ax‖     : %.3e mean, %.3e max\n", diag.mean_resnorm, res_max)
        @printf(" #   │   non-converged solves: %d  (of %d)\n", nonconv, nsolve)
        @printf(" #   │   solve wall time     : %.3f s   (%.1f µs/solve)\n",
                diag.solve_seconds, nsolve > 0 ? diag.solve_seconds / nsolve * 1.0e6 : 0.0)
        @printf(" #   │   final ‖u‖₂           : %.10e\n", unorm)
        println(" #   └────────────────────────────────────────────────────────────────")
        println()
    end

    return nothing
end

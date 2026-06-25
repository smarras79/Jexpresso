#==============================================================================
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
#==============================================================================

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
    if delta == 1
        L_fun!  = inputs[:L_fun]
        build_L = inputs[:build_L]
        upd_L   = inputs[:upd_L]
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

    #==========================================================================
    # Runge-Kutta IMEX
    #==========================================================================
    if method == "RK"
        for key in (:A_RK, :b_RK, :c_RK, :A_RK_tilde, :b_RK_tilde, :c_RK_tilde)
            haskey(coeff, key) || error("IMEX/RK: missing coefficient $(key) in :coeff")
        end
        A_RK       = coeff[:A_RK];       b_RK       = coeff[:b_RK];       c_RK       = coeff[:c_RK]
        A_RK_tilde = coeff[:A_RK_tilde]; b_RK_tilde = coeff[:b_RK_tilde]; c_RK_tilde = coeff[:c_RK_tilde]

        U_stages = Dict{Int,Any}()
        for j = 1:k
            U_stages[j] = KernelAbstractions.zeros(backend, nl_precision, Int64(unkwn))
        end

        while (abs(t_n - inputs[:tend]) > 1.0e-14 && t_n < inputs[:tend])
            for i = 1:k
                time_tilde = t_n + c_RK_tilde[i] * Δt
                λ          = Δt * A_RK_tilde[i, i]

                # Implicit operator for this stage
                L_curr = L_update(U_stages[i], time_tilde, λ)

                # Stage RHS
                rhs = construct_rhs_rk!(rhs_buf, s_j_buf, l_j_buf,
                                        u, U_stages, t_n, i,
                                        A_RK, A_RK_tilde, c_RK)
                if bcs_fun! !== nothing
                    bcs_fun!(rhs, L_curr, time_tilde, params, sem, qp)
                end

                # Nonlinear fixed-point loop: solve L_curr U_i = rhs
                nonl_res = nonl_res_buf
                copyto!(nonl_res, rhs)
                nl_norm_0 = norm(nonl_res)
                nl_norm_k = nl_norm_0
                nl_iter   = 1
                while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                    x = imex_linsolve(L_curr, solver_precision.(nonl_res))
                    U_stages[i] .+= x

                    if delta == 1 && upd_L
                        L_curr = L_update(U_stages[i], time_tilde, λ)
                        bcs_fun! !== nothing && bcs_fun!(rhs, L_curr, time_tilde, params, sem, qp)
                    end

                    mul!(Lu_buf, L_curr, U_stages[i])
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

    #==========================================================================
    # Multistep IMEX
    #==========================================================================
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
        L_curr = L_update(u_next, t_n + Δt, lambda)

        while (abs(t_n - inputs[:tend]) > 1.0e-14 && t_n < inputs[:tend])
            rhs = construct_rhs_multistep!(rhs_buf, s_j_buf, l_j_buf,
                                           u_prev, t_n, alpha, beta, lambda)
            if bcs_fun! !== nothing
                bcs_fun!(rhs, L_curr, t_n + Δt, params, sem, qp)
            end

            nonl_res = nonl_res_buf
            copyto!(nonl_res, rhs)
            nl_norm_0 = norm(nonl_res)
            nl_norm_k = nl_norm_0
            nl_iter   = 1
            while (nl_norm_k > nl_atol && nl_norm_k > nl_rtol * nl_norm_0 && nl_iter < max_nl_iter)
                x = imex_linsolve(L_curr, solver_precision.(nonl_res))
                u_next .+= x

                if delta == 1 && upd_L
                    L_curr = L_update(u_next, t_n + Δt, lambda)
                    bcs_fun! !== nothing && bcs_fun!(rhs, L_curr, t_n + Δt, params, sem, qp)
                end

                mul!(Lu_buf, L_curr, u_next)
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

function imex_time_loop!(inputs, sem, qp, params, u)
    # Time step
    Δt = inputs[:Δt]

    # number of discretization points
    npoint = sem.mesh.npoin
    # number of equations
    neqs = params.neqs
    # total number of unknown
    unkwn = neqs * npoint

    # Initial time
    t_n = inputs[:tinit]

    # Linear solver
    haskey(inputs, :lsolve) ? lsolve = kwargs[:lsolve] : lsolve = nothing
    # Solver parameters
    if lsolve == nothing
        haskey(inputs, :sp) ? solver_parameters = inputs[:sp] : solver_parameters = nothing
        (solver_parameters != nothing && haskey(solver_parameters, :atol)) ? atol = solver_parameters[:atol] : atol = 1.e-06
        (solver_parameters != nothing && haskey(solver_parameters, :rtol)) ? rtol = solver_parameters[:rtol] : rtol = 1.e-10
        (solver_parameters != nothing && haskey(solver_parameters, :restart)) ? restart = solver_parameters[:restart] : restart = true
        (solver_parameters != nothing && haskey(solver_parameters, :memory)) ? memory = solver_parameters[:memory] : memory = 10
        (solver_parameters != nothing && haskey(solver_parameters, :itmax)) ? itmax = solver_parameters[:itmax] : itmax = 100
        (solver_parameters != nothing && haskey(solver_parameters, :verbose)) ? verbose = solver_parameters[:verbose] : verbose = 1
        (solver_parameters != nothing && haskey(solver_parameters, :prec)) ? prec = solver_parameters[:prec] : prec = SmoothedAggregationPreconBuilder()

        haskey(inputs, :solver_precision) ? solver_precision = inputs[:solver_precision] : solver_precision = Float64

        haskey(inputs, :prec_sp) ? prec_sp = inputs[:prec_sp] : prec_sp = Dict(:maxiter      => 1,
                                                                               :abstol       => 1e-8,
                                                                               :precision    => Float64,
                                                                               )
    end

    # Non-linear parameters
    haskey(inputs, :nl_atol) ? nl_atol = inputs[:nl_atol] : nl_atol = 1.e-05
    haskey(inputs, :nl_rtol) ? nl_rtol = inputs[:nl_rtol] : nl_rtol = 1.e-05
    haskey(inputs, :max_nl_iter) ? max_nl_iter = inputs[:max_nl_iter] : max_nl_iter = 10
    haskey(inputs, :nl_precision) ? nl_precision = inputs[:nl_precision] : nl_precision = Float64

    # IMEX
    # With delta = 0, apply explicit method;
    # with delta = 1, apply IMEX
    delta = inputs[:delta]
    # Type of IMEX method
    method = inputs[:method]
    # number of steps for multistep IMEX, or
    # number of stages for ARK
    k = inputs[:k]
    if k < 1
        error("k should be a positive integer")
    end
    # Coefficients of the method
    coeff = inputs[:coeff]

    # Source of the governing equations
    S_fun! = inputs[:S_fun]
    if delta == 1
        # Linear operator representing the fast waves
        L_fun! = inputs[:L_fun]

        # Construction of linear operator L
        build_L = inputs[:build_L]

        # check if L has to be updated
        haskey(inputs, :upd_L) ? upd_L = inputs[:upd_L] : upd_L = false
    end
    # BCs of the solution
    bcs_fun! = inputs[:bcs_fun]

    # coefficients of the method
    if method == "multistep"
        alpha = coeff[:alpha]
        beta = coeff[:beta]
        xi = coeff[:xi]

        lambda = xi * Δt

        # Constructing u_prev
        u_prev = Dict()
        for n_step = 1 : k
            u_prev[n_step] = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
        end
        u_prev[1] .= u

        # warm up
        if k > 1
            haskey(inputs, :Δt_expl) ? Δt_expl = inputs[:Δt_expl] : Δt_expl = Δt / (k - 1)
            for n_step = 2 : k
                rhs = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))

                rhs += u_prev[1]

                s_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
                S_fun!(s_j, u_prev[1], t_n, params)

                # explicit Euler
                rhs += Δt_expl * s_j
 
                # updating solutions
                for i_step = n_step : -1 : 2
                    u_prev[i_step] .= u_prev[i_step - 1]
                end
                u_prev[1] .= rhs
                u .= rhs

                t_n = t_n + Δt_expl
            end
        end
    elseif method == "RK"
        A_RK = coeff[:A_RK]
        b_RK = coeff[:b_RK]
        c_RK = coeff[:c_RK]

        A_RK_tilde = coeff[:A_RK_tilde]
        b_RK_tilde = coeff[:b_RK_tilde]
        c_RK_tilde = coeff[:c_RK_tilde]

        # Constructing U_stages
        U_stages = Dict()
        for j_stages = 1 : k
            U_stages[j_stages] = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
        end
    end

    #------------------------------------------------------------------------
    # Construct rhs
    #------------------------------------------------------------------------
    # Multistep
    function construct_rhs(inputs, params, u_prev, t_n)
        # building rhs
        rhs = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))

        for n_step = 1 : k
            rhs += alpha[n_step] * u_prev[n_step]

            # Evaluating S(u_{n-j}) and L(u_{n-j})
            # s_j and l_j should already be scaled by M_inv
            s_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
            time = t_n - (n_step - 1) * Δt
            S_fun!(s_j, u_prev[n_step], time, params)

            if delta == 1
                l_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
                L_fun!(l_j, u_prev[n_step], time, params)
            end

            if delta == 1
                rhs += lambda * beta[n_step] * (s_j - l_j)
            else
                rhs += lambda * beta[n_step] * s_j
            end
        end

        return rhs
    end

    # Runge-Kutta
    function construct_rhs(inputs, params, u_prev, U_stages, t_n, i)
        # building rhs
        rhs = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))

        rhs += u_prev

        for j_stages = 1 : i - 1
            # Evaluating S(u_{n-j}) and L(u_{n-j})
            # s_j and l_j should already be scaled by M_inv
            s_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
            time = t_n + c_RK[j_stages] * Δt
            S_fun!(s_j, U_stages[j_stages], time, params)

            if delta == 1
                time_tilde = t_n + c_RK_tilde[j_stages] * Δt
                l_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
                L_fun!(l_j, U_stages[j_stages], time, params)
            end

            if delta == 1
                rhs += Δt * (A_RK[i, j_stages] * s_j + (A_RK_tilde[i, j_stages] - A_RK[i, j_stages]) * l_j)
            else
                rhs += Δt * A_RK[i, j_stages] * s_j
            end
        end

        return rhs
    end


    #------------------------------------------------------------------------
    # Construct L
    #------------------------------------------------------------------------
    # Multistep
    function L_update(u, t_n)
        Is = sparse(1.0I, Int64(unkwn), Int64(unkwn))

        if delta == 0
            L = Is
        else
            # build_L needs to have included the scaling by inverse of mass matrix
            L_temp = build_L(u, t_n, params)

            L = Is - lambda * L_temp
        end

        return L
    end

    # Runge-Kutta
    function L_update(u, t_n, a_tilde_ii)
        Is = sparse(1.0I, Int64(unkwn), Int64(unkwn))

        if delta == 0
            L = Is
        else
            # build_L needs to have included the scaling by inverse of mass matrix
            L_temp = build_L(u, t_n, params)

            L = Is - Δt * a_tilde_ii * L_temp
        end

        return L
    end

    #------------------------------------------------------------------------
    # Updating u
    #------------------------------------------------------------------------
    function update_u!(u, u_prev, u_next)
        for n_step = k : -1 : 2
            u_prev[n_step] .= u_prev[n_step - 1]
        end
        u_prev[1] .= u_next
        u .= u_next
    end


    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)

    n_step = 0

    # Multistep
    if method == "multistep"
        # Initialize solution vector
        u_next = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))

        #------------------------------------------------------------------------
        # Building L_curr
        #------------------------------------------------------------------------
        L_curr = L_update(u_next, t_n + Δt)

        #------------------------------------------------------------------------
        # Simulation
        #------------------------------------------------------------------------
        while (abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend])
            println(string("Solving for time ", t_n + Δt))

            rhs = construct_rhs(inputs, params, u_prev, t_n)
            # Apply bcs to rhs and L
            bcs_fun!(rhs, L_curr, t_n + Δt, params, sem, qp)

            # construct vector containing non-linear residual
            nonl_res = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
            nonl_res .= rhs

            nl_norm_0 = norm(nonl_res)
            nl_norm_k = nl_norm_0

            nl_iter = 1

            #------------------------------------------------------------------------
            # Non-linear loop
            #------------------------------------------------------------------------
            while (nl_norm_k > nl_atol && (nl_norm_k > nl_rtol * nl_norm_0) && nl_iter < max_nl_iter)
                println(string("Non-linear solver: iteration ", nl_iter,", non-linear residual ", nl_norm_k))

                @info "Solving linear system..."

                # Solving the linear system
                if lsolve == nothing
                    prec = MyPrecClass.MyPrec(prec_sp[:precision].(solver_precision.(L_curr)),
                                              RugeStubenAMG(), prec_sp)

                    (x, stats) = fgmres(solver_precision.(L_curr), solver_precision.(nonl_res),
                                        memory=memory, N=prec, ldiv=true,
                                        restart=restart, atol=atol, rtol=rtol,
                                        itmax=itmax, verbose = verbose)
                else
                    x = lsolve(solver_precision.(L_curr),
                               solver_precision.(nonl_res))
                end
                @info "Solving linear system... DONE"

                u_next += x

                #------------------------------------------------------------------------
                # Updating L_curr
                #------------------------------------------------------------------------
                if delta == 1 && upd_L
                        L_curr = L_update(u_next, t_n + Δt)

                        # Apply bcs to rhs and L
                        bcs_fun!(rhs, L_curr, t_n + Δt, params, sem, qp)
                end

                nonl_res .= rhs - L_curr * u_next
                nl_norm_k = norm(nonl_res)

                nl_iter += 1
            end
            #------------------------------------------------------------------------
            # End non-linear loop
            #------------------------------------------------------------------------

            # Update solution
            update_u!(u, u_prev, u_next)
        end

        t_n += Δt
        n_step += 1

        write_output(params.SD, u, params.uaux, t_n, n_step,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe)

        # Initialize solution vector
        u_next = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))

        #------------------------------------------------------------------------
        # Updating L_curr
        #------------------------------------------------------------------------
        if delta == 1 && upd_L && abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend]
            L_curr = L_update(u_next, t_n + Δt)
        end
    elseif method == "RK"
        #------------------------------------------------------------------------
        # Simulation
        #------------------------------------------------------------------------
        while (abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend])
            println(string("Solving for time ", t_n + Δt))
            for i_stages = 1 : k
                @info "Stage " i_stages
                time_tilde = t_n + c_RK_tilde[i_stages] * Δt

                #------------------------------------------------------------------------
                # Building L_curr
                #------------------------------------------------------------------------
                L_curr = L_update(U_stages[i_stages], time_tilde, A_RK_tilde[i_stages, i_stages])

                #------------------------------------------------------------------------
                # Building rhs
                #------------------------------------------------------------------------
                rhs = construct_rhs(inputs, params, u, U_stages, t_n, i_stages)
                # Apply bcs to rhs and L
                bcs_fun!(rhs, L_curr, time, params, sem, qp)

                # construct vector containing non-linear residual
                nonl_res = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
                nonl_res .= rhs

                nl_norm_0 = norm(nonl_res)
                nl_norm_k = nl_norm_0

                nl_iter = 1

                #------------------------------------------------------------------------
                # Non-linear loop
                #------------------------------------------------------------------------
                while (nl_norm_k > nl_atol && (nl_norm_k > nl_rtol * nl_norm_0) && nl_iter < max_nl_iter)
                    println(string("Non-linear solver: iteration ", nl_iter,", non-linear residual ", nl_norm_k))

                    @info "Solving linear system..."

                    # Solving the linear system
                    if lsolve == nothing
                        prec = MyPrecClass.MyPrec(prec_sp[:precision].(solver_precision.(L_curr)),
                                                  RugeStubenAMG(), prec_sp)

                        (x, stats) = fgmres(solver_precision.(L_curr), solver_precision.(nonl_res),
                                            memory=memory, N=prec, ldiv=true,
                                            restart=restart, atol=atol, rtol=rtol,
                                            itmax=itmax, verbose = verbose)
                    else
                        x = lsolve(solver_precision.(L_curr),
                                   solver_precision.(nonl_res))
                    end
                    @info "Solving linear system... DONE"

                    U_stages[i_stages] += x

                    #------------------------------------------------------------------------
                    # Updating L_curr
                    #------------------------------------------------------------------------
                    if delta == 1 && upd_L
                            L_curr = L_update(U_stages[i_stages], time, A_RK_tilde[i, i])

                            # Apply bcs to rhs and L
                            bcs_fun!(rhs, L_curr, time, params, sem, qp)
                    end

                    nonl_res .= rhs - L_curr * U_stages[i_stages]
                    nl_norm_k = norm(nonl_res)

                    nl_iter += 1
                end
                #------------------------------------------------------------------------
                # End non-linear loop
                #------------------------------------------------------------------------
            end

            # Update solution
            for i_stages = 1 : k
                s_j = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
                time = t_n + c_RK[i_stages] * Δt
                S_fun!(s_j, U_stages[i_stages], time, params)
                u += Δt * b_RK[i_stages] * s_j
            end

            t_n += Δt
            n_step += 1

            write_output(params.SD, u, params.uaux, t_n, n_step,
                         params.mesh, params.mp,
                         params.connijk_original, params.poin_in_bdy_face_original,
                         params.x_original, params.y_original, params.z_original,
                         inputs[:output_dir], inputs,
                         params.qp.qvars, params.qp.qoutvars,
                         inputs[:outformat];
                         nvar=params.qp.neqs, qexact=params.qp.qe)

            # Initialize stages
            for j_stages = 1 : k
                U_stages[j_stages] = KernelAbstractions.zeros(inputs[:backend], nl_precision, Int64(unkwn))
            end
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return u
end

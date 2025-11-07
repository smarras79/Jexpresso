function imex_time_loop!(inputs, sem, qp, params, u)
    # Time step
    Δt = inputs[:Δt]

    # number of discretization points
    npoint = sem.mesh.npoin
    # number of equations
    neqs = params.neqs
    # total number of unknown
    unkwn = neqs * npoint

    # IMEX
    # With delta = 0, apply explicit method;
    # with delta = 1, apply IMEX
    delta = inputs[:delta]
    # Type of IMEX method
    method = inputs[:method]
    # number of steps for multistep IMEX, or
    # number of stages for ARK
    k = inputs[:k]
    # Coefficients of the method
    coeff = inputs[:coeff]

    # coefficients of the method
    if method == "multistep"
        alpha = coeff[:alpha]
        beta = coeff[:beta]
        xi = coeff[:xi]

        lambda = xi * Δt

        # Constructing u_prev
        u_prev = Dict()
        for n_step = 1 : k
            # initialize u_prev in a better way for multistep
            u_prev[n_step] = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
        end
        u_prev[1] .= u
    end

    # Source of the governing equations
    S_fun! = inputs[:S_fun]
    if delta == 1
        # Linear operator representing the fast waves
        L_fun! = inputs[:L_fun]

        # Construction of linear operator L
        build_L! = inputs[:build_L]

        # check if L has to be updated
        haskey(inputs, :upd_L) ? upd_L = inputs[:upd_L] : upd_L = false
    end
    # BCs of the solution
    bcs_fun! = inputs[:bcs_fun]

    # Linear solver
    haskey(inputs, :lsolve) ? lsolve = kwargs[:lsolve] : lsolve = nothing
    # Solver parameters
    if lsolve == nothing
        haskey(inputs, :sp) ? solver_parameters = inputs[:sp] : solver_parameters = nothing
        (solver_parameters != nothing && haskey(solver_parameters, :atol)) ? atol = solver_parameters[:atol] : atol = 1.e-08
        (solver_parameters != nothing && haskey(solver_parameters, :rtol)) ? rtol = solver_parameters[:rtol] : rtol = 1.e-10
        (solver_parameters != nothing && haskey(solver_parameters, :restart)) ? restart = solver_parameters[:restart] : restart = true
        (solver_parameters != nothing && haskey(solver_parameters, :memory)) ? memory = solver_parameters[:memory] : memory = 10
        (solver_parameters != nothing && haskey(solver_parameters, :verbose)) ? verbose = solver_parameters[:verbose] : verbose = 1
        (solver_parameters != nothing && haskey(solver_parameters, :prec)) ? prec = solver_parameters[:prec] : prec = SmoothedAggregationPreconBuilder()
    end

    # Non-linear parameters
    haskey(inputs, :nl_atol) ? nl_atol = inputs[:nl_atol] : nl_atol = 1.e-05
    haskey(inputs, :nl_rtol) ? nl_rtol = inputs[:nl_rtol] : nl_rtol = 1.e-05
    haskey(inputs, :max_nl_iter) ? max_nl_iter = inputs[:max_nl_iter] : max_nl_iter = 100

    # inverse of mass matrix
    M_inv = diagm(sem.matrix.Minv)

    #------------------------------------------------------------------------
    # Construct rhs
    #------------------------------------------------------------------------
    function construct_rhs(inputs, sem, u_prev, t_n)
        if method == "multistep"
            # building rhs
            rhs = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

            for n_step = 1 : k
                rhs += alpha[n_step] * u_prev[n_step]

                # Evaluating S(u_{n-j}) and L(u_{n-j})
                s_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                time = t_n - n_step * Δt
                # S_fun needs to have included the scaling by inverse of mass matrix
                S_fun!(s_j, u_prev[n_step], time, params)

                if delta == 1
                    l_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                    L_fun!(l_j, u_prev[n_step], time, params)
                    l_j = M_inv * l_j
                end

                if delta == 1
                    rhs += lambda * beta[n_step] * (s_j - l_j)
                else
                    rhs += lambda * beta[n_step] * s_j
                end
            end

        end

        return rhs
    end


    #------------------------------------------------------------------------
    # Construct L
    #------------------------------------------------------------------------
    function L_update(u, t_n)
        if delta == 1
#            L_temp = build_L!(u, t_n + Δt, params)
#            L_temp = sem.matrix.L
#            L = - lambda * M_inv * L_temp
            # L_update needs to have included the scaling by inverse of mass matrix
            L = - lambda * L_temp
        else
            L = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn), Int64(unkwn))
        end
#L_temp = sem.matrix.L
#L = - lambda * L_temp
#@info typeof(L)
L = SparseArrays.sparse(L)
#@info typeof(L)

#        for i = 1 : unkwn
#            L[i, i] += 1.
#        end
        Is = sparse(I, Int64(unkwn), Int64(unkwn))
        L += Is

        return L
    end


    # Initialize solution vector
    u_next = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

    #------------------------------------------------------------------------
    # Updating u
    #------------------------------------------------------------------------
    function update_u!(u, u_prev, u_next)
        for n_step = k : -1 : 2
            u_prev[n_step] .= u_prev[n_step - 1]
        end
        u_prev[1] .= u
        u .= u_next
    end


    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    #------------------------------------------------------------------------
    # Initial time
    #------------------------------------------------------------------------
    t_n = inputs[:tinit]
    n_step = 0

    #------------------------------------------------------------------------
    # Building L_curr
    #------------------------------------------------------------------------
    L_curr = L_update(u_next, t_n + Δt)

    #------------------------------------------------------------------------
    # Simulation
    #------------------------------------------------------------------------
    while (abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend])
        rhs = construct_rhs(inputs, params, u_prev, t_n)
        # construct vector containing non-linear residual
        nonl_res = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
        nonl_res .= rhs

        nl_norm_0 = norm(nonl_res)
        nl_norm_k = nl_norm_0

        nl_iter = 1

        if method == "multistep"
            #------------------------------------------------------------------------
            # Non-linear loop
            #------------------------------------------------------------------------
            while (nl_norm_k > nl_atol && (nl_norm_k / nl_norm_0) > nl_rtol && nl_iter < max_nl_iter)
                # Apply bcs to non-linear residual
#                bcs_fun!(nonl_res, t_n + Δt, params)    # check this, it does not work

                # Solving the linear system
#                sol = solveAx(Lmultiply, rhs, inputs[:lsolver])
#                prob = LinearProblem(Lmultiply, rhs);    
#                sol = solve(prob, inputs[:lsolver])
                if lsolve == nothing
                    prob = LinearProblem(L_curr, nonl_res);
                    prec = SmoothedAggregationPreconBuilder()
                    strategy = KrylovJL_GMRES(
                        FgmresWorkspace(L_curr, nonl_res, memory = memory),
                        restart = restart, verbose = verbose,
                        precs = prec)
                    sol = solve(prob, strategy, atol = atol, rtol = rtol)
                else
                    sol = lsolve(L_curr, nonl_res)
                end

                u_next += sol.u

                #------------------------------------------------------------------------
                # Updating L_curr
                #------------------------------------------------------------------------
                if delta == 1
                    if (upd_L && abs(t_n - inputs[:tend]) > 1.e-14)
                        L_curr = L_update(u_next, t_n + Δt)
                        L_curr = SparseArrays.sparse(L_curr)
                    end
                end

                nonl_res = rhs - L_curr * u_next
                nl_norm_k = norm(nonl_res)

                nl_iter += 1
            end
            #------------------------------------------------------------------------
            # End non-linear loop
            #------------------------------------------------------------------------

            # Update solution
            update_u!(u, u_prev, u_next)

            # Initialize solution vector
            u_next = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
        end

        t_n += Δt
        n_step += 1
@info t_n
@info n_step

        write_output(params.SD, u, params.uaux, t_n, n_step,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe)

        # Initialize solution vector
        u_next = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

        #------------------------------------------------------------------------
        # Updating L_curr
        #------------------------------------------------------------------------
        if delta == 1
            if (upd_L && abs(t_n - inputs[:tend]) > 1.e-14)
                L_curr = L_update(u_next, t_n + Δt)
                L_curr = SparseArrays.sparse(L_curr)
            end
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return u
end

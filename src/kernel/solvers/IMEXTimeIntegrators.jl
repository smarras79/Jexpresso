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
    end

    # Non-linear parameters
    haskey(inputs, :nl_atol) ? nl_atol = inputs[:nl_atol] : nl_atol = 1.e-05
    haskey(inputs, :nl_rtol) ? nl_rtol = inputs[:nl_rtol] : nl_rtol = 1.e-05
    haskey(inputs, :max_nl_iter) ? max_nl_iter = inputs[:max_nl_iter] : max_nl_iter = 10

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
            # initialize u_prev in a better way for multistep
            u_prev[n_step] = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
        end
        u_prev[1] .= u

        # warm up
        if k > 1
            haskey(inputs, :Δt_expl) ? Δt_expl = inputs[:Δt_expl] : Δt_expl = Δt / (k - 1)
            for n_step = 2 : k
                rhs = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

                rhs += u_prev[1]

                s_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
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
    end

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
                # s_j and l_j should already be scaled by M_inv
                s_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                time = t_n - (n_step - 1) * Δt
                S_fun!(s_j, u_prev[n_step], time, params)

                if delta == 1
                    l_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                    L_fun!(l_j, u_prev[n_step], time, params)
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


    # Initialize solution vector
    u_next = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

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

    #------------------------------------------------------------------------
    # Building L_curr
    #------------------------------------------------------------------------
    L_curr = L_update(u_next, t_n + Δt)
    L_solver = L_curr

    #------------------------------------------------------------------------
    # Simulation
    #------------------------------------------------------------------------
    while (abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend])
        println(string("Solving for time ", t_n + Δt))

        rhs = construct_rhs(inputs, params, u_prev, t_n)
        # Apply bcs to rhs and L
        bcs_fun!(rhs, L_solver, t_n + Δt, params, sem, qp)

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
            while (nl_norm_k > nl_atol && (nl_norm_k > nl_rtol * nl_norm_0) && nl_iter < max_nl_iter)
                println(string("Non-linear solver: iteration ", nl_iter,", non-linear residual ", nl_norm_k))

                @info "Solving linear system..."

                # Solving the linear system
                if lsolve == nothing
                    prob = LinearProblem(L_solver, nonl_res);
#                    prob = LinearProblem(L_curr, nonl_res);
                    prec = SmoothedAggregationPreconBuilder()
                    workspace = FgmresWorkspace(L_solver, nonl_res; memory = memory)
#                    workspace = FgmresWorkspace(L_curr, nonl_res; memory = memory)


#function prec!(A, b)
#    boh = SmoothedAggregationPreconBuilder(A)
#    sol = AlgebraicMultigrid._solve(boh, b, maxiter = 1, abstol = 1e-6)
#    x = sol.u
#    return x
#end
#prec = prec!


                    strategy = KrylovJL_GMRES(
                        workspace;
#                        itmax = itmax,
                        restart = restart, verbose = verbose,
                        precs = prec)
                    sol = solve(prob, strategy, atol = atol, rtol = rtol)

#mg = ruge_stuben(L_solver)
#prec = aspreconditioner(mg)
#prec_sp = Dict(
#    :maxiter => 1,
#    :abstol  => 1e-8,
#    )
#prec = MyPrecClass.MyPrec(L_solver, RugeStubenAMG(), prec_sp)

#(x, stats) = fgmres(L_solver, nonl_res, memory=memory, N=prec, ldiv=true,
#                    restart=restart, atol=atol, rtol=rtol,
#                    itmax=itmax, verbose = verbose)
#@info stats.solved
#@info stats.niter
                else
                    sol = lsolve(L_curr, nonl_res)
                end
                @info "Solving linear system... DONE"

                u_next += sol.u
#                u_next += x

                #------------------------------------------------------------------------
                # Updating L_curr
                #------------------------------------------------------------------------
                if delta == 1 && upd_L
                        L_curr = L_update(u_next, t_n + Δt)
                        L_solver = L_curr

                        # Apply bcs to rhs and L
                        bcs_fun!(rhs, L_solver, t_n + Δt, params, sem, qp)
                end

                nonl_res .= rhs - L_solver * u_next
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
        if delta == 1 && upd_L && abs(t_n - inputs[:tend]) > 1.e-14 && t_n < inputs[:tend]
            L_curr = L_update(u_next, t_n + Δt)
            L_solver = L_curr
        end
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return u
end

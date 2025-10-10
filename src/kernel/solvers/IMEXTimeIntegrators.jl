function imex_time_loop!(inputs, sem, qp, params, u)
    # Time step
    Δt = inputs[:Δt]

    # number of discretization points
    npoint = sem.mesh.npoin
    # number of equations
    neqs = params.neqs
    # total number of unknown
    unkwn = neqs * npoint

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

    # RHS of the governing equations
    S_fun = inputs[:S_fun]
    if delta == 1
        # Linear operator representing the fast waves
        L_fun = inputs[:L_fun]
    end

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
    end

    # inverse of mass matrix
    M_inv = diagm(sem.matrix.Minv)
#    # to be changed as function provided by user
#    S = M_inv * sem.matrix.L

    #------------------------------------------------------------------------
    # Construct rhs
    #------------------------------------------------------------------------
    function construct_rhs(inputs, sem, u_prev, t_n)
        if method == "multistep"
            # building rhs
            rhs = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))

            for n_step = 1 : k
                for i = 1 : npoint
                    rhs[i] += alpha[n_step] * u_prev[n_step][i]    # check this expression
                end

                # Evaluating S(u_{n-j}) and L(u_{n-j})
#                s_j = S * u_prev[n_step]    # to be changed as function provided by user
                s_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                time = t_n - n_step * Δt
                S_fun!(s_j, u_prev[:n_step], params, time)
#                s_j = M_inv * s_j    # already done in the application of the source operator

                if delta == 1
                    l_j = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn))
                    L_fun!(l_j, u_prev[n_step], params, time)
                    l_j = M_inv * l_j
                end

                if delta == 1
                    for i = 1 : npoint
                        rhs[i] += lambda * beta[n_step] * (s_j[i] - l_j[i])    # check this expression
                    end
                else
                    for i = 1 : npoint
                        rhs[i] += lambda * beta[n_step] * s_j[i]
                    end
                end
            end

        end

        return rhs
    end


    #------------------------------------------------------------------------
    # Construct L
    #------------------------------------------------------------------------
    function L_curr(u, t_n)
        if delta == 1
            L_temp = L_fun(u, t_n + Δt)
            L = - lambda * M_inv * L_temp
        else
            L = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(unkwn), Int64(unkwn))
        end

        for i = 1 : npoint
            L[i, i] += 1.
        end

        return L
    end


    #------------------------------------------------------------------------
    # Updating u
    #------------------------------------------------------------------------
    function update_u!(u, u_prev, u_sol)
        for n_step = 1 : k - 1
            u_prev[n_step] = u_prev[n_step + 1]
        end
        u_prev[k] = u
        u = u_sol
    end


    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    #
    # Simulation
    #
    t_n = inputs[:tinit]

    while t_n < inputs[:tend]
        rhs = construct_rhs(inputs, params, u_prev, t_n)

        if method == "multistep"
            L = L_curr(u, t_n)

            # Apply bcs to rhs
            apply_boundary_conditions_lin_solve!(sem.matrix.L, 0.0, params.qp.qe,
                                                 params.mesh.x, params.mesh.y, params.mesh.z,
                                                 params.metrics.nx,
                                                 params.metrics.ny,
                                                 params.metrics.nz,
                                                 sem.mesh.npoin, params.mesh.npoin_linear, 
                                                 params.mesh.poin_in_bdy_edge,
                                                 params.mesh.poin_in_bdy_face,
                                                 params.mesh.nedges_bdy,
                                                 params.mesh.nfaces_bdy,
                                                 params.mesh.ngl, params.mesh.ngr,
                                                 params.mesh.nelem_semi_inf,
                                                 params.basis.ψ, params.basis.dψ,
                                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 rhs, 0.0, params.ubdy,
                                                 params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                                 params.mesh.bdy_edge_type,
                                                 params.ω, qp.neqs, params.inputs, params.AD, sem.mesh.SD)

            sol = solveAx(L, rhs, inputs[:lsolver])
            u_sol = sol.u

            update_u!(u, u_prev, u_sol)
        end

        t_n += Δt
    end
    
    println_rank(" # Solving ODE  ................................ DONE"; msg_rank = rank)
    
    return u
end

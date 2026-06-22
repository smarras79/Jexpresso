function solveAx(L, RHS, linear_solver...)
    
    prob = LinearProblem(L, RHS);    
    sol = solve(prob, linear_solver)
    
    return sol
end


#function solveAx_sparse(A, b)
#     return A\b
#end

function standard_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    RHS   = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))
    Mdiag = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))
    
    #-----------------------------------------------------
    # L*q = M*RHS   See algo 12.18 of Giraldo's book
    #-----------------------------------------------------
    for ip =1:sem.mesh.npoin
        RHS[ip] = user_source!(RHS[ip],
                               params.qp.qn[ip],
                               params.qp.qe[ip],
                               sem.mesh.npoin,
                               inputs[:CL], inputs[:SOL_VARS_TYPE];
                               neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                               xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                               ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
    end
    RHS = sem.matrix.M.*RHS

    if inputs[:lsparse] ==  false
        for ip = 1:sem.mesh.npoin
            sem.matrix.L[ip,ip] += inputs[:rconst][1]
        end
    end

    apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                         0.0, params.qp.qe,
                                         params.mesh.coords,
                                         params.metrics.nx,
                                         params.metrics.ny,
                                         params.metrics.nz,
                                         sem.mesh.npoin,
                                         params.mesh.npoin_linear,
                                         params.mesh.poin_in_bdy_edge,
                                         params.mesh.poin_in_bdy_face,
                                         params.mesh.nedges_bdy,
                                         params.mesh.nfaces_bdy,
                                         params.mesh.ngl, params.mesh.ngr,
                                         params.mesh.nelem_semi_inf,
                                         params.basis.ψ, params.basis.dψ,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         RHS, 0.0, params.ubdy,
                                         params.mesh.connijk_lag,
                                         params.mesh.bdy_edge_in_elem,
                                         params.mesh.bdy_edge_type,
                                         params.ω, qp.neqs,
                                         params.inputs, params.AD, sem.mesh.SD)

    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage ..............")))

    #solAxb = @btime solveAx($sem.matrix.L, $RHS, inputs[:ode_solver])
    #sol = solAxb.u
    solAxb = sem.matrix.L \ RHS
    sol = solAxb

    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage .............. DONE")))
    args = (params.SD, sol, params.uaux, 1, 1,
            sem.mesh, nothing,
            nothing, nothing,
            0.0, 0.0, 0.0,
            OUTPUT_DIR, inputs,
            params.qp.qvars,
            params.qp.qoutvars,
            inputs[:outformat])

    # Automatic L2-error check against the exact field qe (e.g. a manufactured
    # solution). No-op when qe carries no exact field.
    print_solution_L2_error(sol, params.qp.qe, sem.matrix.M, sem.mesh.npoin;
                            label="direct solve")

    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe)

    return nothing
end

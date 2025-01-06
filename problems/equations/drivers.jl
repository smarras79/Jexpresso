using HDF5

function driver(nparts,
                distribute,
                inputs::Dict,
                OUTPUT_DIR::String,
                TFloat) 
    comm  = distribute.comm
    rank = MPI.Comm_rank(comm)
    sem = sem_setup(inputs, nparts, distribute)
    
    if (inputs[:backend] != CPU())
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
    end

    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
    # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)
    if inputs[:ladapt] == true
        if rank == 0
            @info "start conformity4ncf_q!"
        end
        @time conformity4ncf_q!(qp.qn, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        @time conformity4ncf_q!(qp.qe, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        MPI.Barrier(comm)
        if rank == 0
            @info "end conformity4ncf_q!"
        end
    end

    if inputs[:ladapt] == true
        amr_freq = inputs[:amr_freq]
        Δt_amr   = amr_freq * inputs[:Δt]
        tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
    else
        tspan = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
    end
    params, u =  params_setup(sem,
                              qp,
                              inputs,
                              OUTPUT_DIR,
                              TFloat,
                              tspan)
    
    if !inputs[:llinsolve]
        #
        # Hyperbolic/parabolic problems that lead to Mdq/dt = RHS
        #
        solution = time_loop!(inputs, params, u)
        
        if (inputs[:ndiagnostics_outputs] > 0)
            write_output(sem.mesh.SD, solution,  sem.mesh,
                         OUTPUT_DIR, inputs,
                         params.qp.qvars,
                         inputs[:outformat];
                         nvar=params.qp.neqs, qexact=params.qp.qe, case="rtb")
        end
        
    else
        #
        # Problems that lead to Ax = b
        #
        RHS = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin), qp.neqs)

        if (inputs[:backend] == CPU())
          
            Minv = diagm(sem.matrix.Minv)
            
            L_temp = Minv * sem.matrix.L
            sem.matrix.L .= L_temp
            
            for ip =1:sem.mesh.npoin
                b = user_source(RHS[ip],
                                params.qp.qn[ip],
                                params.qp.qe[ip],
                                sem.mesh.npoin, inputs[:CL], inputs[:SOL_VARS_TYPE];
                                neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip])
                RHS[ip] = b
            end

            for ip = 1:sem.mesh.npoin
                sem.matrix.L[ip,ip] += inputs[:rconst][1]
            end
            
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
                                                 RHS, 0.0, params.ubdy,
                                                 params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                                 params.mesh.bdy_edge_type,
                                                 params.ω, qp.neqs, params.inputs, params.AD, sem.mesh.SD)
      
         
        else
            k = lin_solve_rhs_gpu_2d!(inputs[:backend])
            k(RHS, qp.qn, qp.qe, sem.mesh.x, sem.mesh.y, qp.neqs; ndrange = sem.mesh.npoin)
            KernelAbstractions.synchronize(inputs[:backend])
            
            Minv = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin), Int64(sem.mesh.npoin))
            k =  diagm_gpu!(inputs[:backend])
            k(Minv , sem.matrix.Minv; ndrange = sem.mesh.npoin)
            L_temp = Minv * sem.matrix.L
            sem.matrix.L .= L_temp

            k = apply_boundary_conditions_gpu_lin_solve!(inputs[:backend])
            k(RHS, sem.matrix.L, sem.mesh.poin_in_bdy_edge, sem.mesh.npoin; ndrange = (sem.mesh.nedges_bdy*sem.mesh.ngl), workgroupsize = (sem.mesh.ngl))
            KernelAbstractions.synchronize(inputs[:backend])
            if ("Laguerre" in sem.mesh.bdy_edge_type)
                k = apply_boundary_conditions_lag_gpu_lin_solve!(inputs[:backend])
                k(RHS, sem.matrix.L, sem.mesh.connijk_lag, sem.mesh.ngl, sem.mesh.ngr, sem.mesh.npoin, sem.mesh.nelem_semi_inf, inputs[:lperiodic_laguerre];
                  ndrange = (sem.mesh.nelem_semi_inf*sem.mesh.ngl,sem.mesh.ngr), workgroupsize = (sem.mesh.ngl,sem.mesh.ngr))
                KernelAbstractions.synchronize(inputs[:backend])
            end
            k = add_to_diag!(inputs[:backend])
            k(sem.matrix.L, TFloat(10.0); ndrange = sem.mesh.npoin)
            KernelAbstractions.synchronize(inputs[:backend])
        end
        
        @time solution = solveAx(sem.matrix.L, RHS, inputs[:ode_solver])

        write_output(sem.mesh.SD, solution,  sem.mesh,
                     OUTPUT_DIR, inputs,
                     params.qp.qvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe, case="none")
        
    end
end

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

    amr_freq = inputs[:amr_freq]
    Δt_amr   = amr_freq * inputs[:Δt]
    tspan    = [inputs[:tinit], inputs[:tinit] + Δt_amr]
    params, u =  params_setup(sem,
                              qp,
                              inputs,
                              OUTPUT_DIR,
                              TFloat,
                              tspan)
    
    if !inputs[:llinsolve]
        
        solution = time_loop!(inputs, params, u)
        
        if (inputs[:ndiagnostics_outputs] > 0)
            write_output(sem.mesh.SD, solution,  sem.mesh,
                         OUTPUT_DIR, inputs,
                         params.qp.qvars,
                         inputs[:outformat];
                         nvar=params.qp.neqs, qexact=params.qp.qe, case="rtb")
        end
        
    else
        
        RHS = KernelAbstractions.zeros(inputs[:backend], TFloat,Int64(sem.mesh.npoin), qp.neqs)
        if (inputs[:backend] == CPU())
            for ip =1:sem.mesh.npoin
                rhs = user_source(RHS[ip],
                                  params.qp.qn[ip],
                                  params.qp.qe[ip],          #ρref
                                  sem.mesh.npoin, inputs[:CL], inputs[:SOL_VARS_TYPE];
                                  neqs=1, x=sem.mesh.x[ip],y=sem.mesh.y[ip])
                RHS[ip] = rhs
            end
       
            Minv = diagm(sem.matrix.Minv)
            L_temp = Minv * sem.matrix.L
            sem.matrix.L .= L_temp 
            apply_boundary_conditions_lin_solve!(sem.matrix.L,RHS,sem.mesh,inputs,sem.mesh.SD)             
            for ip = 1:sem.mesh.npoin
                sem.matrix.L[ip,ip] += 10
            end
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
        
        @time solution = solveAx(-sem.matrix.L, RHS, inputs[:ode_solver]) 

        write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat]; nvar=params.qp.neqs)
        
    end

    
end

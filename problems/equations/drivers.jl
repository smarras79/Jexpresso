using HDF5

function driver(inputs::Dict,        #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    sem = sem_setup(inputs)
    
    if (inputs[:backend] != CPU())
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
    end

    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    params, u =  params_setup(sem,
                              qp,
                              inputs,
                              OUTPUT_DIR,
                              TFloat)
    if !(inputs[:llinsolve])   
        
        solution = time_loop!(inputs, params, u)
        
        if (inputs[:ndiagnostics_outputs] > 0)
            write_output(sem.mesh.SD, solution,  sem.mesh,
                         OUTPUT_DIR, inputs,
                         params.qp.qvars,
                         inputs[:outformat];
                         nvar=params.qp.neqs, qexact=params.qp.qe, case="rtb")
        end
        
    else
        
        RHS = zeros(Float64,sem.mesh.npoin)
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
        
        
        @time solution = solveAx(-sem.matrix.L, RHS, inputs[:ode_solver]) 

        write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat]; nvar=params.qp.neqs)
        
    end

    
end

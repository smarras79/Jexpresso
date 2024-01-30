using HDF5

function driver(inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    sem = sem_setup(inputs)
   
    if !(inputs[:llinsolve])   
  
        qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)
    
        check_length(qp.qn[1,:], qp.neqs+1, "drivers --> initialize.jl")
        
    solution = time_loop!(sem.QT,
                          sem.PT,
                          inputs[:SOL_VARS_TYPE],
                          inputs[:CL],
                          inputs[:AD],
                          sem.mesh, sem.metrics, sem.basis, sem.ω, qp,
                          sem.matrix.M, sem.matrix.Minv,
                          inputs[:Δt],
                          inputs,
                          OUTPUT_DIR,
                          TFloat;fx=sem.fx,fy=sem.fy,fy_lag=sem.fy_lag)

        if (inputs[:ndiagnostics_outputs] > 0)
            write_output(sem.mesh.SD, solution,  sem.mesh, OUTPUT_DIR, inputs, qp.qvars, inputs[:outformat]; nvar=qp.neqs, qexact=qp.qe, case="rtb")
        end
     
    else

       qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)
       
       RHS = zeros(Float64,sem.mesh.npoin)
       for ip =1:sem.mesh.npoin
           rhs = user_source(RHS[ip],
                        qp.qn[ip],
                        qp.qe[ip],          #ρref
                        sem.mesh.npoin, inputs[:CL], inputs[:SOL_VARS_TYPE]; neqs=1, x=sem.mesh.x[ip],y=sem.mesh.y[ip])
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

       write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat]; nvar=qp.neqs)
   
    end

    
end

using HDF5

function driver(inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    sem = sem_setup(inputs)
    
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
    
    #HDF5
    write_output(sem.mesh.SD, solution,  sem.mesh, OUTPUT_DIR, inputs, qp.qvars, HDF5(); nvar=qp.neqs, qexact=qp.qe, case="rtb")
    
    if (inputs[:ndiagnostics_outputs] > 0)
        write_output(sem.mesh.SD, solution,  sem.mesh, OUTPUT_DIR, inputs, qp.qvars, inputs[:outformat]; nvar=qp.neqs, qexact=qp.qe, case="rtb")
    end
    
end

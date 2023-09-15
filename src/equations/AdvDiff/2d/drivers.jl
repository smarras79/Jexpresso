include("./initialize.jl")

function driver(DT::ContGal,       #Space discretization type
                inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 
    
    sem = sem_setup(inputs)
    
    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)
        
    Δt = inputs[:Δt]
    
    solution = time_loop!(sem.QT, sem.PT, sem.mesh, sem.metrics, sem.basis, sem.ω, qp, sem.matrix.M, sem.matrix.De, sem.matrix.Le, Δt, inputs, OUTPUT_DIR, TFloat)
    
    write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat]; nvar=qp.neqs)
    #solution_norms(solution, OUTPUT_DIR, inputs;)
    
end

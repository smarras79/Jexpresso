include("./initialize.jl")

function driver(DT::ContGal,       #Space discretization type
                inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    
    sem = sem_setup(inputs)
    
    #--------------------------------------------------------
    # Initialize q
    #--------------------------------------------------------
    qp = define_q(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, TFloat; neqs=1)

    #Build ∫S(q)dΩ
    RHS = build_rhs_source(sem.mesh.SD, sem.QT, qp.qn, sem.mesh, sem.matrix.M, TFloat)

    #BC
    apply_boundary_conditions!(sem.mesh.SD, zeros(sem.mesh.ngl, sem.mesh.ngl, sem.mesh.nelem), qp.qn, sem.mesh, inputs, sem.QT, sem.metrics, sem.basis.ψ, sem.basis.dψ, sem.ω, 0.0, qp.neqs; sem.matrix.L)
    
    println(" # Solve Lq=RHS ................................")    
    solution = solveAx(sem.matrix.L, RHS, inputs[:ode_solver])
    println(" # Solve Lq=RHS ................................ DONE")
    
    write_output(solution, sem.mesh.SD, sem.mesh, OUTPUT_DIR, inputs, inputs[:outformat]; nvar=qp.neqs)
    
end

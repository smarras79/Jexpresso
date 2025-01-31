function solveAx(L, RHS, linear_solver...)
    
    @info linear_solver
    
    prob = LinearProblem(L, RHS);    
    sol = solve(prob, linear_solver)
    
    return sol
end

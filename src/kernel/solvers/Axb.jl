using LinearSolve
"""
    solveAx(L, RHS, linear_solver... preconditionr)

This function is a wrapper to the LinearSolve.jl library to solve the linear system

    ```math
        Lq = RHS
    ```

"""
function solveAx(L, RHS, linear_solver...)
    
    prob = LinearProblem(L, RHS);    
    @time solution = solve(prob, linear_solver)

    return solution
end

using LinearSolve
using LinearSolve: solve
using SnoopCompile

"""
    solveAx(L, RHS, linear_solver... preconditionr)

This function is a wrapper to the LinearSolve.jl library to solve the linear system

    ```math
        Lq = RHS
    ```

"""
function solveAx(L, RHS, linear_solver...)

    @info linear_solver
    
    prob = LinearProblem(L, RHS);    
    @time sol = solve(prob, linear_solver)

    return sol
end

using LinearSolve
using LinearSolve: solve
using SnoopCompile

function solveAx(L, RHS, linear_solver...)

    @info " Solve Ax = b ................."

    prob = LinearProblem(L, RHS);    
    sol = solve(prob, linear_solver)
    
    @info " Solve Ax = b ................. END"
    return sol
end

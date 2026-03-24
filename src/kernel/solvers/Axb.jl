function solveAx(L, RHS, linear_solver...)
    
    prob = LinearProblem(L, RHS);    
    sol = solve(prob, linear_solver)
    
    return sol
end


#function solveAx_sparse(A, b)
#     return A\b
#end

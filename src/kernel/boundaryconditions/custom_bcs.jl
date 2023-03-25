include("../abstractTypes.jl")

#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
problem_name = ARGS[1]
user_bc_dir = string("../../problems/", problem_name, "/user_bc.jl")
include(user_bc_dir)
#---------------------------------------------------------------------------

function neumann(q,gradq,x,y,t,mesh,metrics,tag)
    
    rhs = user_bc_neumann(q,gradq,x,y,t,tag)
    return rhs
end

function dirichlet!(q,gradq,x,y,t,mesh,metrics,tag)

    q = user_bc_dirichlet!(q,gradq,x,y,t,tag)
    
    return q
end


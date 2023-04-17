include("../abstractTypes.jl")

#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #problem_name
    user_bc_dir = string("../../problems/", ARGS[1], "/user_bc.jl")
elseif (length(ARGS) === 2) #problem_name/problem_case_name
    user_bc_dir = string("../../problems/", ARGS[1], "/", ARGS[2], "/user_bc.jl")
end
include(user_bc_dir)
#---------------------------------------------------------------------------

function neumann(q,gradq,x,y,t, mesh, metrics,tag)
    
    rhs = user_bc_neumann(q,gradq,x,y,t,tag)
    return rhs
end

function neumann(q,gradq,x,t, mesh, metrics)

    rhs = user_bc_neumann(q,gradq,x,t)
    return rhs
end

function dirichlet!(q,gradq,x,y,t, mesh, metrics,tag)

    q = user_bc_dirichlet!(q,gradq,x,y,t,tag)
    
    return q
end

function dirichlet!(q,gradq,x,t, mesh, metrics)

    q = user_bc_dirichlet!(q,gradq,x,t)

    return q
end


include("../abstractTypes.jl")

#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
if (length(ARGS) === 1) #equations
    user_bc_dir = string("../../equations/", ARGS[1], "/user_bc.jl")
elseif (length(ARGS) === 2) #equations/equations_case_name
    user_bc_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_bc.jl")
end
include(user_bc_dir)
#---------------------------------------------------------------------------

function neumann(q, gradq, x, y, t, mesh, metrics, tag)
    
    rhs = user_bc_neumann(q, gradq, x, y, t, tag)
    return rhs
end

function neumann(q, gradq, x, t, mesh, metrics)

    rhs = user_bc_neumann(q, gradq, x, t)
    return rhs
end

function dirichlet!(q, x, y, t, mesh, metrics, tag, qbdy)

    qbdy = user_bc_dirichlet!(q, x, y, t, tag, qbdy)
    
    return qbdy
end

function dirichlet!(q, x, t, mesh, metrics, qbdy)

    qbdy = user_bc_dirichlet!(q, x, t)

    return qbdy
end

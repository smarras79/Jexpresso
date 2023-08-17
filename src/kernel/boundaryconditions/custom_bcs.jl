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

function neumann(q, gradq, x, y, t, mesh, metrics, tag, inputs::Dict)
    
    rhs = user_bc_neumann(q, gradq, x, y, t, tag, inputs::Dict)
    return rhs
end

function neumann(q, gradq, x, t, mesh, metrics, inputs::Dict)

    rhs = user_bc_neumann(q, gradq, x, t, inputs)
    return rhs
end

function dirichlet!(q, gradq, x, y, t, mesh, nx, ny, tag, qbdy, inputs::Dict)

    qbdy = user_bc_dirichlet!(q, gradq, x, y, t, tag, qbdy, inputs::Dict, nx, ny)
    
    return qbdy
end

function dirichlet!(q, gradq, x, t, mesh, metrics, tag, qbdy, inputs::Dict)
    
    qbdy = user_bc_dirichlet!(q, gradq, x, t, tag, qbdy, inputs::Dict)

    return qbdy
end

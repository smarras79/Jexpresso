#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
user_bc_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_bc.jl")    
include(user_bc_dir)
#---------------------------------------------------------------------------
function dirichlet!(qbdy, x, y, t, tag, inputs::Dict)
    user_bc_dirichlet!(qbdy, x, y, t, tag, inputs::Dict)
end

function neumann(q, gradq, x, y, t, tag, inputs::Dict)
    
    rhs = user_bc_neumann(q, gradq, x, y, t, tag, inputs::Dict)
    return rhs
end

function neumann(q, gradq, x, t, inputs::Dict)

    rhs = user_bc_neumann(q, gradq, x, t, inputs)
    return rhs
end

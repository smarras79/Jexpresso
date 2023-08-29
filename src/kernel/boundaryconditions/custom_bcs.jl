#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
user_bc_dir = "../../equations/CompEuler/theta/user_bc.jl"
if (length(ARGS) === 1) #equations
    user_bc_dir = string("../../equations/", ARGS[1], "/user_bc.jl")
elseif (length(ARGS) === 2) #equations/equations_case_name
    user_bc_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_bc.jl")
end
include(user_bc_dir)
#---------------------------------------------------------------------------
function dirichlet!(qbdy, x, y, t, tag, inputs::Dict)
    user_bc_dirichlet!(qbdy, x, y, t, tag, inputs::Dict)
end


function dirichlet!(qbdy, gradq, x, t, tag, inputs::Dict)
    
    user_bc_dirichlet!(qbdy, gradq, x, t, tag, inputs::Dict)

    return qbdy
end

function neumann(q, gradq, x, y, t, tag, inputs::Dict)
    
    rhs = user_bc_neumann(q, gradq, x, y, t, tag, inputs::Dict)
    return rhs
end

function neumann(q, gradq, x, t, inputs::Dict)

    rhs = user_bc_neumann(q, gradq, x, t, inputs)
    return rhs
end

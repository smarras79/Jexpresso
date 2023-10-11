#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
#if (size(ARGS) === 2)
#    user_bc_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_bc.jl")    
#    include(user_bc_dir)
#end
#---------------------------------------------------------------------------

function neumann(q, gradq, x, y, t, tag, inputs::Dict)
    
    rhs = user_bc_neumann(q, gradq, x, y, t, tag, inputs::Dict)
    return rhs
end

function neumann(q, gradq, x, t, inputs::Dict)

    rhs = user_bc_neumann(q, gradq, x, t, inputs)
    return rhs
end

function dirichlet!(q, qbdy, x, y, t, nx, ny, tag,qe,SOL)

   user_bc_dirichlet!(q, x, y, t, tag, qbdy, nx, ny,qe,SOL)
end

function dirichlet!(q, gradq, x, t, mesh, metrics, tag, qbdy, inputs::Dict)
    
    qbdy = user_bc_dirichlet!(q, gradq, x, t, tag, qbdy, inputs::Dict)

    return qbdy
end

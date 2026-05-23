#---------------------------------------------------------------------------
# Fetch problem name to access the user_bc functions
#---------------------------------------------------------------------------
#if (size(ARGS) === 2)
#    user_bc_dir = string("../../equations/", ARGS[1], "/", ARGS[2], "/user_bc.jl")
#    include(user_bc_dir)
#end
#---------------------------------------------------------------------------

function neumann(q::AbstractArray, gradq::AbstractArray, coords::AbstractArray,
                 t::AbstractFloat, tag::AbstractString, inputs)
    rhs = user_bc_neumann(q, gradq, coords, t, tag, inputs)
    return rhs
end

function neumann(q::AbstractArray, gradq::AbstractArray, coords::AbstractArray,
                 t::AbstractFloat, inputs)
    rhs = user_bc_neumann(q, gradq, coords, t, inputs)
    return rhs
end

function dirichlet!(q::AbstractArray, qbdy::AbstractArray, coords::AbstractArray,
                    t::AbstractFloat, nx::AbstractFloat, ny::AbstractFloat,
                    tag::AbstractString, qe::AbstractArray, SOL)
    user_bc_dirichlet!(q, coords, t, tag, qbdy, nx, ny, qe, SOL)
end

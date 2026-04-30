"""
    qibdy is an Array{Floats} of size `nvars`

    src/equations/EQUATIONS_NAME/user_bc.jl contains a set of user-defined boundary conditions functions
    that can be modified as needed.

    The function defined in src/equations/EQUATIONS_NAME/user_bc.jl
    are called by the b.c. functions defined in src/kernel/custom_bcs.jl
    within a boundary-edge loop that detects the "tag" string defined in the user-generated *.msh file.

    For the periodic 1D viscous Burgers case no Dirichlet/Neumann data are needed,
    so the callbacks are empty.
"""
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    nothing
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end


function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs::Dict)

    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    return T(q[1])
end

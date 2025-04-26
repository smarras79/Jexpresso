"""
    qibdy is an Array{Floats} of size `nvars`

    src/equations/EQUATIONS_NAME/user_bc.jl contains a set of user-defined boundary conditions functions
    that can be modified as needed.

    The function defined in src/equations/EQUATIONS_NAME/user_bc.jl 
    are called by the b.c. functions defined in src/kernel/custom_bcs.jl
    within a boundary-edge loop that detects the "tag" string defined in the user-generated *.msh file.

    For example:
    If some domain boundaries of gmsh file mymesh.msh are tagged as "inflow" and "no_slip", then the user
    creating the functions in user_bc.jl must define the behavior of the unknown or its derivatives
    on those boundaries.

    ```math
    if (tag === "inflow")
        qibdy[1] = 3.0
    elseif (tag === "fix_temperature")
        qibdy[2] = 300.0
    end
    return qibdy
    ```
    where  `qibdy[i=1:nvar]` is the value unknown `i`
    
"""
function user_bc_dirichlet!(q, x::AbstractFloat, t::AbstractFloat, tag::String,qbdy::AbstractArray,qe,::TOTAL)

   nothing
end

function user_bc_dirichlet!(q, x::AbstractFloat, t::AbstractFloat, tag::String,qbdy::AbstractArray,qe,::PERT)

    if (tag == "left")
      qbdy[1] = 0.0
    end
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,x,t,lpert)
    T = eltype(q)
    if (x == T(0.0))
        return T(0.0)
    else
        return T(q[1])
    end
end

"""
    qibdy is an Array{Floats} of size `nvars`

    src/problem/PROBLEM_NAME/user_bc.jl contains a set of user-defined boundary conditions functions
    that can be modified as needed.

    The function defined in src/problem/PROBLEM_NAME/user_bc.jl 
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
function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String)
    q[1] = 1.0
    q[2] = 2.5
    q[3] = 0.0 
    return q
end

function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)

    if (x < 0.0)
        q[1] = 1.0
        q[2] = 2.5
        q[3] = 0.0
    else
        nothing
    end
    return q
end

#=function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    q[1] = 0.1
    q[2] = 0.0
    return q
end=#

#=function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    if (x > 1.0)
        q[1] = 2.0
    else
        q[2] = 4.42
    end
    return q
end=#

#=function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    if (x > 1.0)
        if (q[2]/sqrt(9.81*q[1])< 1)
            q[1] = 0.66
        end
    else
        q[2] = 1.53
    end
    return q
end=#

#=function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    if (x > 1.0)
        q[1] = 0.33
    else
        q[2] = 0.18
    end
    return q
end=#

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    flux = zeros(size(q,2),1)
    return flux
end

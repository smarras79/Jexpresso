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

    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ
    
    case = "Sod"
    if case === "Sod"
        
        ρL, uL, vL, pL = 1.000, 0.0, 0.0, 1.0
        ρR, uR, vR, pR = 0.125, 0.0, 0.0, 0.1
        
        xshock_initial = 0.5
        if (x < xshock_initial)
            ρ = ρL
            u = uL
            v = vL
            p = pL
        else
            ρ = ρR
            u = uR
            v = vR
            p = pR
        end
    elseif case === "sound"
        ρ = 1.0
        u = 0.0
        v = 0.0
        p = 1.0
    else
        ρ = 1.0
        u = 0.0
        v = 0.0
        p = 1.0
    end
    ρE = p/(γ - 1.0) + 0.5*ρ*(u*u + v*v)
    q[1] = ρ
    q[2] = ρ*u
    q[3] = ρ*v
    q[4] = ρE
    
    return q
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat)
    flux = zeros(size(q,2),1)
    return flux
end

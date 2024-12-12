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
function user_bc_dirichlet!(q::SubArray{TFloat}, x::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe::SubArray{TFloat},::TOTAL)
    
    PhysConst = PhysicalConst{Float64}()
    γ = 1.4

    case = "sod" # inputs[:case]
    if case == "sod"
        ρL, uL, pL = 1.000, 0.0, 1.0
        ρR, uR, pR = 0.125, 0.0, 0.1
        xshock_initial = 0.5
        if (x < xshock_initial)
            ρ = ρL
            u = uL
            p = pL
        else
            ρ = ρR
            u = uR
            p = pR
        end
                
    elseif case == "sound"
        ρ = 1.0
        u = 0.0
        p = 1.0

    end
    E = p/(γ - 1.0) + 0.5*ρ*u*u
    qbdy[1] = ρ
    qbdy[2] = ρ*u
    qbdy[3] = E
    
    return qbdy
end


function user_bc_dirichlet_gpu(q,qe,x,t,lpert)
    T = eltype(q)
    return T(0.0),T(q[2])
end

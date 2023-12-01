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
function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, npoin_linear, t::AbstractFloat, tag::String,qbdy::AbstractArray,qe::SubArray{Float64},::Any; ip=1)
    
    massflow = 0.579
    ρin = 1.0
    xin = 0.0
    Ain = 1.0 + 2.2*(xin - 1.5)^2
    γ   = 1.4
    Tin = 1.0
    pin = 1.0
    if (tag == "left")
        
        U1 = ρin
        U2 = 2*q[2,2] - q[3,2]
        U3 = pin/(γ - 1.0) + 0.5*U2^2

        qbdy[1] = U1
        qbdy[2] = U2
        qbdy[3] = U3
        
    elseif (tag == "right")
        P = 0.6784

        #Supersonic ouflow:
        #   qbdy[1] = 2*q[2,1] - q[3,1]
     #   qbdy[2] = 2*q[2,2] - q[3,2]
     #   qbdy[3] = 2*q[2,3] - q[3,3]
    end
end

function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, tag::String)

    xl = 0.0
    Al = 1 + 2.2*(xl - 1.5).^2
    qbdy[1] = Al
    qbdy[2] = 2*q[2] - q
    
    return qbdy
end
function user_bc_dirichlet!(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String)

    qbdy[1] = 0.0    
    
    return qbdy
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    
    flux = zeros(size(q,2),1)
    return flux
end

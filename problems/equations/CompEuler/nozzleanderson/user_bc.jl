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
function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, t::AbstractFloat, tag::String,qbdy::AbstractArray,qe::SubArray{Float64},::Any; npoin_linear=1)

    U1in = 0.0
    U2in = 0.0
    U3in = 0.0
    
    U1out = 0.0
    U2out = 0.0
    U3out = 0.0
    
    ip2 = 2 #this is the 2nd point of the linear grid
    ip3 = 3 #this is the 3rd point of the linear grid
    ipN = npoin_linear #last geometric point of the 1D mesh. The H-O node count starts from this one in the first element.

    xin = 0.0
    Ain = 1.0 + 2.2*(xin - 1.5)^2
    Tin = 1.0
    ρin = 1.0
    pin = ρin*Tin
    γ = 1.4
    γm1 = 0.4
    
    if (tag == "left")
        U1in = Ain

        #U2in = 2*q[ip2,2] - q[ip3,2]
        U2in = 2*q[ipN+1,2] - q[ipN+2,2]
        uin  = U2in/U1in
        U3in = U1in*(Tin/γm1 + 0.5*γ*uin*uin)
        
        qbdy[1] = U1in
        #qbdy[2] = U2in
        qbdy[3] = U3in
    end

    if (tag == "right")
        U1out = 2*q[ipN-1,1] - q[ipN-2,1]
        U2out = 2*q[ipN-1,2] - q[ipN-2,2]
        U2out = 2*q[ipN-1,3] - q[ipN-2,3]
        
    #    qbdy[1] = U1out
    #    qbdy[2] = U2out
    #    qbdy[3] = U3out
    end
    
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)
#    if (tag == "free_slip")
      
      if (y<=14950) #(abs(x) < 119500.0 && y<= 14950.0)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
        qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
        qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
      end
      #else 
       # qbdy[2] = 0.0
        #qbdy[3] = 0.0
      #end
      #if (abs(x) > 119500.0 && y < 0.1)
      #  qbdy[2] = 0.0
      #  qbdy[3] = 0.0
      #end
     #@info x,y,nx,ny,qbdy[2],qbdy[3] 
  # return qbdy #, flags
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

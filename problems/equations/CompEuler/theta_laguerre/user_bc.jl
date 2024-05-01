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
function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe::SubArray{Float64},::TOTAL)
    if (tag == "free_slip")
     #if ((x == 5000.0 && y == 0.0) || (x == -5000.0 && y == 0.0) || (x == -5000.0 && y == 10000.0) || (x == 5000.0 && y == 10000.0)) 
     #  a = 1
     #  b = 1
     #  if (x > 0)
     #    a= -1
     #  end
     #  if (y > 0)
     #    b = -1
     #  end
     #  qnl = a*(sqrt(2)/2)*q[2] + b*(sqrt(2)/2)*q[3]
      
    # else
       qnl = nx*q[2] + ny*q[3]
    # end
     qbdy[2] = q[2] - qnl*nx
     qbdy[3] = q[3] - qnl*ny
   end
 #   else
  #    qbdy[2] = 0.0
   # end
   #return qbdy #, flags
    
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)
     #if (y <4995 || abs(x) > 4995)
     if (tag == "free_slip")
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
        qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
        qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
     end

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet(q,x,y,t,nx,ny,qbdy)
    qnl = nx*(q[2]) + ny*(q[3])
    u = q[2] - qnl*nx
    v = q[3] - qnl*ny

    #return Float32(q[1]), Float32(0.0), Float32(0.0), Float32(q[4])
    # return Float32(qbdy[1]), Float32(u), Float32(v), Float32(qbdy[4])
    return Float32(qbdy[1]), Float32(u), Float32(v), Float32(qbdy[4])
    #return Float32(qbdy[1]), Float32(0.0), Float32(0.0), Float32(qbdy[4])
end
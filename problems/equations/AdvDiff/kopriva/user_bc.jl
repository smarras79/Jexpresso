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
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe,::TOTAL)
#    if (tag == "free_slip")
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
    if (y < 0.01)
        qbdy[1] = 0.0#0.025*sinpi(2*30*t/5000.0)
    end
    if (y < 10.0 && abs(x) >9.9)
        qbdy[1] = 0.0
    end  
    # else
 #   else
  #    qbdy[2] = 0.0
   # end
   #return qbdy #, flags
    
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe,::PERT)
#    if (tag == "free_slip")
    if (y < 0.01)
        qbdy[1] = 0.0#0.025*sinpi(2*30*t/5000.0)
    end
    if (y < 10.0 && abs(x) >9.9)
        qbdy[1] = 0.0
    end
#else 
       # qbdy[2] = 0.0
        #qbdy[3] = 0.0
      #end
      #if (abs(x) > 119500.0 && y < 0.1)
      #  qbdy[2] = 0.0
      #  qbdy[3] = 0.0
      #end
     #@info coords,nx,ny,qbdy[2],qbdy[3] 
  # return qbdy #, flags
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,coords,t,nx,ny,qbdy,lpert)
    T = eltype(q)
    if (y < 0.01)
        return T(0.0)#0.025*sinpi(2*30*t/5000.0)
    end
    if (y < 10.0 && abs(x) >9.9)
        return(0.0)
    end
    return T(qbdy[1])
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe,::TOTAL)

#=    if (tag == "free_slip")
        qbdy[2] = 0.0
        qbdy[3] = 0.0
    end

     if (tag == "circle")
        qbdy[2] = 0.0
        qbdy[3] = 0.0

#        @info x y
    end=#
        
    if (tag == "free_slip" || tag == "circle")
        qnl     = nx*q[2] + ny*q[3]
        qbdy[2] = q[2] - qnl*nx
        qbdy[3] = q[3] - qnl*ny
    end
            
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe,::PERT)
#    if (tag == "free_slip")
    
    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    
      #else 
       # qbdy[2] = 0.0
        #qbdy[3] = 0.0
      #end
      #if (abs(x) > 119500.0 && y < 0.1)
      #  qbdy[2] = 0.0
      #  qbdy[3] = 0.0
      #end
     #@info coords,nx,ny,qbdy[2],qbdy[3] 
    
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
    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
        u = (q[2]+qe[2] - qnl*nx) - qe[2]
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
    else
        qnl = nx*(q[2]) + ny*(q[3])
        u = q[2] - qnl*nx
        v = q[3] - qnl*ny
    end
 
    return T(qbdy[1]), T(u), T(v), T(qbdy[4])
end

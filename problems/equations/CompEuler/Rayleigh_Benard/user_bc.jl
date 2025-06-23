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
    
    # else
    qnl = nx*q[2] + ny*q[3]
    # end
    qbdy[2] = 0.0#q[2] - qnl*nx
    qbdy[3] = 0.0#q[3] - qnl*ny
    if (tag == "bottom")
        qbdy[4] = q[1]*301.0
    elseif (tag == "top")
        qbdy[4] = q[1]*299.0
    end
    #if (tag == "lateral")
    #    qbdy[4] = qe[4]
    #end
    #   else
    #    qbdy[2] = 0.0
    # end
    #return qbdy #, flags
    
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe,::PERT)
#    if (tag == "free_slip")
    
    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
    qbdy[2] = 0.0#(q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = 0.0#(q[3]+qe[3] - qnl*ny) - qe[3]
    if (tag == "lateral")
        qdby[4] = 0.0
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

function user_bc_neumann!(F_edge, u, u1, qe, qe1, tag, coords, ::TOTAL)

    if (tag == "bottom")
        F_edge[4] = 0.02*rand()*u[1]
    elseif (tag == "top")
        F_edge[4] = -0.02*rand()*u[1]
    end

end

function user_bc_neumann!(F_edge, u, u1, qe, qe1, tag, coords, τ_f, wθ, ::PERT)

    if (tag == "bottom")
        F_edge[4] = 0.02*rand()*(u[1]+qe[1])
    elseif (tag == "top")
        F_edge[4] = -0.02*rand()*(u[1]+qe[1])
    end
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

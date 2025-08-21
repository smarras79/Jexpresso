function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe::SubArray{Float64},::TOTAL)
    
    #= if(x <= -1.0 || x >= 1.0)
        qbdy[1] = 0.0
    end

    if(y <= -1.0 || y >= 1.0)
        qbdy[1] = 0.0
    end =#
    qbdy[1] = 0.0
    
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)
      
    #= if(x <= -1.0 || x >= 1.0)
        qbdy[1] = 0.0
    end

    if(y <= -1.0 || y >= 1.0)
        qbdy[1] = 0.0
    end =#
    qbdy[1] = 0.0
    
end

# function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
#     flux = zeros(size(q,2),1)
#     return flux
# end

# function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
#     flux = zeros(size(q,2),1)
#     return flux
# end

#=function user_bc_dirichlet_gpu(q,qe,x,y,t,nx,ny,qbdy,lpert)
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
end=#

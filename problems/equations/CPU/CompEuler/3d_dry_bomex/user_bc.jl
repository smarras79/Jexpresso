function user_bc_dirichlet!(q::SubArray{Float64},
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe::SubArray{Float64}, ::TOTAL)

    qnl = nx*q[2] + ny*q[3] + nz*q[4]
    qbdy[2] = (q[2] - qnl*nx) 
    qbdy[3] = (q[3] - qnl*ny) 
    qbdy[4] = (q[4] - qnl*nz) 
    if (z < 0.01) 
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = 0.0
        # @info qe[5], qe[1], qe[5]/qe[1]
        qbdy[5] = 2.0*qe[1] + qe[5]
    end
    
end

function user_bc_dirichlet!(q::SubArray{Float64},
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe::SubArray{Float64}, ::PERT)

    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,x,y,z,t,nx,ny,nz,qbdy,lpert)
    T = eltype(q)
    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        #=if (abs(nx) > T(0.001))
            u = (q[2]+qe[2] - qnl*nx) - qe[2]
        else
            u = qbdy[2]
        end
        if (abs(ny) > T(0.001))
            v = (q[3]+qe[3] - qnl*ny) - qe[3]
        else
            v = qbdy[3]
        end
        if (abs(nz) > T(0.001))
            w = (q[4]+qe[4] - qnl*nz) - qe[4]
        else
            w = qbdy[4]
        end=#
        u = (q[2]+qe[2] - qnl*nx) - qe[2]
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
        w = (q[4]+qe[4] - qnl*nz) - qe[4]
    else
        qnl = nx*(q[2]) + ny*(q[3]) + nz*(q[4])
        #=if (abs(nx) > T(0.001))
            u = q[2] - qnl*nx
        else
            u = qbdy[2]
        end
        if (abs(ny) > T(0.001))
            v = q[3] - qnl*ny
        else
            v = qbdy[3]
        end
        if (abs(nz) > T(0.001))
            w = q[4] - qnl*nz
        else
            w = qbdy[4]
        end=#
        
        u = (q[2] - qnl*nx)
        v = (q[3] - qnl*ny)
        w = (q[4] - qnl*nz)
        if (z < 0.01)
            u = 0.0
            v = 0.0
            w = 0.0
            qbdy[5] = 303.0*q[1]
        end

    end
    return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5])
end

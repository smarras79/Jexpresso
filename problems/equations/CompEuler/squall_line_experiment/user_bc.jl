function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)
    #=if ((z < 0.1 || z > zmax - 10) && ( x < xmin + 10 || x > xmax - 10) ) || (abs(nx) >0.1)
        qbdy[1] = qe[1]
        qbdy[2] = qe[2]
        qbdy[3] = qe[3]
        qbdy[4] = qe[4]
        qbdy[5] = qe[5]
        qbdy[6] = qe[6]
        qbdy[7] = qe[7]
    else=#
        qnl = nx*q[2] + ny*q[3] + nz*q[4]
        qbdy[2] = (q[2] - qnl*nx) #+ qe[2] 
        qbdy[3] = (q[3] - qnl*ny) #+ qe[3]
        qbdy[4] = (q[4] - qnl*nz) #+ qe[4]
    #end
end

function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::PERT)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
        qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
        qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]
        if ((z < 0.1 || z > 24000.0 - 10) && ( x < -60000 + 10 || x > 60000 - 10) ) || (abs(nx) >0.1)
        qbdy[1] = 0.0
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = 0.0
        qbdy[5] = 0.0
        qbdy[6] = 0.0
        qbdy[7] = 0.0
    end
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,x,y,z,t,nx,ny,nz,qbdy,lpert)
    T = eltype(q)
    if ((z < T(0.1) || z > T(24000) - T(10)) && ( x < T(-40000) + T(10) || x > T(40000) - T(10)) ) || (abs(nx) >T(0.1))
        return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
    end
    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        u = (q[2]+qe[2] - qnl*nx) - qe[2]
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
        w = (q[4]+qe[4] - qnl*nz) - qe[4]
        return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5]), T(qbdy[6]), T(qbdy[7])
    else
        qnl = nx*(q[2]) + ny*(q[3]) + nz*(q[4])
        u = (q[2] - qnl*nx)
        v = (q[3] - qnl*ny)
        w = (q[4] - qnl*nz)
        return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5]), T(qbdy[6]), T(qbdy[7])
    end
end

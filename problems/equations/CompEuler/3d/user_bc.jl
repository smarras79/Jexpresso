function user_bc_dirichlet!(q::SubArray{Float64},
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe::SubArray{Float64}, ::TOTAL)

    e = 1.0e-4
    
    xmine = xmin + e
    xmaxe = xmax - e
    ymine = ymin + e
    ymaxe = ymax - e
    zmine = zmin + e
    zmaxe = zmax - e

    if x < xmine || x > xmaxe
        qbdy[2] = 0.0
    end
    
    if y < ymine || y > ymaxe 
        qbdy[3] = 0.0
    end

    if z < zmine || z > zmaxe 
        qbdy[4] = 0.0
    end
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet(q,qe,x,y,z,t,nx,ny,nz,qbdy)

    qnl = nx*(q[2]) + ny*(q[3]) + nz*(q[4])
    if (abs(nx) > Float32(0.001))
        u = q[2] - qnl*nx
    else
        u = qbdy[2]
    end
    if (abs(ny) > Float32(0.001))
        v = q[3] - qnl*ny
    else
        v = qbdy[3]
    end
    if (abs(nz) > Float32(0.001))
        w = q[4] - qnl*nz
    else
        w = qbdy[4]
    end
    return Float32(qbdy[1]), Float32(u), Float32(v), Float32(w), Float32(qbdy[5])
end

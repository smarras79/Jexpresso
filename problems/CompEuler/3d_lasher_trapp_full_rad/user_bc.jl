function user_bc_dirichlet!(q,
                            coords,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)

    qnl = nx*q[2] + ny*q[3] + nz*q[4]
    qbdy[2] = (q[2] - qnl*nx) 
    qbdy[3] = (q[3] - qnl*ny) 
    qbdy[4] = (q[4] - qnl*nz) 
    
end

function user_bc_dirichlet!(q,
                            coords,
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
end

function user_bc_neumann!(F_surf, q, q1, qe, qe1, tag, coords, CL, PhysConst; θ=0, θ1=0, qn0=0, qn1=0)
   
    if (tag == "bottom")
        
        x = coords[1]
        y = coords[2]
        # Gaussian surface heat flux (Lasher-Trapp et al. 2001)
        # Center at domain midpoint (~3 km for a 6x6 km domain); adjust xc/yc/σ to match mesh.
        xc    = 3000.0
        yc    = 3000.0
        σ     = 1000.0   # m; keeps flux < 1% of Q_max at lateral boundaries (3 km away)
        Q_max = 300.0    # W/m² (liquid-ice static energy, no unit conversion needed)
        F_surf[5] = Q_max * exp(-((x - xc)^2 + (y - yc)^2) / (2 * σ^2))
        
    end
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, nz, qbdy, lpert)
    T = eltype(q)
    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
        w = (q[4]+qe[4] - qnl*nz) - qe[4]
    else
        qnl = nx*(q[2]) + ny*(q[3]) + nz*(q[4])
        
        u = (q[2] - qnl*nx)
        v = (q[3] - qnl*ny)
        w = (q[4] - qnl*nz)
        if (coords[3] < 0.01)
            u = 0.0
            v = 0.0
            w = 0.0
        end

    end
    return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5]), T(qbdy[6]), T(qbdy[7])
end

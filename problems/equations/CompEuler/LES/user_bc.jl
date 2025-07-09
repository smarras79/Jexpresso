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
    #if (z < 0.01) 
    #    qbdy[2] = 0.0
    #    qbdy[3] = 0.0
    #    qbdy[4] = 0.0
    #    # e_tot
    #    qbdy[5] = 59169*q[1]
    #    qbdy[6] = 0.0175*q[1]
    #end
    
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

function user_bc_neumann!(F_edge, u, u1, qe, qe1, tag, coords, τ_f, wθ, CL)

    if (tag == "bottom")
        F_edge[4] = 0.02*rand()*(u[1]+qe[1])
    elseif (tag == "top")
        F_edge[4] = -0.02*rand()*(u[1]+qe[1])
    end
end

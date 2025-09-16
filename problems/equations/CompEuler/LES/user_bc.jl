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


function user_bc_neumann!(F_surf, u, u1, qe, qe1, tag, coords, τ_f, wθ, eqtype)

    #if (tag == "wall_model_bottom" || tag == "wall_model_top" || tag == "MOST")
   # if (tag == "MOST")
   #     # Use the pre-computed wall shear stress components
   #     # Apply with correct sign for Neumann BC
   #     F_surf[2] = τ_f[1]  # x-momentum equation
   #     F_surf[3] = τ_f[2]  # y-momentum equation
   #     F_surf[5] = wθ[1]   # θ equation
   #    # @info F_surf[5]
   # end
    if (tag == "bottom")
        F_edge[4] = 0.0
    elseif (tag == "top")
        F_edge[4] = 0.0
    end
end



 

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


function user_bc_neumann!(F_edge, u, u1, qe, qe1, tag, coords, τ_f, wθ, eqtype)
    
    if (tag == "wall_model")
        # Get the wall shear stress magnitude (assuming it's passed via τ_f or accessible)
        #τw_mag = τ_f  # or however you access the stored τ_w value
        τw_mag = sqrt(τ_f[1]*τ_f[1] + τ_f[2]*τ_f[2])
        
        # Get velocity components at the boundary point
        u_vel = u1[2]  # u-velocity
        v_vel = u1[4]  # w-velocity (assuming 3D)
        
        # Compute velocity magnitude in wall-parallel directions
        vel_mag = sqrt(u_vel^2 + v_vel^2)
        
        # Avoid division by zero
        if vel_mag > 1e-12
            # Compute wall shear stress components
            τw_x = -τw_mag * u_vel / vel_mag  # x-component
            τw_y = -τw_mag * v_vel / vel_mag  # z-component
            
            # Apply Neumann BC: F_edge = viscous flux = μ * ∂u/∂n
            # Since τw = μ * ∂u/∂n at the wall
            F_edge[2] = τw_x  # x-momentum equation
            F_edge[3] = τw_y  # z-momentum equation (assuming 3D)
        else
            F_edge[2] = 0.0
            F_edge[3] = 0.0
        end
    end
end

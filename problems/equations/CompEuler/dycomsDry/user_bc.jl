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


function user_bc_neumann!(F_surf, q, q1, qe, qe1, tag, coords, τ_f, wθ, eqtype)

    PhysConst = PhysicalConst{Float32}()

    cp  = PhysConst.cp
    H   = 0.12 #Km/s

    #if (tag == "wall_model_bottom" || tag == "wall_model_top" || tag == "MOST")
    if (tag == "MOST")
        # Use the pre-computed wall shear stress components
        # Apply with correct sign for Neumann BC
        F_surf[2] = -0.2 #τ_f[1]  # x-momentum equation
        F_surf[3] = -0.2 #τ_f[2]  # y-momentum equation
        F_surf[5] = 0.12 #120.0/(q[1]*cp) #150 #ρ*cp*wθ[1]   # θ equation
       # @info F_surf[5]
    end
end

function user_bc_dirichlet!(q,
                            coords,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)

    # CALLED BY DEFAULT.
    #
    # If you don't want Dirichlet to do anything, keep this function empty.
    #
    qnl = nx*q[2] + ny*q[3] + nz*q[4]
    qbdy[2] = (q[2] - qnl*nx) 
    qbdy[3] = (q[3] - qnl*ny)
    qbdy[4] = (q[4] - qnl*nz)

    #=if tag == "top_wall"
        u_geo = 7.0
        α     = 3.0 #deg
        qbdy[2] = u_geo*cospi(α/180.0)
        qbdy[3] = u_geo*sinpi(α/180.0)
    end=#
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


function user_bc_neumann!(F_surf, q, q1, qe, qe1, tag, coords, τ_f, wθ, CL, PhysConst; kwargs...)

    # NOT CALLED BY DEFAULT UNLESS lbdy_flux => true in input
    
    #if (tag == "wall_model_bottom" || tag == "wall_model_top" || tag == "MOST")
    if (tag == "MOST")
        ρ = q[1]
        # Use the pre-computed wall shear stress components
        # Apply with correct sign for Neumann BC
        #if (τ_f[1] != 0.0 || τ_f[2] != 0.0) @info τ_f[1] τ_f[2] end
        
        F_surf[2] = τ_f[1]  # x-momentum equation
        F_surf[3] = τ_f[2]  # y-momentum equation
        F_surf[5] = 0.12   # wθ[1]/(ρ*PhysConst.cp)   # θ equation
   end

    #
    # Example of user defined bulk formulas (in contrast to, e.g., CM_MOST, which is built in its own function called from within rhs.jl
    #
    #=   if (tag == "bottom")
        ρ    = q[1]  #+ qe[1]
        ρ1   = q1[1] #+ qe1[1]
        ρu   = q[2]  #+ qe[2]
        ρu1  = q1[2] #+ qe1[2]
        ρv   = q[3]  #+ qe[3]
        ρv1  = q1[3] #+ qe1[3]

        ρθ   = q[5]
        ρθ1  = q1[5]
        u   = ρu/ρ
        v   = ρv/ρ
        u1  = ρu1/ρ1
        v1  = ρv1/ρ1
        θ = ρθ/ρ
        θ1 = ρθ1/ρ

        u_12  = (u + u1)/2
        v_12  = (v + v1)/2
        θ_12  = (θ + θ1)/2


        cd = 1.1e-3 + 4e-5*sqrt(u_12^2+v_12^2)
        ce = cd

        F_surf[2] = -ρ*cd*u_12*sqrt(u_12^2+v_12^2)
        F_surf[3] = -ρ*cd*v_12*sqrt(u_12^2+v_12^2)
        F_surf[5] = PhysConst.cp*ρ*ce*sqrt(u_12^2+v_12^2)*(θ-θ_12)
    end
    =#

end



 

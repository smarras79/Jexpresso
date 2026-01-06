function user_bc_dirichlet!(q,
                            coords,
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
                            coords,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::PERT)
        PhysConst = PhysicalConst{Float64}()
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
            qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
            qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
            qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]
        #if ((z < 0.1 || z > zmax - 10) && ( x < xmin + 10 || x > xmax - 10) ) || (abs(nx) >0.1)
        #qbdy[4] = 0.0
        #=if (tag == "bottom")
            qbdy[5] = 0.0
            qbdy[6] = 0.0
        end=#
        #=if ((z < 0.1 || z > zmax - 10) && ( x < xmin + 10 || x > xmax - 10) ) || (abs(nx) >0.1)
        qbdy[1] = 0.0
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = 0.0
        qbdy[5] = 0.0
        qbdy[6] = 0.0
        qbdy[7] = 0.0
    end=#
end

function user_bc_neumann!(F_surf, q, q1, qe, qe1, tag, coords, τ_f, wθ, CL, PhysConst; θ=0, θ1 = 0, qn0=0, qn1=0)

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
        ρ    = q[1]  + qe[1]
        ρ1   = q1[1] + qe1[1]
        ρu   = q[2]  + qe[2]
        ρu1  = q1[2] + qe1[2]
        ρv   = q[3]  + qe[3]
        ρv1  = q1[3] + qe1[3]
        ρqt  = q[6]  + qe[6]
        ρqt1 = q1[6] + qe1[6]

        u   = ρu/ρ
        v   = ρv/ρ
        qt  = ρqt/ρ
        u1  = ρu1/ρ1
        v1  = ρv1/ρ1
        qt1 = ρqt1/ρ1

        qv  = qt  - qn0
        qv1 = qt1 - qn1

        u_12  = (u + u1)/2
        v_12  = (v + v1)/2
        θ_12  = (θ + θ1)/2

        qv_12 = (qv + qv1)/2

        cd = 1.1e-3 + 4e-5*sqrt(u_12^2+v_12^2)
        ce = cd

        F_surf[2] = -ρ*cd*u_12*sqrt(u_12^2+v_12^2)
        F_surf[3] = -ρ*cd*v_12*sqrt(u_12^2+v_12^2)
        F_surf[5] = PhysConst.cp*ρ*ce*sqrt(u_12^2+v_12^2)*(θ-θ_12)
        F_surf[6] = ρ*ce*sqrt(u_12^2+v_12^2)*(qv-qv_12)
    end
end

function user_bc_dirichlet_gpu(q,qe,coords,t,nx,ny,nz,qbdy,lpert)
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

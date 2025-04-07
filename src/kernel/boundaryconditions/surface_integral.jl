function compute_surface_integral!(S_face, F_surf, ω, Jac_face, iface, ngl)
    for i = 1:ngl
        for j = 1:ngl

            ωJac = ω[i]*ω[j]*Jac_face[iface,i,j]
            #@info  ωJac*F_surf[i,j,:], ω[i], ω[j], Jac_face[iface,i,j], F_surf[i,j,:]
            S_face[iface,i,j,:] .+= ωJac*F_surf[i,j,:]
        end
    end

end

function DSS_surface_integral!(S_flux, S_face, M_surf_inv, nfaces, ngl, z, zmin, connijk, poin_in_bdy_face, bdy_face_in_elem)

    for iface = 1:nfaces
        if (z[poin_in_bdy_face[iface,3,3]] == zmin)
            for i = 1:ngl
                for j = 1:ngl
                    e = bdy_face_in_elem[iface]
                    ip = poin_in_bdy_face[iface,i,j]
                    ip1 = connijk[e,i,j,2]
                #=if (S_face[iface,i,j,2] > 0.0 || S_face[iface,i,j,2] < 0.0)
                    @info ip, S_face[iface,i,j,:], M_surf_inv[ip]
                end=#
                    S_flux[ip1,:] .+= S_face[iface,i,j,:]*M_surf_inv[ip]
                end
            end
        end
    end

end

function build_surface_mass_matrix(nfaces, npoin, ω, ψ, ngl, Jac_face, poin_in_bdy_face, T, Δx, inputs)

    Mface = zeros(T, ngl^2, ngl^2, nfaces)
    build_mass_matrix!(Mface, NSD_2D(), Inexact(), ψ, ω, nfaces, Jac_face, Δx, ngl-1, ngl-1, T)
    M_surface = zeros(T, npoin) 
    DSS_mass!(M_surface, NSD_2D(), Inexact(), Mface, poin_in_bdy_face, nfaces, npoin, ngl-1, T; llump=inputs[:llump])
    for ip=1:npoin
        if (abs(M_surface[ip]) < 1e-3)
            M_surface[ip] = 1.0
        end
    end
    return M_surface
end

function bulk_surface_flux!(F_surf,q,q1,qe,qe1,θ,θ1,qn,qn1)

    PhysConst = PhysicalConst{Float64}()    
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

    qv  = qt  - qn
    qv1 = qt1 - qn1
    
    u_12  = (u + u1)/2
    v_12  = (v + v1)/2
    θ_12  = (θ + θ1)/2
    qv
    qv_12 = (qv + qv1)/2

    cd = 1.1e-3 + 4e-5*sqrt(u_12^2+v_12^2)
    ce = cd
    
    
    F_surf[2] = ρ*cd*u_12*sqrt(u_12^2+v_12^2)
    F_surf[3] = ρ*cd*v_12*sqrt(u_12^2+v_12^2)
    F_surf[5] = PhysConst.cp*ρ*ce*sqrt(u_12^2+v_12^2)*(θ-θ_12)
    F_surf[6] = ρ*ce*sqrt(u_12^2+v_12^2)*(qv-qv_12)
end


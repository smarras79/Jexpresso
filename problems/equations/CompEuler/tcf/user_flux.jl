function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh, param_set,
                    ::CL, ::TOTAL; neqs=4, ip=1,
                    x=0.0, y=0.0, z=0.0)

    # PhysConst = PhysicalConst{Float64}()
    T      = eltype(q)
    
    ρ      = q[1]
    ρu     = q[2]
    ρv     = q[3]
    ρw     = q[4]
    ρe_tot = q[5]
    
    u     = ρu/ρ
    v     = ρv/ρ
    w     = ρw/ρ
    e_tot = ρe_tot/ρ

    _grav = T(TP.grav(param_set))
    e_pot = _grav * z
    e_kin = T(1 // 2) * (u^2 + v^2 + w^2)
    e_int = e_tot - e_pot - e_kin
   
    Temp = TD.air_temperature(param_set, e_int, 0.0)
    Pressure = TD.air_pressure(param_set,Temp, ρ, 0.0)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρe_tot*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρe_tot*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρe_tot*w
    
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)
    
    @mystop("woring in progress: case user_flux! PERT not implemented for this case!")
end

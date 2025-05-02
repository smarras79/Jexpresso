function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    param_set = create_updated_TD_Parameters(Float64(101325.0))
    T  = eltype(q)
    z  = mesh.z[ip]
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρe_tot = q[5]
    ρqt = q[6]
    ρql = q[7]
    
    e_tot  = ρe_tot/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ
    qt  = ρqt/ρ
    ql  = ρql/ρ

    _grav = T(TP.grav(param_set))
    e_pot = _grav * z
    e_kin = T(1 // 2) * (u^2 + v^2 + w^2)
    e_int = e_tot - e_pot - e_kin
    q_pt = TD.PhasePartition(qt,ql)
    Temp = TD.air_temperature(param_set,e_int,q_pt)
    Pressure = TD.air_pressure(param_set,Temp,ρ,q_pt)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρe_tot*u
    F[6] = ρqt*u
    F[7] = ρql*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρe_tot*v
    G[6] = ρqt*v
    G[7] = ρql*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρe_tot*w
    H[6] = ρqt*w
    H[7] = ρql*w
    
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)

   
    @mystop("woring in progress: case BOMEX user_flux! PERT not implemented!")
end

function user_flux_gpu(q,qe,z,PhysConst, param_set,lpert)
    T = eltype(q)
    if !(lpert)
        ρ  = q[1]
        ρu = q[2]
        ρv = q[3]
        ρw = q[4]
        ρθ = q[5]
        ρqt = q[6]
        ρql = q[7]
        e_tot  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        w  = ρw/ρ
        qt = ρqt/ρ
        ql = ρql/ρ
        _grav = T(TP.grav(param_set))
        e_pot = _grav * z
        e_kin = T(1 // 2) * (u^2 + v^2 + w^2)
        e_int = e_tot - e_pot - e_kin
        q_pt = TD.PhasePartition(qt,ql)
        Temp = TD.air_temperature(param_set,e_int,q_pt)
        Pressure = TD.air_pressure(param_set,Temp,ρ,q_pt)
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρw*u), T(ρθ*u + Pressure*u), T(ρqt*u), T(ρql*u), T(ρv), T(ρu*v), T(ρv*v + Pressure), T(ρw*v), T(ρθ*v + Pressure*v), T(ρqt*v), T(ρql*v), T(ρw), T(ρu*w), T(ρv*w), T(ρw*w + Pressure), T(ρθ*w + Pressure*w), T(ρqt*w), T(ρql*w)
    end
end

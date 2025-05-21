function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh, 
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    T      = eltype(q)
    
    ρ      = q[1]
    ρu     = q[2]
    ρv     = q[3]
    ρw     = q[4]
    ρe_tot = q[5]
    tr     = q[6]
    
    u     = ρu/ρ
    v     = ρv/ρ
    w     = ρw/ρ
    e_tot = ρe_tot/ρ
    
    e_pot = 0.0 #PhysConst.g * z
    e_kin = T(1 // 2) * (u^2 + v^2 + w^2)
    e_int = e_tot - e_pot - e_kin
    
    Temp = e_int/PhysConst.cv
    Pressure = perfectGasLaw_ρTtoP(PhysConst, ρ=ρ, Temp=Temp)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρe_tot*u
    F[6] = tr*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρe_tot*v
    G[6] = tr*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρe_tot*w
    H[6] = tr*w
    
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)
    
    @mystop("woring in progress: case user_flux! PERT not implemented for this case!")
end

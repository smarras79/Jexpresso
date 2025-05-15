function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρhl = q[4]
    ρqt = q[5]
    ρqp = q[6]

    hl  = ρhl/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    qt = ρqt/ρ
    qp = ρqp/ρ

    Pressure = q[end] - qe[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt = qt, qp = qp)
    
    F[1] = ρu
    F[2] = ρu*u + Pressure
    F[3] = ρu*v
    F[4] = ρhl*u
    F[5] = ρqt*u
    F[6] = ρqp*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v + Pressure
    G[4] = ρhl*v
    G[5] = ρqt*v
    G[6] = ρqp*v
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2] + qe[2]
    ρv = q[3] + qe[3]
    ρhl = q[4] + qe[4]
    ρqt = q[5] + qe[5]
    ρqp = q[6] + qe[6]

    hl  = ρhl/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    ρqt = ρqt/ρ
    ρqp = ρqp/ρ

    Pressure = q[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt=qt, qp=qp)
    Pressure = Pressure - qe[end]
    
    F[1] = ρu
    F[2] = ρu*u + Pressure
    F[3] = ρv*u
    F[4] = ρhl*u
    F[5] = ρqt*u
    F[6] = ρqp*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Pressure
    G[4] = ρhl*v
    G[5] = ρqt*v
    G[6] = ρqp*v
    
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    if (lpert)
        ρ  = q[1]+qe[1]
        ρu = q[2]+qe[2]
        ρv = q[3]+qe[3]
        ρw = q[4]+qe[4]
        ρhl = q[5]+qe[5]
        ρqt = q[6]+qe[6]
        ρqp = q[7]+qe[7]
        
        hl  = ρhl/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        w  = ρw/ρ
        qt = ρqt/ρ
        qp = ρqp/ρ

        Pressure = q[end]-qe[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt=qt, qp=qp) - qe[6]
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρw*u), T(ρhl*u), T(ρqt*u), T(ρqp*u), T(ρv), T(ρu*v), T(ρv*v + Pressure), T(ρw*v), T(ρhl*v), T(ρqt*v), T(ρqp*v), T(ρw), T(ρu*w), T(ρv*w), T(ρw*w + Pressure), T(ρhl*w), T(ρqt*w), T(ρqp*w)
    else
        ρ  = q[1]
        ρu = q[2]
        ρv = q[3]
        ρw = q[4]
        ρhl = q[5]
        ρqt = q[6]
        ρqp = q[7]
  
        hl  = ρhl/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        w  = ρw/ρ
        qt = ρqt/ρ
        qp = ρqp/ρ

        Pressure = q[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt=qt, qp=qp)
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρw*u), T(ρhl*u), T(ρqt*u), T(ρqp*u), T(ρv), T(ρu*v), T(ρv*v + Pressure), T(ρw*v), T(ρhl*v), T(ρqt*v), T(ρqp*v), T(ρw), T(ρu*w), T(ρv*w), T(ρw*w + Pressure), T(ρhl*w), T(ρqt*w), T(ρqp*w)
    end
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
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
    Pressure = q[end] - qe[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt = qt, qp = qp)
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρhl*u
    F[6] = ρqt*u
    F[7] = ρqp*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρhl*v
    G[6] = ρqt*v
    G[7] = ρqp*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρhl*w
    H[6] = ρqt*w
    H[7] = ρqp*w
    
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2] + qe[2]
    ρv = q[3] + qe[3]
    ρw = q[4] + qe[4]
    ρhl = q[5] + qe[5]
    ρqt = q[6] + qe[6]
    ρqp = q[7] + qe[7]

    hl  = ρhl/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ
    ρqt = ρqt/ρ
    ρqp = ρqp/ρ

    Pressure = q[end]#perfectGasLaw_ρhltoP(PhysConst, ρ=ρ, hl=hl, qt=qt, qp=qp)
    Pressure = Pressure - qe[end]

    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρv*u
    F[4] = ρw*u
    F[5] = ρhl*u
    F[6] = ρqt*u
    F[7] = ρqp*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure
    G[4] = ρw*v
    G[5] = ρhl*v
    G[6] = ρqt*v
    G[7] = ρqp*v
    
    H[1] = ρw
    H[2] = ρu*w
    H[3] = ρv*w
    H[4] = ρw*w .+ Pressure
    H[5] = ρhl*w
    H[6] = ρqt*w
    H[7] = ρqp*w
    
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

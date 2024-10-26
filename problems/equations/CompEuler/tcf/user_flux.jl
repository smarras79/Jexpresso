function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, H::SubArray{Float64},
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρe = q[5]
    
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ
    
    T = (ρe/ρ - 0.5*sqrt(u*u + v*v + w*w))/PhysConst.cv

    gm1 = PhysConst.γ - 1.0
    #Pressure = gm1*(ρe - 0.5*ρ*(u*u + v*v + w*w))
    Pressure = perfectGasLaw_ρTtoP(PhysConst, ρ=ρ, Temp=T)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = u*(ρe .+ Pressure)

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = v*(ρe .+ Pressure)

    H[1] = ρw
    H[2] = ρw*u
    H[2] = ρw*v
    H[3] = ρw*w .+ Pressure
    H[4] = w*(ρe .+ Pressure)
end

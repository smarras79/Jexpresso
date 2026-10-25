function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh, 
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    T  = eltype(q)
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρθ = q[5]
    
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ
    θ = ρθ/ρ
    
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρθ*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρθ*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρθ*w
    
end

function user_flux!(F, G, H,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)
    
    @mystop("woring in progress: case user_flux! PERT not implemented for this case!")
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    pref::Float64,
                    mesh::St_mesh, ::CL; neqs=4, ip=1)
    
    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    dPress = Press - pref
    
    F[1] = ρu
    F[2] = ρu*u + dPress
    F[3] = ρv*u
    F[4] = ρθ*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + dPress
    G[4] = ρθ*v
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    pref::Float64,
                    mesh::St_mesh, ::NCL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
                
    ρ = q[1]
    u = q[2]
    v = q[3]
    θ = q[4]
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρ*u
    F[2] = u
    F[3] = u
    F[4] = θ
    
    G[1] = ρ*v
    G[2] = v
    G[3] = v
    G[4] = θ
end

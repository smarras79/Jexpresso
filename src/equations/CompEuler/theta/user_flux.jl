function user_flux!(F, G, T, SD::NSD_2D, q::Array, mesh::St_mesh; neqs=4)

    PhysConst = PhysicalConst{Float64}()
    
    #F = zeros(neqs)
    #G = zeros(neqs)
            
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
end


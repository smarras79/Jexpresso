function user_flux(T, SD::NSD_2D, q::Array, mesh::St_mesh; neqs=4)

    PhysConst = PhysicalConst{Float64}()
    
    F = zeros(T, mesh.npoin, neqs)
    G = zeros(T, mesh.npoin, neqs)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρv = q[ip,3]
        ρθ = q[ip,4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        
        Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
        #@info " FLUX USER: " Temp Press ρ
        F[ip,1] = ρu
        F[ip,2] = ρu*u + Press
        F[ip,3] = ρv*u
        F[ip,4] = ρθ*u
        
        G[ip,1] = ρv
        G[ip,2] = ρu*v
        G[ip,3] = ρv*v + Press
        G[ip,4] = ρθ*v
    end
    
    return F, G
end

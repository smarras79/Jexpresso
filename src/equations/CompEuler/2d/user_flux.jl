function user_flux(T, SD::NSD_2D, q::Array, mesh::St_mesh; neqs=4)

    PhysConst = PhysicalConst{Float64}()
    
    F = zeros(T, mesh.npoin, neqs)
    G = zeros(T, mesh.npoin, neqs)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρv = q[ip,3]
        ρE = q[ip,4]
        E  = ρE/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        
        Temp = (E - 0.5*(u*u + v*v))/PhysConst.cv
        Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)
        
        #@info " FLUX USER: " Temp Press ρ
        F[ip,1] = ρu
        F[ip,2] = ρu*u + Press
        F[ip,3] = ρv*u
        F[ip,4] = ρE*u + Press*u
        
        G[ip,1] = ρv
        G[ip,2] = ρu*v
        G[ip,3] = ρv*v + Press
        G[ip,4] = ρE*v + Press*v
    end
    
    return F, G
end

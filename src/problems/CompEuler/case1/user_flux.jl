function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh; neqs=3)

    PhysConst = PhysicalConst{Float64}()
    
    F = zeros(T, mesh.npoin, neqs)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        ρ  = q[ip,1]
        ρu = q[ip,2]
        E  = q[ip,3]
        u  = ρu/ρ
        
        Temp = (E - 0.5*u*u)/PhysConst.cv
        #Press = perfectGasLaw(PhysConst; ρ=ρ, Temp=Temp)

        R = 287.0
        Press = 0.4*(E - 0.5*u*u) #ρ*R*Temp
        
        #@info " FLUX USER: " Temp Press ρ
        F[ip,1] = ρu
        F[ip,2] = ρu*u + Press
        F[ip,3] =  E*u + Press*u
    end
    
    return F
end

function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh; neqs=3)

    PhysConst = PhysicalConst{Float64}()
    
    F = zeros(T, neqs)
    
    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]
    u  = ρu/ρ
        
    Temp = (ρE/ρ - 0.5*u*u)/PhysConst.cv
    Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)
        
    #@info " FLUX USER: " Temp Press ρ
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρE*u + Press*u
    
    return F
end

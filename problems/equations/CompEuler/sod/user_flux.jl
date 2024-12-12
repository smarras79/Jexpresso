function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    
    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]
    E  = ρE/ρ
    u  = ρu/ρ
    
    PhysConst = PhysicalConst{Float64}()
    Temp = (E - 0.5*u*u)/PhysConst.cv
    Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)
    
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρE*u + Press*u
        
end


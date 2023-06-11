function perfectGasLaw_ρTtoP(PhysConst::PhysicalConst; ρ=1.25, Temp=300.0)
    
    return ρ*PhysConst.Rair*Temp #Press
end

function perfectGasLaw_ρPtoT(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    
    return  Press/(ρ*PhysConst.Rair) #Temp
end


function perfectGasLaw_TPtoρ(PhysConst::PhysicalConst; Temp=300.0, Press=100000.0)
    
    return  Press/(ρ*PhysConst.Rair) #ρ
end

function perfectGasLaw_ρθtoP!(Press, PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    
    Press = PhysConst.C0*(ρ*θ)^PhysConst.γ #Press

end


function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    
    return PhysConst.C0*(ρ*θ)^PhysConst.γ #Press

end

function perfectGasLaw_ρPtoθ(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    
    return (1.0/ρ)*(Press/PhysConst.C0)^(1.0/PhysConst.γ) #θ
    
end

function perfectGasLaw_θPtoρ(PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    
    return (1.0/θ)*(Press/PhysConst.C0)^(1.0/PhysConst.γ) #ρ
    
end

function primitive2flux!(F, G, q)

    PhysConst = PhysicalConst{Float64}()

    #p1 = q[1]
    #p2 = ($q[2])./($q[1])
    #p3 = ($q[3])./($q[1])
    #p4 = ($q[4])./($q[1])
    
    prim = SVector{4, Float64}(p1, p2, p3, p4)
    
    Press = 0.0
    perfectGasLaw_ρθtoP!(Press, PhysConst, ρ=prim[1], θ=prim[4])

    F = [q[2] #.*q[2]
         q[2]
         q[2]
         q[2]]
    
    G =  [q[3] #.*q[3]
          q[3]
          q[3]
          q[3]]
    
end

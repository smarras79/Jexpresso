function perfectGasLaw_ρTtoP(PhysConst::PhysicalConst; ρ=1.25, Temp=300.0)
    
    return ρ*PhysConst.Rair*Temp #Press
end

function perfectGasLaw_ρPtoT(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    
    return  Press/(ρ*PhysConst.Rair) #Temp
end


function perfectGasLaw_TPtoρ(PhysConst::PhysicalConst; Temp=300.0, Press=100000.0)
    
    return  Press/(ρ*PhysConst.Rair) #ρ
end

function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    
    return PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
    
end

function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst, ρ::AbstractArray, θ::AbstractArray)
    return PhysConst.C0 .* (ρ .* θ) .^ PhysConst.γ
end

function perfectGasLaw_ρθtoP!(Press::Array{Float64}, PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    
    Press[1] = PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
    
end

function perfectGasLaw_ρPtoθ(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    T = typeof(ρ)    
    return (T(1.0)/ρ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #θ
    
end

function perfectGasLaw_θPtoρ(PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    T = typeof(θ)
    return (T(1.0)/θ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #ρ
    
end

function moistPressure(PhysConst::PhysicalConst; ρ=1.25, Temp=300.0, qv = 0.0)
    T = typeof(Temp)
    #@info ρ*Temp*PhysConst.Rair + ρ*Temp*qv*PhysConst.Rvap, Temp, qv, ρ, PhysConst.Rvap
    return (T(ρ*Temp*PhysConst.Rair + ρ*Temp*qv*PhysConst.Rvap))
end

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

function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst, ρ::AbstractArray, θ::AbstractArray)
    return PhysConst.C0 .* (ρ .* θ) .^ PhysConst.γ
end

function perfectGasLaw_ρθtoP!(Press::Array{Float64}, PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    
    Press[1] = PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
    
end

function perfectGasLaw_ρθqtoP(ρ, θ, qt, ql)
    L_v = 2.5008e6
    kappad = 0.28571428571
    cpv = 1859.0
    cpl = 4181.0
    p0 = 101325.0
    R_d = 287.0024
    molmass_ratio = 1.608079364

    cpd = R_d / kappad
    qd = 1.0 - qt
    qv = qt - ql
    Rm = R_d * (1 + (molmass_ratio - 1) * qt - molmass_ratio * ql)
    cpm = cpd * qd + cpv * qv + cpl * ql
    press = (ρ * R_d * θ)^(cpm / (cpm - Rm)) * (p0)^(-Rm / (cpm - Rm))

    return press
end


function perfectGasLaw_θqPtoT(θ, qt, ql, p)
    L_v = 2.5008e6
    kappad = 0.28571428571
    cpv = 1859.0
    cpl = 4181.0
    p0 = 101325.0
    R_d = 287.0024
    molmass_ratio = 1.608079364

    cpd = R_d / kappad
    qd = 1.0 - qt
    qv = qt - ql
    Rm = R_d * (1 + (molmass_ratio - 1) * qt - molmass_ratio * (qt - qv))
    cpm = cpd * qd + cpv * qv + cpl * ql
    theta = R_d / Rm * θ
    temperature = theta * (p / p0) ^ (Rm / cpm)

    return temperature
end

function perfectGasLaw_ρPtoθ(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    T = typeof(ρ)    
    return (T(1.0)/ρ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #θ
    
end

function perfectGasLaw_θPtoρ(PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    T = typeof(θ)
    return (T(1.0)/θ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #ρ
    
end

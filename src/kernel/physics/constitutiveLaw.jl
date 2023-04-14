function perfectGasLaw(PhysConst::PhysicalConst; ρ=1.25, Temp=300.0, Press=100000.0)
    
    Press = ρ*PhysConst.Rair*Temp

    return Press
end

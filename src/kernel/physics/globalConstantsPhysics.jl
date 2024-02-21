using Parameters

@with_kw struct PhysicalConst{T}

    #Thermodynamic constants at T=300 K
    Rair::T =  287.0 #J/kg.K
    cp::T   = 1005.0 #J/kg.K
    cv::T   =  718.0 #J/kg.K
    γair::T = T(cp/cv)
    γ::T    = T(cp/cv)
    Pr::T   = 0.7
    Prnum::T= 0.1
    pref::T = 100000.0 #Pa
    
    Rovercv::T = T(Rair/cv)
    cpoverR::T = T(cp/Rair)
    cpovercv::T = T(cp/cv)
    cvovercp::T = T(cv/cp)
    C0::T   = T((Rair^γ)/pref^(γ-1.0)) #Rovercv
    
    #Gravity
    g::T = 9.80616 #m/s²
    g2::T= 9.80616*9.80616
    
    #Elasticity
    E::T = 70.0e9                  #Pa
    ν::T = 0.33                    #Poisson's ratio: -dϵ_transverse/dϵ_axial
    λ::T = (E*ν)/((1+ν)*(1-2*ν)) #Lamé parameters λ, μ
    μ::T = E/(2*(1+ν))
    
end

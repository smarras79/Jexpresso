using Parameters

@with_kw struct PhysicalConst{T}

    #Thermodynamic constants at T=300 K
    Rair::T   =  287.0 #J/kg.K
    cp::T     = 1005.0 #J/kg.K
    cv::T     =  718.0 #J/kg.K
    γair::T   = cp/cv
    γ::T      = cp/cv
    Pr::T     = 0.7
    Prnum::T  = 0.1
    pref::T   = 100000.0 #Pa
    ρwater::T = 1000.0   #kg/m3 reference liquid water density
        
    Rovercv = Rair/cv
    cpoverR = cp/Rair
    cpovercv= cp/cv
    cvovercp= cv/cp
    C0::T   = (Rair^γ)/pref^(γ-1.0) #Rovercv
    
    #Gravity
    g::T = 9.80616 #m/s²
    g2::T= 9.80616*9.80616
    
    #Elasticity
    E::T = 70.0e9                  #Pa
    ν::T = 0.33                    #Poisson's ratio: -dϵ_transverse/dϵ_axial
    λ::T = (E*ν)/((1+ν)*(1-2*ν)) #Lamé parameters λ, μ
    μ::T = E/(2*(1+ν))
    
end

using Parameters

@with_kw struct MicrophysicalConst{T}

    #Thermodynamic constants at T=300 K
    Rair::T

    # Constants related to saturation water pressure
    xlv::T    = 2500000.0
    ep2::T    = 0.6217504
    svp1::T   = 0.6112000
    svp2::T   = 17.67000
    svp3::T   = 29.65000
    svpt0::T  = 273.1500
    ρwater::T = 1000.0   #kg/m3 reference liquid water density
    
end


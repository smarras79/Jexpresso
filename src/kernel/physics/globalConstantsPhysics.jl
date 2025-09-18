using Parameters

@with_kw struct PhysicalConst{T}

    #Thermodynamic constants at T=300 K
    Rair::T =  287.0 #J/kg.K
    cp::T   = 1004.0 #J/kg.K
    cv::T   =  718.0 #J/kg.K
    γair::T = T(cp/cv)
    γ::T    = T(cp/cv)
    Pr::T   = 0.7
    Prnum::T= 0.1
    pref::T = 100000.0 #Pa
    Rvap::T = 461.0 #J/kg.K
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

    ## molar masses
    Mol_mass_air   = 28.9647      #g/mol
    Mol_mass_water = 18.02        #g/mol

    # Reference pressure used in potential temperature definition... mainly for BMOEX case, very sensitive to pressure
    potential_temperature_reference_pressure::T = 101325.0 #Pa
end

using Parameters

@with_kw struct MicrophysicalConst{T}
    
    # Constants related to saturation water pressure
    xlv::T        = 2500000.0
    ep2::T        = 0.6217504
    svp1::T       = 0.6112000
    svp2::T       = 17.67000
    svp3::T       = 29.65000
    svpt0::T      = 273.1500
    ρwater::T     = 1000.0   #kg/m3 reference liquid water density
    a_rain::T     = 842      #m^(1-b)s^(-1) Constant in fall speed formula for rain
    a_snow::T     = 4.84     #m^(1-b)s^(-1) Constant in fall speed formula for snow
    a_graupel::T  = 94.5     #m^(1-b)s^(-1) Constant in fall speed formula for graupel
    a_fr::T       = 0.78     #Constant in ventilation factor for rain
    a_fs::T       = 0.65     #Constant in ventilation factor for snow
    a_fg::T       = 0.78     #Constant in ventilation factor for graupel
    b_rain::T     = 0.8      #Exponent in fall speed formula for rain
    b_snow::T     = 0.25     #Exponent in fall speed formula for snow
    b_graupel::T  = 0.5      #Exponent in fall speed formula for graupel
    b_fr::T       = 0.31     #Constant in ventilation factor for rain
    b_fs::T       = 0.44     #Constant in ventilation factor for snow
    b_fg::T       = 0.31     #Constant in ventilation factor for graupel
    C_rain::T     = 1.0      #Rain shape factor
    C_snow::T     = 2/π      #Snow shape factor
    C_graupel::T  = 1.0      #Graupel shape factor
    Da::T         = 2.21e-5  #m^2/s Diffusion coefficient of water vapor at 0C
    Er_c::T       = 1.0      #Collection efficiency of rain for cloud water
    Es_c::T       = 1.0      #Collection efficiency of snow for cloud water
    Eg_c::T       = 1.0      #Collection efficiency of graupel for cloud water
    Er_i::T       = 1.0      #Collection efficiency of rain for cloud ice
    Es_i::T       = 0.1      #Collection efficiency of snow for cloud ice
    Eg_i::T       = 0.1      #Collection efficiency of graupel for cloud ice
    Lc::T         = 2.5104e6 #J/kg Latent heat of vaporization
    Ls::T         = 2.8440e6 #J/kg Latent heat of sublimation
    Lf::T         = 0.3336e6 #J/kg Latent heat of fusion
    N0_rain::T    = 8e6      #m^(-4) Intercept parameter for rain
    N0_snow::T    = 3e6      #m^(-4) Intercept parameter for snow
    N0_graupel::T = 4e6      #m^(-4) Intercept parameter for graupel
    Ka::T         = 2.4e-2   #J m K^(-1)/s Thermal conductivity of air at 0C
    qc0::T        = 1e-3     #kg/kg Threshold cloud water for autoconversion
    qi0::T        = 1e-4     #kg/kg Threshold cloud ice for aggregation
    T0n::T        = 273.16   #K maximum temperature for cloud ice
    T0p::T        = 283.16   #K maximum temperature for snow/graupel
    T0g::T        = 283.16   #K maximum temperature graupel
    T00n::T       = 253.16   #K minimum temperature for cloud water
    T00p::T       = 268.16   #K minimum temperature for rain
    T00g::T       = 223.16   #K minimum temperature for graupel
    VTi::T        = 0.4      #m/s terminal velocity for cloud ice
    α::T          = 0.001    #s^(-1) autoconversion rate
    β::T          = 0.001    #s^(-1) ice aggregation rate
    ρ0::T         = 1.29     #kg/m^3 reference air density
    ρ_rain::T     = 1000     #kg/m^3 rain density
    ρ_snow::T     = 100      #kg/m^3 snow density
    ρ_graupel::T  = 400      #kg/m^3 graupel density
    μ::T          = 1.717e-5 #kg m^(-1)/s Dynamic viscosity of air at 0C
    γ3br::T        = gamma(3+b_rain)
    γ4br::T        = gamma(4+b_rain)
    γ5br::T        = gamma((5+b_rain)/2)
    γ3bs::T        = gamma(3+b_snow)
    γ4bs::T        = gamma(4+b_snow)
    γ5bs::T        = gamma((5+b_snow)/2)
    γ3bg::T        = gamma(3+b_graupel)
    γ4bg::T        = gamma(4+b_graupel)
    γ5bg::T        = gamma((5+b_graupel)/2)
end


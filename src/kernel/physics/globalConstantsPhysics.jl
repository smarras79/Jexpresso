using Parameters

@with_kw struct PhysicalConst{T}

    #Thermodynamic constants at T=300 K
    Rair::T =  287.0 #J/kg.K
    cp::T   = 1005.0 #J/kg.K
    cv::T   =  718.0 #J/kg.K
    γair::T = cp/cv
    Pref::T = 100000.0 #Pa
    
    g::T = 9.80616 #m/s²
    
end

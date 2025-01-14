#### Saturation vapor pressures and mixing ratios
#### Flateau et al (JAM, 1992;1507)

#### Saturation vapor pressure over water
function esatw(T)
    FT = typeof(T)
    a0, a1, a2, a3, a4, a5, a6, a7, a8 = FT(6.11239921), FT(0.443987641), FT(0.142986287e-1),
                                         FT(0.264847430e-3), FT(0.302950461e-5), FT(0.206739458e-7), 
                                         FT(0.640689451e-10), FT(-0.952447341e-13),FT(-0.976195544e-15)
    dT = max(-FT(80),T-FT(273.16))
    esatw = a0 + dT*(a1+dT*(a2+dT*(a3+dT*(a4+dT*(a5+dT*(a6+dT*(a7+dT*a8)))))))
    return FT(esatw)
end

#### Saturation mixing ratio for warm cloud

function qsatw(T,P)
    FT = typeof(T)
    esat = esatw(T)
    qsatw = FT(0.622) * esat/max(esat,P-esat)
    return FT(qsatw)
end

#### Temperature derivative of saturation vapor pressure over water

function dtesatw(T)
    FT = typeof(T)
    a0,a1,a2,a3,a4,a5,a6,a7,a8 = FT(0.443956472), FT(0.285976452e-1), FT(0.794747212e-3), 
                                 FT(0.121167162e-4), FT(0.103167413e-6), FT(0.385208005e-9), 
                                 FT(-0.604119582e-12), FT(-0.792933209e-14), FT(-0.599634321e-17)
    

    dT = max(FT(-80),T-FT(273.16))
    dtesatw = a0 + dT*(a1+dT*(a2+dT*(a3+dT*(a4+dT*(a5+dT*(a6+dT*(a7+dT*a8)))))))
    return FT(dtesatw)
end

#### Temperature derivative of saturation mixing ratio for warm cloud

function dtqsatw(T,P)
    FT = typeof(T)
    dtqsatw = FT(0.622)*dtesatw(T)/P
    return FT(dtqsatw)
end

#### Saturation vapor pressure over ice

function esati(T)
    FT = typeof(T)
    a0,a1,a2,a3,a4,a5,a6,a7,a8 = FT(6.11147274), FT(0.503160820), FT(0.188439774e-1), 
                                 FT(0.420895665e-3), FT(0.615021634e-5),FT(0.602588177e-7), 
                                 FT(0.385852041e-9), FT(0.146898966e-11), FT(0.252751365e-14)
    
    if (T > FT(273.15))
        esati = esatw(T)
    elseif (T > FT(185))
        dT = T -FT(273.16)
        esati = a0 + dT*(a1+dT*(a2+dT*(a3+dT*(a4+dT*(a5+dT*(a6+dT*(a7+dT*a8)))))))
    else
        dT = max(FT(-100), T-FT(273.16))
        esati = FT(0.00763685) + dT*(FT(0.000151069)+dT*FT(7.48215e-07))
    end
    return FT(esati)
end

#### Saturation mixing ratio for ice cloud

function qsati(T,P)
    FT = typeof(T)
    esat = esati(T)
    qsati = FT(0.622)*esat/max(esat,P-esat)
    return FT(qsati)
end

#### Temperature derivative of saturation vapor pressure over ice

function dtesati(T)
    FT = typeof(T)
    a0,a1,a2,a3,a4,a5,a6,a7,a8 = FT(0.503223089), FT(0.377174432e-1),FT(0.126710138e-2), 
                                 FT(0.249065913e-4), FT(0.312668753e-6), FT(0.255653718e-8), 
                                 FT(0.132073448e-10), FT(0.390204672e-13), FT(0.497275778e-16)

    if (T > FT(273.15))
        dtesati = dtesatw(T)
    elseif (T > FT(185))
        dT = T-FT(273.16)
        dtesati = a0 + dT*(a1+dT*(a2+dT*(a3+dT*(a4+dT*(a5+dT*(a6+dT*(a7+dT*a8)))))))
    else
        dT = max(FT(-100),T-FT(273.16))
        dtesati = FT(0.00763685) + dT*(FT(0.000151069)+dT*FT(7.48215e-07))
    end
    return FT(dtesati)
end

#### Temperature derivative of saturation mixing ratio for ice cloud

function dtqsati(T,P)
    FT = typeof(T)
    dtqsati = FT(0.622)*dtesati(T)/P
    return FT(dtqsati)
end

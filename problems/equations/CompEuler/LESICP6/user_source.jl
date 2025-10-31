function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=0.0,
                      y=0.0,
                      z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()
    
    #--------------
    # S(q(x)) = -ρg
    #--------------
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0

    #--------------
    # Coriolis
    #--------------
    lcoriolis = true
    lgeostrophic = true
    
    #--------------
    # SPONGE
    #--------------
    zs = 2250.0
    xr = 0.0
    xl = 0.0
    α  = 0.5
    if (z >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe = α*sinpi(0.5*(z - zs)/(zmax - zs))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = 0.0
    end
    ctop= 1.0*betay_coe

    #if (x >= xr)#nsponge_points * dsy) #&& dbl >= 0.0)
    #    betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    #else
    betaxr_coe = 0.0
    #end
    
    #if (x <= xl)#nsponge_points * dsy) #&& dbl >= 0.0)
    #    betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    #else
    betaxl_coe = 0.0
    #end
    
    cxr = 0.0*betaxr_coe
    cxl = 0.0*betaxl_coe
    cyr = 0.0
    cyl = 0.0
    cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)*(1.0 - cyr)*(1.0 - cyl)

    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    #S[1] -= (cs)*(q[1]-qe[1])
    S[2] -= cs*(q[2]-qe[2])
    S[3] -= cs*(q[3]-qe[3])
    S[4] -= cs*(q[4]-qe[4])
    #S[5] -= cs*(q[5]-qe[5])

    #Coriolis & geostrophic wind    
    if lcoriolis == true
        f = 1.0e-4
        ρ     = q[1]
        u_vel = q[2]
        v_vel = q[3]
        S[2] += ρ * f * v_vel
        S[3] -= ρ * f * u_vel

        if lgeostrophic == true
            U_geo = 10.0
            V_geo = 0.0
            S[2] -= ρ * f * V_geo
            S[3] += ρ * f * U_geo
        end
    end
    
    return  S
end

function user_source!(S, q, qe, npoin, ::CL,::TOTAL; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=24000.0, ngl=5, nely=10,xmin=0.0, xmax=80000.0)
    
    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1] - qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0

    #--------------
    # SPONGE
    #--------------
    zmax = ymax
    if inputs[:lsponge]
        z = y
        zs = inputs[:zsponge]
        xr = 80000#70000.0   # 10 km buffer from right boundary
        xl = 0.0 #10000.0   # 10 km buffer from left boundary
        α  = 0.5

        if (z >= zs)
            betay_coe = α*sinpi(0.5*(z - zs)/(zmax - zs))
        else
            betay_coe = 0.0
        end
        ctop = 1.0*betay_coe

        if (x > xr)
            betaxr_coe = sinpi(0.5*(x - xr)/(xmax - xr))
        else
            betaxr_coe = 0.0
        end

        if (x < xl)
            betaxl_coe = sinpi(0.5*(xl - x)/(xl - xmin))
        else
            betaxl_coe = 0.0
        end

        cxr = 0.25*betaxr_coe
        cxl = 0.25*betaxl_coe
        cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)

        S[2] -= cs*(q[2]-qe[2])
        S[3] -= cs*(q[3]-qe[3])
        S[4] -= cs*(q[4]-qe[4])
    end  #sponge
    
    return S
end 


function user_source!(S, q, qe, npoin, ::CL,::PERT;
                      neqs=1,
                      x=0.0, y=0.0,
                      xmin=0.0, xmax=80000.0,
                      ymin=0.0, ymax=24000.0,
                      ngl=5, nely=10)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0
   
    #--------------
    # SPONGE
    #--------------
    zmax = ymax
    if inputs[:lsponge]
        z = y
        zs = inputs[:zsponge]
        xr = 80000.0#70000.0
        xl = 0.0 #10000.0
        α  = 0.5

        if (z >= zs)
            betay_coe = α*sinpi(0.5*(z - zs)/(zmax - zs))
        else
            betay_coe = 0.0
        end
        ctop = 1.0*betay_coe

        if (x > xr)
            betaxr_coe = sinpi(0.5*(x - xr)/(xmax - xr))
        else
            betaxr_coe = 0.0
        end

        if (x < xl)
            betaxl_coe = sinpi(0.5*(xl - x)/(xl - xmin))
        else
            betaxl_coe = 0.0
        end

        cxr = 0.25*betaxr_coe
        cxl = 0.25*betaxl_coe
        cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)
        
        S[1] -= (cs)*(q[1])
        S[2] -= (cs)*(q[2])
        S[3] -= (cs)*(q[3])
        S[4] -= (cs)*(q[4])
        S[5] -= (cs)*(q[5])
        S[6] -= (cs)*(q[6])
    end  #sponge
   
    return S
end

function user_source_gpu(q,qe,x,y,PhysConst,xmax,xmin,ymax,ymin,lpert)

    T = eltype(x)
    zs = T(15000.0)
    xr = T(70000.0)
    xl = T(10000.0)

    if (y >= zs)
        betay_coe = T(sinpi(T(0.5)*(y-zs)/(ymax-zs))^2)
    else
        betay_coe = T(0.0)
    end
    ctop = betay_coe
    
    if (x > xr)
        betaxr_coe = T(sinpi(T(0.5)*(x-xr)/(xmax-xr))^2)
    else
        betaxr_coe = T(0.0)
    end

    if (x < xl)
        betaxl_coe = T(sinpi(T(0.5)*(xl-x)/(xl-xmin))^2)
    else
        betaxl_coe = T(0.0)
    end

    cxr = T(0.25)*betaxr_coe
    cxl = T(0.25)*betaxl_coe
    ctop = T(0.1)*min(ctop, T(1.0))
    cxr  = min(cxr, T(1.0))
    cxl  = min(cxl, T(1.0))
    cs = T(1.0) - (T(1.0) - ctop)*(T(1.0) - cxr)*(T(1.0) - cxl)

    ρ  = q[1]

    return T(-cs*q[1]), T(-cs*q[2]), T(-cs*q[3]-ρ*PhysConst.g), T(-cs*q[4])
end

function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qe::SubArray{Float64}, npoin, ::CL,::TOTAL; neqs=1,x=0.0, y=0.0, z=0.0, ymin=0.0, zmax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)
   
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1] -qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0
    #### SPONGE

    #
    # clateral
    ys = 18000.0
    xr =  30000.0
    xl = -30000.0
    if (y >= ys)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(y - ys)/(ymax-ys))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = 0.0
    end
    ctop= 1.0*betay_coe
   
    if (x >= xr)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxr_coe = 0.0
    end
   
    if (x <= xl)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxl_coe = 0.0
    end
   
    cxr = 0.0*betaxr_coe
    cxl = 0.0*betaxl_coe
    #@info x,y,cxr,cxl,ctop
    cs = 1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)

    #@info "β x: " ctop,cxr,cxl,cs, ys, x, ymin, ymax, dsy, dbl
      S[1] -= (cs)*(q[1]-qe[1])
      S[2] -= (cs)*(q[2]-qe[2])
      S[3] -= (cs)*(q[3]-qe[3])
      S[4] -= (cs)*(q[4]-qe[4])
    return  S
end 

function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qe::SubArray{Float64}, npoin, ::CL,::PERT; neqs=1,x=0.0, y=0.0, z, zmin=0.0, zmax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0
    #### SPONGE

    #
    # clateral
    nsponge_points = 8

    ys = 18000.0#ymax - 20000.0
    xr = 30000.0
    xl = -30000.0
    
    if (y >= ys)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(z-ys)/(ymax-ys))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
        #betay_coe = 0.9/(1+exp((0.4*ymax-y)/(ymax/18)))
        #betay_coe = 25.0/(1+exp((0.9*ymax-y)/(ymax/15))) ### damps too far down
        #betay_coe = 2.0/(1+exp((0.3*(ymax-15000)-y+15000)/((ymax-15000)/18)))
    else
        betay_coe = 0.0
    end
    #if (abs(x) <=xmin)
      ctop= betay_coe#0.5*betay_coe
    #else
     # ctop = 0.0
    #end 

    if (x >= xr)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxr_coe = 0.0
    end

    if (x <= xl)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxl_coe = 0.0
    end
    
    cxr = 1.0*betaxr_coe#0.25*betaxr_coe
    cxl = 1.0*betaxl_coe#0.25*betaxl_coe
    ctop = 1.0*min(ctop,1)
    cxr  = min(cxr,1)
    cxl  = min(cxl,1)
    cs = 1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)
    
    #@info "β x: " ctop,cxr,cxl,cs, ys, y, x, ymin, ymax, dsy, dbl
    S[1] -= (cs)*(q[1])
    S[2] -= (cs)*(q[2])
    S[3] -= (cs)*q[3]
    S[4] -= (cs)*(q[4])
    S[5] -= (cs)*(q[5])
    #S[6] -= (cs)*(q[6])
    #S[7] -= (cs)*(q[7])
    return  S
end

function user_source_gpu(q,qe,x,y,z,PhysConst,xmax,xmin,ymax,ymin,zmax,zmin,lpert)

    T = eltype(x)
    # distance from the boundary. xs in Restelli's thesis
    ys = T(18000.0)#ymax - 20000.0
    xr = T(30000.0)
    xl = T(-30000.0)

    if (y >= ys)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaz_coe =  T(sinpi(T(0.5)*(y - ys)/(ymax - ys))^2)#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaz_coe = T(0.0)
    end
      ctop= betaz_coe#0.5*betay_coe
    
      if (x > xr)#nsponge_points * dsy) #&& dbl >= 0.0)
          betaxr_coe =  T(sinpi(T(0.5)*(x-xr)/(xmax-xr))^2)#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxr_coe = T(0.0)
    end

    if (x < xl)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe =  T(sinpi(T(0.5)*(xl-x)/(xl-xmin))^2)#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxl_coe = T(0.0)
    end

    cxr = T(1.0)*betaxr_coe#0.25*betaxr_coe
    cxl = T(1.0)*betaxl_coe#0.25*betaxl_coe
    ctop = T(1.0)*min(ctop,1)
    cxr  = min(cxr,1)
    cxl  = min(cxl,1)
    cs = T(1.0) - (T(1.0) -ctop)*(T(1.0)-cxr)*(T(1.0) - cxl)

    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    #return T(0.0), T(0.0), T(0.0), T(-ρ*PhysConst.g), T(0.0), T(0.0), T(0.0) 
    return T(-cs*q[1]), T(-cs*q[2]), T(-cs*q[3]), T(-cs*q[4]-ρ*PhysConst.g), T(-cs*q[5]), T(-cs*q[6]), T(-cs*q[7])
end
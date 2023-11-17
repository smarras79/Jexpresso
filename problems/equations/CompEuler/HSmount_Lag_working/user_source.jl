
function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qe::SubArray{Float64}, npoin, ::CL,::TOTAL; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)
   
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1] -qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    #S[3] = -PhysConst.g
    S[4] = 0.0

    #### SPONGE

    #
    # clateral
    nsponge_points = 8

    # distance from the boundary. xs in Restelli's thesis
    dsy = (ymax - ymin)/(nely*(ngl - 1))# equivalent grid spacing
    dbl = ymax - y
    zs = 14500.0#ymax - 16000.0
    dsx = (xmax - xmin)/(nely*(ngl - 1))# equivalent grid spacing
    dbx = min(xmax - x,x-xmin) 
    xr = 120000.0
    xl = -120000.0
    if (y > zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(y-zs)/(ymax-zs))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = 0.0
    end
    ctop= 0.1*betay_coe
   
    if (x > xr)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxr_coe = 0.0
    end
   
    if (x < xl)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxl_coe = 0.0
    end
   
    cxr = 0.25*betaxr_coe
    cxl = 0.25*betaxl_coe
    #@info x,y,cxr,cxl,ctop
    cs = 1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)

    if (x >= xmin && x <= xmax)     
    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
      S[1] -= (cs)*(q[1]-qe[1])
      S[2] -= (cs)*(q[2]-qe[1]*20.0)
      S[3] -= (cs)*q[3]
      S[4] -= (cs)*(q[4]-qe[4])
    end
    
    return  S
end 

function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qe::SubArray{Float64}, npoin, ::CL,::PERT; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    #S[3] = -PhysConst.g
    S[4] = 0.0

    #### SPONGE

    #
    # clateral
    nsponge_points = 8

    # distance from the boundary. xs in Restelli's thesis
    dsy = (ymax - ymin)/(nely*(ngl - 1))# equivalent grid spacing
    dbl = ymax - y
    zs = 15000.0#ymax - 20000.0
    dsx = (xmax - xmin)/(nely*(ngl - 1))# equivalent grid spacing
    dbx = min(xmax - x,x-xmin)
    xr = 90000.0
    xl = -90000.0
    
    if (y >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(y-zs)/(ymax-zs))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
        #betay_coe = 0.9/(1+exp((0.4*ymax-y)/(ymax/18)))
        #betay_coe = 25.0/(1+exp((0.9*ymax-y)/(ymax/15))) ### damps too far down
        #betay_coe = 10.0/(1+exp((0.67*ymax-y)/(ymax/49)))
    else
        betay_coe = 0.0
    end
    #if (abs(x) <=xmin)
      ctop= betay_coe#0.5*betay_coe
    #else
     # ctop = 0.0
    #end 

    if (x > xr)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxr_coe =  sinpi(0.5*(x-xr)/(xmax-xr))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxr_coe = 0.0
    end

    if (x < xl)#nsponge_points * dsy) #&& dbl >= 0.0)
        betaxl_coe =  sinpi(0.5*(xl-x)/(xl-xmin))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betaxl_coe = 0.0
    end
    
    cxr = 0.05*betaxr_coe#0.25*betaxr_coe
    cxl = 0.05*betaxl_coe#0.25*betaxl_coe
    ctop = 0.5*min(ctop,1)
    cxr  = min(cxr,1)
    cxl  = min(cxl,1)
    cs = 1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)
    
    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    S[1] -= (cs)*(q[1])
    S[2] -= (cs)*(q[2])
    S[3] -= (cs)*q[3]
    S[4] -= (cs)*(q[4])

    return  S
end

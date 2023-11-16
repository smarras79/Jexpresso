
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
   
    cxr = 1.0*betaxr_coe
    cxl = 1.0*betaxl_coe
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
    
    if (x >= 5000.0)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*xmax-x)/(xmax/18)))
    else
        sponge_coe = 0.0
    end
    zt = xmax
    zd = 5000.0
    z = max(x-zd,0.0)
    dgamma  =0.00005
    alpha = 0.85
    sigma = zt/18.0
    fac1 = (alpha*zt-x)/sigma
    #=if (x >=60.0)
      cs = dgamma/(1+exp(fac1))
    else
     cs = 0.0
    end=#
    #cs = 1.0*sinpi(0.5*z/(zt-zd))^2
    cs = min(sponge_coe,1)
    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    S[1] = -(cs)*(q[1])
    S[2] = -(cs)*(q[2])
    return  S
end    

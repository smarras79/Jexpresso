
function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, ρref, npoin, CL; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -77000, xmax =77000)
   
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1] -ρref
    
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
    zs = ymax - 5000.0
    dsx = (xmax - xmin)/(nely*(ngl - 1))# equivalent grid spacing
    dbx = min(xmax - x,x-xmin) 
    xr = 60000.0
    xl = -60000.0
    if (y > zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(y-zs)/(ymax-zs))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = 0.0
    end
    ctop= 0.5*betay_coe
   
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
   
    cxr = 0.1*betaxr_coe
    cxl = 0.1*betaxl_coe
    #@info x,y,cxr,cxl,ctop
    cs = 1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)

    θref = 250.0
    θ0 = 250.0 #K
    T0   = θ0
    p0   = 100000.0

    N    = PhysConst.g/sqrt(PhysConst.cp*T0)
    N2   = N*N
    θ    = θref*exp(N2*y/PhysConst.g)
    #auxi = PhysConst.Rair*θ0
     #            p    = p0*exp(-PhysConst.g*y/auxi)
      #           θ    = θ0*exp(N2*y/PhysConst.g)
    #if (y > 0.1)
      p    = p0*(1.0 + PhysConst.g2*(exp(-y*N2/PhysConst.g) - 1.0)/(PhysConst.cp*θref*N2))^PhysConst.cpoverR
    #else
    #  p = p0
    #end
    #ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p) #kg/m³

    #@info nsponge_points * dsy
     
    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    S[1] -= (cs)*(q[1]-ρref)
    S[2] -= (cs)*(q[2]-ρref*10.0)
    S[3] -= (cs)*q[3]
    S[4] -= (cs)*(q[4]-ρref*θ)
    
    
    return  S
end    

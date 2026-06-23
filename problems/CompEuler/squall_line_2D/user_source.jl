function user_source!(S, q, qe, npoin, ::CL,::TOTAL; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)
    
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1] - qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0

    #### SPONGE
     #--------------
    # SPONGE
    #--------------
    zmax = 24000.0
    if inputs[:lsponge]
        z = y
        @info z
        zs = inputs[:zsponge]
    	xr = 0.0
    	xl = 0.0
    	α  = 0.5
	if (z >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
            betay_coe = α*sinpi(0.5*(z - zs)/(zmax - zs))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
	else
		betay_coe = 0.0
	end
    	ctop= 1.0*betay_coe

    	betaxr_coe = 0.0
        betaxl_coe = 0.0
        
        cxr = 0.0*betaxr_coe
        cxl = 0.0*betaxl_coe
        cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)

        #S[1] -= (cs)*(q[1]-qe[1])
        S[2] -= cs*(q[2]-qe[2])
    	S[3] -= cs*(q[3]-qe[3])
        S[4] -= cs*(q[4]-qe[4])
    end	 #sponge
    
    return  S
end 


function user_source!(S,
                      q,
                      qe, npoin, ::CL,::PERT;
                      neqs=1,
                      x=0.0, y=0.0,
                      xmin=0.0, xmax=1.0,
                      ymin=0.0, ymax=1.0,
                      ngl=5, nely=10)

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
   
    #--------------
    # SPONGE
    #--------------
    zmax = ymax
    if inputs[:lsponge]
        z = y
        zs = inputs[:zsponge]
    	xr = 0.0
    	xl = 0.0
    	α  = 0.5
	if (z >= zs)
            betay_coe = α*sinpi(0.5*(z - zs)/(zmax - zs))#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
	else
		betay_coe = 0.0
	end
    	ctop= 1.0*betay_coe

    	betaxr_coe = 0.0
        betaxl_coe = 0.0
        
        cxr = 0.0*betaxr_coe
        cxl = 0.0*betaxl_coe
        cs  = 1.0 - (1.0 - ctop)*(1.0 - cxr)*(1.0 - cxl)
        
        S[1] -= (cs)*(q[1])
        S[2] -= (cs)*(q[2])
        S[3] -= (cs)*(q[3])
        S[4] -= (cs)*(q[4])
        S[5] -= (cs)*(q[5])
        S[6] -= (cs)*(q[6])
    end	 #sponge
   
    return  S
end

function user_source_gpu(q,qe,x,y,PhysConst,xmax,xmin,ymax,ymin,lpert)

    T = eltype(x)
    # distance from the boundary. xs in Restelli's thesis
    zs = T(15000.0)#ymax - 20000.0
    xr = T(90000.0)
    xl = T(-90000.0)

    if (y >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  T(sinpi(T(0.5)*(y-zs)/(ymax-zs))^2)#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = T(0.0)
    end
      ctop= betay_coe#0.5*betay_coe
    
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

    cxr = T(0.01)*betaxr_coe#0.25*betaxr_coe
    cxl = T(0.01)*betaxl_coe#0.25*betaxl_coe
    ctop = T(0.1)*min(ctop,1)
    cxr  = min(cxr,1)
    cxl  = min(cxl,1)
    cs = T(1.0) - (T(1.0) -ctop)*(T(1.0)-cxr)*(T(1.0) - cxl)

    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    return T(-cs*q[1]), T(-cs*q[2]), T(-cs*q[3]-ρ*PhysConst.g), T(-cs*q[4])
end

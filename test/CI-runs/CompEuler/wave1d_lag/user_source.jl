function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qe::SubArray{Float64}, npoin, ::CL,::TOTAL; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #

    if (x >= 2.49)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
    elseif (x <=-2.49)
        sponge_coe = 2.0/(1+exp(-(0.3*(xmin+2.5)-x-2.5)/((xmax)/18)))
    else
        sponge_coe = 0.0
    end
    #@info xmin,xmax,x,sponge_coe
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

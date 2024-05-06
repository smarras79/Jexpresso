function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ = q[1] #- qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
  
    zs = 8000.0#ymax - 20000.0

    if (y >= zs)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  sinpi(0.5*(abs(x)-zs)/(xmax-zs))^2#1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
        #betay_coe = 0.9/(1+exp((0.4*ymax-y)/(ymax/18)))
        #betay_coe = 25.0/(1+exp((0.9*ymax-y)/(ymax/15))) ### damps too far down
        #betay_coe = 2.0/(1+exp((0.3*(ymax-8000)-y+8000)/((ymax-8000)/18)))
    else
        betay_coe = 0.0
    end

    ctop = 0.0*min(betay_coe,1)
    cs = ctop #1.0 - (1.0 -ctop)*(1.0-cxr)*(1.0 - cxl)

    #@info "β x: " ctop,cxr,cxl,cs, zs, y, x, ymin, ymax, dsy, dbl
    S[1] -= (cs)*(q[1])
    S[2] -= (cs)*(q[2])
    S[3] -= (cs)*q[3]
    S[4] -= (cs)*(q[4])

    return  S
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::NCL,
                      ::AbstractPert; #for NCL() there is no differece between PERT() and TOTAL() in the source
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -g
    #    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0
    
end

function user_source(q,x,y,PhysConst,xmax, xmin, ymax, ymin)


    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    return Float32(0.0), Float32(0.0), Float32(-ρ*PhysConst.g), Float32(0.0)
end

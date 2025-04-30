function user_source(S,
                     q, 
                     qe,
                     npoin::Int64,
                     ::CL, ::TOTAL;
                     neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x))
    #
    L = 2
    alpha = 10
    f   = 0.0 #- (cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y)
    u_e = 0.0 #sin(x/L)*exp(-x/L)*cos(y)
    
    S = f - alpha*u_e
    return S
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    L = 2
    alpha = 10 
    f   = T(0.0) #T(- (cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y))
    u_e = T(0.0) #T(sin(x/L)*exp(-x/L)*cos(y))

    return T(f - alpha*u_e)
end

function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=1.0, y=1.0,
                      xmax=1.0, xmin=0.0,
                      ymax=1.0, ymin=0.0)
  
    
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x))
    #
    #L     = sqrt((xmax-xmin)^2 + (ymax-ymin)^2)
    L     = abs(xmax-xmin)
    
    alpha = 0.0 #1.0
    beta  = 0.0 #1.0
    gamma = 0.0 #1.0
    
    f   = - beta*(cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y)
    
    u_e = gamma*sin(x/L)*exp(-x/L)*cos(y)
    
    return f - alpha*u_e
    
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    L = 2
    alpha = 10 
    f   = T(0.0) #T(- (cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y))
    u_e = T(0.0) #T(sin(x/L)*exp(-x/L)*cos(y))

    return T(f - alpha*u_e)
end

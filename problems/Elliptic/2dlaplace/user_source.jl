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
    L = 2
    alpha = 10
    f   = 0.0 #- (cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y)
    u_e = 0.0 #sin(x/L)*exp(-x/L)*cos(y)
    
    S[1] = f - alpha*u_e
    
end

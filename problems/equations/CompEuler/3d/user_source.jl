function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1)
    
    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0 #Y is the vertical direction in 3D
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0
   
end

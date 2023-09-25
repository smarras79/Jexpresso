function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qref::Float64, npoin::Int64, ::CL; neqs=1)

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

function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, qref::Float64, npoin::Int64, ::NCL; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0
    
end

function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, ρref::Float64, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    dρ  = q[1] - ρref
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -dρ*PhysConst.g
    S[4] = 0.0
    
end

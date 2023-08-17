function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin::Int64; neqs=1)

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
function user_source(T, q::Array, npoin::Int64; neqs=1, x=0.0)

    #
    # This is a fallback function used by rhs
    # if there is no user_source.jl file in the `ARGS[1], ARGS[2],` directory
    # 
    S = @MVector zeros(T, neqs)
    
    #
    # S(q(x)) = 0.0
    #

    S[1] = 0.0
    
    return  S
    
end

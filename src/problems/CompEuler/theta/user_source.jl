function user_source(T, q::Array, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
    S = zeros(T, npoin, neqs)
    
    #
    # S(q(x)) = -ρg
    #
    for ip=1:npoin
        ρ  = q[ip,1]
        
        S[ip,1] = 0.0
        S[ip,2] = 0.0
        S[ip,3] = -ρ*PhysConst.g
        S[ip,4] = 0.0
    end
    
    return  S
    
end

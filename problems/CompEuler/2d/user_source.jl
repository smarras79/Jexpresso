function user_source!(S, q, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{TFloat}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    
end

function user_source(q::Array, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{TFloat}()
    S = zeros(TFloat, neqs)
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    
end

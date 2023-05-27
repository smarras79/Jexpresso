
function user_source(T, q::Array, npoin::Int64; neqs=1, x=0.0)

    PhysConst = PhysicalConst{Float64}()
    S = zeros(T, neqs)
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    
    return  S
    
end

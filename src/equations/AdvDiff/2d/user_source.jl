function user_source(T, q::Array, npoin::Int64; neqs=1, x=0.0)

    S = @MVector zeros(T, neqs)
    
    #
    # S(q(x)) = 0.0
    #

    S[1] = 0.0
    
    return  S
    
end

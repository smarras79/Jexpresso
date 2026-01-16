function user_source!(S::SubArray{TFloat}, q::SubArray{TFloat}, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{TFloat}()
        
    #
    # S(q(x)) = -ρg
    #
    S[:] .= 0.0
    
end

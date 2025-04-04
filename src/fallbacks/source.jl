function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    S[:] .= 0.0
    
end

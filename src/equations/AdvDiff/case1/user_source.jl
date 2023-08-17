function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    S[1] = 0.0
    
end

function user_source!(S, T, q::Array, npoin::Int64, x::Array)
    
    #
    # S(q(x)) = βsin(γx)
    #
    β, γ = 10000, π;

    S[1] = β*sin.(γ*x)
    
end

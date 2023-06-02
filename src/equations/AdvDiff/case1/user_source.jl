function user_source!(S, T, q::Array, npoin::Int64, x::Array)
    
    #
    # S(q(x)) = βsin(γx)
    #
    β, γ = 10000, π;

    S[1] = β*sin.(γ*x)
    
end

function user_source(T, q::Array, npoin::Int64, x::Array)

    S = zeros(T, npoin)
    
    #
    # S(q(x)) = βsin(γx)
    #
    β, γ = 10000, π;

    S = β*sin.(γ*x)
    
    return  S
    
end

function user_flux(T, SD::NSD_1D, q::Array, npoin::Int64)

    F = zeros(T, npoin)
    
    #
    # F(q(x)) = 0.8*q
    #
    F .= 0.8*q[:, 1]
    
    return F
    
end

function user_flux(T, SD::NSD_2D, q::Array, npoin::Int64)

    F = G = zeros(T, npoin)
    
    #
    # F(q(x)) = 0.8*q
    # G(q(x)) = 0.8*q
    #
    F .= 0.8*q[:, 1]
    G .= 0.8*q[:, 1]
    
    return F, G
end

function user_flux(T, SD::NSD_1D, q::Array, npoin::Int64)
    #
    # F(q(x)) = 0.8*q
    #
    F = zeros(T, npoin)
    
    F .= 1.0*q[:,1]
    
    return F
    
end

function user_flux(T, SD::NSD_2D, q::Array, npoin::Int64)
    #
    # F(q(x)) = 0.8*q
    # G(q(x)) = 0.8*q
    #
    F = G = zeros(T, npoin)
    
    F .= 0.8*q[:,1]
    G .= 0.8*q[:,1]
    
    return F, G
end

function user_flux(T, SD::NSD_3D, q::Array, npoin::Int64)

    #
    # F(q(x)) = 0.8*q
    # G(q(x)) = 0.8*q
    # H(q(x)) = 0.0
    #
    F = G = H = zeros(T, npoin)
    
    F .= 0.8*q[:,1]
    G .= 0.8*q[:,1]
    
    return F, G, H
end

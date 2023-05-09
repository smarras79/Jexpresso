function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh)
    #
    # F(q(x)) = 0.8*q
    #
    F = zeros(T, mesh.npoin)
    
    F .= 1.0*q[:,1]
    
    return F
    
end

function user_flux(T, SD::NSD_2D, q::Array, x, y; neqs=1)

    F = @SVector zeros(neqs)
    G = @SVector zeros(neqs)
        
    F = 0.8*q[1]
    G = 0.8*q[1]
    
    return F, G
end


function user_flux_old(T, SD::NSD_2D, q::Array, mesh::St_mesh)
    #
     #F(q(x)) = 0.8*q
     #G(q(x)) = 0.8*q
    #
    F = G = zeros(T, mesh.npoin)
    
   # F .= 0.8*q[:,1]
    #G .= 0.8*q[:,1]

    for ip=1:mesh.npoin
        
        F[ip] = 0.8*q[ip,1]
        G[ip] = 0.8*q[ip,1]
    end
    return F, G
end


function user_flux(T, SD::NSD_3D, q::Array, mesh::St_mesh)

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

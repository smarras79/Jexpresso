function user_flux!(F, G, T, SD::NSD_2D, q::Array, mesh::St_mesh; neqs=1)
    #
    #F(q(x)) = 0.8*q
    #G(q(x)) = 0.8*q
    #    
    F[1] = 0.8*q[1]
    G[1] = 0.8*q[1]
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, mesh::St_mesh; neqs=4)

    #
    #F(q(x)) = 0.8*q
    #G(q(x)) = 0.8*q
    #
    
    F[1] = 0.8*q[1]
    G[1] = 0.8*q[1]
    
end

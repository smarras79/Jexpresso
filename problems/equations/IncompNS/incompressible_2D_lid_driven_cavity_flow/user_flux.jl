# function user_flux!(F::SubArray{TFloat}, G::SubArray{TFloat}, SD::NSD_2D,
#     q,derivative_y,derivative_x,
#     qe::SubArray{TFloat},
#     mesh::St_mesh,
#     ::CL, ::TOTAL; 
#     neqs=1, ip=1)   

#     F[1] = -q[1]*derivative_y
#     G[1] =  q[1]*derivative_x

# end

function user_flux!(F::SubArray{TFloat}, G::SubArray{TFloat}, SD::NSD_2D,
    q,
    qe::SubArray{TFloat},
    mesh::St_mesh,
    ::CL, ::TOTAL; 
    neqs=1, ip=1, ∇ψ = zeros(TFloat,2))   

    F[1] = -q[1]*∇ψ[2]
    G[1] =  q[1]*∇ψ[1]

end

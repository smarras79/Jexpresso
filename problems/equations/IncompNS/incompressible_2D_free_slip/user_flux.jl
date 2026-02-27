function user_flux!(F::SubArray{TFloat}, G::SubArray{TFloat}, SD::NSD_2D,
    q,
    qe::SubArray{TFloat},
    mesh::St_mesh,
    ::CL, ::TOTAL; 
    neqs=1, ip=1, ∇ψ = zeros(TFloat,2))   

    F[1] = -q[1]*∇ψ[2]
    G[1] =  q[1]*∇ψ[1]

end
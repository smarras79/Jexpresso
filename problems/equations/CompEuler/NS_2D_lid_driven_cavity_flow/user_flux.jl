function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
    q,
    qe::SubArray{Float64},
    mesh::St_mesh,
    ::CL, ::TOTAL; 
    neqs=1, ip=1, ∇ψ = zeros(Float64,2))   

    F[1] = -q[1]*∇ψ[2]
    G[1] =  q[1]*∇ψ[1]

end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, x=1)
    γ = 1.4
    γm1 = 0.4

    A  = 1.0 + 2.2*(x - 1.5)^2

    ρ  = q[1]/A
    u  = q[2]/q[1]
    T  = γm1*(q[3]/q[1] - 0.5*γ*u*u)
    p  = ρ*T
    
    F[1] = q[2]
    F[2] = q[2]*u + p*A/γ
    F[3] = q[3]*u + p*A*u
    
end

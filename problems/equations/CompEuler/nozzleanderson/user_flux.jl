function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1, kwargs...)
    γ = 1.4
    γm1 = 0.4

    A  = 1.0 + 2.2*(mesh.x[ip] - 1.5)^2

    ρ  = q[1]/A
    u  = q[2]/q[1]
    T  = γm1*(q[3]/q[1] - 0.5*γ*u*u)
    p  = ρ*T
    
    F[1] = q[2]
    F[2] = q[2]*u + p*A/γ
    F[3] = q[3]*u + p*A*u
    
end

function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    F[1] = q[2]
    F[2] = q[1]
    
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    return T(q[2]), T(q[1])
end

# Max wave-speed hook for the DG numerical flux.
# wave1d is the linear acoustic system u_t + A u_x = 0, A = [0 1; 1 0],
# eigenvalues ±1 ⇒ λ_max = 1.
function user_max_wave_speed(q, qe, SD::NSD_1D, ::TOTAL; neqs=2)
    return 1.0
end

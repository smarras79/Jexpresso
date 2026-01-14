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

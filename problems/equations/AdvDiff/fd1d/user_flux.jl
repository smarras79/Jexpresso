function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1, kwargs...)

    u = 1.0
    F[1] = u*q[1]
    
end

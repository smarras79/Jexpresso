function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)
    
    u = 1.0
    v = 0.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end


function user_flux_gpu(q,qe,PhysConst,lpert)
    
    T = eltype(q) 
    return T(q[1]), T(0.0)

end

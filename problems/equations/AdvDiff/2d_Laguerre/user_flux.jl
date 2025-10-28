function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1, kwargs...)
    
    u = 0.5
    v = 1.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=1, ip=1, kwargs...)
    
    u = 0.5
    v = 1.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end


function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=1, ip=1, kwargs...)
    
    u = 0.5
    v = 1.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    
    T = eltype(q) 
    return T(q[1]*0.5), T(q[1]*1.0)

end

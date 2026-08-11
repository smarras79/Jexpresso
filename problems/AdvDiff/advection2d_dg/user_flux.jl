advection_velocity() = (0.5, 1.0)  # single source: consumed by ALL user_flux! methods,
                                   # user_flux_gpu, and user_max_wave_speed (F0 lesson)

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)
    
    u, v = advection_velocity()
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=1, ip=1)
    
    u, v = advection_velocity()
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end


function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=1, ip=1)
    
    u, v = advection_velocity()
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    u, v = advection_velocity()
    return T(q[1]*u), T(q[1]*v)
end

function user_max_wave_speed(q, qe, SD::NSD_2D, ::TOTAL; neqs=1)
    u, v = advection_velocity()
    return sqrt(u*u + v*v)  # provisional global bound; the 2D surface term
                            # (Phase 5 step 6) decides whether this hook gains
                            # a face-normal argument for |c·n| per face
end

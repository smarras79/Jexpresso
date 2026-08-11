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

# Rusanov λ per face: |c·n|, the wave speed in the face-normal direction.
# (nx, ny) is the face unit normal supplied by surface_rhs_el!(::NSD_2D);
# the defaults keep any normal-free call meaningful (they reduce λ to |u|).
# Settled 2026-08-11 (Q3): the hook carries the normal so that at Euler/SWE,
# where λ = |u·n| + a genuinely needs it, no signature change is required.
function user_max_wave_speed(q, qe, SD::NSD_2D, ::TOTAL; nx::Float64=1.0, ny::Float64=0.0, neqs=1)
    u, v = advection_velocity()
    return abs(u*nx + v*ny)
end

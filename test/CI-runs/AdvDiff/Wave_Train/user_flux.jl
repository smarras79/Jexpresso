function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    U = qe[2]
    H = qe[1]
    u = q[2]
    h = q[1]-qe[1]
    g = 9.81
    F[1] = U*h + H*u
    F[2] = g*h + U*u
    
end

function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)

    U = qe[2]
    H = qe[1]
    u = q[2]
    h = q[1]
    g = 9.81
    F[1] = U*h + H*u
    F[2] = g*h + U*u
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    U = qe[2]
    H = qe[1]
    u = q[2]
    h = q[1]
    g = T(9.81)
    return T(U*h + H*u), T(g*h + U*u)
end

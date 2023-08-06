function user_flux!(F, G, SD::NSD_2D, q, mesh::St_mesh, ip; neqs=1)

    
    prof = exp(-8.0*(mesh.x[ip])^2)
    F[1] = prof*0.5*q[1]
    G[1] = 0.0#0.8*q[1]
    
    return F, G
end

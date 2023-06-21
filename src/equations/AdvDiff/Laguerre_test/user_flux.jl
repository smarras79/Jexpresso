function user_flux(T, SD::NSD_2D, q::Array, mesh::St_mesh; neqs=1,ip)

    F = @MVector zeros(T, neqs)
    G = @MVector zeros(T, neqs)
    
    prof = exp(-8.0*(mesh.x[ip])^2)
    F[1] = prof*0.5*q[1]
    G[1] = 0.0#0.8*q[1]
    
    return F, G
end

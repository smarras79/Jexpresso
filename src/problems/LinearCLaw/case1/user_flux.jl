function user_flux(T, SD::NSD_2D, q::Array, mesh::St_mesh)
    G = zeros(T, mesh.npoin, 3)
    F = zeros(T, mesh.npoin, 3)
    c= 1.0
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        y = mesh.y[ip]
        #x = y = 0.8
        p = q[ip,1]
        u = q[ip,2]
        v = q[ip,3]

        F[ip,1] =  c^2*u
        G[ip,1] =  c^2*v
        F[ip,2] =  p
        G[ip,2] =  0.0
        F[ip,3] =  0.0
        G[ip,3] =  p
    end
    return F, G
end

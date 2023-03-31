include("bathymetry.jl")

function user_flux(T, SD::NSD_2D, q::Array, mesh::St_mesh)
    G = zeros(T, mesh.npoin, 3)
    F = zeros(T, mesh.npoin, 3)
    F1 = zeros(T, mesh.npoin, 3)
    G1 = zeros(T, mesh.npoin, 3)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        y = mesh.y[ip]
        Hb = bathymetry(x,y)
        H = q[ip,1]
        Hu = q[ip,2]
        Hv = q[ip,3]
        u = Hu/H
        v = Hv/H
        F[ip,1]  = Hu
        G[ip,1]  = Hv
        F1[ip,1] = 0.0
        G1[ip,1] = 0.0
        F[ip,2]  = Hb
        G[ip,2]  = 0.0
        F1[ip,2] = (9.81/2) * (H^2 - Hb^2) + Hu * u + Hu * v 
        G1[ip,2] = 0.0
        F[ip,3]  = 0.0
        G[ip,3]  = Hb
        F1[ip,3] = 0.0
        G1[ip,3] = (9.81/2) * (H^2 - Hb^2) + Hv * v + Hv * u
    end
    return F, G, F1, G1
end

function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh)
    F = zeros(T, mesh.npoin, 2)
    F1 = zeros(T, mesh.npoin, 2)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        Hb = bathymetry(x)
        H = q[ip,1]
        Hu = q[ip,2]
        u = Hu/H
        F[ip,1]  = Hu
        F1[ip,1] = 0.0
        F[ip,2]  = Hb
        F1[ip,2] = (9.81/2) * (H^2 - Hb^2) + Hu * u
    end
    return F, F1
end

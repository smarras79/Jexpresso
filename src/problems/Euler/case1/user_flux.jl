function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh; neqs=3)

    PhysConst = PhysicalConst{Float64}()

    neqs = 3
    F = zeros(T, mesh.npoin, neqs)
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρE = q[ip,3]
        u = ρu/ρ
        E = ρE/ρ

        Temp = (E - 0.5*u*u)/PhysConst.cv
        Press = perfectGasLaw(PhysConst; ρ=ρ, Temp=Temp)
        
        F[ip,1] = ρu
        F[ip,2] = ρu*u + Press
        F[ip,3] = ρE*u + u*Press
    end
    return F
end

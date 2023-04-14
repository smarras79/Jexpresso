include("constitutiveLaw.jl")

function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh)

    F = zeros(T, mesh.npoin)
    
    PhysConst = PhysicalConst{Float64}()
    
    for ip=1:mesh.npoin
        x = mesh.x[ip]
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρE = q[ip,3]
        u = ρu/ρ
        E = ρE/ρ

        Temp = (E - 0.5*u*u)/PhysConst.cv
        constitutiveLaw!(ρ, Temp, Press)
        
        F[ip,1] = ρu
        F[ip,2] = ρu*u + Press
        F[ip,3] = ρE*u
    end
    return F, F1
end

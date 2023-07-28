function user_flux!(F::Array{Float64}, G::Array{Float64}, SD::NSD_2D, q::Array, iel, mesh::St_mesh; neqs=4)
    
    PhysConst = PhysicalConst{Float64}()

    for j=1:mesh.ngl, i=1:mesh.ngl
        
        ip = mesh.connijk[i,j,iel]
        
        ρ  = q[ip,1]
        ρu = q[ip,2]
        ρv = q[ip,3]
        ρθ = q[ip,4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        
        Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
        
        F[i,j,1] = ρu
        F[i,j,2] = ρu*u + Press
        F[i,j,3] = ρv*u
        F[i,j,4] = ρθ*u
        
        G[i,j,1] = ρv
        G[i,j,2] = ρu*v
        G[i,j,3] = ρv*v + Press
        G[i,j,4] = ρθ*v
        
    end
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, pref::Float64, mesh::St_mesh; neqs=4)

    PhysConst = PhysicalConst{Float64}()
                
    ρ  = q[1]
    u = q[2]
    v = q[3]
    θ = q[4]
        
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρ*u
    F[2] = u
    F[3] = v
    F[4] = θ
    
    G[1] = ρ*v
    G[2] = u
    G[3] = v
    G[4] = θ
end

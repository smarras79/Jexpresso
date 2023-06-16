function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, mesh::St_mesh; neqs=4)
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    PhysConst = PhysicalConst{Float64}()    
    Pressure = zeros(Float64, 1)
    perfectGasLaw_ρθtoP!(Pressure, PhysConst;  ρ=ρ, θ=θ)
        
    F[1] = ρu
    F[2] = ρu*u .+ Pressure[1]
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure[1]
    G[4] = ρθ*v
    
end

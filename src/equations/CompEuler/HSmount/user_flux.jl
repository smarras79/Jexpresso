function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, Pref, mesh::St_mesh, ip; neqs=4)
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    PhysConst = PhysicalConst{Float64}()    
    Pressure = zeros(Float64, 1)
    #@info mesh.x[ip], mesh.y[ip], ρ,u,v,θ
    perfectGasLaw_ρθtoP!(Pressure, PhysConst;  ρ=ρ, θ=θ)
     
       

    F[1] = ρu
    F[2] = ρu*u .+ Pressure[1] - Pref
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure[1] - Pref
    G[4] = ρθ*v

    #=F[1] = ρ*u
    F[2] = u*u .+ (1/ρ)*Pressure[1]
    F[3] = v*u
    F[4] = θ*u

    G[1] = ρ*v
    G[2] = u*v
    G[3] = v*v .+ (1/ρ)*Pressure[1]
    G[4] = θ*v=#
    
    #@info F,G, mesh.x[ip], mesh.y[ip],ρ,ρu,u 

end

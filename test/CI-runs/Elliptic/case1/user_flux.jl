function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4)

    PhysConst = PhysicalConst{Float64}()

    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0
    #=
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure
    G[4] = ρθ*v
    =#
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4)

    PhysConst = PhysicalConst{Float64}()

    
    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0
end


function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=4)
    
    PhysConst = PhysicalConst{Float64}()
                
    F[1] = 0.0
    F[2] = 0.0
    F[3] = 0.0
    F[4] = 0.0

    G[1] = 0.0 
    G[2] = 0.0
    G[3] = 0.0
    G[4] = 0.0

end

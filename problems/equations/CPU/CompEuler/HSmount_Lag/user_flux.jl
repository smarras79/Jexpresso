function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, qe::SubArray{Float64}, mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    PhysConst = PhysicalConst{Float64}()    
       
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
 
    F[1] = ρu
    F[2] = ρu*u + Press-qe[end]
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press-qe[end]
    G[4] = ρθ*v
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, qe::SubArray{Float64}, mesh::St_mesh, ::CL, ::PERT; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2] + qe[2]
    ρv = q[3] + qe[3]
    ρθ = q[4] + qe[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    PhysConst = PhysicalConst{Float64}()
   
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ) 

    F[1] = ρu
    F[2] = ρu*u + Press-qe[end]
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press-qe[end]
    G[4] = ρθ*v

end


function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D, q::SubArray{Float64}, pref::Float64, mesh::St_mesh, ::NCL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
                
    ρ = q[1]
    u = q[2]
    v = q[3]
    θ = q[4]
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρ*u
    F[2] = u
    F[3] = u
    F[4] = θ
    
    G[1] = ρ*v
    G[2] = v
    G[3] = v
    G[4] = θ
end

function user_flux_gpu(q,qe,PhysConst,lpert)

    T = eltype(q)
    ρ  = q[1]+qe[1]
    ρu = q[2]+qe[2]
    ρv = q[3]+qe[3]
    ρθ = q[4]+qe[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ) - qe[5]
    return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρθ*u), T(ρv),T(ρu*v),T(ρv*v + Pressure),T(ρθ*v)
end

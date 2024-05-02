function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, H::SubArray{Float64},
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρθ = q[5]
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρu*v
    F[4] = ρu*w
    F[5] = ρθ*u

    G[1] = ρv
    G[2] = ρv*u
    G[3] = ρv*v .+ Pressure
    G[4] = ρv*w
    G[5] = ρθ*v
    
    H[1] = ρw
    H[2] = ρw*u
    H[3] = ρw*v
    H[4] = ρw*w .+ Pressure
    H[5] = ρθ*w
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, H::SubArray{Float64},
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρθ = q[5] + qe[5]
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    Pressure = Pressure - qe[end]

    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρv*u
    F[4] = ρw*u
    F[5] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure
    G[4] = ρw*v
    G[5] = ρθ*v
    
    H[1] = ρw
    H[2] = ρu*w
    H[3] = ρv*w
    H[4] = ρw*w .+ Pressure
    H[5] = ρθ*w
    
end

function user_flux(q,PhysConst)

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρθ = q[5]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    w  = ρw/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    return Float32(ρu), Float32(ρu*u + Pressure), Float32(ρv*u), Float32(ρw*u), Float32(ρθ*u), Float32(ρv),Float32(ρu*v),Float32(ρv*v + Pressure),Float32(ρw*v), Float32(ρθ*v), Float32(ρw), Float32(ρu*w), Float32(ρv*w),Float32(ρw*w + Pressure), Float32(ρθ*w)
end

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

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    if (lpert)
        ρ  = q[1]+qe[1]
        ρu = q[2]+qe[2]
        ρv = q[3]+qe[3]
        ρw = q[4]+qe[4]
        ρθ = q[5]+qe[5]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        w  = ρw/ρ
        Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ) - qe[6]
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρw*u), T(ρθ*u), T(ρv), T(ρu*v), T(ρv*v + Pressure), T(ρw*v), T(ρθ*v), T(ρw), T(ρu*w), T(ρv*w), T(ρw*w + Pressure), T(ρθ*w)
    else
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
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρw*u), T(ρθ*u), T(ρv), T(ρu*v), T(ρv*v + Pressure), T(ρw*v), T(ρθ*v), T(ρw), T(ρu*w), T(ρv*w), T(ρw*w + Pressure), T(ρθ*w)
    end
end

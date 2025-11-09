function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=6, ip=1, kwargs...)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    qtr= q[5]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    F[5] = q[5]*u
    F[6] = q[6]*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
    G[5] = q[5]*v
    G[6] = q[6]*v
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=5, ip=1, kwargs...)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4] + qe[4]
    qtr= q[5] + qe[5]
    qtr2=q[6] + qe[6]
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    Press = Press - qe[end]
    
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    F[5] = qtr*u
    F[6] = qtr2*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
    G[5] = qtr*v
    G[6] = qtr2*v
end


function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=5, ip=1, kwargs...)
    
    PhysConst = PhysicalConst{Float64}()
                
    ρ = q[1]
    u = q[2]
    v = q[3]
    θ = q[4]
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρ*u
    F[2] = u
    F[3] = v
    F[4] = θ
    F[5] = q[5]
    F[6] = q[6]
    
    G[1] = ρ*v
    G[2] = u
    G[3] = v
    G[4] = θ
    G[5] = q[5]
    G[6] = q[6]
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    if (lpert)
        ρ  = q[1]+qe[1]
        ρu = q[2]+qe[2]
        ρv = q[3]+qe[3]
        ρθ = q[4]+qe[4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ) - qe[5]
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρθ*u), T(ρv),T(ρu*v),T(ρv*v + Pressure),T(ρθ*v)
    else
        ρ  = q[1]
        ρu = q[2]
        ρv = q[3]
        ρθ = q[4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρθ*u), T(q[5]), T(q[6]), T(ρv),T(ρu*v),T(ρv*v + Pressure),T(ρθ*v), T(q[5]), T(q[6])
    end
end

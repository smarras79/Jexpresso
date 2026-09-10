function user_primitives!(u, qe, uprimitive, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρw = u[4]
    ρE = u[5]
    p  = PhysConst.γm1 * (ρE - 0.5*(ρu^2 + ρv^2 + ρw^2)/ρ)
    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρw/ρ
    uprimitive[5] = p / (ρ * PhysConst.Rair)          # T = p/(ρ·R)
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρw = u[4]
    ρE = u[5]
    p  = PhysConst.γm1 * (ρE - 0.5*(ρu^2 + ρv^2 + ρw^2)/ρ)
    return SVector(ρ, ρu/ρ, ρv/ρ, ρw/ρ, p/(ρ*PhysConst.Rair))
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    PhysConst = PhysicalConst{T}()
    ρ  = u[1]; ρu = u[2]; ρv = u[3]; ρw = u[4]; ρE = u[5]
    p  = PhysConst.γm1 * (ρE - T(0.5)*(ρu^2 + ρv^2 + ρw^2)/ρ)
    return T(ρ), T(ρu/ρ), T(ρv/ρ), T(ρw/ρ), T(p/(ρ*PhysConst.Rair))
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    PhysConst = PhysicalConst{Float64}()
    ρ  = u[1]
    u1 = u[2]/ρ
    u2 = u[3]/ρ
    u3 = u[4]/ρ
    velomagsq = u1*u1 + u2*u2 + u3*u3
    p  = PhysConst.γm1 * (u[5] - 0.5*ρ*velomagsq)
    uout[1] = ρ
    uout[2] = u1
    uout[3] = u2
    uout[4] = u3
    uout[5] = p
end

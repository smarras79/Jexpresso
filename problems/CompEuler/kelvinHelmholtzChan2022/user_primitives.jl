function user_primitives!(u, qe, uprimitive,::TOTAL)

    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    # Slot 4 is ρθ when :energy_equation == "theta", ρE otherwise.
    # SGS_diffusion follows the same convention (see SGS.jl).
    if ENERGY_EQUATION_THETA[]
        ρθ = u[4]
        uprimitive[1] = ρ
        uprimitive[2] = ρu/ρ
        uprimitive[3] = ρv/ρ
        uprimitive[4] = ρθ/ρ                          # θ
    else
        ρE = u[4]
        p  = PhysConst.γm1 * (ρE - 0.5f0 * (ρu^2 + ρv^2)/ρ)
        uprimitive[1] = ρ
        uprimitive[2] = ρu/ρ
        uprimitive[3] = ρv/ρ
        uprimitive[4] = p / (ρ * PhysConst.Rair)      # T
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)

    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    if ENERGY_EQUATION_THETA[]
        ρθ = u[4]
        return SVector(ρ, ρu/ρ, ρv/ρ, ρθ/ρ)
    else
        ρE = u[4]
        p  = PhysConst.γm1 * (ρE - 0.5 * (ρu^2 + ρv^2)/ρ)
        return SVector(ρ, ρu/ρ, ρv/ρ, p / (ρ * PhysConst.Rair))
    end
end


function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = (u[4]+qe[4])/(u[1]+qe[1])-qe[4]/qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()

    uout[1] = u[1]      #ρ
    uout[2] = u[2]/u[1] #u
    uout[3] = u[3]/u[1] #v

    if ENERGY_EQUATION_THETA[]
        θ        = u[4]/u[1]                                          # θ
        uout[4]  = perfectGasLaw_ρθtoP(PhysConst, uout[1], θ)         # P
    else
        velomagsq = uout[2]*uout[2] + uout[3]*uout[3]
        uout[4]   = PhysConst.γm1 * (u[4] - 0.5 * u[1] * velomagsq)   # P
    end

    uout[5] = uout[4]/(PhysConst.Rair*uout[1])                        # T = p/(ρ*R)
end

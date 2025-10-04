function user_primitives!(u,qe,uprimitive,::TOTAL)
        
    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    E  = u[4]

    u  = ρu/ρ
    v  = ρv/ρ
    KE = 0.5 * ρ * (u^2 + v^2)
    p  = PhysConst.γm1 * (E - KE)
    T  = p / (ρ * PhysConst.Rair)

    uprimitive[1] = ρ
    uprimitive[2] = u
    uprimitive[3] = v
    uprimitive[4] = T
    
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
        
    velomagsq = uout[2]*uout[2] + uout[3]*uout[3]
    uout[4] = PhysConst.γm1 * (u[4] - 0.5 * u[1] * velomagsq) #Pressure
    
    uout[5] = uout[4]/(PhysConst.Rair*uout[1]) # T = p/(ρ*R)
    
end

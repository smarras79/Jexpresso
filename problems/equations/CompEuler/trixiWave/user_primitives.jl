function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[1]+qe[2]
end


function user_uout!(ip, ET, uout, u, qe; kwargs...)
    
    PhysConst = PhysicalConst{Float64}()

    ρ       = u[1]
    uout[1] = ρ      # ρ
    uout[2] = u[2]/ρ # u
    E       = u[3]   # E = ρe
    
    velomagsq = uout[2]*uout[2]
    uout[3] = PhysConst.γm1 * (E - 0.5 * ρ * velomagsq) #p   
    uout[4] = uout[3]/(PhysConst.Rair*ρ)                # T = p/(ρ*R)
    
end

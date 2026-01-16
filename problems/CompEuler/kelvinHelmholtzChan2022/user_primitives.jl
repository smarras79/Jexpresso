function user_primitives!(u, qe, uprimitive,::TOTAL)
        
    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    u           = ρu/ρ                            # Velocity u
    v           = ρv/ρ                            # Velocity v
    E_specific  = ρE/ρ                            # Specific total energy
    KE_specific = 0.5 * (u^2 + v^2)               # Specific kinetic energy
    ei_specific = E_specific - KE_specific        # Specific internal energy
    p           = PhysConst.γm1 * ρ * ei_specific # Pressure
    T           = p / (ρ * PhysConst.Rair)        # Temperature
    
    uprimitive[1] = ρ
    uprimitive[2] = u
    uprimitive[3] = v
    uprimitive[4] = T 
end


function user_primitives!(u, qe, uprimitive, ::THETA)
    PhysConst = PhysicalConst{Float64}()
                
    ρ  = u[1] 
    ρu = u[2]
    ρv = u[3]
    ρθ = u[4] 
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    uprimitive[1] = ρ
    uprimitive[2] = u
    uprimitive[3] = v
    uprimitive[4] = θ
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

function user_uout!(ip, ET::THETA, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()
    
    uout[1] = u[1]      #ρ
    uout[2] = u[2]/u[1] #u
    uout[3] = u[3]/u[1] #v
    θ       = u[4]/u[1] #θ

    uout[4] = perfectGasLaw_ρθtoP(PhysConst; ρ=uout[1], θ=θ) #P
   
    uout[5] = uout[4]/(PhysConst.Rair*uout[1]) # T = p/(ρ*R)
   
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    
    PhysConst = PhysicalConst{Float64}()
    
    uout[1] = u[1]      #ρ
    uout[2] = u[2]/u[1] #u
    uout[3] = u[3]/u[1] #v
    uout[4] = u[4]/u[1] #e
        
    velomagsq = uout[2]*uout[2] + uout[3]*uout[3]
    uout[4] = PhysConst.γm1 * (u[4] - 0.5 * u[1] * velomagsq) #Pressure
    
    uout[5] = uout[4]/(PhysConst.Rair*uout[1]) # T = p/(ρ*R)
    
end

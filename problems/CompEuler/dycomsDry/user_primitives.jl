function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = u[4]/(u[1]+qe[1])
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()
    
    uout[1] = u[1]      #ρ
    uout[2] = u[2]/u[1] #u
    uout[3] = u[3]/u[1] #v
    uout[4] = u[4]/u[1] #w
    uout[5] = u[5]/u[1] #θ
    uout[end] = perfectGasLaw_ρθtoP(PhysConst; ρ=uout[1], θ=uout[5]) #P

end

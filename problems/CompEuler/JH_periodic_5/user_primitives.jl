function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = (u[4]+qe[4])/(u[1]+qe[1])-qe[4]/qe[1]
end


function user_uout!(ip, ET, uout, u, qe; mp = mp)
    PhysConst = PhysicalConst{Float64}()

    if ET == TOTAL()
        uout[1] = u[1]
        uout[2] = u[2]/u[1]
        uout[3] = u[3]/u[1]
        uout[4] = u[4]/u[1]
        uout[5] = u[4]/u[1] - qe[4]/qe[1]
        uout[6] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1], θ=u[4]/u[1])
        uout[7] = u[2]/u[1]              # u_total = u (same as uout[2] in TOTAL)
        uout[8] = u[3]/u[1]              # v_total

    elseif ET == PERT()
        uout[1] = u[1] + qe[1]
        uout[2] = u[2]/(u[1] + qe[1])                               # u_pert
        uout[3] = u[3]/(u[1] + qe[1])                               # v_pert
        uout[4] = (u[4] + qe[4])/(u[1] + qe[1])                     # θ total
        uout[5] = (u[4] + qe[4])/(u[1] + qe[1]) - qe[4]/qe[1]      # θ_prime
        uout[6] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1]+qe[1], θ=(u[4]+qe[4])/(u[1]+qe[1]))
        uout[7] = u[2]/(u[1] + qe[1]) + qe[2]/qe[1]                # u_total = u_pert + u_bg
        uout[8] = u[3]/(u[1] + qe[1]) + qe[3]/qe[1]                # v_total = v_pert + v_bg
    end
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
    uprimitive[6] = u[6]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = (u[4]+qe[4])/(u[1]+qe[1])-qe[4]/qe[1]
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
    uprimitive[6] = (u[6]+qe[6])/(u[1]+qe[1])-qe[6]/qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T(u[4]/(u[1]+qe[1])), T((u[5]+qe[5])/(u[1]+qe[1]) - qe[5]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1]), T(u[5]/u[1])
    end
end

function user_uout!(ip, ET, uout, u, ue; mp=mp)

    uout[1] = u[1]
    ρ = uout[1] + ue[1]
    uout[2] = u[2]/ρ
    uout[3] = u[3]/ρ
    uout[4] = u[4]/ρ
    uout[5] = u[5]/ρ
    uout[6] = u[6]/ρ

    if length(uout[:]) > 6
        uout[7]  = mp.Tabs[ip]
        uout[8]  = mp.qn[ip]
        uout[9]  = mp.qc[ip]
        uout[10] = mp.qi[ip]
        uout[11] = mp.qr[ip]
        uout[12] = mp.qs[ip]
        uout[13] = mp.qg[ip]
    end
end

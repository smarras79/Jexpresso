function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::TOTAL)
    uprimitive[1] = u[1]
end

function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::PERT)
    uprimitive[1] = u[1]+qe[1]
end

#=function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end=#

function user_uout!(ip, ET, uout, u, qe; mp=mp, F_data=F_data)

    uout[1] = u[1] #ω
    uout[2] = -F_data[ip,2] #u
    uout[3] = F_data[ip,1] #v
end

function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]
end

function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[1]+qe[2]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return  T(u[1]+qe[1]), T(u[2]+qe[2])
    else
        return T(u[1]),T(u[2])
    end
end

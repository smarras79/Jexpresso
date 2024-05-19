function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::TOTAL)
    uprimitive[1] = u[1]
end

function user_primitives!(u::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::PERT)
    uprimitive[1] = u[1]+qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    return T(u[1])
end

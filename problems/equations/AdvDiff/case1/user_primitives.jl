function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1]+qe[1])
end

function user_uout!(uout, u, qe, ::TOTAL)

    uout[1] = u[1]
end

function user_uout!(uout, u, qe, ::PERT)

    uout[1] = u[1] + qe[1]
    
end

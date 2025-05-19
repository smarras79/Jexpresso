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

function user_uout!(ip, ET, uout, u, qe...)
    uout[1] = u[1] #+qe[1]
end

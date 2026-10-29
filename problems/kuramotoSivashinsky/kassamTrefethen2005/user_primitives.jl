#---------------------------------------------------------------------------------
# One unknown, u. No viscous operator is used (the derivatives are the split
# linear operator), so user_primitives! is a straight copy.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[1] + qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return lpert ? T(u[1] + qe[1]) : T(u[1])
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    uout[1] = u[1]
end

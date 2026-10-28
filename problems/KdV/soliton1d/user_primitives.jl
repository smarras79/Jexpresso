#---------------------------------------------------------------------------------
# The state is already the working set; no viscous operator is used, so
# user_primitives! is a straight copy.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]
    uprimitive[3] = u[3]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[1] + qe[1]
    uprimitive[2] = u[2] + qe[2]
    uprimitive[3] = u[3] + qe[3]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1] + qe[1]), T(u[2] + qe[2]), T(u[3] + qe[3])
    else
        return T(u[1]), T(u[2]), T(u[3])
    end
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

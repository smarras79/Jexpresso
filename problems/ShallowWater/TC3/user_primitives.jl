function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]  # η
    uprimitive[2] = u[2]  # u
    uprimitive[3] = u[3]  # v
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]
    uprimitive[3] = u[3]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1]), T(u[2]), T(u[3])
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]  # η
    uout[2] = u[2]  # u
    uout[3] = u[3]  # v
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

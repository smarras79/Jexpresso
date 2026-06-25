function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[1] + qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1] + qe[1])
end

# NOTE: the second argument is typed (::TOTAL / ::PERT) rather than left
# generic. Jexpresso loads every case's user_primitives.jl into the same
# `Jexpresso` module, and a previously-loaded (or precompiled) multi-variable
# case such as CompEuler/theta or CompEuler/sod1d also defines
# `user_uout!(ip, ::TOTAL, ...)`. A generic `user_uout!(ip, ET, ...)` does NOT
# replace that more specific method, so dispatch would pick the multi-variable
# version and index uout[2..] out of bounds on this 1-variable problem.
# Matching the typed signature makes this definition replace the resident one.
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    uout[1] = u[1] + qe[1]
end

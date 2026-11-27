function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = (u[4]+qe[4])/(u[1]+qe[1])-qe[4]/qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end

# Output the same quantities as the PERT deck (theta_laguerre), so the two can be
# compared field by field: ρ', u (+10, as there), v, θ'. With TOTAL variables u
# holds the total state.
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1] - qe[1]
    uout[2] = u[2]/u[1] + 10
    uout[3] = u[3]/u[1]
    uout[4] = u[4]/u[1] - qe[4]/qe[1]
end

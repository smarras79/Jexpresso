function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    uout[1] = u[1] + qe[1]
end

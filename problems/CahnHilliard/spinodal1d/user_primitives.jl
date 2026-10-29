#---------------------------------------------------------------------------------
# uprimitive[ieq] is the variable whose second derivative the viscous operator
# adds to equation ieq, scaled by :μ[ieq] (user_inputs.jl). The two are
# CROSSED, which is what makes the system Cahn-Hilliard:
#
#     equation 1 (c):   +1     · ∂ₓ²μ        uprimitive[1] = μ
#     equation 2 (μ):   -ε²/τ  · ∂ₓ²c        uprimitive[2] = c
#
# user_uout! is what reaches the output files, and there the state is written
# as it is.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[2]
    uprimitive[2] = u[1]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[2] + qe[2]
    uprimitive[2] = u[1] + qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[2] + qe[2]), T(u[1] + qe[1])
    else
        return T(u[2]), T(u[1])
    end
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
end

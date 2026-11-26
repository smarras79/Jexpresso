#---------------------------------------------------------------------------------
# TC6: state q = [h, hu, hv]. Primitives are [h, u, v]; the output variables are
# [η, u, v] with η = h - H, so TC6's VTK fields carry the same names as TC1-TC4.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2] / u[1]
    uprimitive[3] = u[3] / u[1]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" # ERROR TC6: only SOL_VARS_TYPE => TOTAL() is implemented.")
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1] - TC6_H       # η
    uout[2] = u[2] / u[1]        # u
    uout[3] = u[3] / u[1]        # v
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    error(" # ERROR TC6: only SOL_VARS_TYPE => TOTAL() is implemented.")
end

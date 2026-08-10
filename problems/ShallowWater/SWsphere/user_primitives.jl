#---------------------------------------------------------------------------------
# SWsphere — primitive variables and output variables.
#
# Solution (conservative) state, Marras/Kopera/Giraldo Eq. (8):
#
#   q    = [φ, φu, φv, φw]      φ = gh, (u,v,w) the Cartesian velocity
#
# Primitives — what the viscous operator δν∇²(φu) differentiates, and what the
# artificial-viscosity machinery expects to see:
#
#   uprim = [φ, u, v, w]
#
# Output — what a human reads off the VTK file:
#
#   uout  = [h, u, v, w]
#
# g is a literal here for the same reason ω is one in user_source.jl: these
# callbacks receive no `inputs`. Keep it in step with :galewsky_gravity if the
# deck ever changes it.
#---------------------------------------------------------------------------------

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    φ = u[1]
    uprimitive[1] = φ
    uprimitive[2] = u[2]/φ
    uprimitive[3] = u[3]/φ
    uprimitive[4] = u[4]/φ
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    user_primitives!(u, qe, uprimitive, TOTAL())
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    φ = u[1]
    return T(φ), T(u[2]/φ), T(u[3]/φ), T(u[4]/φ)
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    g = 9.80616
    φ = u[1]
    uout[1] = φ/g        # h  [m]
    uout[2] = u[2]/φ     # u  [m/s]
    uout[3] = u[3]/φ     # v
    uout[4] = u[4]/φ     # w
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end

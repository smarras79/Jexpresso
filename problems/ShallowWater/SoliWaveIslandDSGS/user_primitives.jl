#
# Primitive helpers for the non-linear shallow water solver.
#
# uprimitive is the quantity differentiated by the viscous Laplacian.
# The continuity equation diffuses the depth PERTURBATION H - He
# (He = qe[1] is the lake-at-rest depth defined in initialize.jl): at
# rest the depth itself is cone-shaped (∇²H ≠ 0 over the island), so
# diffusing the full H would pump mass in/out of the wet/dry ring and
# break the discrete equilibrium that user_flux.jl/user_source.jl
# preserve. The momentum components are diffused in conservation form;
# their reference state is zero, so no subtraction is needed.
#
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1] - qe[1]   # H - He (free-surface perturbation)
    uprimitive[2] = u[2]           # Hu
    uprimitive[3] = u[3]           # Hv
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    user_primitives!(u, qe, uprimitive, TOTAL())
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1] - qe[1]), T(u[2]), T(u[3])
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2]
    uout[3] = u[3]
end

#
# Primitive helpers for the non-linear shallow water solver.
#
# uprimitive is the quantity differentiated by the viscous Laplacian.
# Using (H, Hu, Hv) keeps the diffusion in conservation form, which is
# consistent with the residual-based viscosity discussed in the paper.
#
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]   # H
    uprimitive[2] = u[2]   # Hu
    uprimitive[3] = u[3]   # Hv
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    user_primitives!(u, qe, uprimitive, TOTAL())
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1]), T(u[2]), T(u[3])
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

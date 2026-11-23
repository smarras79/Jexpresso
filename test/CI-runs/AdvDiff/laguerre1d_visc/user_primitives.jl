# qe holds the exact solution at tend (comparison only), so it is not added.
function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]
end
function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    return T(u[1]), T(0)
end
function user_uout!(ip, ET, uout, u, qe; kwargs...)
    uout[1] = u[1]
end

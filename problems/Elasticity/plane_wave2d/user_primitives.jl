#---------------------------------------------------------------------------------
# The state is already the working set, so user_primitives! (which feeds the
# viscous operator) is a straight copy. user_uout! is what reaches VTK/HDF5,
# and there the momenta are divided by ρ so the file carries the particle
# velocity — matching the qoutvars declared in initialize.jl.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    for ieq = 1:5
        uprimitive[ieq] = u[ieq]
    end
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    for ieq = 1:5
        uprimitive[ieq] = u[ieq] + qe[ieq]
    end
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    return T(u[1]), T(u[2]), T(u[3]), T(u[4]), T(u[5])
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    P = elastic_properties()
    uout[1] = u[1]/P.ρ     # vx
    uout[2] = u[2]/P.ρ     # vy
    uout[3] = u[3]         # σxx
    uout[4] = u[4]         # σyy
    uout[5] = u[5]         # σxy
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    P = elastic_properties()
    uout[1] = (u[1] + qe[1])/P.ρ
    uout[2] = (u[2] + qe[2])/P.ρ
    uout[3] =  u[3] + qe[3]
    uout[4] =  u[4] + qe[4]
    uout[5] =  u[5] + qe[5]
end

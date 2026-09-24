#---------------------------------------------------------------------------------
# The state is already the working set of the discretisation, so
# user_primitives! (which feeds the viscous operator) is a straight copy.
#
# user_uout! is what reaches HDF5/VTK, and there the momentum densities are
# divided through by their inertia so the file carries v and ω rather than
# ρA v and ρI ω — matching the qoutvars declared in initialize.jl.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]
    uprimitive[3] = u[3]
    uprimitive[4] = u[4]
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    uprimitive[1] = u[1] + qe[1]
    uprimitive[2] = u[2] + qe[2]
    uprimitive[3] = u[3] + qe[3]
    uprimitive[4] = u[4] + qe[4]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1] + qe[1]), T(u[2] + qe[2]), T(u[3] + qe[3]), T(u[4] + qe[4])
    else
        return T(u[1]), T(u[2]), T(u[3]), T(u[4])
    end
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    p = timoshenko_properties()

    uout[1] = u[1]/p.ρA     # v
    uout[2] = u[2]/p.ρI     # ω
    uout[3] = u[3]          # γ
    uout[4] = u[4]          # χ
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    p = timoshenko_properties()

    uout[1] = (u[1] + qe[1])/p.ρA
    uout[2] = (u[2] + qe[2])/p.ρI
    uout[3] =  u[3] + qe[3]
    uout[4] =  u[4] + qe[4]
end

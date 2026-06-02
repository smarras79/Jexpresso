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

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end

# The 6th output (μ_dsgs in initialize.jl) is the Marras DSGS coefficient.
# It is fed through call_user_uout via the μ_dsgs_pnode kwarg — a length-
# npoin Vector that broadcast_dsgs_to_nodes! filled with the per-element
# value of params.μ_dsgs. If the kwarg is absent (e.g. running without
# DSGS) the slot is left as 0.
function user_uout!(ip, ::TOTAL, uout, u, qe; mp=nothing, μ_dsgs_pnode=nothing, kwargs...)

    uout[1] = u[1]
    uout[2] = u[2]/u[1]
    uout[3] = u[3]/u[1]
    uout[4] = u[4]/u[1]

    PhysConst = PhysicalConst{Float64}()
    uout[5] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1], θ=u[4]/u[1])

    if length(uout) >= 6 && μ_dsgs_pnode !== nothing
        uout[6] = μ_dsgs_pnode[ip]
    end
end

function user_uout!(ip, ::PERT, uout, u, qe; mp=nothing, μ_dsgs_pnode=nothing, kwargs...)

    uout[1] = u[1]+qe[1]
    uout[2] = (u[2]+qe[2])/(u[1]+qe[1])
    uout[3] = (u[3]+qe[3])/(u[1]+qe[1])
    uout[4] = (u[4]+qe[4])/(u[1]+qe[1])

    PhysConst = PhysicalConst{Float64}()
    uout[5] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1]+qe[1], θ=(u[4]+qe[4])/(u[1]+qe[1]))

    if length(uout) >= 6 && μ_dsgs_pnode !== nothing
        uout[6] = μ_dsgs_pnode[ip]
    end
end

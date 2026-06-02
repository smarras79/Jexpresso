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

# Output slots 6..8 are the per-equation Marras DSGS coefficients
# pulled from params.μ_dsgs_pnode[ip, :], a npoin × neqs matrix that
# broadcast_dsgs_to_nodes! fills from the per-element coefficients.
# Layout:
#   uout[6] = μ on the x-momentum equation  ( params.μ_dsgs_pnode[ip,2] )
#   uout[7] = μ on the y-momentum equation  ( params.μ_dsgs_pnode[ip,3] )
#   uout[8] = κ on the θ equation, already scaled by Pr/(γ-1)
#                                           ( params.μ_dsgs_pnode[ip,4] )
function user_uout!(ip, ::TOTAL, uout, u, qe; mp=nothing, μ_dsgs_pnode=nothing, kwargs...)

    uout[1] = u[1]
    uout[2] = u[2]/u[1]
    uout[3] = u[3]/u[1]
    uout[4] = u[4]/u[1]

    PhysConst = PhysicalConst{Float64}()
    uout[5] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1], θ=u[4]/u[1])

    if length(uout) >= 8 && μ_dsgs_pnode !== nothing && size(μ_dsgs_pnode, 2) >= 4
        uout[6] = μ_dsgs_pnode[ip, 2]
        uout[7] = μ_dsgs_pnode[ip, 3]
        uout[8] = μ_dsgs_pnode[ip, 4]
    end
end

function user_uout!(ip, ::PERT, uout, u, qe; mp=nothing, μ_dsgs_pnode=nothing, kwargs...)

    uout[1] = u[1]+qe[1]
    uout[2] = (u[2]+qe[2])/(u[1]+qe[1])
    uout[3] = (u[3]+qe[3])/(u[1]+qe[1])
    uout[4] = (u[4]+qe[4])/(u[1]+qe[1])

    PhysConst = PhysicalConst{Float64}()
    uout[5] = perfectGasLaw_ρθtoP(PhysConst, ρ=u[1]+qe[1], θ=(u[4]+qe[4])/(u[1]+qe[1]))

    if length(uout) >= 8 && μ_dsgs_pnode !== nothing && size(μ_dsgs_pnode, 2) >= 4
        uout[6] = μ_dsgs_pnode[ip, 2]
        uout[7] = μ_dsgs_pnode[ip, 3]
        uout[8] = μ_dsgs_pnode[ip, 4]
    end
end

#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping used by the viscous/SGS kernels.
#
# uprimitive[ieq] is the variable whose Laplacian is added to equation ieq by
# _expansion_visc! (per-equation gradient diffusion), so the ordering must
# match the equation ordering (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ):
#
#   1: ρ    (visc multiplier :μ[1] = 0 -> no mass diffusion)
#   2: u    (momentum: full deviatoric stress built by the kernel)
#   3: v    (momentum: full deviatoric stress built by the kernel)
#   4: T    (energy slot: κ∇T heat flux + τ·u viscous work added by the kernel)
#   5: w    (out-of-plane momentum, diffused as a scalar)
#   6-8: B  (magnetic field, scalar diffusion = turbulent resistivity)
#   9: ψ    (GLM field, scalar diffusion)
#
# T = p/ρ is the nondimensional temperature of this nondimensional setup.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)

    ρ = u[1]
    p = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])

    uprimitive[1] = ρ
    uprimitive[2] = u[2]/ρ      # u
    uprimitive[3] = u[3]/ρ      # v
    uprimitive[4] = p/ρ         # T
    uprimitive[5] = u[5]/ρ      # w
    uprimitive[6] = u[6]        # Bx
    uprimitive[7] = u[7]        # By
    uprimitive[8] = u[8]        # Bz
    uprimitive[9] = u[9]        # ψ
end

function user_primitives(u, qe, uprimitive, ::TOTAL)

    ρ = u[1]
    p = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])

    return SVector(ρ, u[2]/ρ, u[3]/ρ, p/ρ, u[5]/ρ, u[6], u[7], u[8], u[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/kelvinHelmholtzChan2022: PERT() solution variables are not supported.")
end

#---------------------------------------------------------------------------------
# Output variables: qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T"]
#---------------------------------------------------------------------------------
function user_uout!(ip, ET, uout, u, qe; kwargs...)

    ρ = u[1]
    p = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])

    uout[1]  = ρ
    uout[2]  = u[2]/ρ   # u
    uout[3]  = u[3]/ρ   # v
    uout[4]  = u[5]/ρ   # w
    uout[5]  = p
    uout[6]  = u[6]     # Bx
    uout[7]  = u[7]     # By
    uout[8]  = u[8]     # Bz
    uout[9]  = u[9]     # ψ
    uout[10] = p/ρ      # T
end

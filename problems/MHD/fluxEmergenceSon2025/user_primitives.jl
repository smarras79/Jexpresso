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
# T = p/ρ is the nondimensional temperature in the units of the paper's
# Table 1 up to the factor γ: T/T₀ = γ p/ρ. The kernels only ever take
# gradients of it, so the factor is immaterial there; user_uout! below
# writes the properly normalized T/T₀.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)

    ρ = u[1]
    p = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])

    uprimitive[1] = ρ
    uprimitive[2] = u[2]/ρ      # u
    uprimitive[3] = u[3]/ρ      # v
    uprimitive[4] = p/ρ         # T (up to γ)
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
    error(" problems/MHD/fluxEmergenceSon2025: PERT() solution variables are not supported.")
end

#---------------------------------------------------------------------------------
# Output variables:
#
#   qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "vA", "β"]
#
# All in the paper's units (ρ/ρ₀, V/C_s, p/p₀, B/B₀, T/T₀):
#   T  = γ p/ρ                     temperature T/T₀ (Table 1: T₀ = m C_s²/(γ k_B))
#   vA = |B|/√ρ                    Alfvén speed V_A/C_s   (paper Fig. 5, row 2)
#   β  = 2p/|B|²                   plasma beta            (paper Fig. 6(e))
# β is capped at 1e6 where B vanishes (the unmagnetized corona at t = 0) so
# that it stays plottable.
#---------------------------------------------------------------------------------
function user_uout!(ip, ET, uout, u, qe; kwargs...)

    ρ  = u[1]
    p  = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    B2 = u[6]*u[6] + u[7]*u[7] + u[8]*u[8]

    uout[1]  = ρ
    uout[2]  = u[2]/ρ   # u
    uout[3]  = u[3]/ρ   # v  (the paper's V_z)
    uout[4]  = u[5]/ρ   # w
    uout[5]  = p
    uout[6]  = u[6]     # Bx
    uout[7]  = u[7]     # By (the paper's B_z)
    uout[8]  = u[8]     # Bz
    uout[9]  = u[9]     # ψ
    uout[10] = γ_mhd*p/ρ             # T/T₀
    uout[11] = sqrt(B2/ρ)            # V_A/C_s
    uout[12] = min(2.0*p/max(B2, 1e-300), 1.0e6)   # plasma β
end

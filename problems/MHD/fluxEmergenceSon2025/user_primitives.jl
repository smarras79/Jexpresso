#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping used by the viscous/SGS kernels.
#
# uprimitive[ieq] is the variable whose Laplacian is added to equation ieq by
# _expansion_visc!. This case runs the DynSGS operator in its
# CONSERVED-VARIABLE form (:dsgs_conserved => true in user_inputs.jl, see
# kernel/physics/SGS.jl): every slot is handed the conserved variable itself
# and receives the same kinematic coefficient, so the dissipation is a
# Laplacian on (ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ).
#
# Why not the physical form (u, v, T with the dynamic ρμ, κ∇T, τ·u) of the
# Orszag-Tang case: the 25× density drop of the chromosphere-corona transition
# is a contact (p continuous, ρ and T jumping) that this mesh under-resolves,
# and keeping it positive needs mass diffusion. Diffusing ρ while the energy
# closure acts on T is thermodynamically inconsistent — mass added to the
# light side carries no energy, p = (γ-1)(E - ½ρv² - ½B²) with γ-1 = 0.05
# went negative within a few τ₀ (measured). With the same Laplacian on ρ and
# on E, E is constant across the isobaric contact and does not diffuse, ρ
# spreads, and p stays what it was. The kinetic energy removed by the ρv
# Laplacian reappears as internal energy because E itself is untouched, the
# right sign; no separate τ·u term is needed.
#
# ... applied to the DEPARTURE from the magnetostatic reference state q_e of
# initialize.jl, ∇·(μ∇(q − q_e)), not to q itself. In a stratified
# atmosphere the conserved variables carry the stratification: ∇²ρ_e = ρ_e/H²
# is not small, and a Laplacian on ρ is a steady upward mass source μρ/H²
# that the residual sensor then feeds on — measured: μ at 70-90% of its cap
# through the whole sheet and the sheet sinking at 0.3 C_s by t = 8 τ₀. On
# q − q_e the operator vanishes at rest (nothing is diffused, the initial
# sheet field is not eroded), acts on a displaced transition-region contact
# exactly as needed for positivity, and on the emerged loop — where q_e is
# negligible — reduces to the plain Laplacian on the conserved variables.
# It stays in divergence form, i.e. conservative. (The reference jump at the
# fixed z_cor contributes −μ∇²ρ_e there once the sensor fires nearby; with
# ρ_e(z_cor) ≈ 10⁻⁷ that is far below the loop densities that then occupy
# the region.)
#---------------------------------------------------------------------------------
# Energy slot: with :dsgs_nazarov_energy the solver sets dsgs_split_energy[]
# and the energy flux is split — slot 4 = δ(E − p/(γ−1)), the non-thermal
# departure, diffused with the ρv slot's ν; slot 11 (= neqs+2) = δ(p/(γ−1)),
# the thermal one, diffused with Dao & Nazarov's κ = ρν/Pr (floor kept).
# See problems/MHD/fluxEmergenceSon2025DSGS/user_primitives.jl.
@inline function fe_energy_split(u)
    p = pressure_mhd(u[1], u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    eth = p/(γ_mhd - 1.0)
    return u[4] - eth, eth
end

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    for ieq = 1:9
        uprimitive[ieq] = u[ieq] - qe[ieq]    # ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ minus the reference state
    end
    if dsgs_split_energy[]
        nth, eth   = fe_energy_split(u)
        nthe, ethe = fe_energy_split(qe)
        uprimitive[4]  = nth - nthe
        uprimitive[11] = eth - ethe
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    e4 = dsgs_split_energy[] ? fe_energy_split(u)[1] - fe_energy_split(qe)[1] : u[4] - qe[4]
    return SVector(u[1] - qe[1], u[2] - qe[2], u[3] - qe[3], e4, u[5] - qe[5],
                   u[6] - qe[6], u[7] - qe[7], u[8] - qe[8], u[9] - qe[9])
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
    uout[11] = sqrt(max(B2/ρ, 0.0))  # V_A/C_s (guarded: a run that went bad must still write its last output)
    uout[12] = min(2.0*p/max(B2, 1e-300), 1.0e6)   # plasma β
end

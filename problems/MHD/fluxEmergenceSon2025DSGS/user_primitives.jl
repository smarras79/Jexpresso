#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping used by the DynSGS-MHD viscous kernel.
#
# uprimitive[ieq] is the variable whose gradient the dissipation of equation
# ieq acts on (kernel/operators/rhs.jl, _expansion_visc!). This case runs the
# operator in the CONSERVED form of the sibling case (:dsgs_conserved) and,
# in addition, in its REFERENCE-WEIGHTED form (:dsgs_ref_weight, see
# kernel/physics/SGS.jl): for the fluid slots
#
#     uprimitive[1:5]  = (q − q_e)/ρ_e,        coefficient  μ·ρ_e,
#
# i.e. the dissipation is  ∇·(μ ρ_e ∇((q − q_e)/ρ_e))  — the weighted
# diffusion of the RELATIVE departure from the magnetostatic reference
# state q_e of initialize.jl, with the reference density ρ_e(z) handed to
# the kernel in the spare slot neqs+1 of uprimitive. For the magnetic slots
#
#     uprimitive[6:9]  = q − q_e,               coefficient  μ,
#
# as in the sibling (B_e is a sheet field, not a density-like weight).
#
# Why this form. Both ∇·(μ∇(q − q_e)) (the sibling) and this one vanish at
# rest and are conservative. They differ in what they do across the 25×
# reference jump of the chromosphere–corona transition once the sensor
# fires there. The absolute form diffuses ρ − ρ_e: a departure of −10% on
# the dense side is −2×10⁻⁹ in absolute terms, more than the entire density
# of the light side (8×10⁻¹⁰), so diffusing it across the contact drives the
# light side negative — that is the mechanism by which the sibling's run,
# with μ already at the wave-speed cap, still evacuated the coronal foot of
# the transition region at t ≈ 14 τ₀ and needed its floors. The relative
# form diffuses r = ρ/ρ_e with the positive weight ρ_e: r is pulled toward
# its neighbours' values (0.9 in the example), never past them, so ρ = r ρ_e
# stays positive wherever the neighbours are — the discrete maximum
# principle of a weighted Laplacian (up to the non-M-matrix entries of the
# high-order stiffness matrix, which is the same caveat as for any SEM
# artificial viscosity). Momentum and energy take the same weight so that
# the fluid slots are diffused as one scaled state (ρ, ρv, E − E_e)/ρ_e and
# the pressure of a diffused node stays a convex combination of admissible
# ones; the kinetic energy removed from ρv reappears as internal energy
# because E is diffused, not p (as in the sibling).
#
# In the emerged loop, where q_e is negligible against q (ρ_e ≈ 10⁻⁸ in the
# corona against 10⁻⁵ in the loop), (q − q_e)/ρ_e ≈ q/ρ_e and the operator
# is a Laplacian on q with the coefficient μ, i.e. exactly the sibling's
# conserved-variable Laplacian, so the shock capturing in the loop is the
# same DynSGS as before.
#---------------------------------------------------------------------------------
# Energy slot. E = ½ρ|v|² + ½|B|² (+ ½ψ²) + p/(γ−1), and the Laplacian on E
# contains ∇·(νB·∇B) and the kinetic analogue: exactly the conservative
# energy fluxes of the ν∇B and ν∇(ρv) Laplacians of the other slots, which
# is what makes the single-ν conserved form energy-consistent (the Joule and
# viscous heating are what is left). With :dsgs_nazarov_energy the solver
# sets dsgs_split_energy[] (kernel/physics/SGS.jl) and the energy flux is
# split: slot 4 gets the NON-THERMAL departure, diffused by the kernel with
# the ρv slot's ν, and the spare slot 11 (= neqs+2) the THERMAL one
# δ(p/(γ−1)), diffused with the E-slot coefficient max(γ(γ−1)/Pr·ν_res,
# ν_floor), i.e. Dao & Nazarov's κ = ρν/Pr with the background floor kept.
# Without it slot 4 is δE and slot 11 is unused. (Scaling the whole slot
# instead — the first implementation — let ν∇B spread the sheet's field
# while its magnetic energy stayed put, and cut the floor that damps the
# node-to-node mode: the sheet core overheated, a temperature sawtooth grew
# across the corona and the emergence stalled.)
@inline function fe_energy_split(u)
    p = pressure_mhd(u[1], u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    eth = p/(γ_mhd - 1.0)
    return u[4] - eth, eth      # (non-thermal, thermal) energy
end

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    ρe  = qe[1]
    iρe = 1.0/ρe
    uprimitive[1] = (u[1] - qe[1])*iρe            # ρ  : relative departure
    uprimitive[2] = (u[2] - qe[2])*iρe            # ρu
    uprimitive[3] = (u[3] - qe[3])*iρe            # ρv
    if dsgs_split_energy[]
        nth, eth   = fe_energy_split(u)
        nthe, ethe = fe_energy_split(qe)
        uprimitive[4]  = (nth - nthe)*iρe         # E: non-thermal part (coefficient ν)
        uprimitive[11] = (eth - ethe)*iρe         #    thermal part (Nazarov's κ)
    else
        uprimitive[4]  = (u[4] - qe[4])*iρe       # E
    end
    uprimitive[5] = (u[5] - qe[5])*iρe            # ρw
    for ieq = 6:9
        uprimitive[ieq] = u[ieq] - qe[ieq]         # Bx, By, Bz, ψ: absolute departure
    end
    uprimitive[10] = ρe                            # the weight, read by _expansion_visc!
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    iρe = 1.0/qe[1]
    e4  = dsgs_split_energy[] ? (fe_energy_split(u)[1] - fe_energy_split(qe)[1])*iρe : (u[4] - qe[4])*iρe
    return SVector((u[1] - qe[1])*iρe, (u[2] - qe[2])*iρe, (u[3] - qe[3])*iρe, e4, (u[5] - qe[5])*iρe,
                   u[6] - qe[6], u[7] - qe[7], u[8] - qe[8], u[9] - qe[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/fluxEmergenceSon2025DSGS: PERT() solution variables are not supported.")
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

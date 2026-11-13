#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping of the DynSGS-MHD viscous kernel: the same
# reference-weighted, conserved-form operator as
# problems/MHD/fluxEmergenceSon2025DSGS (see its user_primitives.jl and README),
#
#     ∂_t q_i = … + ∇·(μ ρ_e ∇((q_i − q_e,i)/ρ_e)),    i = ρ, ρu, ρv, E, ρw,
#     ∂_t q_i = … + ∇·(μ ∇(q_i − q_e,i)),              i = Bx, By, Bz, ψ,
#
# zero at rest, conservative, and a maximum principle on ρ/ρ_e across the
# 10⁸-fold reference density drop between the photosphere and the corona. q_e is
# the magnetohydrostatic reference of initialize.jl (atmosphere + potential
# coronal field). ρ_e goes to the kernel in the spare slot neqs+1 of uprimitive.
# The explicit ν and η of the paper (:mhd_dissipation, :mhd_anomalous_eta) are
# added by the kernel to the coefficient of the ρu, ρv, ρw and B slots.
#
# Energy slot: E as a whole (:dsgs_nazarov_energy => false in the deck). At
# γ = 5/3 the split of the energy flux of the parent case changes the thermal
# coefficient by γ(γ−1)/Pr = 1.11 only, and it is incompatible with the explicit
# ν, η (it diffuses the non-thermal energy with the ρv slot's coefficient). The
# split path is kept for completeness.
#---------------------------------------------------------------------------------
@inline function lk10_energy_split(u)
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
        nth, eth   = lk10_energy_split(u)
        nthe, ethe = lk10_energy_split(qe)
        uprimitive[4]  = (nth - nthe)*iρe
        uprimitive[11] = (eth - ethe)*iρe
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
    e4  = dsgs_split_energy[] ? (lk10_energy_split(u)[1] - lk10_energy_split(qe)[1])*iρe : (u[4] - qe[4])*iρe
    return SVector((u[1] - qe[1])*iρe, (u[2] - qe[2])*iρe, (u[3] - qe[3])*iρe, e4, (u[5] - qe[5])*iρe,
                   u[6] - qe[6], u[7] - qe[7], u[8] - qe[8], u[9] - qe[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/fluxEmergenceLeake2010: PERT() solution variables are not supported.")
end

#---------------------------------------------------------------------------------
# Output variables (initialize.jl):
#
#   ρ, u, v, w, p, Bx, By, Bz, ψ   code units (ρ_ph, v0, p_ph, B0 = 1030.6 G)
#   T             p/ρ, in T_ph = 5700 K
#   β             2p/|B|² (capped at 1e8 where B vanishes)
#   log10rho_cgs  log₁₀ of ρ in g cm⁻³: the colour scale of the paper's density figures
#   Bz_G          the shear field B_z in G: the colour scale of its B_z figures
#---------------------------------------------------------------------------------
function user_uout!(ip, ET, uout, u, qe; kwargs...)
    ρ  = u[1]
    p  = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    B2 = u[6]*u[6] + u[7]*u[7] + u[8]*u[8]
    uout[1]  = ρ
    uout[2]  = u[2]/ρ
    uout[3]  = u[3]/ρ
    uout[4]  = u[5]/ρ
    uout[5]  = p
    uout[6]  = u[6]
    uout[7]  = u[7]
    uout[8]  = u[8]
    uout[9]  = u[9]
    uout[10] = p/ρ
    uout[11] = min(2.0*p/max(B2, 1e-300), 1.0e8)
    uout[12] = log10(max(ρ*LK10_RHO0, 1e-300))
    uout[13] = u[8]*LK10_BU
end

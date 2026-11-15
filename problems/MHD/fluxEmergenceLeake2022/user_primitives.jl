#---------------------------------------------------------------------------------
# DynSGS in conserved form (:dsgs_conserved): the Laplacian acts on the departure
# of the conserved variables from the reference state qe (atmosphere + dipole),
# so the stratification and the potential dipole are not diffused, as in
# problems/MHD/fluxEmergenceSon2025 and fluxEmergenceLeake2014.
#---------------------------------------------------------------------------------
#
# With :dsgs_ref_weight (kernel/physics/SGS.jl, rhs.jl) the slots ρ, ρu, ρv, E, ρw are
# diffused as the RELATIVE departure (q − q_e)/ρ_e with the coefficient μ ρ_e (the weight
# goes in the spare slot 10): ∇·(μ ρ_e ∇((q − q_e)/ρ_e)) instead of ∇·(μ ∇(q − q_e)). Both
# vanish at rest; only the relative form obeys a maximum principle across the 10³-10⁷ density
# drop of the chromosphere and transition region, where a departure of the dense side,
# diffused upward, can exceed the whole local state (README). B and ψ keep the absolute form.
#
# With :dsgs_nazarov_energy (the split energy flux, SGS.jl dsgs_split_energy) slot 4 holds the
# NON-THERMAL departure δ(½ρ|v|² + ½|B|² + ½ψ²), diffused with the ρu slot's coefficient so that
# the kinetic and magnetic energy fluxes of E are those of the ρv and B Laplacians, and the
# spare slot 11 the THERMAL one δ(p/(γ−1)), diffused with slot 4's coefficient (both relative
# to ρ_e with :dsgs_ref_weight). Without the split, E with the weighting diffuses the magnetic
# energy differently from B, and in the low-β transition region the pressure, the small
# remainder E − ½ρ|v|² − ½|B|², takes the difference.
@inline function lk22_energy_split(u)
    eth = pressure_mhd(u[1], u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])/(γ_mhd - 1.0)
    return u[4] - eth, eth      # (non-thermal, thermal) energy
end

function user_primitives!(u, qe, uprimitive, ::TOTAL)
    w = dsgs_ref_weight[] ? 1.0/qe[1] : 1.0
    for ieq = 1:5
        uprimitive[ieq] = (u[ieq] - qe[ieq])*w
    end
    for ieq = 6:9
        uprimitive[ieq] = u[ieq] - qe[ieq]
    end
    if dsgs_split_energy[]
        nth, eth   = lk22_energy_split(u)
        nthe, ethe = lk22_energy_split(qe)
        uprimitive[4]  = (nth - nthe)*w   # E: non-thermal part (the ρu slot's coefficient)
        uprimitive[11] = (eth - ethe)*w   #    thermal part (slot 4's coefficient)
    end
    dsgs_ref_weight[] && (uprimitive[10] = qe[1])   # the weight, read by _expansion_visc!
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    w  = dsgs_ref_weight[] ? 1.0/qe[1] : 1.0
    e4 = dsgs_split_energy[] ? (lk22_energy_split(u)[1] - lk22_energy_split(qe)[1])*w : (u[4] - qe[4])*w
    return SVector((u[1] - qe[1])*w, (u[2] - qe[2])*w, (u[3] - qe[3])*w, e4, (u[5] - qe[5])*w,
                   u[6] - qe[6], u[7] - qe[7], u[8] - qe[8], u[9] - qe[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/fluxEmergenceLeake2022: PERT() solution variables are not supported.")
end

# Output: ρ, u, v, w, p, Bx, By, Bz, ψ, T = m_f p/ρ (in T₀; m_f = :lk22_mf), plasma β = 2p/|B|²,
# |B|, log10 ρ (the atmosphere spans twelve decades).
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
    uout[10] = lk22_mf[]*p/ρ
    uout[11] = min(2.0*p/max(B2, 1e-300), 1.0e12)
    uout[12] = sqrt(B2)
    uout[13] = log10(max(ρ, 1e-300))
end

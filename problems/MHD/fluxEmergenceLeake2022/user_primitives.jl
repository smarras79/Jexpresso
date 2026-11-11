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
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    if dsgs_ref_weight[]
        iρe = 1.0/qe[1]
        for ieq = 1:5
            uprimitive[ieq] = (u[ieq] - qe[ieq])*iρe
        end
        for ieq = 6:9
            uprimitive[ieq] = u[ieq] - qe[ieq]
        end
        uprimitive[10] = qe[1]            # the weight, read by _expansion_visc!
    else
        for ieq = 1:9
            uprimitive[ieq] = u[ieq] - qe[ieq]
        end
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    w = dsgs_ref_weight[] ? 1.0/qe[1] : 1.0
    return SVector((u[1] - qe[1])*w, (u[2] - qe[2])*w, (u[3] - qe[3])*w, (u[4] - qe[4])*w, (u[5] - qe[5])*w,
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

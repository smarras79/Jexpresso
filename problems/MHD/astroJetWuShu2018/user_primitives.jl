#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping used by the DynSGS-MHD viscous kernel.
#
# uprimitive[ieq] is the variable whose gradient the dissipation of equation ieq
# acts on (kernel/operators/rhs.jl, _expansion_visc!). This case runs DynSGS in
# its CONSERVED form (:dsgs_conserved => true), exactly as the MHD shock tube
# problems/MHD/brioWu1d does, so the primitives ARE the conserved variables and
# the single kinematic ν from the residual is applied to every slot:
#
#     ∂_t q + ∇·F(q) = ∇·(ν ∇q)   on (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
#
# This is the Lax-Friedrichs-type viscous regularization that the residual
# method reduces to at its cap μ_max = C_max Δ (|v| + c_f), written on the
# conserved variables. WHY THAT FORM AND NOT THE PHYSICAL ONE (ν on ρ, ρν on u,
# κ on T, η on B — what problems/MHD/orszagTangBormanis2024 uses):
#
#   * The beam is a 10:1 density contact AND a 5.6-parts-per-million pressure
#     residue of the total energy (see the note in user_flux.jl). Diffusing T
#     and ρ separately does not keep p = (γ-1)(ρE - ½ρ|v|² - ½|B|²) admissible
#     across such a jump — it is the mechanism that cost
#     problems/MHD/fluxEmergenceSon2025 its positivity at the 25× contact of
#     the solar transition region. Diffusing the conserved state instead makes
#     the diffused node a convex combination of its neighbours' states, which
#     is where any positivity argument for an artificial-viscosity
#     regularization comes from (Guermond & Popov's viscous regularization of
#     the Euler system; Dao & Nazarov 2022 use the same form for Brio-Wu).
#   * It is in divergence form on every slot, so mass, momentum, total energy
#     and B are conserved to machine precision and the SHOCK SPEEDS — the thing
#     this benchmark is looked at for — are not altered by the stabilization.
#   * The magnetic and kinetic energy that the B and ρv slots remove reappears
#     in E automatically: the ν∇E Laplacian carries ∇·(νB·∇B) and the kinetic
#     analogue, which is why rhs.jl drops the separate τ·u viscous-work term in
#     this mode (add_tau_u is false under :dsgs_conserved).
#
# ONE CAVEAT on the 2D kernel, for the record. _expansion_visc! always treats
# slots 2 and 3 as momentum and builds a deviatoric stress from their
# gradients, so those two slots get ∇·τ(ρv) rather than ν∇²(ρu), ν∇²(ρv) — for
# a normal shock that is (4/3)ν instead of ν on the normal momentum component.
# It is still a divergence-form, dissipative operator and it is the established
# behaviour of the conserved-form MHD path here
# (problems/MHD/fluxEmergenceSon2025DSGS); the 4/3 could be taken out of :μ[2],
# :μ[3] if an exact Lax-Friedrichs form were ever wanted.
#
# :dsgs_nazarov_energy is deliberately left OFF. It would split the energy
# primitive into a non-thermal part at ν and a thermal part at κ = ρν/Pr, which
# is the right thing for a stratified atmosphere that has to conduct heat but
# breaks the single-ν convex-combination property above — and this case has
# nothing to gain from it (no stratification, no reference state).
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    for ieq = 1:9
        uprimitive[ieq] = u[ieq]
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    return SVector(u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/astroJetWuShu2018: PERT() solution variables are not supported.")
end

#---------------------------------------------------------------------------------
# Output variables:
#
#   qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T",
#               "log10rho", "log10p", "beta", "Mach"]
#
# log10rho and log10p are what the papers plot for this problem (Wu & Shu show
# schlieren images of log₁₀ρ and log₁₀p at t = 0.002, as do Balsara 2012 and the
# OEDG paper, the latter at t = 0.001, 0.0015 and 0.002), so they are written
# directly rather than left to a ParaView calculator. They are floored at 1e-300 so that a node that
# has gone inadmissible still writes a finite number instead of taking the whole
# output file down with a DomainError — if log10p saturates at -300 anywhere,
# the run has lost positivity there and the snapshot is a diagnostic, not a
# result. :lschlieren => true in the deck adds the two schlieren fields from
# |∇ρ| that those figures actually are.
#
#   T    = p/ρ        the nondimensional temperature of this nondimensional setup
#   beta = 2p/|B|²    plasma beta (10⁻² in the ambient medium at t = 0)
#   Mach = |v|/c      sonic Mach number, c = sqrt(γp/ρ) (800 in the beam)
#---------------------------------------------------------------------------------
function user_uout!(ip, ET, uout, u, qe; kwargs...)

    ρ  = u[1]
    p  = pressure_mhd(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8], u[9])
    B2 = u[6]*u[6] + u[7]*u[7] + u[8]*u[8]
    v2 = (u[2]*u[2] + u[3]*u[3] + u[5]*u[5])/(ρ*ρ)
    c2 = γ_mhd*p/ρ

    uout[1]  = ρ
    uout[2]  = u[2]/ρ   # u
    uout[3]  = u[3]/ρ   # v
    uout[4]  = u[5]/ρ   # w
    uout[5]  = p
    uout[6]  = u[6]     # Bx
    uout[7]  = u[7]     # By
    uout[8]  = u[8]     # Bz
    uout[9]  = u[9]     # ψ
    uout[10] = p/ρ                                   # T
    uout[11] = log10(max(ρ, 1.0e-300))               # log10 ρ
    uout[12] = log10(max(p, 1.0e-300))               # log10 p
    uout[13] = min(2.0*p/max(B2, 1.0e-300), 1.0e6)   # plasma β
    uout[14] = sqrt(max(v2, 0.0))/sqrt(max(c2, 1.0e-300))   # sonic Mach number
end

#---------------------------------------------------------------------------------
# Primitive variables handed to the viscous assembly, and the output map.
#
# uprimitive is what _expansion_visc! differentiates: slot i of the viscous
# flux is visc_coeff[i]·∇uprimitive[i] (plus the τ·u work term on the energy
# slot). The DynSGS coefficients produced by compute_dsgs_viscosity! for the
# total-energy system are those of Nazarov & Hoffman eq. (3.3)/(3.7)
#
#     mass      β ∇ρ        with β = μ/‖ρ‖∞,K
#     momentum  μ ε(u)
#     energy    κ ∇T        with κ = P/(γ-1)·μ
#
# so the primitive set must be (ρ, u, v, T) with T the paper's temperature
#
#     T = E/ρ - |u|²/2 = p/((γ-1)ρ),
#
# i.e. the specific internal energy — the paper scales cv = 1, and in SI
# that quantity is cv·T_physical. Putting the physical temperature in K here
# instead would silently rescale κ by cv ≈ 718 and over-diffuse the energy
# equation by that factor.
#
# user_uout! is independent of this: it writes the human-readable output
# fields named by qoutvars in initialize.jl, where T IS the physical
# temperature in K.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)

    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)   # T = p/((γ-1)ρ)
end

function user_primitives!(u, qe, uprimitive, ::PERT)

    ρ  = u[1] + qe[1]
    ρu = u[2] + qe[2]
    ρv = u[3] + qe[3]
    ρE = u[4] + qe[4]

    uprimitive[1] = ρ
    uprimitive[2] = ρu/ρ
    uprimitive[3] = ρv/ρ
    uprimitive[4] = ρE/ρ - 0.5*(ρu*ρu + ρv*ρv)/(ρ*ρ)
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    return T(ρ), T(ρu/ρ), T(ρv/ρ), T(ρE/ρ - T(0.5)*(ρu*ρu + ρv*ρv)/(ρ*ρ))
end

#
# Output fields, in the order of qoutvars = ["ρ", "u", "v", "p", "T"].
#
function user_uout!(ip, ET, uout, u, qe; kwargs...)

    PhysConst = PhysicalConst{Float64}()

    ρ  = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]

    uout[1] = ρ                                              # ρ
    uout[2] = ρu/ρ                                           # u
    uout[3] = ρv/ρ                                           # v
    uout[4] = PhysConst.γm1*(ρE - 0.5*(ρu*ρu + ρv*ρv)/ρ)     # p
    uout[5] = uout[4]/(PhysConst.Rair*ρ)                     # T [K]
end

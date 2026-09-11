#---------------------------------------------------------------------------------
# Conserved -> "primitive" mapping used by the DynSGS-MHD viscous kernel.
#
# This case runs DynSGS in its CONSERVED form (:dsgs_conserved => true): the
# 1D viscous loop applies one scalar Laplacian per slot, so with the
# conserved variables themselves as primitives and one kinematic ν on every
# slot the regularization is ∇·(ν∇q) on (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz) —
# in divergence form, hence exactly conservative (the shock speeds of a
# shock tube depend on it), and the form in which the magnetic and kinetic
# energy removed from the B and ρv slots is accounted for in E (the ν∇E
# Laplacian carries ∇·(νB·∇B) and the kinetic analogue; see
# problems/MHD/fluxEmergenceSon2025DSGS/user_primitives.jl). It is the
# Lax-Friedrichs-type viscosity Dao & Nazarov's residual method reduces to
# at its cap, on the conserved variables.
#---------------------------------------------------------------------------------
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    for ieq = 1:8
        uprimitive[ieq] = u[ieq]
    end
end

function user_primitives(u, qe, uprimitive, ::TOTAL)
    return SVector(u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8])
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    error(" problems/MHD/brioWu1d: PERT() solution variables are not supported.")
end

#---------------------------------------------------------------------------------
# Output variables: qoutvars = ["ρ", "u", "v", "p", "By"]
#---------------------------------------------------------------------------------
function user_uout!(ip, ET, uout, u, qe; kwargs...)
    ρ = u[1]
    uout[1] = ρ
    uout[2] = u[2]/ρ
    uout[3] = u[3]/ρ
    uout[4] = pressure_mhd1d(ρ, u[2], u[3], u[5], u[4], u[6], u[7], u[8])
    uout[5] = u[7]
end

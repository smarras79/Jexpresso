#---------------------------------------------------------------------------------
# Inviscid fluxes of the 2D compressible Euler equations in total-energy
# form, q = (ρ, ρu, ρv, ρE):
#
#   F = ( ρu,  ρu² + p,  ρuv,      (ρE + p)u )
#   G = ( ρv,  ρuv,      ρv² + p,  (ρE + p)v )
#
# with p = (γ-1)(ρE - ½ρ|u|²).
#
# Pressure is taken straight from the conservative state rather than
# through a temperature round-trip: across a shock the round-trip loses
# nothing analytically, but the direct form is one operation and cannot
# drift from the equation of state the DynSGS residual is built on.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]

    u = ρu/ρ
    v = ρv/ρ

    Pressure = PhysConst.γm1*(ρE - 0.5*(ρu*u + ρv*v))

    F[1] = ρu
    F[2] = ρu*u + Pressure
    F[3] = ρv*u
    F[4] = (ρE + Pressure)*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Pressure
    G[4] = (ρE + Pressure)*v
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρE = q[4]

    u = ρu/ρ
    v = ρv/ρ

    Pressure = PhysConst.γm1*(ρE - T(0.5)*(ρu*u + ρv*v))

    return T(ρu), T(ρu*u + Pressure), T(ρv*u), T((ρE + Pressure)*u),
           T(ρv), T(ρu*v), T(ρv*v + Pressure), T((ρE + Pressure)*v)
end

#---------------------------------------------------------------------------------
# KINETIC-ENERGY / ENTROPY-PRESERVING FLUX DIFFERENCING (:lkep => true).
#
# With :lkep the RHS is not assembled from the pointwise flux above. Instead
# _expansion_inviscid_KEP! (rhs.jl) builds it from SYMMETRIC TWO-POINT volume
# fluxes, and what this routine supplies per node is the auxiliary state those
# means need. Copied verbatim from CompEuler/kelvinHelmholtzChan2022, which is
# where this path is exercised for the 2D total-energy system.
#
# WHY IT IS HERE. A collocation CG/SEM integrates the nonlinear flux with the
# same LGL rule it interpolates on, so the flux is ALIASED: energy is
# transferred into the grid-scale modes, and CG has nothing to take it back
# out. At Mach 7 that matters far more than at Mach 3, because
#
#     p = (γ-1)(ρE - ½ρ|u|²)
#
# is a difference of two nearly equal numbers. A relative error in ρE or ρu
# comes out of that subtraction amplified by
#
#     (γ-1)·ρE/p = 1 + γ(γ-1)M²/2
#
# which is 3.5 at Mach 3 and 14.7 at Mach 7 — it grows like M². So the SAME
# aliasing error buys four times the pressure error here, the pressure drives
# the momentum flux, and the loop closes. Flux differencing with a KEP/EC
# two-point flux removes the aliasing-driven transfer by construction, which
# is exactly the term this failure is made of.
#
# :volume_flux picks the two-point flux: ranocha() (entropy conservative,
# the default), kennedy_gruber(), or central_euler() (plain central — the
# flux-differencing form of the standard scheme, useful as the control).
#---------------------------------------------------------------------------------
function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::central_euler)
    aux[1] = q[1]
    aux[2] = q[2]
    aux[3] = q[3]
    aux[4] = q[4]
end

function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::kennedy_gruber)

    PhysConst = PhysicalConst{Float64}()

    rho   = q[1]
    rho_u = q[2]
    rho_v = q[3]
    rho_e = q[4]

    u = rho_u/rho
    v = rho_v/rho

    p = PhysConst.γm1*(rho_e - (0.5*rho_u*u + 0.5*v*rho_v))

    aux[1] = rho
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = rho_e/rho
end

function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::ranocha)

    PhysConst = PhysicalConst{Float64}()

    rho   = q[1]
    rho_u = q[2]
    rho_v = q[3]
    rho_e = q[4]

    u = rho_u/rho
    v = rho_v/rho

    p = PhysConst.γm1*(rho_e - (0.5*rho_u*u + 0.5*v*rho_v))

    aux[1] = rho
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = rho_e
    aux[6] = log(rho)
    aux[7] = log(p)
end

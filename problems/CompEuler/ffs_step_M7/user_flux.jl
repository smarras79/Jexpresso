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

#---------------------------------------------------------------------------------
# The two-point volume fluxes themselves.
#
# _expansion_inviscid_KEP! (rhs.jl:2403) calls flux_turbo(aux_l, aux_r,
# volume_flux_type) for every node pair in the element. Like user_fluxaux!,
# flux_turbo is a PER-CASE function — the kernel declares the call, the case
# supplies the methods — so a deck that sets :lkep => true without these gets
#
#     UndefVarError: `flux_turbo` not defined
#
# from inside the first RHS evaluation. Copied verbatim from
# CompEuler/kelvinHelmholtzChan2022, which is where this path is exercised for
# the 2D total-energy system.
#
# Three volume fluxes, matching the three user_fluxaux! methods above:
#
#   ranocha()         entropy conservative (Ranocha 2018). The default here,
#                     and the one the aliasing argument in user_inputs.jl is
#                     aimed at. Needs aux slots 1-7, log(ρ) and log(p)
#                     included; the logarithmic means are inlined with the
#                     series expansion for nearly-equal arguments, so there is
#                     no 0/0 in the smooth free stream.
#   kennedy_gruber()  kinetic-energy preserving only, aux slots 1-5.
#   central_euler()   the plain central flux written in flux-differencing
#                     form, aux slots 1-4. The CONTROL: if a case behaves the
#                     same under this as under ranocha(), the two-point
#                     machinery is not what is helping. It needs flux(q,
#                     ::central_euler) below, which is the same pointwise flux
#                     as user_flux! above but returning SVectors.
#---------------------------------------------------------------------------------

@inline function flux(q, ::central_euler)
    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρe = q[4]

    e  = ρe/ρ
    u  = ρu/ρ
    v  = ρv/ρ

    γ   = PhysConst.γ
    γm1 = γ - 1.0
    
    velomagsq = (u*u + v*v)
    ke        = 0.5*ρ*velomagsq
    Pressure  = γm1*(ρe - ke)
    
    f1 = ρu
    f2 = ρu*u .+ Pressure
    f3 = ρv*u
    f4 = u*(ke + γ*Pressure/γm1)

    g1 = ρv
    g2 = ρu*v
    g3 = ρv*v .+ Pressure
    g4 = v*(ke + γ*Pressure/γm1)

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_turbo(u_ll, u_rr, ::ranocha)
    PhysConst = PhysicalConst{Float64}()
	rho_ll, v1_ll, v2_ll, p_ll, rho_e_ll, log_rho_ll, log_p_ll = u_ll
	rho_rr, v1_rr, v2_rr, p_rr, rho_e_rr, log_rho_rr, log_p_rr = u_rr
	    x1 = rho_ll
            log_x1 = log_rho_ll
            y1 = rho_rr
            log_y1 = log_rho_rr
            x1_plus_y1 = x1 + y1
            y1_minus_x1 = y1 - x1
            z1 = y1_minus_x1^2 / x1_plus_y1^2
            special_path1 = x1_plus_y1 / (2 + z1 * (2 / 3 + z1 * (2 / 5 + 2 / 7 * z1)))
            regular_path1 = y1_minus_x1 / (log_y1 - log_x1)
            rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

            # algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
            # in exact arithmetic since
            #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
            #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
            # inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
            x2 = rho_ll * p_rr
            log_x2 = log_rho_ll + log_p_rr
            y2 = rho_rr * p_ll
            log_y2 = log_rho_rr + log_p_ll
            x2_plus_y2 = x2 + y2
            y2_minus_x2 = y2 - x2
            z2 = y2_minus_x2^2 / x2_plus_y2^2
            special_path2 = (2 + z2 * (2 / 3 + z2 * (2 / 5 + 2 / 7 * z2))) / x2_plus_y2
            regular_path2 = (log_y2 - log_x2) / y2_minus_x2
            inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

            v1_avg = 0.5 * (v1_ll + v1_rr)
            v2_avg = 0.5 * (v2_ll + v2_rr)
            p_avg = 0.5 * (p_ll + p_rr)
            velocity_square_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr)
	    gamma = PhysConst.cp/PhysConst.cv
            # calculate fluxes depending on cartesian orientation
            f1 = rho_mean * v1_avg
            f2 = f1 * v1_avg + p_avg
            f3 = f1 * v2_avg
            f4 = f1 * (velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) + 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

            g1 = rho_mean * v2_avg
            g2 = g1 * v1_avg 
	    g3 = g1 * v2_avg + p_avg
            g4 = g1 * (velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) + 0.5 * (p_ll * v2_rr + p_rr * v2_ll)
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_turbo(u_ll, u_rr, ::kennedy_gruber)
    PhysConst = PhysicalConst{Float64}()
    rho_ll, v1_ll, v2_ll, p_ll, e_ll = u_ll
    rho_rr, v1_rr, v2_rr, p_rr, e_rr = u_rr

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (e_ll + e_rr)

    # Calculate fluxes depending on orientation
    
        f1 = rho_avg * v1_avg
        f2 = rho_avg * v1_avg * v1_avg + p_avg
        f3 = rho_avg * v1_avg * v2_avg
        f4 = (rho_avg * e_avg + p_avg) * v1_avg
    
        g1 = rho_avg * v2_avg
        g2 = rho_avg * v2_avg * v1_avg
        g3 = rho_avg * v2_avg * v2_avg + p_avg
        g4 = (rho_avg * e_avg + p_avg) * v2_avg

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_turbo(u_ll, u_rr, sol_type::central_euler)
	return 0.5f0 .* (flux(u_ll,sol_type) .+ flux(u_rr, sol_type))
end

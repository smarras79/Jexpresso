"""
    Minimal Surface Fluxes Implementation
    Based on Monin-Obukhov Similarity Theory

    This is a self-contained implementation of surface flux calculations
    following the formulations used in ClimateMachineAtmos.jl/SurfaceFluxes.jl
    but written for maximum portability and simplicity.

    References:
    - Businger et al. (1971): Flux-Profile Relationships in the Atmospheric Surface Layer
    - Dyer (1974): A review of flux-profile relationships
    - Monin & Obukhov (1954): Basic laws of turbulent mixing in the surface layer

    To run it in stand-alone mode, simply do:

    Julia> include("./path/to/this/file/CM_MOST.jl')
    Julia> MinimalSurfaceFluxes.example_with_plots()


    """

#module MinimalSurfaceFluxes
#=
try
    using Plots
    gr()
catch
    error("Plots.jl not available. Install with: using Pkg; Pkg.add(\"Plots\")")
end
=#

# Physical constants from PhysicalConst (defined in globalConstantsPhysics.jl)

"""
    qsat(T, p; over_water=true)

Saturation specific humidity [kg/kg] at temperature T [K] and pressure p [Pa].

Uses Tetens/Bolton formula for saturation vapor pressure.

# Arguments
- `T`: Temperature [K]
- `p`: Pressure [Pa]
- `over_water`: If true, use formula for water; if false, use formula for ice

# Example
```julia
T_sfc = 300.0      # Surface temperature [K]
p_sfc = 101325.0   # Surface pressure [Pa]
PhysConst = PhysicalConst{Float64}()
q_s = PhysConst.salt_factor * qsat(T_sfc, p_sfc)  # Surface specific humidity over ocean
```
"""
function qsat(T, p, PhysConst; over_water=true)
    # Tetens form, ula for saturation vapor pressure [Pa]
    if over_water
        # Valid for T > 253 K (over liquid water)
        e_s = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
    else
        # Over ice (for cold surfaces)
        e_s = 611.2 * exp(21.87 * (T - 273.15) / (T - 7.66))
    end
    # Convert to specific humidity: q = ε * e / (p - (1-ε)*e)
    ε = PhysConst.ε_ratio
    return ε * e_s / (p - (1 - ε) * e_s)
end

# Universal function parameters (Businger-Dyer)
const a_m = 16.0              # momentum stability parameter (unstable)
const a_h = 16.0              # heat stability parameter (unstable) 
const b_m = 5.0               # momentum stability parameter (stable)
const b_h = 5.0               # heat stability parameter (stable)

#==============================================================================
 MOST GUARD RAILS

 Businger-Dyer MOST has two singular limits, and an LES surface layer visits
 both of them at some node on almost every step:

   |u| -> 0   at a convergence line between thermals. u* -> 0, Q_H -> 0 with
              it, and L ~ u*^2 -> 0, so zeta = z/L -> -inf. Nothing in the
              formulation bounds the stability correction that comes back.
   Q_H -> 0   near-neutral. The old code returned a POSITIVE L = 1e6 here
              regardless of the sign of the flux, i.e. it reported a STABLE
              surface layer in the middle of a convective one. That is a sign
              discontinuity in the middle of the map a 20-step undamped Picard
              iteration is iterating, with no convergence flag and no fallback.

 Three settings fix both, and they are Refs rather than deck arguments because
 CM_MOST! is called per boundary node per RHS evaluation and a kwargs call
 there allocates a NamedTuple every time (the same reason the _surface_scales
 helpers below are positional-only). params_setup.jl sets them once from
 :most_u_min / :most_zeta_min / :most_zeta_max.

   MOST_U_MIN     Gustiness floor on the wind speed entering the DRAG
                  COEFFICIENT [m/s]. Not a floor on the stress: CM_MOST!
                  divides the stress direction by the same floored speed, so
                  tau stays linear in u_i and goes to zero with it, in a
                  well-defined direction, instead of turning into a
                  finite-magnitude random walk.
   MOST_ZETA_*    Bounds on z/L. The correction is evaluated at the bound
                  rather than extrapolated past it; [-5, 2] is the usual
                  range over which Businger-Dyer was fitted.

 The iteration itself is recast in 1/L rather than L. That is not cosmetic:
 1/L -> 0 IS the neutral limit and it is approached continuously from both
 sides, so the near-neutral case needs no special branch and cannot flip sign;
 zeta = z/L and zeta0 = z0/L are then z*invL and z0*invL, which stay mutually
 consistent when the bound binds (clamping zeta alone would leave the
 roughness term computed from an L the rest of the correction no longer uses).
==============================================================================#
const MOST_U_MIN    = Ref(0.1)     # :most_u_min    [m/s]
const MOST_ZETA_MIN = Ref(-5.0)    # :most_zeta_min
const MOST_ZETA_MAX = Ref(2.0)     # :most_zeta_max
const MOST_MAXITER  = Ref(20)
# Under-relaxation on 1/L. Undamped Picard on Businger-Dyer oscillates on the
# unstable side; 0.7 converges in a handful of steps over the whole clamped
# range without slowing the neutral case noticeably.
const MOST_RELAX    = Ref(0.7)
# |L| past which every z of interest is neutral. Only the sign matters here,
# and getting it right is the whole point -- see the header.
const MOST_L_BIG    = 1.0e6

@inline function most_clamp_invL(z_ref, invL)
    (isfinite(invL) && z_ref > 0.0) || return 0.0
    return clamp(invL, MOST_ZETA_MIN[] / z_ref, MOST_ZETA_MAX[] / z_ref)
end

# 1/L = -kappa*g*Q_H / (rho*cp*u*^3*theta).  Reciprocal of obukhov_length
# below, including its rho, and returns the neutral 0 rather than an Inf or a
# NaN for every degenerate input.
@inline function most_inv_L(u_star, theta_ref, Q_H, PhysConst, rho)
    den = rho * PhysConst.cp * u_star^3 * theta_ref
    (isfinite(den) && den > 0.0) || return 0.0
    invL = -PhysConst.karman * PhysConst.g * Q_H / den
    return isfinite(invL) ? invL : 0.0
end

#------------------------------------------------------------------------------
# The shared fixed point behind every _surface_scales_* below.
#
# wtheta_imposed = NaN  -> classic MOST: theta* is diagnosed from the air/surface
#                          theta difference.
# wtheta_imposed finite -> PRESCRIBED SURFACE HEAT FLUX: theta* = -wtheta/u*, so
#                          Q_H = rho*cp*wtheta is fixed and L is built from the
#                          flux the model actually applies. Used when the deck
#                          sets :user_heatflux, where the alternative is a drag
#                          stability-corrected by a flux that is then thrown
#                          away (BCs.jl discards MOST's own wtheta when
#                          delta_hf = 1) and diagnosed from a free surface node
#                          that nothing pins.
#
# Returns (u_star, theta_star, invL, denom_h, converged). denom_h is handed
# back so the moist caller can build q* on the same stability correction
# instead of recomputing psi_h.
#------------------------------------------------------------------------------
function _most_iterate(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h,
                       PhysConst, rho, wtheta_imposed)
    κ  = PhysConst.karman
    cp = PhysConst.cp

    lnzm = log(z_ref / z0_m)
    lnzh = log(z_ref / z0_h)
    if !(isfinite(lnzm) && lnzm > 0.0 && isfinite(lnzh) && lnzh > 0.0 &&
         isfinite(theta_ref) && theta_ref > 0.0 && rho > 0.0)
        # z_ref below the roughness length, a non-positive theta, or a bad rho:
        # there is no surface layer to model. Zero fluxes beat a NaN that
        # reaches the RHS.
        return (0.0, 0.0, 0.0, lnzh, false)
    end

    u_eff = max(u_ref, MOST_U_MIN[])
    lflux = isfinite(wtheta_imposed)

    u_star     = κ * u_eff / lnzm
    theta_star = lflux ? -wtheta_imposed / u_star :
                          κ * (theta_ref - theta_s) / lnzh
    invL = most_clamp_invL(z_ref,
             most_inv_L(u_star, theta_ref, -rho * cp * u_star * theta_star, PhysConst, rho))

    w      = MOST_RELAX[]
    dh     = lnzh
    conv   = false
    for _ = 1:MOST_MAXITER[]
        ζ   = z_ref * invL
        ζ0m = z0_m * invL
        ζ0h = z0_h * invL

        dm = lnzm - psi_m(ζ) + psi_m(ζ0m)
        dh = lnzh - psi_h(ζ) + psi_h(ζ0h)
        # Both are comfortably positive across the clamped zeta range (the
        # unstable branch tends to (1/4)log(z/z0), not to zero). This is a net
        # against a pathological roughness ratio, not a physical bound.
        dm = max(dm, 0.1 * lnzm)
        dh = max(dh, 0.1 * lnzh)

        u_star_new     = κ * u_eff / dm
        theta_star_new = lflux ? -wtheta_imposed / u_star_new :
                                  κ * (theta_ref - theta_s) / dh
        invL_raw = most_inv_L(u_star_new, theta_ref,
                              -rho * cp * u_star_new * theta_star_new, PhysConst, rho)
        invL_new = most_clamp_invL(z_ref, invL + w * (invL_raw - invL))

        # Converge on zeta, which is what the correction is a function of --
        # not on L, whose relative change is meaningless as it passes through
        # the neutral limit.
        conv = abs((invL_new - invL) * z_ref) < 1.0e-5
        u_star     = u_star_new
        theta_star = theta_star_new
        invL       = invL_new
        conv && break
    end

    if !(isfinite(u_star) && isfinite(theta_star))
        # Fall back to the neutral log law rather than propagating a NaN.
        u_star     = κ * u_eff / lnzm
        theta_star = lflux ? -wtheta_imposed / u_star :
                              κ * (theta_ref - theta_s) / lnzh
        invL = 0.0
        conv = false
    end
    return (u_star, theta_star, invL, dh, conv)
end

# Internal positional-only helpers — no keyword arguments to avoid NamedTuple heap allocation
# per boundary-point call.  The kwargs version of _surface_scales below delegates here.

function _surface_scales_dry(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho)
    return _surface_scales_dry(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho, NaN)
end

function _surface_scales_dry(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho,
                             wtheta_imposed)
    u_star, theta_star, _, _, _ =
        _most_iterate(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho, wtheta_imposed)
    return u_star, theta_star
end

function _surface_scales_moist(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho, q_ref, q_s)
    return _surface_scales_moist(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho,
                                 q_ref, q_s, NaN)
end

function _surface_scales_moist(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho,
                               q_ref, q_s, wtheta_imposed)
    κ = PhysConst.karman
    u_star, theta_star, _, dh, _ =
        _most_iterate(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho, wtheta_imposed)
    # Same stability correction as heat (Businger-Dyer), so q* rides on the
    # denominator the iteration converged on rather than on a fresh psi_h call.
    q_star = dh > 0.0 ? κ * (q_ref - q_s) / dh : 0.0
    return u_star, theta_star, q_star
end

# Backward-compatible kwargs wrapper (not on the hot path — BCs.jl calls the positional overloads below).
function _surface_scales(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst;
                                  q_ref=nothing, q_s=nothing, rho=1.225, max_iter=20, tol=1e-4)
    if !isnothing(q_ref) && !isnothing(q_s)
        return _surface_scales_moist(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho, q_ref, q_s)
    else
        u_star, theta_star = _surface_scales_dry(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst, rho)
        return u_star, theta_star, zero(u_ref)
    end
end

# --- Hot-path CM_MOST! overloads: all positional, no kwargs ---

# The stress direction. |u| is floored by the SAME MOST_U_MIN the drag
# coefficient sees, so tau_i = -rho*u*^2*u_i/max(|u|,u_min) stays linear in u_i
# and vanishes with it. The old form divided by |u| + 2.22e-16, which is a
# guard against a division by zero and not against anything physical: with a
# wind floor in u* it would have handed a finite stress an arbitrary direction
# at every stagnation node.
@inline function _most_stress!(τ_f, ρ, u_star, u_ref, v_ref, w_ref, u_magnitude)
    u_eff = max(u_magnitude, MOST_U_MIN[])
    if !(u_eff > 0.0) || !isfinite(u_star)
        # Only reachable with :most_u_min => 0 at a node whose tangential wind
        # is exactly zero. There is no direction to put a stress in.
        τ_f[1] = 0.0; τ_f[2] = 0.0; τ_f[3] = 0.0
        return
    end
    τ_magnitude = ρ * u_star * u_star
    τ_f[1] = -τ_magnitude * (u_ref / u_eff)
    τ_f[2] = -τ_magnitude * (v_ref / u_eff)
    τ_f[3] = -τ_magnitude * (w_ref / u_eff)
    return
end

# Dry (no moisture): z0_m and z0_h are positional
function CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst, z0_m, z0_h)
    CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst, z0_m, z0_h, NaN)
end

# Dry, with a PRESCRIBED surface kinematic heat flux (K m/s). theta*, and so
# the stability correction on the drag, is built from that flux instead of from
# the air/surface theta difference; wθ[1] comes back equal to it by
# construction. See the _most_iterate header.
function CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst,
                  z0_m, z0_h, wtheta_imposed)
    u_magnitude = sqrt(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref)
    u_star, theta_star = _surface_scales_dry(u_magnitude, theta_ref, z_ref, theta_s,
                                             z0_m, z0_h, PhysConst, ρ, wtheta_imposed)
    _most_stress!(τ_f, ρ, u_star, u_ref, v_ref, w_ref, u_magnitude)
    wθ[1] = -u_star * theta_star
end


# Backward-compatible kwargs wrapper (kept for external callers; NOT the hot path).
function CM_MOST!(τ_f, wθ, wq, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst;
                  q_ref=nothing, q_s=nothing, z0_m=0.1, z0_h=0.01)
    if !isnothing(q_ref) && !isnothing(q_s)
        CM_MOST!(τ_f, wθ, wq, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst, q_ref, q_s, z0_m, z0_h)
    else
        CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst, z0_m, z0_h)
    end
end

# Backward compatible version without wq argument (dry, default roughness lengths).
function CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst)
    CM_MOST!(τ_f, wθ, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref, PhysConst, 0.1, 0.01, NaN)
end

function CM_MOST!(τ_f, wθ, wqv, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref,
                  PhysConst, qv_in, qv_sfc, z0_m, z0_h)

    CM_MOST!(τ_f, wθ, wqv, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref,
             PhysConst, qv_in, qv_sfc, z0_m, z0_h, NaN)
end

function CM_MOST!(τ_f, wθ, wqv, ρ, u_ref, v_ref, w_ref, theta_ref, theta_s, z_ref,
                  PhysConst, qv_in, qv_sfc, z0_m, z0_h, wtheta_imposed)

    u_magnitude = sqrt(u_ref*u_ref + v_ref*v_ref + w_ref*w_ref)
    # One fixed point for all three scales, on the same stability correction.
    # This used to run the dry solve for tau and wtheta and then a SECOND,
    # independent surface_conditions() solve for the moisture flux, so C_H came
    # from a different converged L than the one the drag used.
    u_star, theta_star, q_star =
        _surface_scales_moist(u_magnitude, theta_ref, z_ref, theta_s, z0_m, z0_h,
                              PhysConst, ρ, qv_in, qv_sfc, wtheta_imposed)
    _most_stress!(τ_f, ρ, u_star, u_ref, v_ref, w_ref, u_magnitude)
    wθ[1]  = -u_star * theta_star
    wqv[1] = -u_star * q_star

    return
end


"""
    Universal stability functions φ_m and φ_h following Businger-Dyer (1971)

    For unstable conditions (ζ < 0):
    - φ_m(ζ) = (1 - a_m * ζ)^(-1/4)
    - φ_h(ζ) = (1 - a_h * ζ)^(-1/2)

    For stable conditions (ζ >= 0):
    - φ_m(ζ) = 1 + b_m * ζ
    - φ_h(ζ) = 1 + b_h * ζ
    """
function phi_m(zeta)
    if zeta < 0
        return (1 - a_m * zeta)^(-0.25)
    else
        return 1 + b_m * zeta
    end
end

function phi_h(zeta)
    if zeta < 0
        return (1 - a_h * zeta)^(-0.5)
    else
        return 1 + b_h * zeta
    end
end

"""
    Integrated stability functions ψ_m and ψ_h (diabatic correction functions)

    These are the integrals of the universal functions, used in profile calculations.
    """
function psi_m(zeta)
    if zeta < 0
        # Unstable conditions
        x = (1 - a_m * zeta)^0.25
        return 2 * log((1 + x) / 2) + log((1 + x^2) / 2) - 2 * atan(x) + π/2
    else
        # Stable conditions
        return -b_m * zeta
    end
end

function psi_h(zeta)
    if zeta < 0
        # Unstable conditions  
        y = (1 - a_h * zeta)^0.5
        return 2 * log((1 + y) / 2)
    else
        # Stable conditions
        return -b_h * zeta
    end
end

"""
    Calculate Obukhov length L from surface fluxes

    L = -rho * cp * u_star³ * T_ref / (κ * g * Q_H)

    where:
    - u_star: friction velocity [m/s]
    - T_ref: reference potential temperature [K]
    - Q_H: surface sensible heat flux [W/m²]
    - rho: TOTAL air density at the reference height [kg/m³]

    THE rho IS NOT OPTIONAL AND WAS MISSING. Q_H is handed in as a DYNAMIC flux
    in W/m² (-rho*cp*u_star*theta_star), so the rho*cp that turns the kinematic
    u_star*theta_star into it has to appear in the numerator too. Without it
    this returned L/rho, i.e. |L| about 20% small at sea level, and every
    zeta = z/L built from it about 20% too large -- a systematic overstatement
    of the stability correction that grows with altitude as rho falls.

    Callers must pass the TOTAL density. Under :SOL_VARS_TYPE => PERT() that is
    uaux[ip,1] + qe[ip,1], not uaux[ip,1]; BCs.jl already reconstructs it that
    way before calling CM_MOST!.

    SIGN CONVENTION, and the second thing that was wrong here: Q_H > 0 is an
    UPWARD (surface-heating) flux, which is UNSTABLE, which is L < 0. The
    near-neutral guard used to return +1e6 for either sign, i.e. it reported a
    stable surface layer in the middle of a convective one -- and it was
    reached whenever u_star got small, which in an LES surface layer is at
    some node on essentially every step. The guard now carries the sign of the
    flux through.
    """
function obukhov_length(u_star, T_ref, Q_H, PhysConst, rho)
    # L < 0 for Q_H > 0. Q_H == 0 is exactly neutral, where the sign is
    # immaterial and either bound is as good as the other.
    Lsign = Q_H > 0.0 ? -1.0 : 1.0
    (abs(Q_H) > 1e-6 && rho > 0.0 && u_star > 0.0 && isfinite(T_ref)) ||
        return Lsign * MOST_L_BIG
    L = -rho * PhysConst.cp * u_star^3 * T_ref / (PhysConst.karman * PhysConst.g * Q_H)
    isfinite(L) || return Lsign * MOST_L_BIG
    abs(L) > MOST_L_BIG && return Lsign * MOST_L_BIG
    return L
end

"""
    Wind profile based on Monin-Obukhov similarity theory

    u(z) = (u_star/κ) * [ln(z/z0_m) - ψ_m(z/L) + ψ_m(z0_m/L)]

    where:
    - z: height above surface [m]
    - z0_m: momentum roughness length [m]
    - L: Obukhov length [m]
    - u_star: friction velocity [m/s]
    """
function wind_profile(z, z0_m, L, u_star, PhysConst)
    zeta = z / L
    zeta0 = z0_m / L
    return (u_star / PhysConst.karman) * (log(z / z0_m) - psi_m(zeta) + psi_m(zeta0))
end

"""
    Temperature profile based on Monin-Obukhov similarity theory

    θ(z) = θ_s + (θ_star/κ) * [ln(z/z0_h) - ψ_h(z/L) + ψ_h(z0_h/L)]

    where:
    - z: height above surface [m]
    - z0_h: thermal roughness length [m] 
    - L: Obukhov length [m]
    - θ_s: surface potential temperature [K]
    - θ_star: temperature scale [K]
    """
function temperature_profile(z, z0_h, L, theta_s, theta_star, PhysConst)
    zeta = z / L
    zeta0 = z0_h / L
    return theta_s + (theta_star / PhysConst.karman) * (log(z / z0_h) - psi_h(zeta) + psi_h(zeta0))
end

"""
    Calculate friction velocity from wind speed and stability

    This solves iteratively:
    u = (u_star/κ) * [ln(z/z0_m) - ψ_m(z/L)]

    Given u, z, z0_m, and L, find u_star
    """
function friction_velocity_from_wind(u, z, z0_m, L, PhysConst; max_iter=20, tol=1e-6)
    zeta = z / L
    zeta0 = z0_m / L
    κ = PhysConst.karman

    # Initial guess
    u_star = κ * u / log(z / z0_m)

    for i in 1:max_iter
        # Calculate wind from current u_star
        u_calc = (u_star / κ) * (log(z / z0_m) - psi_m(zeta) + psi_m(zeta0))
        
        # Check convergence
        error = abs(u_calc - u)
        if error < tol
            break
        end
        
        # Update u_star (simple iteration)
        u_star *= u / u_calc
    end
    
    return u_star
end

"""
    Calculate temperature scale from temperature difference and stability

    θ_star = -Q_H / (ρ * cp * u_star)

    where Q_H is the sensible heat flux [W/m²]
    """
function temperature_scale(Q_H, rho, u_star, PhysConst)
    return -Q_H / (rho * PhysConst.cp * u_star)
end

"""
    Main surface flux calculation routine

    Given atmospheric conditions at measurement height, calculate surface fluxes
    and similarity parameters.

    Inputs:
    - u_ref: wind speed at reference height [m/s]
    - theta_ref: potential temperature at reference height [K]
    - z_ref: reference height [m]
    - theta_s: surface potential temperature [K]
    - z0_m: momentum roughness length [m]
    - z0_h: thermal roughness length [m]
    - q_ref: specific humidity at reference height [kg/kg] (optional, default nothing)
    - q_s: surface specific humidity [kg/kg] (optional, default nothing)
    - rho: air density [kg/m³] (optional, default 1.225)

    Returns NamedTuple with:
    - u_star: friction velocity [m/s]
    - theta_star: temperature scale [K]
    - q_star: humidity scale [kg/kg] (if humidity provided)
    - L: Obukhov length [m]
    - Q_H: sensible heat flux [W/m²]
    - Q_E: latent heat flux [W/m²] (if humidity provided)
    - zeta: stability parameter z_ref/L
    - C_D: drag coefficient
    - C_H: heat transfer coefficient
    - C_E: moisture transfer coefficient
    """
function surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst;
                            q_ref=nothing, q_s=nothing, rho=1.225, max_iter=20, tol=1e-4)

    κ  = PhysConst.karman
    cp = PhysConst.cp

    # Check if humidity is provided
    compute_lhf = !isnothing(q_ref) && !isnothing(q_s)

    # Same gustiness floor and zeta bounds as the hot path (see the MOST GUARD
    # RAILS block): this function reports the transfer coefficients for the
    # fluxes CM_MOST! applies, so it has to be limited the same way or the two
    # disagree at exactly the low-wind nodes where the limits bind.
    u_ref = max(u_ref, MOST_U_MIN[])

    # Initial guess assuming neutral conditions
    u_star = κ * u_ref / log(z_ref / z0_m)
    theta_star = κ * (theta_ref - theta_s) / log(z_ref / z0_h)
    q_star = compute_lhf ? κ * (q_ref - q_s) / log(z_ref / z0_h) : 0.0

    # Initial heat flux guess
    Q_H = -rho * cp * u_star * theta_star
    L = obukhov_length(u_star, theta_ref, Q_H, PhysConst, rho)

    # Iterative solution
    for iter in 1:max_iter
        # Update stability parameter. Bounding zeta and re-deriving L from the
        # bound keeps zeta and zeta0 mutually consistent -- clamping zeta alone
        # would leave the roughness terms on an L the rest no longer uses.
        zeta = clamp(z_ref / L, MOST_ZETA_MIN[], MOST_ZETA_MAX[])
        L    = zeta == 0.0 ? L : z_ref / zeta
        zeta0_m = z0_m / L
        zeta0_h = z0_h / L

        # Calculate new friction velocity
        u_star_new = κ * u_ref / (log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m))

        # Calculate new temperature scale
        theta_star_new = κ * (theta_ref - theta_s) / (log(z_ref / z0_h) - psi_h(zeta) + psi_h(zeta0_h))

        # Calculate new humidity scale (uses same stability function as heat)
        if compute_lhf
            q_star_new = κ * (q_ref - q_s) / (log(z_ref / z0_h) - psi_h(zeta) + psi_h(zeta0_h))
        else
            q_star_new = 0.0
        end

        # Update heat flux and Obukhov length
        Q_H_new = -rho * cp * u_star_new * theta_star_new
        L_new = obukhov_length(u_star_new, theta_ref, Q_H_new, PhysConst, rho)

        # Check convergence
        error = abs(L_new - L) / max(abs(L), abs(L_new))

        if error < tol
            u_star = u_star_new
            theta_star = theta_star_new
            q_star = q_star_new
            Q_H = Q_H_new
            L = L_new
            break
        end

        # Update for next iteration
        u_star = u_star_new
        theta_star = theta_star_new
        q_star = q_star_new
        Q_H = Q_H_new
        L = L_new
    end

    # Calculate transfer coefficients
    zeta = clamp(z_ref / L, MOST_ZETA_MIN[], MOST_ZETA_MAX[])
    L    = zeta == 0.0 ? L : z_ref / zeta
    zeta0_m = z0_m / L
    zeta0_h = z0_h / L

    C_D = (κ / (log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m)))^2
    C_H = κ^2 / ((log(z_ref / z0_m) - psi_m(zeta) + psi_m(zeta0_m)) *
        (log(z_ref / z0_h) - psi_h(zeta) + psi_h(zeta0_h)))
    C_E = C_H  # Moisture transfer coefficient (same as heat for Businger-Dyer)

    # Calculate latent heat flux: Q_E = -ρ * Lv * u* * q*
    # Positive Q_E means upward flux (evaporation)
    Q_E = compute_lhf ? -rho * PhysConst.Lc * u_star * q_star : 0.0

    return (
        u_star = u_star,
        theta_star = theta_star,
        q_star = q_star,
        L = L,
        Q_H = Q_H,
        Q_E = Q_E,
        zeta = zeta,
        C_D = C_D,
        C_H = C_H,
        C_E = C_E
    )
end

"""
    Calculate momentum flux (wind stress) from friction velocity

    τ = ρ * u_star²

    Returns momentum flux in [N/m²] or [Pa]
    """
function momentum_flux(u_star, rho=1.225)
    return rho * u_star^2
end


#=
"""
    Plot vertical profiles of wind speed and potential temperature

    Requires Plots.jl package. Install with: using Pkg; Pkg.add("Plots")

    Parameters:
    - result: output from surface_conditions()
    - z0_m, z0_h: roughness lengths [m]
    - theta_s: surface potential temperature [K]
    - z_max: maximum height for plots [m] (default: 100)
    - n_points: number of points in profiles (default: 100)
    """
function plot_profiles(result, z0_m, z0_h, theta_s; z_max=100.0, n_points=100)
    
    # Create height array (logarithmic spacing for better resolution near surface)
    z_min = max(z0_m, z0_h) * 1.1  # Start just above roughness length
    heights = exp.(range(log(z_min), log(z_max), length=n_points))
    
    # Calculate profiles
    wind_speeds = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temperatures = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Create plots
    p1 = plot(wind_speeds, heights, 
              xlabel="Wind Speed [m/s]", 
              ylabel="Height [m]",
              title="Wind Speed Profile\n(ζ = $(round(result.zeta, digits=3)), L = $(round(result.L, digits=1)) m)",
              linewidth=2,
              color=:blue,
              grid=true,
              legend=false)
    
    p2 = plot(temperatures, heights,
              xlabel="Potential Temperature [K]",
              ylabel="Height [m]", 
              title="Temperature Profile\n(Q_H = $(round(result.Q_H, digits=1)) W/m²)",
              linewidth=2,
              color=:red,
              grid=true,
              legend=false)
    
    # Add reference level markers if within range
    ref_height = 10.0  # Assuming 10m reference from example
    if z_min <= ref_height <= z_max
        u_ref = wind_profile(ref_height, z0_m, result.L, result.u_star)
        theta_ref = temperature_profile(ref_height, z0_h, result.L, theta_s, result.theta_star)
        
        scatter!(p1, [u_ref], [ref_height], 
                 markersize=6, color=:blue, markershape=:circle,
                 label="Reference ($(ref_height)m)")
        
        scatter!(p2, [theta_ref], [ref_height],
                 markersize=6, color=:red, markershape=:circle, 
                 label="Reference ($(ref_height)m)")
    end
    
    # Combine plots
    plot(p1, p2, layout=(1, 2), size=(800, 400))
end

"""
    Plot comparison between neutral and stratified profiles

    Shows how atmospheric stability affects the profiles compared to neutral conditions
    """
function plot_stability_comparison(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; z_max=100.0)
    
    # Height array  
    z_min = max(z0_m, z0_h) * 1.1
    heights = exp.(range(log(z_min), log(z_max), length=100))
    
    # Stratified profiles (actual)
    wind_stratified = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temp_stratified = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Neutral profiles for comparison (L → ∞)
    L_neutral = 1e6  # Very large L approximates neutral conditions
    wind_neutral = [wind_profile(z, z0_m, L_neutral, result.u_star) for z in heights]
    temp_neutral = [temperature_profile(z, z0_h, L_neutral, theta_s, result.theta_star) for z in heights]
    
    # Wind speed comparison
    p1 = plot(wind_neutral, heights,
              label="Neutral", 
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Wind Speed [m/s]",
              ylabel="Height [m]",
              title="Wind Profiles: Stratified vs Neutral")
    
    plot!(p1, wind_stratified, heights,
          label="Stratified (ζ=$(round(result.zeta, digits=3)))",
          linewidth=2,
          color=:blue)
    
    # Temperature comparison  
    p2 = plot(temp_neutral, heights,
              label="Neutral",
              linestyle=:dash, 
              linewidth=2,
              color=:gray,
              xlabel="Potential Temperature [K]",
              ylabel="Height [m]",
              title="Temperature Profiles: Stratified vs Neutral")
    
    plot!(p2, temp_stratified, heights,
          label="Stratified (Q_H=$(round(result.Q_H, digits=1)) W/m²)",
          linewidth=2,
          color=:red)
    
    # Add reference point
    scatter!(p1, [u_ref], [z_ref], markersize=6, color=:blue, label="Reference")
    scatter!(p2, [theta_ref], [z_ref], markersize=6, color=:red, label="Reference")
    
    plot(p1, p2, layout=(1, 2), size=(900, 400))
end

"""
    Plot universal functions φ_m and φ_h vs stability parameter ζ
    """
function plot_universal_functions(zeta_range=(-2.0, 2.0); n_points=200)
    
    zeta_vals = range(zeta_range[1], zeta_range[2], length=n_points)
    phi_m_vals = [phi_m(ζ) for ζ in zeta_vals]
    phi_h_vals = [phi_h(ζ) for ζ in zeta_vals]
    
    p1 = plot(zeta_vals, phi_m_vals,
              xlabel="ζ = z/L",
              ylabel="φ_m(ζ)",
              title="Momentum Universal Function",
              linewidth=2,
              color=:blue,
              label="φ_m (Businger-Dyer)",
              grid=true)
    
    p2 = plot(zeta_vals, phi_h_vals,
              xlabel="ζ = z/L", 
              ylabel="φ_h(ζ)",
              title="Heat Universal Function",
              linewidth=2,
              color=:red,
              label="φ_h (Businger-Dyer)",
              grid=true)
    
    # Add neutral line
    hline!(p1, [1.0], linestyle=:dash, color=:gray, label="Neutral")
    hline!(p2, [1.0], linestyle=:dash, color=:gray, label="Neutral")
    
    # Add stability regime labels
    annotate!(p1, [(-1, 3)], text("Unstable\n(Convective)", 10))
    annotate!(p1, [(1, 6)], text("Stable", 10))
    annotate!(p2, [(-1, 2)], text("Unstable\n(Convective)", 10)) 
    annotate!(p2, [(1, 6)], text("Stable", 10))
    
    plot(p1, p2, layout=(1, 2), size=(800, 400))
end

"""
    Create comprehensive plot matrix showing all aspects of surface flux analysis

    Returns a matrix layout with 6 subplots:
    - Top row: Wind profile, Temperature profile  
    - Middle row: Wind comparison (stratified vs neutral), Temperature comparison
    - Bottom row: Universal functions φ_m(ζ), φ_h(ζ)
    """
function create_comprehensive_plots(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; z_max=50.0)
    
    # Height array for profiles
    z_min = max(z0_m, z0_h) * 1.1
    heights = exp.(range(log(z_min), log(z_max), length=100))
    
    # Calculate all profiles
    wind_stratified = [wind_profile(z, z0_m, result.L, result.u_star) for z in heights]
    temp_stratified = [temperature_profile(z, z0_h, result.L, theta_s, result.theta_star) for z in heights]
    
    # Neutral profiles for comparison
    L_neutral = 1e6
    wind_neutral = [wind_profile(z, z0_m, L_neutral, result.u_star) for z in heights]
    temp_neutral = [temperature_profile(z, z0_h, L_neutral, theta_s, result.theta_star) for z in heights]
    
    # Universal function data
    zeta_vals = range(-2.0, 2.0, length=200)
    phi_m_vals = [phi_m(ζ) for ζ in zeta_vals]
    phi_h_vals = [phi_h(ζ) for ζ in zeta_vals]
    
    # Create individual plots
    
    # Plot 1: Wind profile
    p1 = plot(wind_stratified, heights,
              xlabel="Wind Speed [m/s]",
              ylabel="Height [m]", 
              title="Wind Speed Profile",
              linewidth=2.5,
              color=:blue,
              grid=true,
              legend=false,
              titlefontsize=11)
    scatter!(p1, [u_ref], [z_ref], markersize=5, color=:blue, markershape=:circle)
    annotate!(p1, [(u_ref+0.3, z_ref)], text("Ref", 8))
    
    # Plot 2: Temperature profile  
    p2 = plot(temp_stratified, heights,
              xlabel="Temperature [K]",
              ylabel="Height [m]",
              title="Temperature Profile", 
              linewidth=2.5,
              color=:red,
              grid=true,
              legend=false,
              titlefontsize=11)
    scatter!(p2, [theta_ref], [z_ref], markersize=5, color=:red, markershape=:circle)
    annotate!(p2, [(theta_ref-0.15, z_ref)], text("Ref", 8))
    
    # Plot 3: Wind comparison
    p3 = plot(wind_neutral, heights,
              label="Neutral",
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Wind Speed [m/s]", 
              ylabel="Height [m]",
              title="Wind: Stratified vs Neutral",
              titlefontsize=11,
              legendfontsize=8)
    plot!(p3, wind_stratified, heights,
          label="Stratified",
          linewidth=2.5,
          color=:blue)
    scatter!(p3, [u_ref], [z_ref], markersize=4, color=:blue, label="")
    
    # Plot 4: Temperature comparison
    p4 = plot(temp_neutral, heights,
              label="Neutral", 
              linestyle=:dash,
              linewidth=2,
              color=:gray,
              xlabel="Temperature [K]",
              ylabel="Height [m]",
              title="Temperature: Stratified vs Neutral", 
              titlefontsize=11,
              legendfontsize=8)
    plot!(p4, temp_stratified, heights,
          label="Stratified",
          linewidth=2.5,
          color=:red)
    scatter!(p4, [theta_ref], [z_ref], markersize=4, color=:red, label="")
    
    # Plot 5: Universal function φ_m
    p5 = plot(zeta_vals, phi_m_vals,
              xlabel="ζ = z/L",
              ylabel="φ_m(ζ)",
              title="Momentum Universal Function",
              linewidth=2.5,
              color=:blue,
              grid=true,
              legend=false,
              titlefontsize=11)
    hline!(p5, [1.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p5, [0.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p5, [result.zeta], linestyle=:dash, color=:blue, alpha=0.8, linewidth=2)
    
    # Plot 6: Universal function φ_h  
    p6 = plot(zeta_vals, phi_h_vals,
              xlabel="ζ = z/L",
              ylabel="φ_h(ζ)", 
              title="Heat Universal Function",
              linewidth=2.5,
              color=:red,
              grid=true,
              legend=false,
              titlefontsize=11)
    hline!(p6, [1.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p6, [0.0], linestyle=:dot, color=:gray, alpha=0.7)
    vline!(p6, [result.zeta], linestyle=:dash, color=:red, alpha=0.8, linewidth=2)
    
    # Combine into matrix layout
    plot(p1, p2, p3, p4, p5, p6,
         layout=(3, 2),
         size=(800, 900),
         margin=4Plots.mm,
         plot_title="Surface Flux Analysis: MOST Theory\nζ=$(round(result.zeta,digits=3)), L=$(round(result.L,digits=1))m, Q_H=$(round(result.Q_H,digits=1))W/m²",
         plot_titlefontsize=14)
end

"""
    Save comprehensive surface flux analysis to PDF

    Creates a complete analysis including:
    - Calculated parameters summary 
    - 6-panel plot matrix with all key visualizations
    - Saves to specified filename

    Parameters:
    - result: output from surface_conditions()
    - filename: output PDF filename (default: "surface_flux_analysis.pdf")
    - Additional parameters for plotting
    """
function save_analysis_to_pdf(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref; 
                              filename="surface_flux_analysis.pdf", z_max=50.0)
    println("Creating comprehensive surface flux analysis...")
    
    # Create the comprehensive plot matrix
    p = create_comprehensive_plots(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref, z_max=z_max)
    
    # Save to PDF
    try
        savefig(p, filename)
        println("Analysis saved to: $filename")
        
        # Also display the plot
        display(p)
        
        return p
        
    catch e
        println("Error saving PDF: $e")
        println("Note: PDF saving requires a working LaTeX installation or try PNG format instead")
        
        # Fallback to PNG
        png_filename = replace(filename, ".pdf" => ".png")
        savefig(p, png_filename)
        println("Saved as PNG instead: $png_filename")
        
        return p
    end
end

"""
    Enhanced example with comprehensive plotting and PDF output
    """
function example_with_plots()
    println("=== Minimal Surface Fluxes: Comprehensive Analysis ===")
    
    # Example atmospheric conditions
    u_ref = 10.0      # wind speed at 10m [m/s]
    theta_ref = 298.0 # 285.0  # potential temperature at 10m [K] 
    theta_s = 302.0   # 288.0    # surface potential temperature [K]
    z_ref = 10.0      # reference height [m]
    z0_m = 0.1        # momentum roughness length [m]
    z0_h = 0.01       # thermal roughness length [m]
    
    # Calculate surface conditions
    result = surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h)
    
    println("\n=== CALCULATED PARAMETERS ===")
    println("  Friction velocity u* = $(round(result.u_star, digits=4)) m/s")
    println("  Temperature scale θ* = $(round(result.theta_star, digits=4)) K") 
    println("  Obukhov length L = $(round(result.L, digits=1)) m")
    println("  Sensible heat flux = $(round(result.Q_H, digits=1)) W/m²")
    println("  Stability parameter ζ = $(round(result.zeta, digits=4))")
    println("  Drag coefficient CD = $(round(result.C_D * 1000, digits=2)) × 10⁻³")
    println("  Heat transfer coeff CH = $(round(result.C_H * 1000, digits=2)) × 10⁻³")
    
    # Calculate momentum flux
    tau = momentum_flux(result.u_star)
    println("  Momentum flux τ = $(round(tau, digits=4)) N/m²")
    
    # Determine stability regime
    if result.zeta < -0.1
        stability = "Moderately Unstable (Convective)"
    elseif result.zeta < 0
        stability = "Weakly Unstable" 
    elseif result.zeta < 0.1
        stability = "Near Neutral"
    else
        stability = "Stable"
    end
    println("  Atmospheric stability: $stability")
    
    println("\n=== GENERATING COMPREHENSIVE ANALYSIS ===")
    
    # Create and save comprehensive analysis
    try
        p = save_analysis_to_pdf(result, z0_m, z0_h, theta_s, u_ref, theta_ref, z_ref,
                                 filename="MOST_surface_flux_analysis.pdf")
        
        println("\n=== ANALYSIS COMPLETE ===")
        println("✓ Surface flux calculations completed")
        println("✓ Comprehensive 6-panel visualization created")
        println("✓ Results saved to PDF")
        println("✓ Physical consistency verified")
        
    catch e
        println("Plotting failed: $e")
        println("Install Plots.jl with: using Pkg; Pkg.add(\"Plots\")")
    end
end

#end # module MinimalSurfaceFluxes
=#


# Run example if this file is executed directly
#MinimalSurfaceFluxes.example_with_plots()

"""
    Example usage and validation function
    """
function example_usage()
    println("=== Minimal Surface Fluxes Example ===")
    
    # Example atmospheric conditions
    u_ref     = 10.0  # wind speed at 10m [m/s]
    theta_ref = 298.0 # 285.0  # potential temperature at 10m [K] 
    theta_s   = 302.0 # 288.0  # surface potential temperature [K]
    z_ref     = 10.0  # reference height [m]
    z0_m      = 0.1   # momentum roughness length [m]
    z0_h      = 0.01  # thermal roughness length [m]
    PhysConst = PhysicalConst{Float64}()

    # Calculate surface conditions
    result = surface_conditions(u_ref, theta_ref, z_ref, theta_s, z0_m, z0_h, PhysConst)
    
    println("Results:")
    println("  Friction velocity u* = $(round(result.u_star, digits=4)) m/s")
    println("  Temperature scale θ* = $(round(result.theta_star, digits=4)) K") 
    println("  Obukhov length L = $(round(result.L, digits=1)) m")
    println("  Sensible heat flux = $(round(result.Q_H, digits=1)) W/m²")
    println("  Stability parameter ζ = $(round(result.zeta, digits=4))")
    println("  Drag coefficient CD = $(round(result.C_D * 1000, digits=2)) × 10⁻³")
    println("  Heat transfer coeff CH = $(round(result.C_H * 1000, digits=2)) × 10⁻³")
    
    # Calculate momentum flux
    tau = momentum_flux(result.u_star)
    println("  Momentum flux τ = $(round(tau, digits=4)) N/m²")
    
    # Test profile functions
    heights = [2.0, 5.0, 10.0, 20.0, 50.0]
    println("\nWind profile:")
    for z in heights
        u = wind_profile(z, z0_m, result.L, result.u_star, PhysConst)
        println("  u($(z)m) = $(round(u, digits=2)) m/s")
    end

    println("\nTemperature profile:")
    for z in heights
        theta = temperature_profile(z, z0_h, result.L, theta_s, result.theta_star, PhysConst)
        println("  θ($(z)m) = $(round(theta, digits=2)) K")
    end
end

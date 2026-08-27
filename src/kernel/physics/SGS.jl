#----------------------------------------------------------------------
# SHARED SGS HELPERS
#----------------------------------------------------------------------
#
# Buoyancy stability function for the eddy viscosity (Lilly 1962;
# Deardorff 1980; Moeng 1984):
#
#     f(Ri) = sqrt(max(0, 1 - Ri/Pr_t)),   Ri = N²/S²
#
# so that
#
#     ν_t = (C_s Δ)² |S| f(Ri) = (C_s Δ)² sqrt(S² - N²/Pr_t).
#
# Stable (Ri > 0) reduces mixing and shuts it off at Ri = Pr_t; unstable
# (Ri < 0) enhances it. The right-hand form shows the unstable branch is
# well behaved without a cap: as |S| → 0 it tends to (C_s Δ)² sqrt(-N²/Pr_t),
# it does not blow up.
#
# This replaces the earlier branch that used the Monin-Obukhov SURFACE-LAYER
# function min(sqrt(1 - 16 Ri), 3) on the unstable side. The 16 there belongs
# to the φ_m profile relationship, not to an SGS closure; against Lilly's
# coefficient of 1/Pr_t ≈ 1.43 it over-predicts the eddy viscosity by up to
# the 3x cap throughout a convective mixed layer, i.e. exactly where a CBL
# LES is most sensitive to spurious dissipation.
#
@inline function sgs_stability_function(Ri, Pr_t)
    return sqrt(max(0.0, 1.0 - Ri/Pr_t))
end

#
# SGS mixing length with the standard near-wall limit
#
#     1/ℓ² = 1/(C_s Δ)² + 1/(κ z)²      (Mason & Thomson 1992 with n = 2)
#
# Away from the wall ℓ → C_s Δ; approaching it ℓ → κ z, which is what the
# surface layer actually supports. Without this limit the first few nodes of
# a wall-modelled LES carry an eddy viscosity set by the grid rather than by
# the distance to the wall, over-mixing the surface layer.
#
# Returns ℓ², so callers use μ_turb = ρ ℓ² |S| f_Ri in place of
# ρ C_s² Δ² |S| f_Ri.
#
# NOTE: `z` must be the distance to the WALL. For a flat lower boundary that
# is the height above it, which is what the callers pass. On a terrain-
# following mesh the height above terrain is not currently available at this
# point, so :lwall_damping must stay false for warped grids.
#
@inline function sgs_mixing_length2(C_s2, Δ2, z, karman, lwall_damping)
    CsΔ2 = C_s2 * Δ2
    (lwall_damping && z > 0.0) || return CsΔ2
    κz2 = (karman * z)^2
    return CsΔ2 * κz2 / (CsΔ2 + κz2)
end

#----------------------------------------------------------------------
# SMAGORINSKY
#----------------------------------------------------------------------
#
# LEGACY per-node path. When :visc_model is SMAG()/VREM(), params.sgs is a
# concrete AbstractSGSModel and the RHS goes through compute_sgs_cache! plus
# the cache-reading SGS_diffusion further down instead of these methods. They
# are kept for the sgs === nothing dispatches only, and they read C_s from
# PhysConst directly, so they do NOT honour an inputs[:C_s] override — put new
# closure work in compute_sgs_cache!, not here.
#
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::SMAG, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_temperature = (ieq == 4)
    
    #
    # Neutral/unstable: Pr_t ≈ 0.7 - 0.85
    # Stable:           Pr_t ≈ 1.0 - 2.0 (usually handled with Richardson corrections)
    # Very unstable:    Pr_t ≈ 1/3
    #
    Pr_t  = PhysConst.Pr_t   # Turbulent Prandtl number
    Sc_t  = PhysConst.Sc_t   # Turbulent Schmidt number for other scalars
    μ_mol = PhysConst.μ_mol  # Molecular viscosity [Pa·s]
    κ_mol = PhysConst.κ_mol  # Molecular thermal diffusivity [m²/s]
    C_s   = PhysConst.C_s    # Smagorinsky constant
    cp    = PhysConst.cp
    C_s2  = C_s*C_s
    # EOS-consistent cp = γ·R/(γ-1). Using PhysConst.cp directly is unsafe
    # for non-dimensional setups where Rair is rescaled but cp is left at its
    # SI value, which would over-estimate k_eff by orders of magnitude.
    cp    = PhysConst.γ * PhysConst.Rair / PhysConst.γm1


    # Smagorinsky
    # Strain rate tensor (symmetric part of velocity gradient)
    S11 = u11
    S22 = u22
    S12 = 0.5 * (u12 + u21)
    S21 = S12

    # Strain rate magnitude
    # |S| = sqrt(2 * S_ij * S_ij)
    S_ij_S_ij = S11*S11 + S22*S22 + 2.0*S12*S12
    Sij       = sqrt(2.0 * S_ij_S_ij)

    # Turbulent viscosity (same for all equations)
    μ_turb = ρ * C_s2 * Δ2 * Sij
    if is_u_momentum || is_v_momentum

        return (μ_mol + μ_turb) * visc_coeffieq[ieq] # effective viscosity

    elseif is_temperature
        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (ρ*κ_mol + κ_turb) * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end


#=============================================================================
 NOT CALLED. NOTHING IN THE TREE REACHES THIS METHOD.

 The 3D viscous path goes through the CACHE form instead. _expansion_visc!
 (rhs.jl, NSD_3D, ContGal) calls compute_sgs_cache! ONCE per element and then,
 inside the equation loop, the seven-argument reader

     SGS_diffusion(visc_coeffieq, ieq, rho, ip, sgs, ltheta_eqn, SD)

 which is defined far below in this file and just looks up sgs.mu_turb[ip].
 The strain rate, the buoyancy correction and the FILTER WIDTH Delta^2 all
 live in compute_sgs_cache! now -- so an @info added here never prints, and a
 change made here has no effect on any 3D run.

 The 2D long form directly above IS live: the NSD_2D _expansion_visc! still
 passes the velocity gradients and Delta2 explicitly.

 AND IT HAS DRIFTED, which is the reason for a banner rather than a deletion.
 This copy has no near-wall length limit -- no lwall_damping, no karman, no
 l = min(C_s*Delta, kappa*z) -- while compute_sgs_cache! does. Routing 3D back
 through here would silently drop the wall damping and inflate nu_t in the
 first element above the surface, which is exactly where an LES with a wall
 model can least afford it.

 To instrument the filter width, put the @info in compute_sgs_cache! (the
 SGS_SMAG method), whose Delta2 argument is the square of
 mesh.Delta_elem_filter[iel]/nop. Or read it off the mesh driver's
 "LES filter width" line, which prints the global range at startup and needs
 no code change at all.
=============================================================================#
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u33,
                               u12, u21,
                               u13, u31,
                               u23, u32,
                               θ_ref,
                               dθdz,
                               PhysConst, Δ2,
                               ::SMAG, ::NSD_3D;
                               ltheta_eqn=true,
                               lrichardson=false)
    
    PhysConst = PhysicalConst{Float64}()
    C_s   = PhysConst.C_s       # Smagorinsky constant
    Pr_t  = PhysConst.Pr_t      # Turbulent Prandtl number
    Sc_t  = PhysConst.Sc_t      # Turbulent Schmidt number
    μ_mol = PhysConst.μ_mol     # Molecular viscosity [Pa·s]
    κ_mol = PhysConst.κ_mol     # Molecular thermal diffusivity [m²/s]
    g     = PhysConst.g
    cp    = PhysConst.cp
    C_s2  = C_s*C_s
    
    # Equation type identification
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_w_momentum  = (ieq == 4)
    is_temperature = (ieq == 5)
    
    # ===== 3D Strain Rate Tensor =====
    # Diagonal components
    S11 = u11  # ∂u/∂x
    S22 = u22  # ∂v/∂y
    S33 = u33  # ∂w/∂z
    
    # Off-diagonal components (symmetrized)
    S12 = 0.5 * (u12 + u21)  # 0.5*(∂u/∂y + ∂v/∂x)
    S13 = 0.5 * (u13 + u31)  # 0.5*(∂u/∂z + ∂w/∂x)
    S23 = 0.5 * (u23 + u32)  # 0.5*(∂v/∂z + ∂w/∂y)
    
    # Strain rate magnitude squared (for Richardson number)
    # S² = 2*S_ij*S_ij
    S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
    Sij2      = 2.0 * S_ij_S_ij  # S² = 2*S_ij*S_ij
    Sij       = sqrt(Sij2)         # |S|
    
    # ===== Richardson Number Correction =====
    # Applied to all equations (momentum, temperature, tracers) via shared f_Ri
    f_Ri = 1.0  # Default: no correction
    
    if ltheta_eqn && lrichardson
        
        # Buoyancy frequency squared: N² = (g/θ) * dθ/dz
        # Positive N² indicates stable stratification
        # Note: dθdz should be the actual vertical derivative (not just computational)
        N2 = abs(θ_ref) > 1.0f-12 ? (g / θ_ref) * dθdz : 0.0
        
        # Richardson number: Ri = N²/S²
        # Ri > 0: stable stratification (suppresses turbulence)
        # Ri < 0: unstable stratification (enhances turbulence)
        # Ri >= Pr_t: turbulence completely suppressed
        Ri = (Sij2 > 1.0f-12) ? N2 / Sij2 : 0.0
        
        # Stability function for Richardson correction
        # Various formulations exist in literature
        f_Ri = sgs_stability_function(Ri, Pr_t)
    elseif lrichardson && !ltheta_eqn
        # ===== Moist Richardson Number Logic =====
        # Note: In this mode, the caller has pre-calculated:
        # θ_ref  => T_abs (Absolute Temperature in Kelvin)
        # dθdz   => dhl_eff_dz = [1/(cp*(1+γ)) * dhl/dz] - [T_abs * dqn/dz]
        # 
        # This effective gradient accounts for:
        # 1. Latent heat release via the (1+γ) moist adjustment factor.
        # 2. Hydrometeor loading (weight of liquid/ice) via the dqndz term.

        # Buoyancy frequency squared using the moist-effective gradient: 
        # N²m = (g / T_abs) * dhl_eff_dz
        # Units: [m/s²] / [K] * [K/m] = [s⁻²]
        N2 = abs(θ_ref) > 1.0f-12 ? (g / θ_ref) * dθdz : 0.0
        
        # Richardson number: Ratio of buoyancy resistance to shear production
        # Ri = N²m / S²
        Ri = (Sij2 > 1.0f-12) ? N2 / Sij2 : 0.0
        
        # Stability function for Richardson correction (Smagorinsky scaling)
        f_Ri = sgs_stability_function(Ri, Pr_t)
    end
    
    # Turbulent viscosity with Richardson correction
    # μ_turb = ρ * (C_s * Δ)² * |S| * f(Ri)
    μ_turb = ρ * C_s2 * Δ2 * Sij * f_Ri
    
    # ===== Return appropriate coefficient based on equation type =====
    if is_u_momentum || is_v_momentum || is_w_momentum
        # Momentum equations use effective viscosity
        return (μ_mol + μ_turb) * visc_coeffieq[ieq]
        
    elseif is_temperature
        # Temperature equation uses effective thermal diffusivity
        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t
        if ltheta_eqn
            # Potential temperature equation
            return κ_turb * visc_coeffieq[ieq]
        else
            # Internal energy or enthalpy equation
            return (ρ*κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
        
    else
        # Other scalar equations (species, TKE, etc.)
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end
end

#----------------------------------------------------------------------
# VREMAN
#----------------------------------------------------------------------
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::VREM, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    
    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_temperature = (ieq == 4)
    
    #
    # Neutral/unstable: Pr_t ≈ 0.7 - 0.85
    # Stable:           Pr_t ≈ 1.0 - 2.0 (usually handled with Richardson corrections)
    # Very unstable:    Pr_t ≈ 1/3
    #
    Pr_t       = PhysConst.Pr_t   # Turbulent Prandtl number
    Sc_t       = PhysConst.Sc_t   # Turbulent Schmidt number for other scalars
    μ_mol      = PhysConst.μ_mol  # Molecular viscosity [Pa·s]
    κ_mol      = PhysConst.κ_mol  # Molecular thermal diffusivity [m²/s]
    C_s        = PhysConst.C_s    # Smagorinsky constant
    C_s2       = C_s*C_s
    # EOS-consistent cp = γ·R/(γ-1); see 2D-SMAG comment for rationale.
    cp         = PhysConst.γ * PhysConst.Rair / PhysConst.γm1
    C_vrem     = 2.5 * C_s2  # Vreman coefficient
    eps_vreman = eps(1.0)    # Safety epsilon

    # Vreman β tensor
    β11 = Δ2 * (u11*u11 + u12*u12)
    β12 = Δ2 * (u11*u21 + u12*u22)
    β22 = Δ2 * (u21*u21 + u22*u22)

    B_β = β11*β22 - β12*β12
    
    # Frobenius norm squared of velocity gradient
    u_ij_u_ij = u11*u11 + u12*u12 + u21*u21 + u22*u22

    
    # Vreman eddy viscosity with safety checks
    if u_ij_u_ij > eps_vreman && B_β > 0.0
        μ_turb = ρ * C_vrem * sqrt(B_β / u_ij_u_ij)
    else
        μ_turb = 0.0
    end
    
    if is_u_momentum || is_v_momentum
        return (μ_mol + μ_turb) * visc_coeffieq[ieq] # effective viscosity

    elseif  is_temperature # Assuming potential temperature equation is at index 4

        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            # Total-energy / enthalpy form: thermal conductivity k_eff = cp * μ_eff / Pr_t
            return cp * (μ_mol + μ_turb) / Pr_t * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end



#=============================================================================
 NOT CALLED. NOTHING IN THE TREE REACHES THIS METHOD.

 The 3D viscous path goes through the CACHE form instead. _expansion_visc!
 (rhs.jl, NSD_3D, ContGal) calls compute_sgs_cache! ONCE per element and then,
 inside the equation loop, the seven-argument reader

     SGS_diffusion(visc_coeffieq, ieq, rho, ip, sgs, ltheta_eqn, SD)

 which is defined far below in this file and just looks up sgs.mu_turb[ip].
 The strain rate, the buoyancy correction and the FILTER WIDTH Delta^2 all
 live in compute_sgs_cache! now -- so an @info added here never prints, and a
 change made here has no effect on any 3D run.

 The 2D long form directly above IS live: the NSD_2D _expansion_visc! still
 passes the velocity gradients and Delta2 explicitly.

 AND IT HAS DRIFTED, which is the reason for a banner rather than a deletion.
 This copy has no near-wall length limit -- no lwall_damping, no karman, no
 l = min(C_s*Delta, kappa*z) -- while compute_sgs_cache! does. Routing 3D back
 through here would silently drop the wall damping and inflate nu_t in the
 first element above the surface, which is exactly where an LES with a wall
 model can least afford it.

 To instrument the filter width, put the @info in compute_sgs_cache! (the
 SGS_VREM method), whose Delta2 argument is the square of
 mesh.Delta_elem_filter[iel]/nop. Or read it off the mesh driver's
 "LES filter width" line, which prints the global range at startup and needs
 no code change at all.
=============================================================================#
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u12, u13,
                               u21, u22, u23,
                               u31, u32, u33,
                               θ_ref, dθdz,
                               PhysConst, Δ2,
                               ::VREM, ::NSD_3D;
                               ltheta_eqn=true,
                               lrichardson=false)

    # NOTE: this row-wise parameter list (u11, u12, u13, u21, …) does
    # NOT match the diagonals-then-off-diagonals order that the call
    # sites in src/kernel/operators/rhs.jl pass — which is what the
    # SMAG 3D signature above expects. So `u12` here actually receives
    # `dvdy`, `u13` receives `dwdz`, etc.
    #
    # Because the user's `e95cb259` reference (last commit where
    # CompEuler/3d worked) had exactly this mismatch and was treated
    # as the correct numerical baseline, we keep the parameter order
    # as-is. The Vreman β below is consequently *not* the textbook
    # formula but the de-facto formula this code has been calibrated
    # against. If the textbook β is what you want, reorder the
    # parameter list to `u11, u22, u33, u12, u21, u13, u31, u23, u32`
    # (matching SMAG) — the body uses u_ij as ∂u_i/∂x_j and will then
    # compute the standard β_ij = Σ_m Δ² u_im u_jm.

    is_u_momentum  = (ieq == 2)
    is_v_momentum  = (ieq == 3)
    is_w_momentum  = (ieq == 4)
    is_temperature = (ieq == 5)

    Pr_t       = PhysConst.Pr_t   # Turbulent Prandtl number
    Sc_t       = PhysConst.Sc_t   # Turbulent Schmidt number for other scalars
    μ_mol      = PhysConst.μ_mol  # Molecular viscosity [Pa·s]
    κ_mol      = PhysConst.κ_mol  # Molecular thermal diffusivity [m²/s]
    g          = PhysConst.g         # Gravitational acceleration (m/s²)
    C_s        = PhysConst.C_s    # Smagorinsky constant
    C_s2       = C_s*C_s
    cp         = PhysConst.cp
    C_vrem     = 2.5 * C_s2  # Vreman coefficient
    eps_vreman = eps(1.0)    # Safety epsilon
    
    # Vreman β tensor (3D)
    # β_ij = Δ_m^2 * u_im * u_jm (sum over m=1,2,3)
    β11 = Δ2 * (u11*u11 + u12*u12 + u13*u13)
    β12 = Δ2 * (u11*u21 + u12*u22 + u13*u23)
    β13 = Δ2 * (u11*u31 + u12*u32 + u13*u33)
    β22 = Δ2 * (u21*u21 + u22*u22 + u23*u23)
    β23 = Δ2 * (u21*u31 + u22*u32 + u23*u33)
    β33 = Δ2 * (u31*u31 + u32*u32 + u33*u33)
    
    # B_β for 3D
    B_β = β11*β22 + β11*β33 + β22*β33 - (β12*β12 + β13*β13 + β23*β23)
    
    # Frobenius norm squared of 3x3 velocity gradient tensor
    u_ij_u_ij =
        u11*u11 + u12*u12 + u13*u13 +
        u21*u21 + u22*u22 + u23*u23 +
        u31*u31 + u32*u32 + u33*u33

    
    f_Ri = 1.0
    if ltheta_eqn && lrichardson
        
        # Strain rate tensor (symmetric part of velocity gradient)
        S11 = u11
        S22 = u22
        S33 = u33
        S12 = 0.5 * (u12 + u21)
        S13 = 0.5 * (u13 + u31)
        S23 = 0.5 * (u23 + u32)
        
        # Strain rate magnitude
	# |S| = sqrt(2 * S_ij * S_ij)
        S_ij_S_ij  = S11^2 + S22^2 + S33^2 + 2.0*(S12^2 + S13^2 + S23^2)
        Sij2        = 2.0 * S_ij_S_ij

        # Buoyancy frequency squared: N² = (g/θ) * dθ/dz
        # Note: assuming z is vertical (modify if different coordinate system)
        N2 = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdz : 0.0
        
        # Richardson number: Ri = N²/S²
        Ri = (Sij2 > 1e-12) ? N2 / Sij2 : 0.0
        
        # Stability function for Richardson correction
        # Various formulations exist; using a smooth transition
        f_Ri = sgs_stability_function(Ri, Pr_t)
    end
    
    # Vreman eddy viscosity with safety checks
    if u_ij_u_ij > eps_vreman && B_β > 0.0
        μ_turb = ρ * C_vrem * sqrt(B_β / u_ij_u_ij) * f_Ri
    else
        μ_turb = 0.0
    end
    
    if is_u_momentum || is_v_momentum
        return (μ_mol + μ_turb) * visc_coeffieq[ieq] # effective viscosity
    elseif  is_temperature # Assuming potential temperature equation is at index 4

        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return cp * (ρ*κ_mol + κ_turb) * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end

#----------------------------------------------------------------------
# DYNAMIC SGS (Marras et al., 2015, JCP 301:77-101) — residual-based
# artificial viscosity, parameter free.
#
# Per-element coefficient:
#     μ_res = C1 · Δ² · max_i ‖R_i‖∞,Ω / ‖q_i − ⟨q_i⟩‖∞,Ω
#     μ_max = C2 · Δ · max(|u| + c)
#     μ_dsgs[iel] = max(0, min(μ_res, μ_max))
# where R_i is the STRONG-form BDF2 residual of conservation law i —
# (3qⁿ − 4qⁿ⁻¹ + qⁿ⁻²)/(2Δt) − M⁻¹·RHS. Since jexpresso assembles RHS
# in weak form (post-DSS, pre-mass-matrix division), the rhs argument
# is multiplied by Minv[ip] inline before the BDF2 minus, which makes
# μ_dsgs dimensionally a kinematic viscosity (m²/s) regardless of SD.
#
# Both numerators and denominators are L∞ norms over a region larger than
# one element — the rank's subdomain by default, the whole domain under
# :ldsgs_global_norms (see _dsgs_norm_scope below) — so the coefficient
# cannot be inlined into the (k,l) loop the way
# SMAG/VREM are — it is precomputed once per RHS call into the
# pre-allocated μ_dsgs[1:nelem] buffer. SGS_diffusion(::DSGS, ::SD)
# is the standard per-quadrature-point accessor — the caller updates
# visc_coeffieq with the current element's μ_dsgs[iel] before
# entering the (k,l) loop, so this just returns it.
#----------------------------------------------------------------------

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
                               ::DSGS, ::NSD_1D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return visc_coeffieq[ieq]

end

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
                               ::DSGS, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return visc_coeffieq[ieq]

end

# The no-`inputs` counterpart. The generic 2D _expansion_visc! calls
# SGS_diffusion with two different argument lists — the momentum and
# scalar branches pass (…, PhysConst, Δ2, VT, SD), only the τ·u
# viscous-work term passes `inputs` before VT — and until this method
# existed the 2D DSGS() path raised MethodError on its first call
# (ieq = 1, the "other scalars" branch at rhs.jl:2207), i.e.
# problems/CompEuler/theta_dsgs could not start at all.
#
# NOTE this makes that case *run*; it does not make it correct. Two
# other defects on the same path are untouched and deliberate to leave
# alone (see DSGS.md §6): the residual there omits M⁻¹, and the two
# momentum slots are still zeroed by a leftover diagnostic block, so
# only the ρθ equation is actually stabilized.
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::DSGS, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return visc_coeffieq[ieq]

end

# Same accessor for the MHD variant: compute_dsgs_viscosity!(::DSGS_MHD)
# has already packed the per-element, per-equation coefficient into
# visc_coeffieq, so the assembly loop just reads it back.
#
# TWO methods are needed because the 2D _expansion_visc! calls
# SGS_diffusion with two different argument lists: the momentum/scalar
# branches pass (…, PhysConst, Δ2, VT, SD) while the τ·u viscous-work
# term in the total-energy branch passes (…, PhysConst, Δ2, inputs, VT, SD).
# (The ::DSGS, ::NSD_2D pair above only defines the `inputs` form, so the
# Euler-θ 2D DSGS path MethodErrors on the first call — it has evidently
# never been exercised. Not touched here.)
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::DSGS_MHD, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return visc_coeffieq[ieq]

end

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
                               ::DSGS_MHD, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return visc_coeffieq[ieq]

end

# ================================================================================
# _dsgs_norm_scope — how far the DynSGS normalising scales reach
#
# Every implementation below normalises the element residual by a mean ⟨q_i⟩
# and an L∞ spread ‖q_i − ⟨q_i⟩‖. Marras eq. (9) and Nazarov & Hoffman
# eq. (3.5) write both over the whole domain Ω. Under MPI that is a
# collective, and it is on the critical path of every rank on EVERY RHS call
# — five times per step under CarpenterKennedy2N54, times two or three
# reductions each.
#
# Default: RANK-LOCAL (`lglobal_norms = false`). These two quantities only set
# the SCALE the residual indicator is measured against; what the model needs
# from them is the order of magnitude of the solution's variation, and a
# partition of a connected domain resolves that as well as the whole domain
# does. μ is bounded by min(μ_res, μ_max) either way, so the flow solution
# differs only at the level of the usual round-off divergence. No
# communication at all.
#
# Opt-in: the paper's domain norms, with
#
#     :ldsgs_global_norms => true      # in user_inputs.jl
#
# threaded down from rhs.jl. Use it when you want μ reproducible across rank
# counts — a regression test that compares fields bit-for-bit between a
# 1-rank and an N-rank run — or when a subdomain genuinely cannot see the
# solution's scale (a partition that lies entirely inside a uniform region
# while the interesting structure lives on another rank). Costs 2-3
# Allreduce per RHS call.
#
# The communicator is Jexpresso's own get_mpi_comm(), NOT MPI.COMM_WORLD:
# under MPMD coupling COMM_WORLD also carries Alya's ranks, which never call
# into DynSGS, so a collective on it deadlocks every Jexpresso rank.
# ================================================================================

# ---------------- 1D --------------------------------------------------
#
# Conservation form q = (ρ, ρu, ρE) on a 1D LGL mesh. The signature is
# a hand-typed function barrier (concrete arrays, no params.* lookups)
# so Julia can specialize and the inner loop is allocation-free.
#
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS, ::NSD_1D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 wres::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δx::AbstractVector{TT},
                                 nelem::Int, ngl::Int;
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    # 1D CompEuler in total-energy form q = (ρ, ρu, ρE). Marras's
    # unified formula gives ONE residual-based coefficient per element;
    # for visualisation parity with the 2D version we replicate it into
    # every column of μ_dsgs[iel, :] so the caller / VTU sees per-
    # equation slots even when they are identical. The user-supplied
    # inputs[:μ] vector enters as a per-equation multiplicative
    # factor so the user can scale the DSGS contribution down (or off)
    # equation by equation.

    invnp = one(TT)/(nelem*ngl)
    γ     = TT(1.4)
    C1    = TT(1.0)
    C2    = TT(0.5)
    eps   = Base.eps(TT)
    neqs  = size(μ_dsgs, 2)

    # qe is accepted for forward compatibility with the 2D signature
    # but the 1D test cases (case1, sod1d) have qe ≈ 0 so subtracting
    # it would not change the denominators meaningfully.

    # --- Pass 1: averages of q (see _dsgs_norm_scope) -------------------
    ρ_avg  = zero(TT); ρu_avg = zero(TT); ρE_avg = zero(TT)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            ρ_avg  += q[ip,1]
            ρu_avg += q[ip,2]
            ρE_avg += q[ip,3]
        end
    end
    if lglobal_norms
        sums = TT[ρ_avg, ρu_avg, ρE_avg, TT(nelem*ngl)]
        MPI.Allreduce!(sums, MPI.SUM, get_mpi_comm())
        invnp_g = one(TT)/max(sums[4], one(TT))
        ρ_avg = sums[1]*invnp_g; ρu_avg = sums[2]*invnp_g; ρE_avg = sums[3]*invnp_g
    else
        ρ_avg  *= invnp
        ρu_avg *= invnp
        ρE_avg *= invnp
    end

    # --- Pass 2: L∞ norms of |q - ⟨q⟩| ---------------------------------
    denom1 = zero(TT); denom2 = zero(TT); denom3 = zero(TT)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = max(denom1, abs(q[ip,1] - ρ_avg))
            denom2 = max(denom2, abs(q[ip,2] - ρu_avg))
            denom3 = max(denom3, abs(q[ip,3] - ρE_avg))
        end
    end
    if lglobal_norms
        norms = TT[denom1, denom2, denom3]
        MPI.Allreduce!(norms, MPI.MAX, get_mpi_comm())
        denom1 = norms[1]; denom2 = norms[2]; denom3 = norms[3]
    end
    denom1 += eps; denom2 += eps; denom3 += eps

    # --- Pass 3: per-element loop --------------------------------------
    inv2Δt = one(TT)/(2*Δt)
    @inbounds for ie = 1:nelem
        Δ = Δx[ie]/ngl

        n1   = zero(TT); n2 = zero(TT); n3 = zero(TT)
        uTmx = zero(TT)
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            # w = 0 on a domain boundary, where the residual is the boundary
            # condition rather than a discretisation error. See dsgs_wres in
            # params_setup.jl.
            w  = wres[ip]
            R1 = w*abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 = w*abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 = w*abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
            n1 = max(n1, R1); n2 = max(n2, R2); n3 = max(n3, R3)

            ρl = q[ip,1]
            ul = q[ip,2]/ρl
            el = q[ip,3]/ρl
            # Specific internal energy, then the sound speed. For a perfect
            # gas p = (γ-1)ρ·e_int, so a² = γp/ρ = γ(γ-1)·e_int. The (γ-1)
            # was previously missing, which inflated the wave-speed cap by
            # 1/sqrt(γ-1) ≈ 1.58 at γ = 1.4 and let μ_res govern more often
            # than the Marras bound intends.
            eint = max(el - TT(0.5)*ul*ul, zero(TT))
            uTmx = max(uTmx, abs(ul) + sqrt(γ*(γ - one(TT))*eint))
        end

        μ_res = C1*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3)
        μ_max = C2*Δ*uTmx
        μ     = max(zero(TT), min(μ_max, μ_res))

        # Same coefficient on every equation (1D E-form, Marras eq. 10),
        # scaled per equation by the user-supplied inputs[:μ] vector.
        for ieq = 1:neqs
            μ_dsgs[ie, ieq] = visc_coeff[ieq] * μ
        end
    end

    return nothing
end

# ---------------- 2D --------------------------------------------------
#
# Conservation form q = (ρ, ρu, ρv, ρθ) for the Euler-θ system. Δ is
# min(Δx, Δy)/(N+1) (Marras et al. eq. 8), and c is built from the
# perfect-gas-law for θ:  p = C0·(ρθ)^γ ⇒ c² = γp/ρ. Same
# function-barrier discipline as the 1D variant — no params accesses,
# no struct constructions, no allocations.
#
# `ltheta` selects which system slot 4 belongs to, and is passed down
# from inputs[:energy_equation] by the rhs.jl call site:
#
#   ltheta = true  (default, :energy_equation => "theta")
#       q = (ρ, ρu, ρv, ρθ), the Marras et al. (2015) Euler-θ form
#       implemented in the body below.
#
#   ltheta = false (:energy_equation => "energy")
#       q = (ρ, ρu, ρv, ρE), the total-energy form of Nazarov &
#       Hoffman, Int. J. Numer. Meth. Fluids 71 (2013) 339-357,
#       eq. (3.4)-(3.7) — see _dsgs_2d_energy! below. This is the
#       variant to use for shock capturing: across a shock ρθ is not
#       conserved, so the θ system cannot carry the right shock speed
#       in the first place.
#
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 wres::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 PhysConst::PhysicalConst{TT},
                                 Pr::TT,
                                 nelem::Int, ngl::Int;
                                 ltheta::Bool=true,
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    if !ltheta
        # `wres` MUST be forwarded. It is the boundary mask (1 interior, 0 on a
        # domain boundary -- dsgs_wres in params_setup.jl) that keeps the
        # residual from reading the boundary condition instead of a
        # discretisation error, and _dsgs_2d_energy! takes it in position 7.
        # Omitting it does not shift the meaning of the later arguments
        # harmlessly: it is a MethodError at the first call, which is where
        # CompEuler/ffs_step died.
        _dsgs_2d_energy!(μ_dsgs, q, q1, q2, rhs, Minv, wres, visc_coeff,
                         Δt, connijk, Δelem, PhysConst, Pr, nelem, ngl,
                         lglobal_norms)
        return nothing
    end

    # Marras et al. (JCP 2015) eq. (8-10), implemented exactly as in
    # the lineage from fp/mymaster — the version that was already
    # known to run the rising-bubble case to completion. Residual is
    # the weak-form rhs[ip, i] directly (post-DSS, pre-mass-matrix
    # division); attempts to "correct" it with M⁻¹·rhs (the strong-
    # form residual) shrink the residual by ~10³ on 2D atmospheric
    # meshes and effectively turn DSGS off, which is not what the
    # algorithm was designed for in this lineage.
    #
    #     μ_res|e = C1 · Δ² · max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω
    #     μ_max|e = C2 · Δ · (|u| + c)_∞,e
    #     μ|e     = max(0, min(μ_max, μ_res))
    #
    # Per-equation split (Marras eq. 10), with the user-supplied
    # inputs[:μ] multiplier on each slot:
    #     μ_dsgs[iel, 1] = 0                              (no mass diffusion)
    #     μ_dsgs[iel, 2] = visc_coeff[2] · μ              (ρu)
    #     μ_dsgs[iel, 3] = visc_coeff[3] · μ              (ρv)
    #     μ_dsgs[iel, 4] = visc_coeff[4] · Pr/(γ-1) · μ   (ρθ)
    #
    # Minv and qe stay in the function-barrier signature so the rhs.jl
    # call site doesn't have to change, but they are unused here.

    invnp = one(TT)/(nelem*ngl*ngl)
    γ     = PhysConst.γ
    C0    = PhysConst.C0
    C1    = TT(1.0)
    C2    = TT(0.5)
    γm1   = γ - one(TT)
    eps   = TT(1.0e-16)

    # --- Pass 1: averages of (ρ, ρu, ρv, ρθ) — see _dsgs_norm_scope -----
    ρ_avg  = zero(TT); ρu_avg = zero(TT)
    ρv_avg = zero(TT); ρθ_avg = zero(TT)
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                ρ_avg  += q[ip,1]
                ρu_avg += q[ip,2]
                ρv_avg += q[ip,3]
                ρθ_avg += q[ip,4]
            end
        end
    end
    if lglobal_norms
        sums = TT[ρ_avg, ρu_avg, ρv_avg, ρθ_avg, TT(nelem*ngl*ngl)]
        MPI.Allreduce!(sums, MPI.SUM, get_mpi_comm())
        invnp_g = one(TT)/max(sums[5], one(TT))
        ρ_avg  = sums[1]*invnp_g; ρu_avg = sums[2]*invnp_g
        ρv_avg = sums[3]*invnp_g; ρθ_avg = sums[4]*invnp_g
    else
        ρ_avg  *= invnp; ρu_avg *= invnp
        ρv_avg *= invnp; ρθ_avg *= invnp
    end

    # --- Pass 2: L∞ norms of |q - ⟨q⟩| ---------------------------------
    denom1 = zero(TT); denom2 = zero(TT)
    denom3 = zero(TT); denom4 = zero(TT)
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                denom1 = max(denom1, abs(q[ip,1] - ρ_avg))
                denom2 = max(denom2, abs(q[ip,2] - ρu_avg))
                denom3 = max(denom3, abs(q[ip,3] - ρv_avg))
                denom4 = max(denom4, abs(q[ip,4] - ρθ_avg))
            end
        end
    end
    if lglobal_norms
        norms = TT[denom1, denom2, denom3, denom4]
        MPI.Allreduce!(norms, MPI.MAX, get_mpi_comm())
        denom1 = norms[1]; denom2 = norms[2]
        denom3 = norms[3]; denom4 = norms[4]
    end
    # Machine-zero floor on every denominator (Marras eq. 9 prescribes
    # ‖q − ⟨q⟩‖∞,Ω in the denominator; we add eps to guarantee a finite
    # ratio even before any spatial variation has developed).
    denom1 += eps; denom2 += eps
    denom3 += eps; denom4 += eps

    # The momentum slots need a slightly larger physical-scale floor:
    # at t = 0 the fluid is at rest globally, so ‖ρu − ⟨ρu⟩‖∞,Ω and
    # ‖ρv − ⟨ρv⟩‖∞,Ω literally start at zero. With only machine eps to
    # absorb that, the R/denom ratio runs away and caps μ at the
    # wave-speed bound C2·Δ·(|u|+c) before any flow has developed,
    # which on this case is enough to push ρθ past zero in the very
    # first RK substage. The floor is a tiny fraction (1e-3) of the
    # natural momentum scale ρ_avg·c_avg — large enough to keep the
    # cold-start ratio bounded, small enough to vanish once actual
    # momentum perturbations have grown above it.
    θ_avg  = ρθ_avg/max(abs(ρ_avg), eps)
    p_avg  = C0*(max(ρ_avg*θ_avg, zero(TT)))^γ
    c_avg  = sqrt(max(γ*p_avg/max(abs(ρ_avg), eps), zero(TT)))
    mom_floor = TT(1.0e-3) * abs(ρ_avg) * c_avg
    denom2 = max(denom2, mom_floor)
    denom3 = max(denom3, mom_floor)

    # --- Pass 3: per-element residual L∞, μ_max bound, μ_dsgs[ie] ------
    DSGS_MONITOR[] && fill!(DSGS_DBG, 0.0)
    @inbounds for ie = 1:nelem
        # Marras's element size: min(Δx, Δy)/(N+1). Δelem[ie] is the
        # min corner-to-corner distance in the element; ngl = N+1.
        Δ = Δelem[ie]/ngl

        n1   = zero(TT); n2 = zero(TT)
        n3   = zero(TT); n4 = zero(TT)
        uTmx = zero(TT)
        ρ_el = zero(TT)

        for j = 1:ngl
            @simd for i = 1:ngl
                ip = connijk[ie,i,j,1]
                Mi = Minv[ip]

                # Strong-form residual. rhs[] here is the DSS-assembled
                # WEAK-form RHS (rhs! divides by the mass matrix later), so
                # it must be multiplied by M⁻¹ for the difference to be
                # dimensionally meaningful: ∂q/∂t has units q/time, the raw
                # weak RHS has units (mass matrix)·q/time. The 1D path has
                # always done this; the 2D path did not.
                # w = 0 on a domain boundary; see dsgs_wres in params_setup.jl.
                w  = wres[ip]
                R1 = w*abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) - Mi*rhs[ip,1])
                R2 = w*abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) - Mi*rhs[ip,2])
                R3 = w*abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) - Mi*rhs[ip,3])
                R4 = w*abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])/(2*Δt) - Mi*rhs[ip,4])
                n1 = max(n1, R1); n2 = max(n2, R2)
                n3 = max(n3, R3); n4 = max(n4, R4)
                if DSGS_MONITOR[]
                    _dsgs_dbg!(R1/denom1, ip, 1, R1); _dsgs_dbg!(R2/denom2, ip, 2, R2)
                    _dsgs_dbg!(R3/denom3, ip, 3, R3); _dsgs_dbg!(R4/denom4, ip, 4, R4)
                end

                ρl = q[ip,1]
                ul = q[ip,2]/ρl
                vl = q[ip,3]/ρl
                θl = q[ip,4]/ρl
                # Equation of state p = C0·(ρθ)^γ  ⇒  c² = γp/ρ
                pl  = C0 * (ρl*θl)^γ
                c_l = sqrt(max(γ*pl/ρl, zero(TT)))
                uTmx = max(uTmx, sqrt(ul*ul + vl*vl) + c_l)
                ρ_el += ρl
            end
        end
        ρ_el /= TT(ngl*ngl)

        μ_res = C1*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3, n4/denom4)
        μ_max = C2*Δ*uTmx
        μ     = max(zero(TT), min(μ_max, μ_res))   # kinematic, m²/s

        # μ above is KINEMATIC. _expansion_visc! applies visc_coeff·∇²(prim)
        # and user_primitives! hands this system (ρ, u, v, θ), so momentum
        # and θ both need the DYNAMIC coefficient ρ̄·μ.
        μ_dyn = ρ_el*μ

        # Per-equation split (Marras eq. 10), scaled by the user-supplied
        # inputs[:μ] multiplier so the case can be run with DSGS off
        # (visc_coeff = [0,…]).
        μ_dsgs[ie,1] = zero(TT)                             # ρ : no mass diffusion
        μ_dsgs[ie,2] = visc_coeff[2] * μ_dyn                # ρu (eq. 10a)
        μ_dsgs[ie,3] = visc_coeff[3] * μ_dyn                # ρv (eq. 10a)
        μ_dsgs[ie,4] = visc_coeff[4] * (Pr/γm1) * μ_dyn     # ρθ (eq. 10b)
    end

    return nothing
end

# ================================================================================
# Residual-based artificial viscosity (DynSGS) — 2D, compressible Euler in
# TOTAL-ENERGY form q = (ρ, ρu, ρv, ρE).
#
#   M. Nazarov, J. Hoffman, "Residual-based artificial viscosity for
#   simulation of turbulent compressible flow using adaptive finite element
#   methods", Int. J. Numer. Meth. Fluids 71 (2013) 339-357.
#
# This is the shock-capturing variant. The Euler-θ version above transports
# ρθ, which is an entropy variable: it is constant across a contact but NOT
# conserved across a shock, so no amount of stabilization makes that system
# produce the right shock speed. Slot 4 here is ρE, the conserved total
# energy, and the viscosity is built from the residual of that system.
#
# Per element K, with a constant Δt and the BDF2 stencil over the three
# stored states (qⁿ, qⁿ⁻¹, qⁿ⁻²)  —  paper eq. (3.4):
#
#     R_ρ = (3ρⁿ − 4ρⁿ⁻¹ + ρⁿ⁻²)/(2Δt) + ∇·(ρu)
#     R_m = (3mⁿ − 4mⁿ⁻¹ + mⁿ⁻²)/(2Δt) + ∇·(m⊗u + pI)
#     R_E = (3Eⁿ − 4Eⁿ⁻¹ + Eⁿ⁻²)/(2Δt) + ∇·((E + p)u)
#
# The divergence terms are read off the assembled inviscid RHS: `rhs` is
# the DSS-assembled WEAK-form residual (rhs! divides by the mass matrix
# later), so it is multiplied by M⁻¹ here to get ∂q/∂t units — same
# convention as the 1D and MHD implementations in this file.
#
# Then eq. (3.5)-(3.7):
#
#     μ₁|K   = C1·h_K²·‖ρ−ρ̄‖_{∞,Ω}·max( ‖R_ρ‖_{∞,K}/‖ρ−ρ̄‖_{∞,Ω},
#                                        ‖R_m‖_{∞,K}/‖m−m̄‖_{∞,Ω},
#                                        ‖R_E‖_{∞,K}/‖E−Ē‖_{∞,Ω} )
#     μ_max|K = C2·h_K·‖ρ‖_{∞,K}·‖ |u| + √(γT) ‖_{∞,K}
#     μ|K     = min(μ_max|K, μ₁|K)
#     κ|K     = P/(γ−1)·μ|K            (heat conduction, on ∇T)
#     β|K     = μ|K/‖ρ‖_{∞,K}          (density diffusion, on ∇ρ)
#
# with C1 = 1, C2 = 0.5 and P ≈ 0.1 the artificial Prandtl number
# (inputs[:Pr]). NOTE that the leading ‖ρ−ρ̄‖_{∞,Ω} factor in μ₁ and the
# ‖ρ‖_{∞,K} factor in μ_max make μ a DYNAMIC viscosity, which is what
# _expansion_visc! wants for the momentum slots — so, unlike the θ path
# above, there is no separate ρ̄_el multiplication at the end.
#
# T is the paper's temperature, T = E/ρ − |u|²/2, i.e. the specific
# internal energy in the paper's cv = 1 scaling (p = (γ−1)ρT, eq. 2.3).
# It stays dimensionally consistent in SI: T = p/((γ−1)ρ) = cv·T_physical.
# The case's user_primitives! must therefore put THAT quantity in slot 4
# for the κ·∇T flux to match eq. (3.3) — see
# problems/CompEuler/ffs_step/user_primitives.jl.
#
# Unlike the θ path, ρ is NOT left undiffused: eq. (3.3) carries β∇ρ in
# the mass flux, and for shock capturing it is what keeps the density
# jump from ringing. The user's inputs[:μ][1] multiplier scales it and
# can switch it off with 0.0.
#
# ⟨q⟩ and ‖q−⟨q⟩‖ are rank-local unless :ldsgs_global_norms is set — see
# _dsgs_norm_scope above. In the default (rank-local) mode everything in this
# routine is allocation-free, same discipline as the other implementations
# here; the global mode allocates the two small reduction buffers, once per
# RHS call and not per node.
# ================================================================================
function _dsgs_2d_energy!(μ_dsgs::AbstractMatrix{TT},
                          q::AbstractMatrix{TT},
                          q1::AbstractMatrix{TT},
                          q2::AbstractMatrix{TT},
                          rhs::AbstractMatrix{TT},
                          Minv::AbstractVector{TT},
                          wres::AbstractVector{TT},
                          visc_coeff::AbstractVector{TT},
                          Δt::TT,
                          connijk::AbstractArray{TI,4},
                          Δelem::AbstractVector{TT},
                          PhysConst::PhysicalConst{TT},
                          Pr::TT,
                          nelem::Int, ngl::Int,
                          lglobal_norms::Bool) where {TT<:AbstractFloat, TI<:Integer}

    γ    = PhysConst.γ
    γm1  = γ - one(TT)
    C1   = TT(1.0)
    C2   = TT(0.5)
    eps  = TT(1.0e-16)

    # --- Pass 1: rank-local means ⟨ρ⟩, ⟨ρu⟩, ⟨ρv⟩, ⟨ρE⟩ ----------------
    ρ_avg = zero(TT); ρu_avg = zero(TT)
    ρv_avg = zero(TT); ρE_avg = zero(TT)
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                ρ_avg  += q[ip,1]
                ρu_avg += q[ip,2]
                ρv_avg += q[ip,3]
                ρE_avg += q[ip,4]
            end
        end
    end
    if lglobal_norms
        sums = TT[ρ_avg, ρu_avg, ρv_avg, ρE_avg, TT(nelem*ngl*ngl)]
        MPI.Allreduce!(sums, MPI.SUM, get_mpi_comm())
        inv_npts = one(TT)/max(sums[5], one(TT))
        ρ_avg  = sums[1]*inv_npts; ρu_avg = sums[2]*inv_npts
        ρv_avg = sums[3]*inv_npts; ρE_avg = sums[4]*inv_npts
    else
        inv_npts = one(TT)/max(TT(nelem*ngl*ngl), one(TT))
        ρ_avg  *= inv_npts; ρu_avg *= inv_npts
        ρv_avg *= inv_npts; ρE_avg *= inv_npts
    end

    # --- Pass 2: L∞ of |q − ⟨q⟩| ---------------------------------------
    #
    # The momentum norm is the one the paper writes, ‖m − m̄‖_{∞,Ω} on the
    # momentum VECTOR, not two independent per-component norms.
    dρ = zero(TT); dm = zero(TT); dE = zero(TT)
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip  = connijk[ie,i,j,1]
                du  = q[ip,2] - ρu_avg
                dv  = q[ip,3] - ρv_avg
                dρ  = max(dρ, abs(q[ip,1] - ρ_avg))
                dm  = max(dm, sqrt(du*du + dv*dv))
                dE  = max(dE, abs(q[ip,4] - ρE_avg))
            end
        end
    end
    if lglobal_norms
        norms = TT[dρ, dm, dE]
        MPI.Allreduce!(norms, MPI.MAX, get_mpi_comm())
        dρ = norms[1]; dm = norms[2]; dE = norms[3]
    end

    # Physical-scale floors. A uniform free stream — which is exactly the
    # t = 0 state of a shock-tube or a supersonic-inflow problem — has
    # ‖q−⟨q⟩‖_{∞,Ω} = 0 identically, and R/eps would then blow the ratio
    # up and pin μ at the μ_max cap over the whole domain before any
    # flow structure exists. Each denominator is floored at a small
    # fraction (1e-3) of that field's natural scale, built from the mean
    # state; the floor vanishes from the picture as soon as real
    # perturbations grow past it.
    ρ_ref = max(abs(ρ_avg), eps)
    p_avg = γm1*max(ρE_avg - TT(0.5)*(ρu_avg*ρu_avg + ρv_avg*ρv_avg)/ρ_ref, zero(TT))
    c_avg = sqrt(max(γ*p_avg/ρ_ref, eps))
    rel   = TT(1.0e-3)
    dρ = max(dρ, rel*ρ_ref)              + eps
    dm = max(dm, rel*ρ_ref*c_avg)        + eps
    dE = max(dE, rel*ρ_ref*c_avg*c_avg)  + eps

    # --- Pass 3: per-element residual L∞, wave-speed cap, split --------
    inv2Δt = one(TT)/(2*Δt)
    @inbounds for ie = 1:nelem

        # Marras's element length scale: min edge / (N+1). Δelem[ie] is
        # the min corner-to-corner distance in the element, ngl = N+1.
        h = Δelem[ie]/ngl

        ratio = zero(TT)   # max_i ‖R_i‖_{∞,K}/‖q_i − ⟨q_i⟩‖_{∞,Ω}
        wmax  = zero(TT)   # ‖ |u| + √(γT) ‖_{∞,K}
        ρmax  = zero(TT)   # ‖ρ‖_{∞,K}

        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                Mi = Minv[ip]

                # w = 0 on a domain boundary; see dsgs_wres in params_setup.jl.
                w   = wres[ip]
                Rρ  = w*abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
                Rmu = (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
                Rmv = (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
                Rm  = w*sqrt(Rmu*Rmu + Rmv*Rmv)
                RE  = w*abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])*inv2Δt - Mi*rhs[ip,4])

                ratio = max(ratio, Rρ/dρ, Rm/dm, RE/dE)

                ρl = max(q[ip,1], eps)
                ul = q[ip,2]/ρl
                vl = q[ip,3]/ρl
                # T = E/ρ − |u|²/2 (paper eq. 3.2, cv = 1 scaling), clamped
                # at zero so a transient negative internal energy in an
                # under-resolved cell cannot produce a NaN wave speed.
                Tl = max(q[ip,4]/ρl - TT(0.5)*(ul*ul + vl*vl), zero(TT))
                wmax = max(wmax, sqrt(ul*ul + vl*vl) + sqrt(γ*Tl))
                ρmax = max(ρmax, ρl)
            end
        end

        # eq. (3.5)-(3.7). Both branches carry a density, so μ is DYNAMIC.
        μ_res = C1*h*h*dρ*ratio
        μ_cap = C2*h*ρmax*wmax
        μ     = max(zero(TT), min(μ_cap, μ_res))

        μ_dsgs[ie,1] = visc_coeff[1] * μ/max(ρmax, eps)   # β on ∇ρ
        μ_dsgs[ie,2] = visc_coeff[2] * μ                  # μ on ∇u
        μ_dsgs[ie,3] = visc_coeff[3] * μ                  # μ on ∇v
        μ_dsgs[ie,4] = visc_coeff[4] * (Pr/γm1) * μ       # κ on ∇T
    end

    return nothing
end

# ================================================================================
# Marras-Nazarov Dynamic SGS (DynSGS) — 2D, ideal GLM-MHD (nine fields)
#
#   S. Marras, M. Nazarov, F. X. Giraldo, "Stabilized high-order Galerkin
#   methods based on a parameter-free dynamic SGS model for LES",
#   J. Comput. Phys. 301 (2015) 77-101.
#   M. Nazarov, J. Hoffman, "Residual-based artificial viscosity for
#   simulation of turbulent compressible flow using adaptive FE methods",
#   Int. J. Numer. Meth. Fluids 71 (2013) 339-357.
#
# The model is parameter-free in the sense that the eddy viscosity is set
# by the local residual of the governing equations rather than by a tuned
# constant: where the discrete solution satisfies the PDE the residual is
# small and so is the viscosity; at shocks and under-resolved features the
# residual spikes and viscosity appears exactly there. That is the whole
# point of using it here instead of Smagorinsky, whose  ρ Cs² Δ² |S|  is
# blind to whether the flow is resolved and therefore has to be scaled up
# globally (8x on this grid) to survive the shocks — which then over-damps
# the smooth 90% of the domain.
#
#     μ_res|e = C1 · Δ² · max_i ( ‖R_i‖_{∞,e} / ‖q_i − ⟨q_i⟩‖_{∞,Ω} )
#     μ_max|e = C2 · Δ  · (‖v‖ + c_f)_{∞,e}
#     μ|e     = max(0, min(μ_max, μ_res))
#
# with the BDF2 residual of equation i
#
#     R_i = (3qⁿ_i − 4qⁿ⁻¹_i + qⁿ⁻²_i)/(2Δt) − M⁻¹·RHS_i
#
# Notes specific to this implementation:
#
#  *  M⁻¹·RHS, not RHS.  R must have units of q/time for the ratio
#     R/‖q−⟨q⟩‖ to be a frequency and μ_res to come out as m²/s. params.RHS
#     inside viscous_rhs_el! is the DSS-assembled *weak-form* inviscid RHS
#     (the mass-matrix division happens later in rhs!), so it is multiplied
#     by Minv here — matching the 1D implementation.
#
#  *  Step-cadenced history.  params.qp.qnm1/qnm2 are advanced on every RK
#     *stage*, so they are stage snapshots, not states one Δt apart, and a
#     BDF2 stencil built on them does not approximate ∂q/∂t. DynSGS-MHD
#     therefore carries its own pair (params.dsgs_qnm1/qnm2), advanced once
#     per time step by rhs!.
#
#  *  Residual set.  The max runs over the eight genuine conservation laws
#     (ρ, ρu, ρv, E, ρw, Bx, By, Bz) and EXCLUDES ψ: the GLM field is a
#     numerical constraint carrier, not a conserved quantity, and its
#     residual is dominated by the Dedner damping source rather than by any
#     under-resolution of the flow.
#
#  *  Denominator floors.  ‖q_i − ⟨q_i⟩‖_{∞,Ω} is zero for any field that
#     is uniform (all of them at t=0) or identically zero (ρw and Bz for
#     Orszag-Tang, for all time). Each denominator is floored at a small
#     fraction of that field's natural physical scale, built from the
#     domain-mean state, so a degenerate field contributes 0/floor = 0
#     instead of 0/eps = garbage.
#
#  *  Units.  μ as defined above is KINEMATIC (m²/s). _expansion_visc!
#     applies visc_coeff·∇²(primitive) to each equation, so the momentum
#     and energy slots — whose primitives are u, v, w and T — need the
#     DYNAMIC coefficient ρ̄·μ, while the magnetic slots take μ directly
#     as a turbulent resistivity (∇²B already carries the right units).
#     ρ̄ is the element-mean density.
#
# Per-equation split (Marras eq. 10, adapted to the GLM-MHD field set):
#     [1] ρ  : 0                       — mass stays conservative
#     [2,3,5] ρu, ρv, ρw : ρ̄·μ
#     [4] E  : ρ̄·μ·γ/((γ−1)·Pr_t)      — see below
#     [6,7,8] B          : μ           — turbulent resistivity
#     [9] ψ              : μ
# each scaled by the user's inputs[:μ][ieq] multiplier.
#
# The energy factor: user_primitives! hands slot 4 the temperature
# T = p/ρ (= R·T_phys), and the physical flux is ∇·(k∇T_phys) with
# k = μ_dyn·cp/Pr_t. Rewriting in terms of T gives the coefficient
# k/R = μ_dyn·γ/((γ−1)·Pr_t), since cp = γR/(γ−1).
#
# ⟨q⟩ and ‖q−⟨q⟩‖ are rank-local unless :ldsgs_global_norms is set — see
# _dsgs_norm_scope above. `comm` is what the global mode reduces over.
# ================================================================================
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS_MHD, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 wres::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 avg::AbstractVector{TT},
                                 denom::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 γ::TT, Pr_t::TT, C1::TT, C2::TT,
                                 comm,
                                 nelem::Int, ngl::Int;
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)          # residual max excludes the ψ slot
    γm1  = γ - one(TT)
    eps  = TT(1.0e-16)

    # --- Pass 1: rank-local means ⟨q_i⟩ --------------------------------
    @inbounds for ieq = 1:neqs
        avg[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                for ieq = 1:neqs
                    avg[ieq] += q[ip,ieq]
                end
            end
        end
    end
    inv_npts = one(TT)/max(TT(nelem*ngl*ngl), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(nelem*ngl*ngl), MPI.SUM, comm)
        MPI.Allreduce!(avg, MPI.SUM, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end

    # --- Pass 2: L∞ of |q_i − ⟨q_i⟩| -----------------------------------
    @inbounds for ieq = 1:neqs
        denom[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem
        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                for ieq = 1:neqs
                    denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
                end
            end
        end
    end
    if lglobal_norms
        MPI.Allreduce!(denom, MPI.MAX, comm)
    end

    # Physical-scale floors (see header). Built from the mean state:
    #   ρ̄, the mean sound-ish speed c̄, and the mean field strength.
    ρ_avg = max(abs(avg[1]), eps)
    p_avg = γm1*max(avg[4] - TT(0.5)*(avg[2]*avg[2] + avg[3]*avg[3] + avg[5]*avg[5])/ρ_avg
                    - TT(0.5)*(avg[6]*avg[6] + avg[7]*avg[7] + avg[8]*avg[8]), zero(TT))
    c_avg = sqrt(max(γ*p_avg/ρ_avg, eps))
    rel   = TT(1.0e-3)
    @inbounds begin
        denom[1] = max(denom[1], rel*ρ_avg)                 # ρ
        mom_fl   = rel*ρ_avg*c_avg
        denom[2] = max(denom[2], mom_fl)                    # ρu
        denom[3] = max(denom[3], mom_fl)                    # ρv
        denom[4] = max(denom[4], rel*ρ_avg*c_avg*c_avg)     # E
        if neqs >= 5; denom[5] = max(denom[5], mom_fl); end # ρw
        b_fl = rel*sqrt(ρ_avg)*c_avg
        for ieq = 6:min(neqs,8)
            denom[ieq] = max(denom[ieq], b_fl)              # B
        end
        for ieq = 1:neqs
            denom[ieq] += eps
        end
    end

    # --- Pass 3: per-element residual L∞, wave-speed cap, split --------
    inv2Δt = one(TT)/(2*Δt)
    @inbounds for ie = 1:nelem

        # Marras's element length scale: min edge / (N+1).
        Δ = Δelem[ie]/ngl

        ratio = zero(TT)      # max_i ‖R_i‖∞,e / denom_i
        wmax  = zero(TT)      # (‖v‖ + c_f)∞,e
        ρ_el  = zero(TT)      # element-mean density

        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                Mi = Minv[ip]

                for ieq = 1:NRES
                    # wres[ip] = 0 on a domain boundary; see params_setup.jl.
                    R = wres[ip]*abs((3*q[ip,ieq] - 4*q1[ip,ieq] + q2[ip,ieq])*inv2Δt - Mi*rhs[ip,ieq])
                    r = R/denom[ieq]
                    ratio = max(ratio, r)
                end

                ρl = max(q[ip,1], eps)
                ul = q[ip,2]/ρl
                vl = q[ip,3]/ρl
                wl = (neqs >= 5) ? q[ip,5]/ρl : zero(TT)
                B2 = q[ip,6]*q[ip,6] + q[ip,7]*q[ip,7] + q[ip,8]*q[ip,8]
                ψl = (neqs >= 9) ? q[ip,9] : zero(TT)
                pl = γm1*(q[ip,4] - TT(0.5)*ρl*(ul*ul + vl*vl + wl*wl)
                          - TT(0.5)*B2 - TT(0.5)*ψl*ψl)
                # c_f ≤ sqrt(a² + b²): the fast magnetosonic speed bounded
                # over all propagation directions (a² = γp/ρ, b² = |B|²/ρ).
                cf = sqrt(max(γ*pl/ρl + B2/ρl, zero(TT)))
                wmax  = max(wmax, sqrt(ul*ul + vl*vl + wl*wl) + cf)
                ρ_el += ρl
            end
        end
        ρ_el /= TT(ngl*ngl)

        μ_res = C1*Δ*Δ*ratio
        μ_max = C2*Δ*wmax
        μ     = max(zero(TT), min(μ_max, μ_res))    # kinematic, m²/s

        μ_dyn = ρ_el*μ                              # dynamic, for u/v/w/T

        μ_dsgs[ie,1] = zero(TT)                                    # ρ
        μ_dsgs[ie,2] = visc_coeff[2]*μ_dyn                         # ρu
        μ_dsgs[ie,3] = visc_coeff[3]*μ_dyn                         # ρv
        μ_dsgs[ie,4] = visc_coeff[4]*μ_dyn*γ/(γm1*Pr_t)            # E
        if neqs >= 5
            μ_dsgs[ie,5] = visc_coeff[5]*μ_dyn                     # ρw
        end
        for ieq = 6:min(neqs,8)
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μ                     # B (resistivity)
        end
        if neqs >= 9
            μ_dsgs[ie,9] = visc_coeff[9]*μ                         # ψ
        end
    end

    return nothing
end

# The BDF2 numerator, shared by both 3D branches.
#
# 3qⁿ - 4qⁿ⁻¹ + qⁿ⁻² is algebraically identical to 3(qⁿ - qⁿ⁻¹) - (qⁿ⁻¹ - qⁿ⁻²)
# and is not the same rounding: the direct form builds 3q and 4q¹ and cancels
# them, so its error is eps*|q| taken on the TOTAL state, which on an
# atmospheric column is ρθ ~ 360 rather than the perturbation ~ 1. Divided by
# 2Δt and by a floored denominator and multiplied by Δ², that noise is a
# spurious eddy viscosity on a fluid at rest. Differencing consecutive states
# gives an exact zero on a steady state instead.
@inline dsgs_bdf2(q, q1, q2, ip, ieq) =
    3*(q[ip,ieq] - q1[ip,ieq]) - (q1[ip,ieq] - q2[ip,ieq])

# ================================================================================
# DynSGS — 3D. TWO EQUATION SETS, TWO DIFFERENT MODELS.
#
# `:energy_equation` decides which, exactly as it does in 2D, and the two are
# not variants of each other -- they differ in the equation of state, in the
# wave speed, in what the normalising norms are taken on, in whether the
# coefficient comes out kinematic or dynamic, and in whether the mass equation
# is diffused at all. Getting the branch wrong is not a small error: the
# theta-form closure p = C0(ρθ)^γ evaluated on ρE ~ 2.5e5 J/m³ returns a
# "pressure" of order 1e14 Pa and a cap to match.
#
#   :energy_equation => "theta"    q = (ρ, ρu, ρv, ρw, ρθ [, ρqₓ …])
#                                  _dsgs_3d_theta!  — Marras, Nazarov &
#                                  Giraldo, JCP 301 (2015), eqs. (8)-(10).
#
#   :energy_equation => "energy"   q = (ρ, ρu, ρv, ρw, ρE [, ρqₓ …])
#                                  _dsgs_3d_energy! — Nazarov & Hoffman,
#                                  IJNMF 71 (2013), eqs. (3.4)-(3.7). The 3D
#                                  sibling of _dsgs_2d_energy! above, which is
#                                  what CompEuler/ffs_step runs.
#
# WHAT THE TWO SHARE. The residual itself (BDF2 on the step-cadenced history
# against M⁻¹·RHS), the per-element L∞, the min(cap, residual) structure, the
# filter width Δelem_filter/nop, and the decision to take the normalising
# norms on q − qe rather than on the total field (see the header of
# _dsgs_3d_theta! for why that matters on a stratified column, and note it
# degrades to the total field automatically when qe is uniform or zero, which
# is the shock-capturing case).
#
# WHAT THEY DO NOT SHARE
#
#                        theta form                energy form
#   EOS                  p = C₀(ρθ)^γ             p = (γ-1)(E - ½ρ|v|²)
#   wave speed           √(γp/ρ)                  √(γT), T = E/ρ - ½|v|²
#   momentum norm        per component            on the VECTOR ‖m - m̄‖
#   μ_res                C₁Δ²·max(Rᵢ/dᵢ)          C₁Δ²·‖ρ-ρ̄‖·max(Rᵢ/dᵢ)
#   μ_cap                C₂Δ(‖v‖+c)               C₂Δ‖ρ‖(‖v‖+c)
#   coefficient is       KINEMATIC → ×ρ̄_el        already DYNAMIC
#   mass equation        undiffused (Marras 10)   carries β∇ρ (N&H 3.3)
#   energy scale floor   1e-3·ρ̄·θ̄                 1e-3·ρ̄·c̄²
#
# The mass-diffusion difference is not a detail on the shock-capturing side:
# ffs_step's own deck records `:μ[1] = 0` as the single most destabilising
# change measured on that case.
#
# WHERE THE 3D SPLIT DIFFERS FROM THE 2D ONE, AND WHY IT HAS TO
# -------------------------------------------------------------
# In 2D the per-equation coefficients go straight into visc_coeff_dsgs, so
# each slot can carry its own formula. In 3D the coefficient reaches the RHS
# as ONE per-node scalar, sgs.μ_turb, and the split is whatever
# SGS_diffusion(::AbstractSGSModel, ::NSD_3D) makes of it — which is the
# standard LES one, Pr_t on the energy slot and Sc_t on the others. That is
# the price of going through the closure struct, and the closure struct is
# what buys the full stress tensor, :implicit_vdiff and the LES statistics
# (see SGS_DSGS in sgsStructs.jl).
#
# For the theta form the two agree: Pr_t IS the right split for a potential
# temperature. For the energy form they do not:
#
#   Nazarov & Hoffman eq. (3.7)      SGS_diffusion, 3D
#   κ = P/(γ-1)·μ,  P ≈ 0.1          κ = μ/Pr_t
#   β = μ/‖ρ‖_{∞,K}                  μ/Sc_t
#
# so a 3D shock-capturing run that wants the paper's coefficients sets them
# through the per-equation mask:
#
#   :μ[5] = Pr_t·P/(γ-1) = 0.7·0.1/0.4 = 0.175      to get κ = P/(γ-1)·μ
#   :μ[1] = Sc_t/ρ_typical                          to approximate β
#
# The second is only exact for a uniform ρ, which is why it is written here
# rather than applied silently: on a strongly compressed field β and μ/Sc_t
# differ by the local density ratio. If that matters for a case, the honest
# fix is a second per-node coefficient array, not a mask.
# ================================================================================
function compute_dsgs_viscosity!(sgs::SGS_DSGS{TT},
                                 ::NSD_3D,
                                 μ_dsgs::AbstractMatrix{TT},
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 wres::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem_filter::AbstractVector{TT},
                                 Δfallback::TT,
                                 nop::Int,
                                 PhysConst::PhysicalConst{TT},
                                 comm,
                                 nelem::Int, ngl::Int, neqs::Int;
                                 lpert::Bool=false,
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    if sgs.ltheta_eqn
        _dsgs_3d_theta!(sgs, μ_dsgs, q, q1, q2, qe, rhs, Minv, wres, visc_coeff,
                        Δt, connijk, Δelem_filter, Δfallback, nop, PhysConst,
                        comm, nelem, ngl, neqs; lpert = lpert,
                        lglobal_norms = lglobal_norms)
    else
        _dsgs_3d_energy!(sgs, μ_dsgs, q, q1, q2, qe, rhs, Minv, wres, visc_coeff,
                         Δt, connijk, Δelem_filter, Δfallback, nop, PhysConst,
                         comm, nelem, ngl, neqs; lpert = lpert,
                         lglobal_norms = lglobal_norms)
    end
    return nothing
end

# ================================================================================
# DynSGS — 3D, compressible Euler in POTENTIAL-TEMPERATURE form
#          q = (ρ, ρu, ρv, ρw, ρθ [, ρq_x ...])
#
#   S. Marras, M. Nazarov, F. X. Giraldo, "Stabilized high-order Galerkin
#   methods based on a parameter-free dynamic SGS model for LES",
#   J. Comput. Phys. 301 (2015) 77-101, eqs. (8)-(10).
#
# This is the 3D sibling of compute_dsgs_viscosity!(::DSGS, ::NSD_2D) above.
# Three things are genuinely different, and each of them is a correctness
# requirement rather than a preference.
#
# 1. IT WRITES INTO A CLOSURE STRUCT, NOT INTO visc_coeff_dsgs.
#
#    The 2D path packs one number per element into visc_coeff_dsgs[ieq] and
#    lets the scalar `visc_coeff[ieq]*grad(primitive)` kernel apply it. The 3D
#    viscous kernel is not that: it is the full stress tensor
#    tau_ij = mu(du_i/dx_j + du_j/dx_i) - (2/3) mu div(u) delta_ij, and it is
#    reached only through the `sgs::AbstractSGSModel` method of
#    _expansion_visc!, which reads SGS_diffusion(..., sgs, ...). Filling
#    sgs.mu_turb therefore gives DynSGS the correct 3D stress tensor for free,
#    and with it a truthful viscous row in the CFL table, a working
#    :implicit_vdiff, and correct subfilter stresses in the LES statistics --
#    all three read the same struct. See the note on SGS_DSGS in sgsStructs.jl.
#
# 2. THE NORMALISING SCALES ARE TAKEN ON THE PERTURBATION, NOT THE TOTAL FIELD.
#
#    Eq. (9) divides each residual by ‖q_i - <q_i>‖_inf over the domain. On the
#    2D rising bubble the background is nearly uniform and it makes little
#    difference which state that norm is taken on. On a 5 km stratified column
#    it decides whether the model does anything at all:
#
#        field    ‖q - <q>‖_inf on the TOTAL state   turbulent scale
#        rho      ~3.5e-1 kg/m^3  (hydrostatic)      ~1e-3
#        rho*th   ~6.5e+1 K kg/m^3 (the sounding)    ~5e-1
#
#    i.e. the hydrostatic background is 100-300x the turbulence the sensor is
#    supposed to see, so a denominator built on the total state divides the
#    residual by the wrong number by two orders of magnitude and switches
#    DynSGS off exactly where an LES needs it. The norms here are therefore
#    taken on q - qe, the departure from the reference sounding, which is the
#    quantity whose spread is the flow's own scale. The BDF2 numerator is
#    unaffected: 3qe - 4qe + qe = 0, so subtracting a time-independent
#    reference state changes nothing there.
#
# 3. THE FILTER WIDTH IS THE ONE THE REST OF THE 3D LES MACHINERY USES.
#
#    Delta = Delta_elem_filter[iel]/nop, i.e. exactly what SMAG and VREM get
#    from compute_element_size_driver (:les_filter_width, default :max) and
#    what les_statistics.jl reports against. The 2D path uses Delta_elem/ngl,
#    the SHORTEST corner-to-corner distance -- on a 160 x 160 x 40 m element
#    that is 40/5 = 8 m against 160/4 = 40 m here, a factor 25 in mu (Delta
#    enters squared). Using two different filter widths for the closure and for
#    the diagnostic that reports it is how a run ends up unable to explain its
#    own output.
#
# DENOMINATOR FLOORS. Every field of this problem starts horizontally uniform,
# so ‖q'_i - <q'_i>‖_inf is exactly zero at t = 0 for every slot and the ratio
# R_i/denom_i would be garbage rather than large. Each denominator is floored at
# 1e-3 of that field's natural scale, built from the domain-mean TOTAL state:
#
#     rho                1e-3 * rhobar
#     rho u, rho v, rho w 1e-3 * rhobar * cbar
#     rho theta          1e-3 * rhobar * thetabar
#     other scalars      1e-3 * <|q_i|>
#
# A degenerate field then contributes 0/floor = 0 rather than 0/eps.
#
# THE CAP IS LOOSE IN A LOW-MACH FLOW, AND THAT IS NOT A BUG BUT IT IS WORTH
# KNOWING. mu_max = C2*Delta*(‖v‖+c) is the first-order-upwind viscosity, and in
# an atmosphere c ~ 340 m/s against ‖v‖ ~ 10, so at Delta = 40 m it is
# 0.5*40*350 = 7000 m^2/s -- two to three orders above anything a PBL closure
# produces. min(mu_max, mu_res) therefore never binds here and mu_res governs
# alone, which is the regime the model is designed for; the cap is doing its
# job as a safety net for the shock case it was written for, not as a limiter
# on this one. Lower :dsgs_C2 if you want it to bite.
#
# MPI. <q_i> and ‖q_i - <q_i>‖ are DOMAIN norms by definition, so with
# :ldsgs_global_norms they are Allreduced -- two small collectives per RHS
# call. Rank-local (the default) costs no communication but makes the eddy
# viscosity depend on the partitioning, which is fine for a production run and
# not fine for a rank-count study.
# ================================================================================
function _dsgs_3d_theta!(sgs::SGS_DSGS{TT},
                                 μ_dsgs::AbstractMatrix{TT},
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 wres::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem_filter::AbstractVector{TT},
                                 Δfallback::TT,
                                 nop::Int,
                                 PhysConst::PhysicalConst{TT},
                                 comm,
                                 nelem::Int, ngl::Int, neqs::Int;
                                 lpert::Bool=false,
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    avg   = sgs.avg
    denom = sgs.denom
    scale = sgs.scale
    ν_el  = sgs.ν_el
    μ_el  = sgs.μ_el

    C1  = sgs.C1
    C2  = sgs.C2
    γ   = PhysConst.γ
    C0  = PhysConst.C0
    eps = TT(1.0e-16)

    inv2Δt = one(TT)/(TT(2)*Δt)
    npts   = nelem*ngl*ngl*ngl

    # `Δelem_filter` is per element; `Δfallback` (mesh.Δeffective_l) is the one
    # number for the whole domain kept for meshes built before that array
    # existed. Same test as _viscous_rhs_el_3d!, so the two cannot drift.
    luse_local_Δ = length(Δelem_filter) == nelem && nop > 0

    #--- Pass 1: means of the perturbation, and of |total| for the floors ----
    @inbounds for ieq = 1:neqs
        avg[ieq]   = zero(TT)
        scale[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,k]
        for ieq = 1:neqs
            # ⟨q⟩ and ‖q-⟨q⟩‖ are taken on the STORED PROGNOSTIC VARIABLE, with
            # no qe arithmetic -- eq. (9) to the letter, and the same as every
            # other path in this file. `scale` carries the mean TOTAL state
            # instead, because the FLOORS below are physical scales (ρ̄, c̄) and
            # those are not defined on a perturbation.
            avg[ieq]   += q[ip,ieq]
            scale[ieq] += lpert ? q[ip,ieq] + qe[ip,ieq] : q[ip,ieq]
        end
    end
    if lglobal_norms
        MPI.Allreduce!(avg,   MPI.SUM, comm)
        MPI.Allreduce!(scale, MPI.SUM, comm)
        npts_g = MPI.Allreduce(npts, MPI.SUM, comm)
        invn   = one(TT)/max(TT(npts_g), one(TT))
    else
        invn = one(TT)/max(TT(npts), one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq]   *= invn
        scale[ieq] *= invn
    end

    #--- Pass 2: L∞ of |q' - <q'>| ------------------------------------------
    @inbounds for ieq = 1:neqs
        denom[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,k]
        for ieq = 1:neqs
            denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
        end
    end
    lglobal_norms && MPI.Allreduce!(denom, MPI.MAX, comm)

    #--- Denominator floors, built from the mean TOTAL state ----------------
    # `raw` keeps the measured norm and `floors` the floor, so the gate below
    # can tell "this field has structure" from "this field is still flat".
    raw    = sgs.raw
    floors = sgs.floors
    @inbounds for ieq = 1:neqs
        raw[ieq] = denom[ieq]
    end
    ρ̄  = max(abs(scale[1]), eps)
    θ̄  = neqs >= 5 ? max(abs(scale[5]), eps)/ρ̄ : TT(300)
    p̄  = C0*(max(ρ̄*θ̄, zero(TT)))^γ
    c̄  = sqrt(max(γ*p̄/ρ̄, zero(TT)))
    FLOOR = TT(1.0e-3)
    @inbounds begin
        floors[1] = FLOOR*ρ̄
        for ieq = 2:min(neqs,4); floors[ieq] = FLOOR*ρ̄*c̄; end
        neqs >= 5 && (floors[5] = FLOOR*ρ̄*θ̄)
        for ieq = 6:neqs; floors[ieq] = FLOOR*max(abs(scale[ieq]), eps); end
        denom[1] = max(denom[1] + eps, FLOOR*ρ̄)
        for ieq = 2:min(neqs,4)
            denom[ieq] = max(denom[ieq] + eps, FLOOR*ρ̄*c̄)
        end
        if neqs >= 5
            denom[5] = max(denom[5] + eps, FLOOR*ρ̄*θ̄)
        end
        for ieq = 6:neqs
            denom[ieq] = max(denom[ieq] + eps, FLOOR*max(abs(scale[ieq]), eps))
        end
    end

    #--- A DENOMINATOR THAT IS STILL ITS FLOOR MEANS THAT FIELD HAS NO SCALE
    #    YET, AND ITS RESIDUAL MUST NOT DRIVE THE COEFFICIENT.
    #
    # The floors exist so that a field with ‖q-⟨q⟩‖ = 0 gives 0/floor = 0
    # instead of 0/eps. They are physical scales, and for the momenta that
    # scale is ACOUSTIC: 1e-3·ρ̄·c̄ = 0.25 kg/(m²s) here, i.e. |u| = 0.33 m/s in
    # an atmosphere whose flow speeds are 10 m/s. A bubble starting from rest
    # therefore spends its first seconds dividing a real vertical-momentum
    # tendency (buoyancy, ρgΔθ/θ ≈ 0.065) by a floor 40x below the scale the
    # flow will actually reach, and the ratio comes out O(1) 1/s -- which at
    # Δ² = 62500 m² is μ ≈ 4e4, i.e. straight into the C₂Δ(‖v‖+c) cap, before
    # any flow exists at all. Measured: that is what killed rtb3d_dsgs at
    # t = 1.5 s, on eq. 4, at interior nodes, with denom[4] sitting exactly on
    # its floor of 0.246.
    #
    # The fix is not a bigger floor -- that is tuning. It is to use the floor
    # as a GATE rather than as a value: while a denominator is still the floor,
    # that equation has no resolved scale for the residual to be measured
    # against, so it is excluded from the max. Setting it to Inf does that
    # branch-free, and makes it visible in the monitor. Each equation rejoins
    # by itself the moment its own field develops structure above the floor
    # (for the momenta here, within a few seconds).
    #
    # If every denominator is at its floor -- the t = 0 state of any case that
    # starts uniform -- the ratio is 0 and the model is off, which is correct:
    # there is nothing yet to stabilise.
    #
    # NOT applied to the 1D or 2D paths. 1D carries no floors at all, and the
    # 2D ones are what CompEuler/sod1d, theta_dsgs and ffs_step were tuned
    # against; the same change is very likely right for them too, but that is a
    # measurement to make on those cases, not one to make blind.
    @inbounds for ieq = 1:neqs
        raw[ieq] <= floors[ieq] && (denom[ieq] = TT(Inf))
    end

    #--- Pass 3: per-element residual L∞, wave-speed cap, coefficient -------
    DSGS_MONITOR[] && fill!(DSGS_DBG, 0.0)
    @inbounds for ie = 1:nelem
        Δ  = luse_local_Δ ? Δelem_filter[ie]/nop : Δfallback
        Δ2 = Δ*Δ

        ratio = zero(TT)     # max_i ‖R_i‖_{∞,e} / denom_i
        uTmx  = zero(TT)     # (‖v‖ + c)_{∞,e}
        ρ_el  = zero(TT)

        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,k]
            Mi = Minv[ip]

            # Strong-form residual. params.RHS holds the DSS-assembled WEAK
            # form at this point in _build_rhs! (the mass-matrix division
            # happens at the end), so it MUST be premultiplied by M^-1 for the
            # subtraction against dq/dt to be dimensionally a residual at all.
            #
            # THE BDF2 STENCIL IS EVALUATED AS 3(q - q1) - (q1 - q2), WHICH IS
            # THE SAME NUMBER AND NOT THE SAME ROUNDING.
            #
            #     3q - 4q1 + q2  =  3(q - q1) - (q1 - q2)
            #
            # exactly, but the left form forms 3q and 4q1 first and then
            # cancels them, so its rounding error is eps*|q| -- taken on the
            # TOTAL state, which on this problem is rho*theta ~ 360 and not the
            # perturbation ~ 1. Divided by 2*dt and by a FLOORED denominator
            # (1e-3*rho*theta) and multiplied by Delta^2 = 1600 m^2, that noise
            # comes out as a spurious eddy viscosity on a state that is exactly
            # steady:
            #
            #     3*1.2 - 4*1.2 + 1.2 = -4.4e-16    (not 0)
            #     -> nu ~ 6e-10 m^2/s on a fluid at rest
            #
            # Harmless in magnitude, and still wrong in kind: this model's
            # entire claim is that an exact solution gets ZERO viscosity, and a
            # sensor with a floor under it should not be reporting the floating
            # point representation of the background state. The right form
            # differences CONSECUTIVE states, whose difference is O(dq), so a
            # steady state gives an exact zero and the noise on a moving one
            # scales with the change rather than with the state. Measured on
            # the same case: 6e-10 -> 0.
            for ieq = 1:neqs
                dq = TT(3)*(q[ip,ieq] - q1[ip,ieq]) - (q1[ip,ieq] - q2[ip,ieq])
                R  = wres[ip]*abs(dq*inv2Δt - Mi*rhs[ip,ieq])
                r  = R/denom[ieq]
                ratio = max(ratio, r)
                # Debug record, opt-in. The branch is on a Ref{Bool} that does
                # not change during a run, so it predicts perfectly.
                DSGS_MONITOR[] && _dsgs_dbg!(r, ip, ieq, R)
            end

            ρl = lpert ? q[ip,1] + qe[ip,1] : q[ip,1]
            ρl = ρl > zero(TT) ? ρl : eps
            ul = q[ip,2]/ρl
            vl = q[ip,3]/ρl
            wl = q[ip,4]/ρl
            θl = lpert ? (q[ip,5] + qe[ip,5])/ρl : q[ip,5]/ρl
            pl = C0*(max(ρl*θl, zero(TT)))^γ
            cl = sqrt(max(γ*pl/ρl, zero(TT)))
            uTmx = max(uTmx, sqrt(ul*ul + vl*vl + wl*wl) + cl)
            ρ_el += ρl
        end
        ρ_el /= TT(ngl*ngl*ngl)

        μ_res = C1*Δ2*ratio          # kinematic, m²/s
        μ_max = C2*Δ*uTmx            # kinematic, m²/s
        ν     = max(zero(TT), min(μ_max, μ_res))

        ν_el[ie] = ν
        # DYNAMIC from here on: SGS_diffusion(::AbstractSGSModel, ::NSD_3D)
        # returns μ_turb straight into the stress tensor, whose primitives are
        # u, v, w and θ — so the coefficient it wants is ρ·ν, exactly as for
        # Smagorinsky (which builds μ_turb = ρ ℓ² |S|).
        μ_el[ie] = ρ_el*ν

        # The per-equation μ_dsgs matrix is DIAGNOSTIC ONLY here: it is what
        # broadcast_dsgs_to_nodes! writes to VTK as mu_dsgs_<var>. What the RHS
        # actually applies comes from μ_el via SGS_diffusion, so these entries
        # mirror that split rather than defining it —
        #     mass     visc_coeff[1]·μ/Sc_t   (0 with the usual :μ[1] = 0)
        #     momentum visc_coeff[ieq]·μ
        #     θ        visc_coeff[5]·μ/Pr_t
        # exactly the three branches of SGS_diffusion.
        μd = μ_el[ie]
        μ_dsgs[ie,1] = visc_coeff[1]*μd/sgs.Sc_t
        for ieq = 2:min(neqs,4)
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μd
        end
        if neqs >= 5
            μ_dsgs[ie,5] = visc_coeff[5]*μd/sgs.Pr_t
        end
        for ieq = 6:neqs
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μd/sgs.Sc_t
        end
    end

    return nothing
end

# ================================================================================
# DynSGS — 3D, compressible Euler/Navier-Stokes in TOTAL-ENERGY form
#          q = (ρ, ρu, ρv, ρw, ρE [, ρq_x ...])
#
#   M. Nazarov, J. Hoffman, "Residual-based artificial viscosity for
#   simulation of turbulent compressible flow using adaptive finite element
#   methods", Int. J. Numer. Meth. Fluids 71 (2013) 339-357, eqs. (3.4)-(3.7).
#
# The 3D sibling of _dsgs_2d_energy! above -- the path CompEuler/ffs_step
# runs -- and it is transcribed from it term for term rather than adapted from
# the theta-form 3D path, because the two models differ in more than a
# dimension. What changes going from 2D to 3D is only that the momentum vector
# has three components; every formula is the 2D one.
#
#     R_ρ = (3ρⁿ - 4ρⁿ⁻¹ + ρⁿ⁻²)/(2Δt) - M⁻¹RHS_ρ
#     R_m = the same on (ρu, ρv, ρw), taken as a VECTOR norm
#     R_E = the same on ρE
#
#     μ|K   = min( C₂·h·‖ρ‖_{∞,K}·‖ ‖v‖+√(γT) ‖_{∞,K} ,
#                  C₁·h²·‖ρ-ρ̄‖_{∞,Ω}·max_i ‖R_i‖_{∞,K}/‖q_i-⟨q_i⟩‖_{∞,Ω} )
#
# with T = E/ρ - ½‖v‖² the paper's temperature in its cv = 1 scaling
# (p = (γ-1)ρT, eq. 2.3), so c = √(γT) is the sound speed and stays
# dimensionally consistent in SI.
#
# NOTE μ IS ALREADY DYNAMIC HERE. Both branches of the min carry a density --
# the leading ‖ρ-ρ̄‖ in μ_res and the ‖ρ‖_{∞,K} in μ_cap -- so unlike the
# theta path there is NO ρ̄_el multiplication at the end. Applying one would
# be a factor of ρ, which at 1.2 kg/m³ is quiet and at the 5x compression
# behind a Mach-3 bow shock is not.
#
# THE MOMENTUM NORM IS ON THE VECTOR, not per component: ‖m - m̄‖_{∞,Ω} and
# ‖R_m‖_{∞,K}, exactly as the paper writes them and as _dsgs_2d_energy! does.
# Per-component norms would make the sensor depend on the orientation of the
# mesh relative to the flow, which for a shock at an angle to the grid is
# precisely the case that matters.
#
# NORMALISATION ON q - qe. Shared with the theta path, and it costs nothing
# here: a shock-capturing case has a uniform or zero reference state, so
# q - qe IS q and this reduces to the paper's expression. It is kept because a
# 3D total-energy ATMOSPHERE (:energy_equation => "energy" with a sounding,
# which is what CompEuler/giga_les_MOST is) has the same 100-300x hydrostatic
# background the theta path has to divide out.
#
# The per-equation split is SGS_diffusion's, not eq. (3.7)'s -- see the
# dispatcher above for the two mask values that recover the paper's.
# ================================================================================
function _dsgs_3d_energy!(sgs::SGS_DSGS{TT},
                          μ_dsgs::AbstractMatrix{TT},
                          q::AbstractMatrix{TT},
                          q1::AbstractMatrix{TT},
                          q2::AbstractMatrix{TT},
                          qe::AbstractMatrix{TT},
                          rhs::AbstractMatrix{TT},
                          Minv::AbstractVector{TT},
                          wres::AbstractVector{TT},
                          visc_coeff::AbstractVector{TT},
                          Δt::TT,
                          connijk::AbstractArray{TI,4},
                          Δelem_filter::AbstractVector{TT},
                          Δfallback::TT,
                          nop::Int,
                          PhysConst::PhysicalConst{TT},
                          comm,
                          nelem::Int, ngl::Int, neqs::Int;
                          lpert::Bool=false,
                          lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    avg   = sgs.avg
    denom = sgs.denom
    scale = sgs.scale
    ν_el  = sgs.ν_el
    μ_el  = sgs.μ_el

    C1  = sgs.C1
    C2  = sgs.C2
    γ   = PhysConst.γ
    γm1 = γ - one(TT)
    eps = TT(1.0e-16)

    inv2Δt = one(TT)/(TT(2)*Δt)
    npts   = nelem*ngl*ngl*ngl

    luse_local_Δ = length(Δelem_filter) == nelem && nop > 0

    #--- Pass 1: means of the perturbation, and signed means of the total ----
    @inbounds for ieq = 1:neqs
        avg[ieq]   = zero(TT)
        scale[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,k]
        for ieq = 1:neqs
            # See _dsgs_3d_theta!: the norms are on the stored prognostic
            # variable, `scale` on the total, because the floors are physical.
            avg[ieq]   += q[ip,ieq]
            scale[ieq] += lpert ? q[ip,ieq] + qe[ip,ieq] : q[ip,ieq]
        end
    end
    if lglobal_norms
        MPI.Allreduce!(avg,   MPI.SUM, comm)
        MPI.Allreduce!(scale, MPI.SUM, comm)
        npts_g = MPI.Allreduce(npts, MPI.SUM, comm)
        invn   = one(TT)/max(TT(npts_g), one(TT))
    else
        invn = one(TT)/max(TT(npts), one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq]   *= invn
        scale[ieq] *= invn
    end

    #--- Pass 2: L∞ norms. Slot 2 carries the MOMENTUM VECTOR norm for all
    #    three components; slots 3 and 4 are left at zero and never read, so
    #    the ratio below cannot pick up a per-component norm by accident.
    @inbounds for ieq = 1:neqs
        denom[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, k = 1:ngl, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,k]
        du = q[ip,2] - avg[2]
        dv = q[ip,3] - avg[3]
        dw = q[ip,4] - avg[4]
        denom[2] = max(denom[2], sqrt(du*du + dv*dv + dw*dw))
        for ieq in (1, 5)
            ieq <= neqs || continue
            denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
        end
        for ieq = 6:neqs
            denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
        end
    end
    lglobal_norms && MPI.Allreduce!(denom, MPI.MAX, comm)

    #--- Floors, from the mean TOTAL state -----------------------------------
    #
    # A uniform free stream -- the t = 0 state of a shock tube or a supersonic
    # inflow -- has ‖q - ⟨q⟩‖ = 0 identically, and R/eps would pin μ at the cap
    # over the whole domain before any structure exists. Same 1e-3 fractions
    # as _dsgs_2d_energy!.
    ρ̄  = max(abs(scale[1]), eps)
    Ē  = neqs >= 5 ? scale[5]/ρ̄ : zero(TT)
    KĒ = TT(0.5)*(scale[2]*scale[2] + scale[3]*scale[3] + scale[4]*scale[4])/(ρ̄*ρ̄)
    T̄  = max(Ē - KĒ, eps)
    c̄  = sqrt(max(γ*T̄, eps))
    FLOOR = TT(1.0e-3)
    raw    = sgs.raw
    floors = sgs.floors
    @inbounds for ieq = 1:neqs
        raw[ieq] = denom[ieq]
    end
    @inbounds begin
        floors[1] = FLOOR*ρ̄
        floors[2] = FLOOR*ρ̄*c̄
        floors[3] = floors[2]; floors[4] = floors[2]
        neqs >= 5 && (floors[5] = FLOOR*ρ̄*c̄*c̄)
        for ieq = 6:neqs; floors[ieq] = FLOOR*max(abs(scale[ieq]), eps); end
        denom[1] = max(denom[1], FLOOR*ρ̄)        + eps
        denom[2] = max(denom[2], FLOOR*ρ̄*c̄)      + eps
        if neqs >= 5
            denom[5] = max(denom[5], FLOOR*ρ̄*c̄*c̄) + eps
        end
        for ieq = 6:neqs
            denom[ieq] = max(denom[ieq] + eps, FLOOR*max(abs(scale[ieq]), eps))
        end
    end
    # Floor-as-gate; see the identical block in _dsgs_3d_theta! for why.
    @inbounds for ieq = 1:neqs
        ieq in (3, 4) && continue          # the momentum norm lives in slot 2
        raw[ieq] <= floors[ieq] && (denom[ieq] = TT(Inf))
    end
    dρ = denom[1] === TT(Inf) ? max(raw[1], floors[1]) : denom[1]

    #--- Pass 3: per element -------------------------------------------------
    DSGS_MONITOR[] && fill!(DSGS_DBG, 0.0)
    @inbounds for ie = 1:nelem
        Δ  = luse_local_Δ ? Δelem_filter[ie]/nop : Δfallback
        Δ2 = Δ*Δ

        ratio = zero(TT)     # max_i ‖R_i‖_{∞,e} / denom_i
        wmax  = zero(TT)     # ‖ ‖v‖ + √(γT) ‖_{∞,e}
        ρmax  = zero(TT)     # ‖ρ‖_{∞,e}
        ρ_el  = zero(TT)

        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,k]
            Mi = Minv[ip]

            # dsgs_bdf2 is 3q - 4q1 + q2 written as 3(q-q1) - (q1-q2): same
            # number, eps*|dq| of cancellation instead of eps*|q|. A named
            # helper rather than a closure over `ip` -- a closure capturing a
            # loop variable boxes it, and this is the innermost loop of the
            # model.
            # w = 0 on a domain boundary; see dsgs_wres in params_setup.jl.
            w   = wres[ip]
            Rρ  = w*abs(dsgs_bdf2(q, q1, q2, ip, 1)*inv2Δt - Mi*rhs[ip,1])
            Rmu = dsgs_bdf2(q, q1, q2, ip, 2)*inv2Δt - Mi*rhs[ip,2]
            Rmv = dsgs_bdf2(q, q1, q2, ip, 3)*inv2Δt - Mi*rhs[ip,3]
            Rmw = dsgs_bdf2(q, q1, q2, ip, 4)*inv2Δt - Mi*rhs[ip,4]
            Rm  = w*sqrt(Rmu*Rmu + Rmv*Rmv + Rmw*Rmw)
            ratio = max(ratio, Rρ/denom[1], Rm/denom[2])
            DSGS_MONITOR[] && (_dsgs_dbg!(Rρ/denom[1], ip, 1, Rρ);
                               _dsgs_dbg!(Rm/denom[2], ip, 2, Rm))
            if neqs >= 5
                RE = w*abs(dsgs_bdf2(q, q1, q2, ip, 5)*inv2Δt - Mi*rhs[ip,5])
                ratio = max(ratio, RE/denom[5])
                DSGS_MONITOR[] && _dsgs_dbg!(RE/denom[5], ip, 5, RE)
            end
            for ieq = 6:neqs
                Rq = w*abs(dsgs_bdf2(q, q1, q2, ip, ieq)*inv2Δt - Mi*rhs[ip,ieq])
                ratio = max(ratio, Rq/denom[ieq])
                DSGS_MONITOR[] && _dsgs_dbg!(Rq/denom[ieq], ip, ieq, Rq)
            end

            ρl = lpert ? q[ip,1] + qe[ip,1] : q[ip,1]
            ρl = ρl > zero(TT) ? ρl : eps
            ul = q[ip,2]/ρl
            vl = q[ip,3]/ρl
            wl = q[ip,4]/ρl
            El = neqs >= 5 ? (lpert ? (q[ip,5] + qe[ip,5])/ρl : q[ip,5]/ρl) : zero(TT)
            # T = E/ρ - ½‖v‖², clamped at zero so a transient negative internal
            # energy in an under-resolved cell cannot produce a NaN wave speed.
            Tl = max(El - TT(0.5)*(ul*ul + vl*vl + wl*wl), zero(TT))
            wmax = max(wmax, sqrt(ul*ul + vl*vl + wl*wl) + sqrt(γ*Tl))
            ρmax = max(ρmax, ρl)
            ρ_el += ρl
        end
        ρ_el /= TT(ngl*ngl*ngl)

        # eq. (3.5)-(3.6). BOTH branches carry a density, so this is DYNAMIC
        # already and must NOT be multiplied by ρ̄_el.
        μ_res = C1*Δ2*dρ*ratio
        μ_cap = C2*Δ*ρmax*wmax
        μ     = max(zero(TT), min(μ_cap, μ_res))

        μ_el[ie] = μ
        ν_el[ie] = μ/max(ρ_el, eps)      # kinematic, for reporting only

        # Diagnostic slots for VTK, mirroring the split SGS_diffusion applies.
        # Slot 1 is NOT zero here: eq. (3.3) carries beta*grad(rho) in the mass
        # flux and ffs_step's deck records :μ[1] = 0 as the single most
        # destabilising change measured on that case.
        μ_dsgs[ie,1] = visc_coeff[1]*μ/sgs.Sc_t
        for ieq = 2:min(neqs,4)
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μ
        end
        if neqs >= 5
            μ_dsgs[ie,5] = visc_coeff[5]*μ/sgs.Pr_t
        end
        for ieq = 6:neqs
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μ/sgs.Sc_t
        end
    end

    return nothing
end


# ================================================================================
# compute_sgs_cache!(::SGS_DSGS, ..., ::NSD_3D)
#
# The per-element half of DynSGS-3D. compute_dsgs_viscosity! has already run
# for the whole mesh and left this element's coefficient in sgs.μ_el[iel]; all
# that is left per element is
#
#   * broadcast μ_el[iel] onto the element's nodes in μ_turb, which is what
#     SGS_diffusion reads inside the ieq loop that follows this call;
#   * fill Sij, which DynSGS itself does not need but les_statistics.jl does
#     (_fill_sgs_cached!, the branch taken whenever params.sgs is a closure);
#   * optionally ADD the Smagorinsky viscosity — see :dsgs_add_smagorinsky on
#     the struct for why that switch exists.
#
# N² and f_Ri are filled when :lrichardson is on so the diagnostic fields mean
# the same thing they do under SMAG. They multiply the SMAGORINSKY part only:
# a residual is not a strain rate and a buoyancy stability function has no
# business rescaling it — the residual already knows whether the flow is doing
# anything.
# ================================================================================
function compute_sgs_cache!(sgs::SGS_DSGS,
                            uprimitive,
                            mp, uaux,
                            ngl, dψ,
                            dξdx, dξdy, dξdz,
                            dηdx, dηdy, dηdz,
                            dζdx, dζdy, dζdz,
                            connijk, iel, Δ2,
                            micro, ::NSD_3D)

    lrichardson = sgs.lrichardson
    ladd_smag   = sgs.ladd_smagorinsky
    g       = sgs.g
    Pr_t    = sgs.Pr_t
    C_s2    = sgs.C_s2
    karman  = sgs.karman
    lwall_damping = sgs.lwall_damping
    zwall   = sgs.zwall
    μ_e     = sgs.μ_el[iel]

    for m = 1:ngl, l = 1:ngl, k = 1:ngl
        ip = connijk[iel, k, l, m]

        dudξ = 0.0; dudη = 0.0; dudζ = 0.0
        dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
        dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
        dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0

        for ii = 1:ngl
            dudξ += dψ[ii,k] * uprimitive[ii,l,m,2]
            dudη += dψ[ii,l] * uprimitive[k,ii,m,2]
            dudζ += dψ[ii,m] * uprimitive[k,l,ii,2]
            dvdξ += dψ[ii,k] * uprimitive[ii,l,m,3]
            dvdη += dψ[ii,l] * uprimitive[k,ii,m,3]
            dvdζ += dψ[ii,m] * uprimitive[k,l,ii,3]
            dwdξ += dψ[ii,k] * uprimitive[ii,l,m,4]
            dwdη += dψ[ii,l] * uprimitive[k,ii,m,4]
            dwdζ += dψ[ii,m] * uprimitive[k,l,ii,4]
            dθdξ += dψ[ii,k] * uprimitive[ii,l,m,5]
            dθdη += dψ[ii,l] * uprimitive[k,ii,m,5]
            dθdζ += dψ[ii,m] * uprimitive[k,l,ii,5]
        end

        dξdx_klm = dξdx[iel,k,l,m];  dξdy_klm = dξdy[iel,k,l,m];  dξdz_klm = dξdz[iel,k,l,m]
        dηdx_klm = dηdx[iel,k,l,m];  dηdy_klm = dηdy[iel,k,l,m];  dηdz_klm = dηdz[iel,k,l,m]
        dζdx_klm = dζdx[iel,k,l,m];  dζdy_klm = dζdy[iel,k,l,m];  dζdz_klm = dζdz[iel,k,l,m]

        dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
        dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
        dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
        dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
        dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
        dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
        dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm
        dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm
        dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

        S11 = dudx
        S22 = dvdy
        S33 = dwdz
        S12 = 0.5*(dudy + dvdx)
        S13 = 0.5*(dudz + dwdx)
        S23 = 0.5*(dvdz + dwdy)

        sgs.S11[ip] = S11;  sgs.S22[ip] = S22;  sgs.S33[ip] = S33
        sgs.S12[ip] = S12;  sgs.S13[ip] = S13;  sgs.S23[ip] = S23

        S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
        Sij2_val  = 2.0 * S_ij_S_ij

        N2_val = 0.0
        f_Ri_val = 1.0
        if lrichardson
            θ_ref  = uprimitive[k,l,m,5]
            dθdz   = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm
            N2_val = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdz : 0.0
            Ri     = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = sgs_stability_function(Ri, Pr_t)
        end
        sgs.N2[ip]   = N2_val
        sgs.f_Ri[ip] = f_Ri_val

        # The DynSGS coefficient is constant over the element by construction.
        # Shared (DSS) nodes are overwritten by whichever element writes last,
        # which is harmless: this call happens immediately before the ieq loop
        # for THIS element, so every value the viscous kernel reads inside it
        # is this element's own. Downstream readers (the CFL table, the LES
        # statistics, :implicit_vdiff) see the same last-writer-wins field the
        # VTK mu_dsgs_* outputs already use.
        μ = μ_e
        if ladd_smag
            ρ  = uprimitive[k,l,m,1]
            ℓ2 = sgs_mixing_length2(C_s2, Δ2, zwall[ip], karman, lwall_damping)
            μ += ρ * ℓ2 * sqrt(Sij2_val) * f_Ri_val
        end
        sgs.μ_turb[ip] = μ
    end
    return
end

# Helper: expand the per-element, per-equation μ_dsgs[1:nelem,1:neqs]
# onto every node so the per-equation coefficients can be written to
# PNG / VTU like any other field. Shared (DSS) nodes get the value of
# the last element they belong to — that's fine for visualization.
function broadcast_dsgs_to_nodes!(μ_dsgs_pnode::AbstractMatrix{TT},
                                  μ_dsgs::AbstractMatrix{TT},
                                  connijk::AbstractArray{TI,4},
                                  nelem::Int, ngl::Int,
                                  SD::AbstractSpaceDimensions) where {TT,TI}
    neqs = size(μ_dsgs, 2)
    if SD === NSD_1D()
        @inbounds for ie = 1:nelem
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                for ieq = 1:neqs
                    μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, ieq]
                end
            end
        end
    elseif SD === NSD_2D()
        @inbounds for ie = 1:nelem
            for j = 1:ngl
                for i = 1:ngl
                    ip = connijk[ie,i,j,1]
                    for ieq = 1:neqs
                        μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, ieq]
                    end
                end
            end
        end
    elseif SD === NSD_3D()
        @inbounds for ie = 1:nelem
            for k = 1:ngl
                for j = 1:ngl
                    for i = 1:ngl
                        ip = connijk[ie,i,j,k]
                        for ieq = 1:neqs
                            μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, ieq]
                        end
                    end
                end
            end
        end
    end
    return nothing
end
# ================================================================================
# Cache-reading SGS_diffusion — NSD_3D
# Called inside the ieq loop after compute_sgs_cache! has run for the element.
# Reads pre-computed μ_turb[ip] from the sgs struct; no Sij recomputation.
# Dispatches on AbstractSGSModel so one method covers SMAG and VREM.
# ================================================================================
@inline function SGS_diffusion(visc_coeffieq, ieq, ρ, ip,
                                sgs::AbstractSGSModel,
                                ltheta_eqn, ::NSD_3D)
    μ_turb = sgs.μ_turb[ip]
    Pr_t   = sgs.Pr_t
    Sc_t   = sgs.Sc_t
    μ_mol  = sgs.μ_mol
    κ_mol  = sgs.κ_mol

    if ieq == 2 || ieq == 3 || ieq == 4  # momentum
        return (μ_mol + μ_turb) * visc_coeffieq[ieq]
    elseif ieq == 5                        # temperature / energy
        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t
        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (ρ*κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
    else                                   # other scalars (moisture, species)
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end
end

# ================================================================================
# compute_sgs_cache!
# One pass over all GLL points of element iel — fills sgs cache arrays.
# Called once per element before the ieq loop in viscous_rhs_el!, replacing
# the redundant per-equation Sij recomputation.
# ================================================================================

function compute_sgs_cache!(sgs::SGS_SMAG,
                             uprimitive,
                             mp, uaux,
                             ngl, dψ,
                             dξdx, dξdy, dξdz,
                             dηdx, dηdy, dηdz,
                             dζdx, dζdy, dζdz,
                             connijk, iel, Δ2,
                             micro, ::NSD_3D)

    lrichardson = sgs.lrichardson
    g       = sgs.g
    cp      = sgs.cp
    Lc      = sgs.Lc
    Ls      = sgs.Ls
    Rvap    = sgs.Rvap
    Rair    = sgs.Rair
    ε_ratio = sgs.ε_ratio
    # Pr_t, not Ri_crit, is the cutoff: sgs_stability_function follows Lilly,
    # where mixing shuts off at Ri = Pr_t. sgs.Ri_crit is retained on the
    # struct for other consumers but no longer drives the eddy viscosity.
    Pr_t    = sgs.Pr_t
    C_s2    = sgs.C_s2
    karman  = sgs.karman
    lwall_damping = sgs.lwall_damping
    zwall   = sgs.zwall

    for m = 1:ngl, l = 1:ngl, k = 1:ngl
        ip = connijk[iel, k, l, m]

        dudξ = 0.0; dudη = 0.0; dudζ = 0.0
        dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
        dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
        dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0
        dTdξ = 0.0; dTdη = 0.0; dTdζ = 0.0
        dqndξ = 0.0; dqndη = 0.0; dqndζ = 0.0

        for ii = 1:ngl
            dudξ  += dψ[ii,k] * uprimitive[ii,l,m,2]
            dudη  += dψ[ii,l] * uprimitive[k,ii,m,2]
            dudζ  += dψ[ii,m] * uprimitive[k,l,ii,2]
            dvdξ  += dψ[ii,k] * uprimitive[ii,l,m,3]
            dvdη  += dψ[ii,l] * uprimitive[k,ii,m,3]
            dvdζ  += dψ[ii,m] * uprimitive[k,l,ii,3]
            dwdξ  += dψ[ii,k] * uprimitive[ii,l,m,4]
            dwdη  += dψ[ii,l] * uprimitive[k,ii,m,4]
            dwdζ  += dψ[ii,m] * uprimitive[k,l,ii,4]
            dθdξ  += dψ[ii,k] * uprimitive[ii,l,m,5]
            dθdη  += dψ[ii,l] * uprimitive[k,ii,m,5]
            dθdζ  += dψ[ii,m] * uprimitive[k,l,ii,5]
            if micro > 1
                ip_ii = connijk[iel,ii,l,m]
                ip_il = connijk[iel,k,ii,m]
                ip_im = connijk[iel,k,l,ii]
                dTdξ  += dψ[ii,k] * mp.Tabs[ip_ii]
                dTdη  += dψ[ii,l] * mp.Tabs[ip_il]
                dTdζ  += dψ[ii,m] * mp.Tabs[ip_im]
                dqndξ += dψ[ii,k] * mp.qn[ip_ii]
                dqndη += dψ[ii,l] * mp.qn[ip_il]
                dqndζ += dψ[ii,m] * mp.qn[ip_im]
            end
        end

        dξdx_klm = dξdx[iel,k,l,m];  dξdy_klm = dξdy[iel,k,l,m];  dξdz_klm = dξdz[iel,k,l,m]
        dηdx_klm = dηdx[iel,k,l,m];  dηdy_klm = dηdy[iel,k,l,m];  dηdz_klm = dηdz[iel,k,l,m]
        dζdx_klm = dζdx[iel,k,l,m];  dζdy_klm = dζdy[iel,k,l,m];  dζdz_klm = dζdz[iel,k,l,m]

        dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
        dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
        dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
        dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
        dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
        dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
        dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm
        dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm
        dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

        S11 = dudx
        S22 = dvdy
        S33 = dwdz
        S12 = 0.5*(dudy + dvdx)
        S13 = 0.5*(dudz + dwdx)
        S23 = 0.5*(dvdz + dwdy)

        S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
        Sij2_val  = 2.0 * S_ij_S_ij
        Sij_val   = sqrt(Sij2_val)

        sgs.S11[ip] = S11;  sgs.S22[ip] = S22;  sgs.S33[ip] = S33
        sgs.S12[ip] = S12;  sgs.S13[ip] = S13;  sgs.S23[ip] = S23

        # N² — dry or moist (Shi et al. 2019 eqs. 17–22)
        N2_val = 0.0
        if lrichardson
            if micro == 1
                θ_ref  = uprimitive[k,l,m,5]
                dθdz   = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm
                N2_val = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdz : 0.0
            else
                T_ref = mp.Tabs[ip]
                p_ref = uaux[ip, end]
                dTdz  = dTdξ*dξdz_klm + dTdη*dηdz_klm + dTdζ*dζdz_klm
                dqndz = dqndξ*dξdz_klm + dqndη*dηdz_klm + dqndζ*dζdz_klm

                # eq. (21): phase fraction β
                β     = T_ref >= 273.15 ? 1.0 :
                        T_ref >  233.15 ? (T_ref - 233.15)/40.0 : 0.0

                qs_w  = qsatw(T_ref, p_ref)
                qs_i  = qsati(T_ref, p_ref)
                qs_bl = β * qs_w + (1.0 - β) * qs_i  # eq. (20)

                if mp.qn[ip] > qs_bl  # eq. (22): saturated
                    dqsdT = β * dtqsatw(T_ref, p_ref) + (1.0 - β) * dtqsati(T_ref, p_ref)
                    Γ_m_w = (g/cp) * (1.0 + Lc*qs_w/(Rair*T_ref)) /
                                     (1.0 + Lc^2*qs_w/(cp*Rvap*T_ref^2))
                    Γ_m_i = (g/cp) * (1.0 + Ls*qs_i/(Rair*T_ref)) /
                                     (1.0 + Ls^2*qs_i/(cp*Rvap*T_ref^2))
                    Γ_m   = β * Γ_m_w + (1.0 - β) * Γ_m_i
                    N2_val = (g/T_ref) * (dTdz + Γ_m) *
                             (1.0 + T_ref/(ε_ratio + qs_bl) * dqsdT) -
                             g/(1.0 + mp.qn[ip]) * dqndz
                else  # subsaturated: dry N² using T
                    N2_val = (g/T_ref) * (dTdz + g/cp)
                end
            end
        end
        sgs.N2[ip] = N2_val

        # Richardson stability function
        f_Ri_val = 1.0
        if lrichardson
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = sgs_stability_function(Ri, Pr_t)
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,m,1]
        # ℓ² = (C_sΔ)² away from the wall, → (κz)² approaching it.
        ℓ2 = sgs_mixing_length2(C_s2, Δ2, zwall[ip], karman, lwall_damping)
        sgs.μ_turb[ip] = ρ * ℓ2 * Sij_val * f_Ri_val
    end
    return
end

function compute_sgs_cache!(sgs::SGS_VREM,
                             uprimitive,
                             mp, uaux,
                             ngl, dψ,
                             dξdx, dξdy, dξdz,
                             dηdx, dηdy, dηdz,
                             dζdx, dζdy, dζdz,
                             connijk, iel, Δ2,
                             micro, ::NSD_3D)

    lrichardson = sgs.lrichardson
    g       = sgs.g
    cp      = sgs.cp
    Lc      = sgs.Lc
    Ls      = sgs.Ls
    Rvap    = sgs.Rvap
    Rair    = sgs.Rair
    ε_ratio = sgs.ε_ratio
    # Pr_t, not Ri_crit, is the cutoff: sgs_stability_function follows Lilly,
    # where mixing shuts off at Ri = Pr_t. sgs.Ri_crit is retained on the
    # struct for other consumers but no longer drives the eddy viscosity.
    Pr_t    = sgs.Pr_t
    C_vrem  = sgs.C_vrem
    eps_v   = eps(1.0)
    karman        = sgs.karman
    lwall_damping = sgs.lwall_damping
    zwall         = sgs.zwall

    for m = 1:ngl, l = 1:ngl, k = 1:ngl
        ip = connijk[iel, k, l, m]

        dudξ = 0.0; dudη = 0.0; dudζ = 0.0
        dvdξ = 0.0; dvdη = 0.0; dvdζ = 0.0
        dwdξ = 0.0; dwdη = 0.0; dwdζ = 0.0
        dθdξ = 0.0; dθdη = 0.0; dθdζ = 0.0
        dTdξ = 0.0; dTdη = 0.0; dTdζ = 0.0
        dqndξ = 0.0; dqndη = 0.0; dqndζ = 0.0

        for ii = 1:ngl
            dudξ  += dψ[ii,k] * uprimitive[ii,l,m,2]
            dudη  += dψ[ii,l] * uprimitive[k,ii,m,2]
            dudζ  += dψ[ii,m] * uprimitive[k,l,ii,2]
            dvdξ  += dψ[ii,k] * uprimitive[ii,l,m,3]
            dvdη  += dψ[ii,l] * uprimitive[k,ii,m,3]
            dvdζ  += dψ[ii,m] * uprimitive[k,l,ii,3]
            dwdξ  += dψ[ii,k] * uprimitive[ii,l,m,4]
            dwdη  += dψ[ii,l] * uprimitive[k,ii,m,4]
            dwdζ  += dψ[ii,m] * uprimitive[k,l,ii,4]
            dθdξ  += dψ[ii,k] * uprimitive[ii,l,m,5]
            dθdη  += dψ[ii,l] * uprimitive[k,ii,m,5]
            dθdζ  += dψ[ii,m] * uprimitive[k,l,ii,5]
            if micro > 1
                ip_ii = connijk[iel,ii,l,m]
                ip_il = connijk[iel,k,ii,m]
                ip_im = connijk[iel,k,l,ii]
                dTdξ  += dψ[ii,k] * mp.Tabs[ip_ii]
                dTdη  += dψ[ii,l] * mp.Tabs[ip_il]
                dTdζ  += dψ[ii,m] * mp.Tabs[ip_im]
                dqndξ += dψ[ii,k] * mp.qn[ip_ii]
                dqndη += dψ[ii,l] * mp.qn[ip_il]
                dqndζ += dψ[ii,m] * mp.qn[ip_im]
            end
        end

        dξdx_klm = dξdx[iel,k,l,m];  dξdy_klm = dξdy[iel,k,l,m];  dξdz_klm = dξdz[iel,k,l,m]
        dηdx_klm = dηdx[iel,k,l,m];  dηdy_klm = dηdy[iel,k,l,m];  dηdz_klm = dηdz[iel,k,l,m]
        dζdx_klm = dζdx[iel,k,l,m];  dζdy_klm = dζdy[iel,k,l,m];  dζdz_klm = dζdz[iel,k,l,m]

        dudx = dudξ*dξdx_klm + dudη*dηdx_klm + dudζ*dζdx_klm
        dudy = dudξ*dξdy_klm + dudη*dηdy_klm + dudζ*dζdy_klm
        dudz = dudξ*dξdz_klm + dudη*dηdz_klm + dudζ*dζdz_klm
        dvdx = dvdξ*dξdx_klm + dvdη*dηdx_klm + dvdζ*dζdx_klm
        dvdy = dvdξ*dξdy_klm + dvdη*dηdy_klm + dvdζ*dζdy_klm
        dvdz = dvdξ*dξdz_klm + dvdη*dηdz_klm + dvdζ*dζdz_klm
        dwdx = dwdξ*dξdx_klm + dwdη*dηdx_klm + dwdζ*dζdx_klm
        dwdy = dwdξ*dξdy_klm + dwdη*dηdy_klm + dwdζ*dζdy_klm
        dwdz = dwdξ*dξdz_klm + dwdη*dηdz_klm + dwdζ*dζdz_klm

        # Vreman β tensor (uses full velocity gradient, not symmetrized)
        β11 = Δ2*(dudx*dudx + dudy*dudy + dudz*dudz)
        β12 = Δ2*(dudx*dvdx + dudy*dvdy + dudz*dvdz)
        β13 = Δ2*(dudx*dwdx + dudy*dwdy + dudz*dwdz)
        β22 = Δ2*(dvdx*dvdx + dvdy*dvdy + dvdz*dvdz)
        β23 = Δ2*(dvdx*dwdx + dvdy*dwdy + dvdz*dwdz)
        β33 = Δ2*(dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)
        B_β = β11*β22 + β11*β33 + β22*β33 - (β12*β12 + β13*β13 + β23*β23)
        u_ij_u_ij = dudx*dudx + dudy*dudy + dudz*dudz +
                    dvdx*dvdx + dvdy*dvdy + dvdz*dvdz +
                    dwdx*dwdx + dwdy*dwdy + dwdz*dwdz

        # N² (same logic as SGS_SMAG)
        N2_val = 0.0
        if lrichardson
            if micro == 1
                θ_ref  = uprimitive[k,l,m,5]
                dθdz   = dθdξ*dξdz_klm + dθdη*dηdz_klm + dθdζ*dζdz_klm
                N2_val = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdz : 0.0
            else
                T_ref = mp.Tabs[ip]
                p_ref = uaux[ip, end]
                dTdz  = dTdξ*dξdz_klm + dTdη*dηdz_klm + dTdζ*dζdz_klm
                dqndz = dqndξ*dξdz_klm + dqndη*dηdz_klm + dqndζ*dζdz_klm

                β     = T_ref >= 273.15 ? 1.0 :
                        T_ref >  233.15 ? (T_ref - 233.15)/40.0 : 0.0
                qs_w  = qsatw(T_ref, p_ref)
                qs_i  = qsati(T_ref, p_ref)
                qs_bl = β * qs_w + (1.0 - β) * qs_i

                if mp.qn[ip] > qs_bl
                    dqsdT = β * dtqsatw(T_ref, p_ref) + (1.0 - β) * dtqsati(T_ref, p_ref)
                    Γ_m_w = (g/cp) * (1.0 + Lc*qs_w/(Rair*T_ref)) /
                                     (1.0 + Lc^2*qs_w/(cp*Rvap*T_ref^2))
                    Γ_m_i = (g/cp) * (1.0 + Ls*qs_i/(Rair*T_ref)) /
                                     (1.0 + Ls^2*qs_i/(cp*Rvap*T_ref^2))
                    Γ_m   = β * Γ_m_w + (1.0 - β) * Γ_m_i
                    N2_val = (g/T_ref) * (dTdz + Γ_m) *
                             (1.0 + T_ref/(ε_ratio + qs_bl) * dqsdT) -
                             g/(1.0 + mp.qn[ip]) * dqndz
                else
                    N2_val = (g/T_ref) * (dTdz + g/cp)
                end
            end
        end
        sgs.N2[ip] = N2_val

        # Strain rate, computed and stored UNCONDITIONALLY (SGS_SMAG does the
        # same). Vreman's own eddy viscosity does not use it -- that comes from
        # B_β/‖∇u‖² above -- but les_statistics reads these six out of the
        # cache, and that is what lets it report the stress built from the SAME
        # μ_turb the model applied, near-wall limiter and Richardson factor and
        # filter width included, rather than a recomputed one that sees none of
        # them.
        S11 = dudx;  S22 = dvdy;  S33 = dwdz
        S12 = 0.5*(dudy + dvdx)
        S13 = 0.5*(dudz + dwdx)
        S23 = 0.5*(dvdz + dwdy)
        sgs.S11[ip] = S11;  sgs.S22[ip] = S22;  sgs.S33[ip] = S33
        sgs.S12[ip] = S12;  sgs.S13[ip] = S13;  sgs.S23[ip] = S23

        f_Ri_val = 1.0
        if lrichardson
            S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
            Sij2_val  = 2.0 * S_ij_S_ij
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = sgs_stability_function(Ri, Pr_t)
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,m,1]
        # NEAR-WALL LIMIT, applied exactly where Smagorinsky applies it.
        #
        # Vreman's β_ij carries Δ² (see the β11.. lines above), so B_β goes as
        # Δ⁴ and sqrt(B_β/‖∇u‖²) as Δ²: the constant multiplying Δ² here is
        # C_vrem, filling precisely the role C_s² fills in Smagorinsky. So the
        # limiter is the same function on the same footing --
        # ℓ² = C_vrem Δ² (κz)² / (C_vrem Δ² + (κz)²) -- and dividing it by Δ²
        # puts it back in the place C_vrem occupied.
        #
        # With lwall_damping false, sgs_mixing_length2 returns C_vrem*Δ² and
        # ℓ2/Δ2 is C_vrem to the last bit: this is a no-op unless asked for.
        ℓ2 = sgs_mixing_length2(C_vrem, Δ2, zwall[ip], karman, lwall_damping)
        C_eff = Δ2 > 0.0 ? ℓ2 / Δ2 : C_vrem
        μ_base = (u_ij_u_ij > eps_v && B_β > 0.0) ?
                 ρ * C_eff * sqrt(B_β / u_ij_u_ij) : 0.0
        sgs.μ_turb[ip] = μ_base * f_Ri_val
    end
    return
end

# ================================================================================
# Cache-filling compute_sgs_cache! — NSD_2D
# Mirrors the 3D versions above but drops the ζ/w terms. y is the
# vertical coordinate in 2D (see compute_vertical_derivative_q!), so
# the buoyancy gradient uses dξdy/dηdy in place of dξdz/dηdz.
# ================================================================================
function compute_sgs_cache!(sgs::SGS_SMAG,
                             uprimitive,
                             mp, uaux,
                             ngl, dψ,
                             dξdx, dξdy,
                             dηdx, dηdy,
                             connijk, iel, Δ2,
                             micro, ::NSD_2D)

    lrichardson = sgs.lrichardson
    g       = sgs.g
    cp      = sgs.cp
    Lc      = sgs.Lc
    Ls      = sgs.Ls
    Rvap    = sgs.Rvap
    Rair    = sgs.Rair
    ε_ratio = sgs.ε_ratio
    # Pr_t, not Ri_crit, is the cutoff: sgs_stability_function follows Lilly,
    # where mixing shuts off at Ri = Pr_t. sgs.Ri_crit is retained on the
    # struct for other consumers but no longer drives the eddy viscosity.
    Pr_t    = sgs.Pr_t
    C_s2    = sgs.C_s2
    karman  = sgs.karman
    lwall_damping = sgs.lwall_damping
    zwall   = sgs.zwall

    for l = 1:ngl, k = 1:ngl
        ip = connijk[iel, k, l]

        dudξ = 0.0; dudη = 0.0
        dvdξ = 0.0; dvdη = 0.0
        dθdξ = 0.0; dθdη = 0.0
        dTdξ = 0.0; dTdη = 0.0
        dqndξ = 0.0; dqndη = 0.0

        for ii = 1:ngl
            dudξ  += dψ[ii,k] * uprimitive[ii,l,2]
            dudη  += dψ[ii,l] * uprimitive[k,ii,2]
            dvdξ  += dψ[ii,k] * uprimitive[ii,l,3]
            dvdη  += dψ[ii,l] * uprimitive[k,ii,3]
            dθdξ  += dψ[ii,k] * uprimitive[ii,l,4]
            dθdη  += dψ[ii,l] * uprimitive[k,ii,4]
            if micro > 1
                ip_ii = connijk[iel,ii,l]
                ip_il = connijk[iel,k,ii]
                dTdξ  += dψ[ii,k] * mp.Tabs[ip_ii]
                dTdη  += dψ[ii,l] * mp.Tabs[ip_il]
                dqndξ += dψ[ii,k] * mp.qn[ip_ii]
                dqndη += dψ[ii,l] * mp.qn[ip_il]
            end
        end

        dξdx_kl = dξdx[iel,k,l];  dξdy_kl = dξdy[iel,k,l]
        dηdx_kl = dηdx[iel,k,l];  dηdy_kl = dηdy[iel,k,l]

        dudx = dudξ*dξdx_kl + dudη*dηdx_kl
        dudy = dudξ*dξdy_kl + dudη*dηdy_kl
        dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl
        dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl

        S11 = dudx
        S22 = dvdy
        S12 = 0.5*(dudy + dvdx)

        S_ij_S_ij = S11*S11 + S22*S22 + 2.0*S12*S12
        Sij2_val  = 2.0 * S_ij_S_ij
        Sij_val   = sqrt(Sij2_val)

        sgs.S11[ip] = S11;  sgs.S22[ip] = S22;  sgs.S12[ip] = S12

        # N² — dry or moist (Shi et al. 2019 eqs. 17–22); y is vertical in 2D.
        N2_val = 0.0
        if lrichardson
            if micro == 1
                θ_ref  = uprimitive[k,l,4]
                dθdy   = dθdξ*dξdy_kl + dθdη*dηdy_kl
                N2_val = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdy : 0.0
            else
                T_ref = mp.Tabs[ip]
                p_ref = uaux[ip, end]
                dTdy  = dTdξ*dξdy_kl + dTdη*dηdy_kl
                dqndy = dqndξ*dξdy_kl + dqndη*dηdy_kl

                β     = T_ref >= 273.15 ? 1.0 :
                        T_ref >  233.15 ? (T_ref - 233.15)/40.0 : 0.0

                qs_w  = qsatw(T_ref, p_ref)
                qs_i  = qsati(T_ref, p_ref)
                qs_bl = β * qs_w + (1.0 - β) * qs_i

                if mp.qn[ip] > qs_bl
                    dqsdT = β * dtqsatw(T_ref, p_ref) + (1.0 - β) * dtqsati(T_ref, p_ref)
                    Γ_m_w = (g/cp) * (1.0 + Lc*qs_w/(Rair*T_ref)) /
                                     (1.0 + Lc^2*qs_w/(cp*Rvap*T_ref^2))
                    Γ_m_i = (g/cp) * (1.0 + Ls*qs_i/(Rair*T_ref)) /
                                     (1.0 + Ls^2*qs_i/(cp*Rvap*T_ref^2))
                    Γ_m   = β * Γ_m_w + (1.0 - β) * Γ_m_i
                    N2_val = (g/T_ref) * (dTdy + Γ_m) *
                             (1.0 + T_ref/(ε_ratio + qs_bl) * dqsdT) -
                             g/(1.0 + mp.qn[ip]) * dqndy
                else
                    N2_val = (g/T_ref) * (dTdy + g/cp)
                end
            end
        end
        sgs.N2[ip] = N2_val

        f_Ri_val = 1.0
        if lrichardson
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = sgs_stability_function(Ri, Pr_t)
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,1]
        # ℓ² = (C_sΔ)² away from the wall, → (κz)² approaching it.
        ℓ2 = sgs_mixing_length2(C_s2, Δ2, zwall[ip], karman, lwall_damping)
        sgs.μ_turb[ip] = ρ * ℓ2 * Sij_val * f_Ri_val
    end
    return
end

function compute_sgs_cache!(sgs::SGS_VREM,
                             uprimitive,
                             mp, uaux,
                             ngl, dψ,
                             dξdx, dξdy,
                             dηdx, dηdy,
                             connijk, iel, Δ2,
                             micro, ::NSD_2D)

    lrichardson = sgs.lrichardson
    g       = sgs.g
    cp      = sgs.cp
    Lc      = sgs.Lc
    Ls      = sgs.Ls
    Rvap    = sgs.Rvap
    Rair    = sgs.Rair
    ε_ratio = sgs.ε_ratio
    # Pr_t, not Ri_crit, is the cutoff: sgs_stability_function follows Lilly,
    # where mixing shuts off at Ri = Pr_t. sgs.Ri_crit is retained on the
    # struct for other consumers but no longer drives the eddy viscosity.
    Pr_t    = sgs.Pr_t
    C_vrem  = sgs.C_vrem
    eps_v   = eps(1.0)
    karman        = sgs.karman
    lwall_damping = sgs.lwall_damping
    zwall         = sgs.zwall

    for l = 1:ngl, k = 1:ngl
        ip = connijk[iel, k, l]

        dudξ = 0.0; dudη = 0.0
        dvdξ = 0.0; dvdη = 0.0
        dθdξ = 0.0; dθdη = 0.0
        dTdξ = 0.0; dTdη = 0.0
        dqndξ = 0.0; dqndη = 0.0

        for ii = 1:ngl
            dudξ  += dψ[ii,k] * uprimitive[ii,l,2]
            dudη  += dψ[ii,l] * uprimitive[k,ii,2]
            dvdξ  += dψ[ii,k] * uprimitive[ii,l,3]
            dvdη  += dψ[ii,l] * uprimitive[k,ii,3]
            dθdξ  += dψ[ii,k] * uprimitive[ii,l,4]
            dθdη  += dψ[ii,l] * uprimitive[k,ii,4]
            if micro > 1
                ip_ii = connijk[iel,ii,l]
                ip_il = connijk[iel,k,ii]
                dTdξ  += dψ[ii,k] * mp.Tabs[ip_ii]
                dTdη  += dψ[ii,l] * mp.Tabs[ip_il]
                dqndξ += dψ[ii,k] * mp.qn[ip_ii]
                dqndη += dψ[ii,l] * mp.qn[ip_il]
            end
        end

        dξdx_kl = dξdx[iel,k,l];  dξdy_kl = dξdy[iel,k,l]
        dηdx_kl = dηdx[iel,k,l];  dηdy_kl = dηdy[iel,k,l]

        dudx = dudξ*dξdx_kl + dudη*dηdx_kl
        dudy = dudξ*dξdy_kl + dudη*dηdy_kl
        dvdx = dvdξ*dξdx_kl + dvdη*dηdx_kl
        dvdy = dvdξ*dξdy_kl + dvdη*dηdy_kl

        β11 = Δ2*(dudx*dudx + dudy*dudy)
        β12 = Δ2*(dudx*dvdx + dudy*dvdy)
        β22 = Δ2*(dvdx*dvdx + dvdy*dvdy)
        B_β = β11*β22 - β12*β12
        u_ij_u_ij = dudx*dudx + dudy*dudy + dvdx*dvdx + dvdy*dvdy

        # N² — dry or moist; y is vertical in 2D.
        N2_val = 0.0
        if lrichardson
            if micro == 1
                θ_ref  = uprimitive[k,l,4]
                dθdy   = dθdξ*dξdy_kl + dθdη*dηdy_kl
                N2_val = abs(θ_ref) > 1e-12 ? (g / θ_ref) * dθdy : 0.0
            else
                T_ref = mp.Tabs[ip]
                p_ref = uaux[ip, end]
                dTdy  = dTdξ*dξdy_kl + dTdη*dηdy_kl
                dqndy = dqndξ*dξdy_kl + dqndη*dηdy_kl

                β     = T_ref >= 273.15 ? 1.0 :
                        T_ref >  233.15 ? (T_ref - 233.15)/40.0 : 0.0

                qs_w  = qsatw(T_ref, p_ref)
                qs_i  = qsati(T_ref, p_ref)
                qs_bl = β * qs_w + (1.0 - β) * qs_i

                if mp.qn[ip] > qs_bl
                    dqsdT = β * dtqsatw(T_ref, p_ref) + (1.0 - β) * dtqsati(T_ref, p_ref)
                    Γ_m_w = (g/cp) * (1.0 + Lc*qs_w/(Rair*T_ref)) /
                                     (1.0 + Lc^2*qs_w/(cp*Rvap*T_ref^2))
                    Γ_m_i = (g/cp) * (1.0 + Ls*qs_i/(Rair*T_ref)) /
                                     (1.0 + Ls^2*qs_i/(cp*Rvap*T_ref^2))
                    Γ_m   = β * Γ_m_w + (1.0 - β) * Γ_m_i
                    N2_val = (g/T_ref) * (dTdy + Γ_m) *
                             (1.0 + T_ref/(ε_ratio + qs_bl) * dqsdT) -
                             g/(1.0 + mp.qn[ip]) * dqndy
                else
                    N2_val = (g/T_ref) * (dTdy + g/cp)
                end
            end
        end
        sgs.N2[ip] = N2_val

        # Same as the 3D method above; 2D stores the three components its
        # SGS_SMAG counterpart stores, leaving S33/S13/S23 at zero.
        S11 = dudx;  S22 = dvdy;  S12 = 0.5*(dudy + dvdx)
        sgs.S11[ip] = S11;  sgs.S22[ip] = S22;  sgs.S12[ip] = S12

        f_Ri_val = 1.0
        if lrichardson
            Sij2_val = 2.0*(S11*S11 + S22*S22 + 2.0*S12*S12)
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = sgs_stability_function(Ri, Pr_t)
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,1]
        # NEAR-WALL LIMIT, applied exactly where Smagorinsky applies it.
        #
        # Vreman's β_ij carries Δ² (see the β11.. lines above), so B_β goes as
        # Δ⁴ and sqrt(B_β/‖∇u‖²) as Δ²: the constant multiplying Δ² here is
        # C_vrem, filling precisely the role C_s² fills in Smagorinsky. So the
        # limiter is the same function on the same footing --
        # ℓ² = C_vrem Δ² (κz)² / (C_vrem Δ² + (κz)²) -- and dividing it by Δ²
        # puts it back in the place C_vrem occupied.
        #
        # With lwall_damping false, sgs_mixing_length2 returns C_vrem*Δ² and
        # ℓ2/Δ2 is C_vrem to the last bit: this is a no-op unless asked for.
        ℓ2 = sgs_mixing_length2(C_vrem, Δ2, zwall[ip], karman, lwall_damping)
        C_eff = Δ2 > 0.0 ? ℓ2 / Δ2 : C_vrem
        μ_base = (u_ij_u_ij > eps_v && B_β > 0.0) ?
                 ρ * C_eff * sqrt(B_β / u_ij_u_ij) : 0.0
        sgs.μ_turb[ip] = μ_base * f_Ri_val
    end
    return
end

# ================================================================================
# Cache-reading SGS_diffusion — NSD_2D
# Mirrors the NSD_3D cache reader above; ieq layout for 2D is
# (ρ, ρu, ρv, ρθ/energy, [scalars...]).
# ================================================================================
@inline function SGS_diffusion(visc_coeffieq, ieq, ρ, ip,
                                sgs::AbstractSGSModel,
                                ltheta_eqn, ::NSD_2D)
    μ_turb = sgs.μ_turb[ip]
    Pr_t   = sgs.Pr_t
    Sc_t   = sgs.Sc_t
    μ_mol  = sgs.μ_mol
    κ_mol  = sgs.κ_mol

    if ieq == 2 || ieq == 3                # momentum (u, v)
        return (μ_mol + μ_turb) * visc_coeffieq[ieq]
    elseif ieq == 4                         # temperature / energy
        # DYNAMIC, not kinematic: the prognostic variable is ρθ (ρhl, ρqx),
        # so the diffusive flux assembled by _expansion_visc! as κ_eff·∇θ must
        # carry the ρ of ρ·ν_t/Pr_t·∇θ. Dividing μ_turb by ρ here cancelled it
        # and left every scalar flux short by a factor ρ (≈14% at the surface,
        # and inconsistent with the momentum branch above, which is dynamic).
        κ_turb = μ_turb / Pr_t
        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (ρ*κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
    else                                      # other scalars (moisture, species)
        κ_turb_scalar = μ_turb / Sc_t   # dynamic, see the θ branch above
        return (ρ*κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end
end

# ================================================================================
# Helpers for LES statistics: compute S_ij and μ_turb at a single GLL point.
# Returns (μ_turb, S11, S22, S33, S12, S13, S23, S_ij_S_ij).
# Called only at statistics output time (not on the hot RHS path).
# ================================================================================
@inline function compute_sij_and_mu_turb(ρ,
                                          dudx, dudy, dudz,
                                          dvdx, dvdy, dvdz,
                                          dwdx, dwdy, dwdz,
                                          PhysConst, Δ2, ::SMAG; C_s = PhysConst.C_s)
    C_s2 = C_s * C_s
    S11  = dudx;  S22 = dvdy;  S33 = dwdz
    S12  = 0.5 * (dudy + dvdx)
    S13  = 0.5 * (dudz + dwdx)
    S23  = 0.5 * (dvdz + dwdy)
    S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2*(S12*S12 + S13*S13 + S23*S23)
    μ_turb = ρ * C_s2 * Δ2 * sqrt(2 * S_ij_S_ij)
    return μ_turb, S11, S22, S33, S12, S13, S23, S_ij_S_ij
end

@inline function compute_sij_and_mu_turb(ρ,
                                          dudx, dudy, dudz,
                                          dvdx, dvdy, dvdz,
                                          dwdx, dwdy, dwdz,
                                          PhysConst, Δ2, ::VREM; C_s = PhysConst.C_s)
    C_s2   = C_s * C_s
    C_vrem = 2.5 * C_s2
    eps_v  = eps(1.0)
    β11 = Δ2 * (dudx*dudx + dudy*dudy + dudz*dudz)
    β12 = Δ2 * (dudx*dvdx + dudy*dvdy + dudz*dvdz)
    β13 = Δ2 * (dudx*dwdx + dudy*dwdy + dudz*dwdz)
    β22 = Δ2 * (dvdx*dvdx + dvdy*dvdy + dvdz*dvdz)
    β23 = Δ2 * (dvdx*dwdx + dvdy*dwdy + dvdz*dwdz)
    β33 = Δ2 * (dwdx*dwdx + dwdy*dwdy + dwdz*dwdz)
    B_β = β11*β22 + β11*β33 + β22*β33 - (β12*β12 + β13*β13 + β23*β23)
    u_ij_u_ij = dudx*dudx + dudy*dudy + dudz*dudz +
                dvdx*dvdx + dvdy*dvdy + dvdz*dvdz +
                dwdx*dwdx + dwdy*dwdy + dwdz*dwdz
    μ_turb = (u_ij_u_ij > eps_v && B_β > 0.0) ?
             ρ * C_vrem * sqrt(B_β / u_ij_u_ij) : 0.0
    S11  = dudx;  S22 = dvdy;  S33 = dwdz
    S12  = 0.5 * (dudy + dvdx)
    S13  = 0.5 * (dudz + dwdx)
    S23  = 0.5 * (dvdz + dwdy)
    S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2*(S12*S12 + S13*S13 + S23*S23)
    return μ_turb, S11, S22, S33, S12, S13, S23, S_ij_S_ij
end

@inline function compute_sij_and_mu_turb(ρ,
                                          dudx, dudy, dudz,
                                          dvdx, dvdy, dvdz,
                                          dwdx, dwdy, dwdz,
                                          PhysConst, Δ2, ::Any; C_s = PhysConst.C_s)
    return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

# ================================================================================
# DynSGS monitor -- JEXPRESSO_DSGS_MONITOR=1
#
# One line per time step with what the model actually produced, plus WHERE the
# worst residual came from. "μ is too big" is not actionable; "μ is too big
# because equation 5's residual at (x, z) = (500, 7500) is 40x the tendency
# there" is, and it is the difference between tuning a constant and finding a
# defect.
#
# What to read:
#   nu       max/mean of the coefficient over the elements. Compare against
#            C₂Δ(‖v‖+c), which in a low-Mach atmosphere is O(1e4) m²/s -- if
#            nu is at that, DynSGS has degenerated to first-order upwind and
#            the run's fate is decided by whether Δt respects THAT scheme's
#            CFL rather than by anything the sensor did.
#   den      the normalising denominators of eq. (9). A coefficient two orders
#            too large usually comes from here, not from the residual.
#   worst    the ratio R_i/den_i that won, the equation it came from, and the
#            coordinates of the node. A worst node sitting ON A BOUNDARY is the
#            signature of the residual picking up a strongly-imposed boundary
#            condition rather than a discretisation error.
#
# Works for every path: the record is written by whichever
# compute_dsgs_viscosity! ran.
# ================================================================================
const DSGS_MONITOR = Ref(false)

# ================================================================================
# WHICH RESIDUAL, AND WHY THE DEFAULT IS THE ONE THAT IS NOT LITERALLY BDF2
#
# `:dsgs_residual => :tendency` (default) | `:strict`
#
# The history pair is advanced once per step. WHERE that advance sits relative
# to the read decides what the sensor measures, and the two choices are not a
# matter of taste -- they measure different things, and the difference was
# established by measurement on CompEuler/sod1d, not by argument.
#
#   :strict     roll AFTER the read, rebuild only on the stage where uaux is
#               qⁿ. The stencil is then literally
#                   (3qⁿ - 4qⁿ⁻¹ + qⁿ⁻²)/2Δt  ≈  ∂q/∂t
#               and R = |∂q/∂t - M⁻¹·RHS|. But `params.RHS` at that point in
#               _build_rhs! holds the INVISCID flux divergence only -- the
#               viscous term is added a few lines later -- while ∂q/∂t is the
#               derivative of the solution the viscous term helped produce. So
#                   R = |the artificial viscosity's own contribution|
#               and the model is reading itself. That is Nazarov & Hoffman
#               eq. (3.4) taken literally (R of the EXACT PDE, on a solution
#               that satisfies PDE + artificial viscosity = 0), and on a
#               weak-form Galerkin RHS it is a contraction: mu decays instead
#               of locking onto the feature.
#
#   :tendency   roll BEFORE the read, rebuild every stage -- what this code has
#               always done. `dsgs_qnm2` is then the CURRENT state, the stencil
#               collapses to -(qⁿ - qⁿ⁻¹)/2Δt, and
#                   R ≈ 1.5·|∂q/∂t|  at the first stage of each step.
#               It is a TENDENCY sensor, not the strict residual. It fires
#               wherever the solution is changing fast, which for a shock, a
#               contact or a steepening front is where the stabilisation is
#               needed, and it does not read its own output.
#
# MEASURED, sod1d at t = 0.2, 100 elements at nop = 4 (element 85 is the shock,
# element 1 is x = 0):
#
#     :tendency  mu max 1.28e-3  at element 85   <- the shock
#     :strict    mu max 1.19e-4  at element  1   <- 10x smaller, and at the
#                                                   boundary, with visible
#                                                   oscillations in rho and rhoE
#
# So the default is :tendency, and the name says what it is rather than
# implying a residual it is not. :strict is kept because the question of how to
# get a true residual out of a weak-form RHS is worth returning to -- the
# answer is probably to compare against RHS_inviscid + RHS_visc from the
# PREVIOUS step, which is the consistency error of the scheme rather than
# either of these -- and because having both makes the difference measurable
# on any case in one run.
# ================================================================================
const DSGS_STRICT = Ref(false)          # :dsgs_residual => :strict
const DSGS_NOMASK = Ref(false)          # JEXPRESSO_DSGS_NOMASK=1, diagnostic

# (max ratio, node, equation, raw residual). A module-level buffer rather than
# a field on SGS_DSGS because the 1D and 2D paths have no closure struct and
# the monitor has to work for them too.
const DSGS_DBG = zeros(Float64, 4)

@inline function _dsgs_dbg!(r, ip, ieq, R)
    if r > DSGS_DBG[1]
        DSGS_DBG[1] = r;   DSGS_DBG[2] = Float64(ip)
        DSGS_DBG[3] = Float64(ieq); DSGS_DBG[4] = R
    end
    return nothing
end

function dsgs_monitor(params, time)
    DSGS_MONITOR[] || return nothing
    MPI.Comm_rank(get_mpi_comm()) == 0 || return nothing
    nelem = Int(params.mesh.nelem)
    nelem > 0 || return nothing

    sgs = params.sgs
    if sgs isa SGS_DSGS
        νmax = maximum(sgs.ν_el); νmean = sum(sgs.ν_el)/nelem; iemax = argmax(sgs.ν_el)
        den  = sgs.denom
    else
        # 1D/2D: the per-equation matrix is all there is. Slot 2 (a momentum
        # slot on every system that uses this) is the representative one.
        col  = size(params.μ_dsgs, 2) >= 2 ? 2 : 1
        v    = @view params.μ_dsgs[1:nelem, col]
        νmax = maximum(v); νmean = sum(v)/nelem; iemax = argmax(v)
        den  = Float64[]
    end

    nop  = Int(params.mesh.nop)
    Δmax = length(params.mesh.Δelem_filter) == nelem && nop > 0 ?
           maximum(params.mesh.Δelem_filter)/nop : Float64(params.mesh.Δeffective_l)

    @printf(" #  DSGS t=%.3f  nu max %.4g (el %d) mean %.4g  Δ %.4g%s\n",
            time, νmax, iemax, νmean, Δmax,
            isempty(den) ? "" : "  den " *
                string([round(d; sigdigits = 3) for d in den]))

    ip = Int(DSGS_DBG[2])
    if ip > 0
        c  = params.mesh.coords
        nsd = size(c, 1)
        pos = nsd >= 3 ? @sprintf("(%.0f, %.0f, %.0f)", c[1,ip], c[2,ip], c[3,ip]) :
                         @sprintf("(%.0f, %.0f)",       c[1,ip], c[2,ip])
        @printf(" #        worst ratio %.4g  eq %d  R %.4g  at ip %d %s%s\n",
                DSGS_DBG[1], Int(DSGS_DBG[3]), DSGS_DBG[4], ip, pos,
                _dsgs_on_boundary(params.mesh, ip) ? "  <-- ON A BOUNDARY" : "")
    end
    return nothing
end

# Is this node on a domain boundary? Answered from the mesh's own boundary
# node lists rather than from the coordinates, so it stays right on a warped
# or partially-periodic mesh.
function _dsgs_on_boundary(mesh, ip::Int)
    pf = mesh.poin_in_bdy_face
    if length(pf) > 0
        for k in eachindex(pf)
            pf[k] == ip && return true
        end
    end
    pe = mesh.poin_in_bdy_edge
    if length(pe) > 0
        for k in eachindex(pe)
            pe[k] == ip && return true
        end
    end
    return false
end

#----------------------------------------------------------------------
# SMAGORINSKY
#----------------------------------------------------------------------
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
        κ_turb = μ_turb / (ρ * Pr_t)

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (κ_mol + κ_turb) * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end


#
#
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
    Ri_crit = PhysConst.Ri_crit # Critical Richardson number (typically 0.25)
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
        # Ri > Ri_crit: turbulence completely suppressed
        Ri = (Sij2 > 1.0f-12) ? N2 / Sij2 : 0.0
        
        # Stability function for Richardson correction
        # Various formulations exist in literature
        f_Ri = if Ri >= Ri_crit
            # Stable stratification above critical Richardson number
            # Turbulence is completely suppressed
            0.0
            
        elseif Ri >= 0.0
            # Stable but sub-critical: reduce mixing
            # Smooth transition to zero at Ri_crit
            # Common formulation: f(Ri) = (1 - Ri/Ri_crit)²
            ratio = Ri / Ri_crit
            (1.0 - ratio) * (1.0 - ratio)
            
        else
            # Unstable stratification (Ri < 0): enhance mixing
            # Various formulations:
            # - sqrt(1 - 16*Ri): from Monin-Obukhov similarity
            # - (1 - 16*Ri)^(1/4): alternative formulation
            # Cap at maximum enhancement factor (e.g., 3x)
            min(sqrt(1.0 - 16.0*Ri), 3.0)
        end
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
        f_Ri = if Ri >= Ri_crit
            # Laminar regime: Stratification is strong enough to kill turbulence
            0.0
            
        elseif Ri >= 0.0
            # Stable regime: Turbulence is present but suppressed by buoyancy
            # Using the quadratic suppression: (1 - Ri/Ri_crit)²
            ratio = Ri / Ri_crit
            (1.0 - ratio) * (1.0 - ratio)
            
        else
            # Unstable regime (Ri < 0): Buoyancy enhances turbulent mixing
            # Enhancement factor capped at 3.0 to maintain numerical stability
            min(sqrt(1.0 - 16.0*Ri), 3.0)
        end
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
        κ_turb = μ_turb / (ρ * Pr_t)
        if ltheta_eqn
            # Potential temperature equation
            return κ_turb * visc_coeffieq[ieq]
        else
            # Internal energy or enthalpy equation
            return (κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
        
    else
        # Other scalar equations (species, TKE, etc.)
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
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

        κ_turb = μ_turb / (ρ * Pr_t)

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            # Total-energy / enthalpy form: thermal conductivity k_eff = cp * μ_eff / Pr_t
            return cp * (μ_mol + μ_turb) / Pr_t * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end



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
    Ri_crit    = PhysConst.Ri_crit   # Critical Richardson number
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
        f_Ri = if Ri >= Ri_crit
            # Stable stratification suppresses turbulence
            0.0
        elseif Ri >= 0.0
            # Stable but sub-critical: reduce mixing
            (1.0 - Ri/Ri_crit)^2
        else
            # Unstable stratification: enhance mixing
            min(sqrt(1.0 - 16.0*Ri), 3.0)  # Cap at 3x base mixing
        end
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

        κ_turb = μ_turb / (ρ * Pr_t)

        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return cp * (κ_mol + κ_turb) * visc_coeffieq[ieq]
        end

    else
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
    end

end

#----------------------------------------------------------------------
# DYNAMIC SGS (Marras et al., 2015, JCP 301:77-101) — residual-based
# artificial viscosity, parameter free.
#
# Per-element coefficient:
#     μ_res = CR · Δ² · max_i ‖R_i‖∞,Ω / ‖q_i − ⟨q_i⟩‖∞,Ω
#     μ_max = Cmax · Δ · max(|u| + c)
#     μ_dsgs[iel] = max(0, min(μ_res, μ_max))
# where R_i is the element-wise strong residual of conservation law i,
# ∂ₜq_i − rhs_el[K,i]/m_i^K (the stage-consistent time derivative minus the
# element's own weak RHS divided by its lumped mass entry; see
# _dsgs_nodal_residual_1d! for why the assembled RHS cannot be used), which
# makes μ_dsgs dimensionally a kinematic viscosity (m²/s) regardless of SD.
#
# Both numerators and denominators are L∞ norms over a region larger than
# one element — the rank's subdomain by default, the whole domain under
# :dsgs_norms => "domain" (see _dsgs_norm_scope below) — so the coefficient
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
#
# :dsgs_nodal_rho (see compute_dsgs_viscosity!(::DSGS_MHD)): when set, the
# momentum and energy slots of μ_dsgs hold the KINEMATIC coefficient and the
# dynamic one is formed here with the density of the quadrature point, i.e.
# the viscous flux is ∇·(ρ μ ∇u) with the local ρ instead of ∇·(ρ̄ μ ∇u)
# with the element mean. rhs.jl sets the Ref from the inputs before each
# DynSGS assembly.
#
const dsgs_nodal_rho = Ref{Bool}(false)

# :dsgs_ref_weight (fluxEmergenceSon2025DSGS): the coefficient of slots 1-5
# is multiplied, at the quadrature point, by the weight the case's
# user_primitives! stores in the spare slot neqs+1 of uprimitive (there:
# the reference density ρ_e), and the primitives of those slots are the
# conserved perturbations divided by that weight. The operator is then
# ∇·(μ ρ_e ∇((q − q_e)/ρ_e)): zero at rest like the q − q_e form, but a
# weighted diffusion of the RELATIVE departure, which obeys a maximum
# principle across a reference jump — the q − q_e form does not (a −10%
# departure of the 25× denser chromosphere side of the solar transition
# region, diffused into the coronal side, exceeds the whole coronal
# density). Read by _expansion_visc!(…, ::ContGal, NSD_2D); rhs.jl sets it
# from the inputs before each DynSGS assembly.
const dsgs_ref_weight = Ref{Bool}(false)

# Split energy flux of a conserved-form DynSGS-MHD case with
# :dsgs_nazarov_energy (set by rhs.jl from the inputs). The case's
# user_primitives! then hands the E slot the NON-THERMAL departure
# δ(½ρ|v|² + ½|B|² + ½ψ²) (relative to ρ_e where :dsgs_ref_weight) in slot 4
# and the THERMAL one δ(p/(γ−1)) in the spare slot neqs+2, and
# _expansion_visc! diffuses the first with the ρv-slot coefficient ν — so
# that the magnetic and kinetic energy fluxes of the E equation stay those
# of the B and ρv Laplacians (measured: scaling the whole slot let ν∇B
# spread the flux sheet while its magnetic energy stayed put, the sheet
# core overheated and the emergence stalled) — and the second with the
# E-slot coefficient max(γ(γ−1)/Pr_t·ν_res, ν_floor): Dao & Nazarov's
# κ = ρν/Pr on the residual part, the full Cmin floor kept (measured: with
# the floor cut 19× too, a node-to-node temperature mode grew across the
# corona within 10 τ₀).
const dsgs_split_energy = Ref{Bool}(false)

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::DSGS_MHD, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return (dsgs_nodal_rho[] && 2 <= ieq <= 5) ? ρ*visc_coeffieq[ieq] : visc_coeffieq[ieq]

end

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
                               ::DSGS_MHD, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)

    return (dsgs_nodal_rho[] && 2 <= ieq <= 5) ? ρ*visc_coeffieq[ieq] : visc_coeffieq[ieq]

end

# Shallow water (DSGS_SW): the element's (or, in the nodal form, the node's)
# kinematic coefficient of the slot, nothing else.
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               ::DSGS_SW, ::NSD_2D;
                               ltheta_eqn=true,
                               lrichardson=false)
    return visc_coeffieq[ieq]
end

@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
                               ::DSGS_SW, ::NSD_2D;
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
# Default: RANK-LOCAL (`lglobal_norms = false`, :dsgs_norms => "rank" in
# mod_inputs.jl; every kernel below and every call site in rhs.jl takes the
# same flag). These two quantities only set
# the SCALE the residual indicator is measured against; what the model needs
# from them is the order of magnitude of the solution's variation, and a
# partition of a connected domain resolves that as well as the whole domain
# does. μ is bounded by min(μ_res, μ_max) either way, so the flow solution
# differs only at the level of the usual round-off divergence. No
# communication at all.
#
# Opt-in: the paper's domain norms, with
#
#     :dsgs_norms => "domain"          # in user_inputs.jl
#
# (params_setup.jl turns it into params.dsgs_global_norms, the Bool the
# call sites in rhs.jl thread down; "element", DSGS_MHD only, normalizes
# per element: params.dsgs_local_norms). Use it when you want μ reproducible across rank
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
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 wt::NTuple{3,TT},
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
    CR    = TT(1.0)
    Cmax    = TT(0.5)
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
    @inbounds for ie = 1:nelem
        Δ = Δx[ie]/ngl

        n1   = zero(TT); n2 = zero(TT); n3 = zero(TT)
        uTmx = zero(TT)
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            imK = one(TT)/(ω[i]*Je[ie,i])   # element lumped mass at the node

            R1 = abs((wt[1]*q[ip,1] + wt[2]*q1[ip,1] + wt[3]*q2[ip,1]) - imK*rhs_el[ie,i,1])
            R2 = abs((wt[1]*q[ip,2] + wt[2]*q1[ip,2] + wt[3]*q2[ip,2]) - imK*rhs_el[ie,i,2])
            R3 = abs((wt[1]*q[ip,3] + wt[2]*q1[ip,3] + wt[3]*q2[ip,3]) - imK*rhs_el[ie,i,3])
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

        μ_res = CR*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3)
        μ_max = Cmax*Δ*uTmx
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
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 wt::NTuple{3,TT},
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 PhysConst::PhysicalConst{TT},
                                 Pr::TT,
                                 nelem::Int, ngl::Int;
                                 ltheta::Bool=true,
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    if !ltheta
        _dsgs_2d_energy!(μ_dsgs, q, q1, q2, rhs_el, ω, Je, visc_coeff,
                         wt, connijk, Δelem, PhysConst, Pr, nelem, ngl,
                         lglobal_norms)
        return nothing
    end

    # Marras et al. (JCP 2015) eq. (8-10). The residual is the
    # ELEMENT-WISE strong residual, ∂ₜq_i − rhs_el[K,i]/m_i^K (see
    # _dsgs_nodal_residual_1d!): the earlier lineage found that the
    # assembled M⁻¹·RHS "shrinks the residual by ~10³ and turns DSGS off"
    # — because with a lumped mass matrix it is the time-integration
    # error, not a residual — and used the un-divided weak RHS instead;
    # the element residual has the units and the meaning both lacked.
    #
    #     μ_res|e = CR · Δ² · max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω
    #     μ_max|e = Cmax · Δ · (|u| + c)_∞,e
    #     μ|e     = max(0, min(μ_max, μ_res))
    #
    # Per-equation split (Marras eq. 10), with the user-supplied
    # inputs[:μ] multiplier on each slot:
    #     μ_dsgs[iel, 1] = 0                              (no mass diffusion)
    #     μ_dsgs[iel, 2] = visc_coeff[2] · μ              (ρu)
    #     μ_dsgs[iel, 3] = visc_coeff[3] · μ              (ρv)
    #     μ_dsgs[iel, 4] = visc_coeff[4] · Pr/(γ-1) · μ   (ρθ)
    #
    # qe stays in the function-barrier signature so the rhs.jl
    # call site doesn't have to change, but they are unused here.

    neqs  = size(μ_dsgs, 2)
    invnp = one(TT)/(nelem*ngl*ngl)
    γ     = PhysConst.γ
    C0    = PhysConst.C0
    CR    = TT(1.0)
    Cmax    = TT(0.5)
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
    # wave-speed bound Cmax·Δ·(|u|+c) before any flow has developed,
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
                imK = one(TT)/(ω[i]*ω[j]*Je[ie,i,j])   # element lumped mass at the node

                # Strong-form residual. rhs[] here is the DSS-assembled
                # WEAK-form RHS (rhs! divides by the mass matrix later), so
                # it must be multiplied by M⁻¹ for the difference to be
                # dimensionally meaningful: ∂q/∂t has units q/time, the raw
                # weak RHS has units (mass matrix)·q/time. The 1D path has
                # always done this; the 2D path did not.
                R1 = abs((wt[1]*q[ip,1] + wt[2]*q1[ip,1] + wt[3]*q2[ip,1]) - imK*rhs_el[ie,i,j,1])
                R2 = abs((wt[1]*q[ip,2] + wt[2]*q1[ip,2] + wt[3]*q2[ip,2]) - imK*rhs_el[ie,i,j,2])
                R3 = abs((wt[1]*q[ip,3] + wt[2]*q1[ip,3] + wt[3]*q2[ip,3]) - imK*rhs_el[ie,i,j,3])
                R4 = abs((wt[1]*q[ip,4] + wt[2]*q1[ip,4] + wt[3]*q2[ip,4]) - imK*rhs_el[ie,i,j,4])
                n1 = max(n1, R1); n2 = max(n2, R2)
                n3 = max(n3, R3); n4 = max(n4, R4)

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

        μ_res = CR*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3, n4/denom4)
        μ_max = Cmax*Δ*uTmx
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
        # Passive tracers (slots 5..neqs, e.g. CompEuler/thetaTracers): the
        # case transports them as un-weighted scalars (∂ₜq + ∇·(q u) = 0)
        # and hands the scalar itself to _expansion_visc!, so the
        # diffusivity is the KINEMATIC ν (m²/s), as the scalar branch of
        # the Smagorinsky model (μ_turb/(ρ Sc_t)). Left unfilled before
        # this, the tracers ran without any stabilization.
        for ieq = 5:neqs
            μ_dsgs[ie,ieq] = visc_coeff[ieq] * μ
        end
    end

    return nothing
end

# ---------------- 3D --------------------------------------------------
#
# Conservation form q = (ρ, ρu, ρv, ρw, ρθ) for the Euler-θ system in
# three dimensions — the same model as the 2D kernel above, node loop and
# element scale extended to the third direction. Written because
# `:visc_model => DSGS()` on a 3D case (problems/CompEuler/3d, the LES
# cases) had no kernel to dispatch to at all: params.sgs is `nothing` for
# every model but Smagorinsky and Vreman, so the 3D viscous assembly fell
# through to its constant-coefficient branch and ran the deck's :μ as a
# plain Laplacian coefficient in m²/s — a DynSGS run that was silently not
# DynSGS, and, with :μ[1] ≠ 0, a mass diffusion the model never asks for.
#
# Same structure as the 2D θ kernel:
#
#     ν_K = max(0, min(C_max Δ λ_K, C_R Δ² R_K)),  Δ = Δ_K/(k+1)
#
# with the element residual of DSGS.md §1.2 and the slot split of Marras
# et al. eq. (10). The primitives handed to the viscous operator are
# (ρ, u, v, w, θ), so the momentum and θ slots carry the DYNAMIC ρ̄ν.
#
# The energy form (:energy_equation => "energy", q with ρE in the last
# slot) has no 3D kernel: `ltheta = false` raises rather than silently
# building a θ-form coefficient from a total-energy state.
#
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS, ::NSD_3D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 wt::NTuple{3,TT},
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 PhysConst::PhysicalConst{TT},
                                 Pr::TT,
                                 nelem::Int, ngl::Int;
                                 ltheta::Bool=true,
                                 lglobal_norms::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    ltheta || error(" compute_dsgs_viscosity!(::DSGS, ::NSD_3D): only the θ form is implemented in 3D.\n" *
                    "   Set :energy_equation => \"theta\", or use :visc_model => SMAG() / VREM() / AV().")

    neqs  = size(μ_dsgs, 2)
    invnp = one(TT)/(nelem*ngl*ngl*ngl)
    γ     = PhysConst.γ
    C0    = PhysConst.C0
    CR    = TT(1.0)
    Cmax  = TT(0.5)
    γm1   = γ - one(TT)
    eps   = TT(1.0e-16)

    # --- Pass 1: averages of (ρ, ρu, ρv, ρw, ρθ) — see _dsgs_norm_scope --
    ρ_avg  = zero(TT); ρu_avg = zero(TT); ρv_avg = zero(TT)
    ρw_avg = zero(TT); ρθ_avg = zero(TT)
    @inbounds for ie = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,k]
            ρ_avg  += q[ip,1]
            ρu_avg += q[ip,2]
            ρv_avg += q[ip,3]
            ρw_avg += q[ip,4]
            ρθ_avg += q[ip,5]
        end
    end
    if lglobal_norms
        sums = TT[ρ_avg, ρu_avg, ρv_avg, ρw_avg, ρθ_avg, TT(nelem*ngl*ngl*ngl)]
        MPI.Allreduce!(sums, MPI.SUM, get_mpi_comm())
        invnp_g = one(TT)/max(sums[6], one(TT))
        ρ_avg  = sums[1]*invnp_g; ρu_avg = sums[2]*invnp_g; ρv_avg = sums[3]*invnp_g
        ρw_avg = sums[4]*invnp_g; ρθ_avg = sums[5]*invnp_g
    else
        ρ_avg  *= invnp; ρu_avg *= invnp; ρv_avg *= invnp
        ρw_avg *= invnp; ρθ_avg *= invnp
    end

    # --- Pass 2: L∞ norms of |q - ⟨q⟩| ---------------------------------
    denom1 = zero(TT); denom2 = zero(TT); denom3 = zero(TT)
    denom4 = zero(TT); denom5 = zero(TT)
    @inbounds for ie = 1:nelem
        for k = 1:ngl, j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,k]
            denom1 = max(denom1, abs(q[ip,1] - ρ_avg))
            denom2 = max(denom2, abs(q[ip,2] - ρu_avg))
            denom3 = max(denom3, abs(q[ip,3] - ρv_avg))
            denom4 = max(denom4, abs(q[ip,4] - ρw_avg))
            denom5 = max(denom5, abs(q[ip,5] - ρθ_avg))
        end
    end
    if lglobal_norms
        norms = TT[denom1, denom2, denom3, denom4, denom5]
        MPI.Allreduce!(norms, MPI.MAX, get_mpi_comm())
        denom1 = norms[1]; denom2 = norms[2]; denom3 = norms[3]
        denom4 = norms[4]; denom5 = norms[5]
    end
    denom1 += eps; denom2 += eps; denom3 += eps
    denom4 += eps; denom5 += eps

    # Momentum floor, exactly as in 2D: the atmosphere starts globally at
    # rest, so the three momentum spreads start at zero and only `eps`
    # would separate the ratio from infinity — which pins ν at the
    # wave-speed cap on the very first stage.
    θ_avg  = ρθ_avg/max(abs(ρ_avg), eps)
    p_avg  = C0*(max(ρ_avg*θ_avg, zero(TT)))^γ
    c_avg  = sqrt(max(γ*p_avg/max(abs(ρ_avg), eps), zero(TT)))
    mom_floor = TT(1.0e-3) * abs(ρ_avg) * c_avg
    denom2 = max(denom2, mom_floor)
    denom3 = max(denom3, mom_floor)
    denom4 = max(denom4, mom_floor)

    # --- Pass 3: per-element residual L∞, μ_max bound, μ_dsgs[ie] ------
    @inbounds for ie = 1:nelem
        Δ = Δelem[ie]/ngl

        n1 = zero(TT); n2 = zero(TT); n3 = zero(TT)
        n4 = zero(TT); n5 = zero(TT)
        uTmx = zero(TT)
        ρ_el = zero(TT)

        for k = 1:ngl, j = 1:ngl
            @simd for i = 1:ngl
                ip  = connijk[ie,i,j,k]
                imK = one(TT)/(ω[i]*ω[j]*ω[k]*Je[ie,i,j,k])  # element lumped mass at the node

                R1 = abs((wt[1]*q[ip,1] + wt[2]*q1[ip,1] + wt[3]*q2[ip,1]) - imK*rhs_el[ie,i,j,k,1])
                R2 = abs((wt[1]*q[ip,2] + wt[2]*q1[ip,2] + wt[3]*q2[ip,2]) - imK*rhs_el[ie,i,j,k,2])
                R3 = abs((wt[1]*q[ip,3] + wt[2]*q1[ip,3] + wt[3]*q2[ip,3]) - imK*rhs_el[ie,i,j,k,3])
                R4 = abs((wt[1]*q[ip,4] + wt[2]*q1[ip,4] + wt[3]*q2[ip,4]) - imK*rhs_el[ie,i,j,k,4])
                R5 = abs((wt[1]*q[ip,5] + wt[2]*q1[ip,5] + wt[3]*q2[ip,5]) - imK*rhs_el[ie,i,j,k,5])
                n1 = max(n1, R1); n2 = max(n2, R2); n3 = max(n3, R3)
                n4 = max(n4, R4); n5 = max(n5, R5)

                ρl = q[ip,1]
                ul = q[ip,2]/ρl
                vl = q[ip,3]/ρl
                wl = q[ip,4]/ρl
                θl = q[ip,5]/ρl
                # p = C0·(ρθ)^γ  ⇒  c² = γp/ρ. ρθ is clamped at 0 so that a
                # solution already going negative is reported by the flux,
                # which says which equation broke, and not by a DomainError
                # raised inside the viscosity kernel.
                pl  = C0 * max(ρl*θl, zero(TT))^γ
                c_l = sqrt(max(γ*pl/ρl, zero(TT)))
                uTmx = max(uTmx, sqrt(ul*ul + vl*vl + wl*wl) + c_l)
                ρ_el += ρl
            end
        end
        ρ_el /= TT(ngl*ngl*ngl)

        μ_res = CR*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3, n4/denom4, n5/denom5)
        μ_max = Cmax*Δ*uTmx
        μ     = max(zero(TT), min(μ_max, μ_res))   # kinematic, m²/s
        μ_dyn = ρ_el*μ

        μ_dsgs[ie,1] = zero(TT)                             # ρ : no mass diffusion
        μ_dsgs[ie,2] = visc_coeff[2] * μ_dyn                # ρu (eq. 10a)
        μ_dsgs[ie,3] = visc_coeff[3] * μ_dyn                # ρv (eq. 10a)
        μ_dsgs[ie,4] = visc_coeff[4] * μ_dyn                # ρw (eq. 10a)
        μ_dsgs[ie,5] = visc_coeff[5] * (Pr/γm1) * μ_dyn     # ρθ (eq. 10b)
        for ieq = 6:neqs                                    # passive tracers: kinematic ν
            μ_dsgs[ie,ieq] = visc_coeff[ieq] * μ
        end
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
#     μ₁|K   = CR·h_K²·‖ρ−ρ̄‖_{∞,Ω}·max( ‖R_ρ‖_{∞,K}/‖ρ−ρ̄‖_{∞,Ω},
#                                        ‖R_m‖_{∞,K}/‖m−m̄‖_{∞,Ω},
#                                        ‖R_E‖_{∞,K}/‖E−Ē‖_{∞,Ω} )
#     μ_max|K = Cmax·h_K·‖ρ‖_{∞,K}·‖ |u| + √(γT) ‖_{∞,K}
#     μ|K     = min(μ_max|K, μ₁|K)
#     κ|K     = P/(γ−1)·μ|K            (heat conduction, on ∇T)
#     β|K     = μ|K/‖ρ‖_{∞,K}          (density diffusion, on ∇ρ)
#
# with CR = 1, Cmax = 0.5 and P ≈ 0.1 the artificial Prandtl number
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
# ⟨q⟩ and ‖q−⟨q⟩‖ are rank-local unless :dsgs_norms => "domain" — see
# _dsgs_norm_scope above. In the default (rank-local) mode everything in this
# routine is allocation-free, same discipline as the other implementations
# here; the global mode allocates the two small reduction buffers, once per
# RHS call and not per node.
# ================================================================================
function _dsgs_2d_energy!(μ_dsgs::AbstractMatrix{TT},
                          q::AbstractMatrix{TT},
                          q1::AbstractMatrix{TT},
                          q2::AbstractMatrix{TT},
                          rhs_el::AbstractArray{TT},
                          ω::AbstractVector{TT},
                          Je::AbstractArray{TT},
                          visc_coeff::AbstractVector{TT},
                          wt::NTuple{3,TT},
                          connijk::AbstractArray{TI,4},
                          Δelem::AbstractVector{TT},
                          PhysConst::PhysicalConst{TT},
                          Pr::TT,
                          nelem::Int, ngl::Int,
                          lglobal_norms::Bool) where {TT<:AbstractFloat, TI<:Integer}

    γ    = PhysConst.γ
    γm1  = γ - one(TT)
    CR   = TT(1.0)
    Cmax   = TT(0.5)
    neqs = size(μ_dsgs, 2)
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
                imK = one(TT)/(ω[i]*ω[j]*Je[ie,i,j])   # element lumped mass at the node

                Rρ  = abs((wt[1]*q[ip,1] + wt[2]*q1[ip,1] + wt[3]*q2[ip,1]) - imK*rhs_el[ie,i,j,1])
                Rmu = (wt[1]*q[ip,2] + wt[2]*q1[ip,2] + wt[3]*q2[ip,2]) - imK*rhs_el[ie,i,j,2]
                Rmv = (wt[1]*q[ip,3] + wt[2]*q1[ip,3] + wt[3]*q2[ip,3]) - imK*rhs_el[ie,i,j,3]
                Rm  = sqrt(Rmu*Rmu + Rmv*Rmv)
                RE  = abs((wt[1]*q[ip,4] + wt[2]*q1[ip,4] + wt[3]*q2[ip,4]) - imK*rhs_el[ie,i,j,4])

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
        μ_res = CR*h*h*dρ*ratio
        μ_cap = Cmax*h*ρmax*wmax
        μ     = max(zero(TT), min(μ_cap, μ_res))

        μ_dsgs[ie,1] = visc_coeff[1] * μ/max(ρmax, eps)   # β on ∇ρ
        μ_dsgs[ie,2] = visc_coeff[2] * μ                  # μ on ∇u
        μ_dsgs[ie,3] = visc_coeff[3] * μ                  # μ on ∇v
        μ_dsgs[ie,4] = visc_coeff[4] * (Pr/γm1) * μ       # κ on ∇T
        for ieq = 5:neqs                                  # passive tracers: kinematic ν
            μ_dsgs[ie,ieq] = visc_coeff[ieq] * μ/max(ρmax, eps)
        end
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
#     μ_res|e = CR · Δ² · max_i ( ‖R_i‖_{∞,e} / ‖q_i − ⟨q_i⟩‖_{∞,Ω} )
#     μ_max|e = Cmax · Δ  · (‖v‖ + c_f)_{∞,e}
#     μ|e     = max(0, min(μ_max, μ_res))
#
# with the BDF2 residual of equation i
#
#     R_i = ∂ₜq_i − M⁻¹·RHS_i,   ∂ₜq_i = wt₁ q_i + wt₂ q1_i + wt₃ q2_i
#
#     with (q1, q2, wt) the stage-consistent stencil chosen by rhs.jl
#     (_dsgs_stencil): at the first stage of a step the BDF2 of
#     (qⁿ, qⁿ⁻¹, qⁿ⁻²), at a later stage t = tⁿ + τ the second-order
#     three-point derivative at τ through (q(τ), qⁿ, qⁿ⁻¹).
#
# Notes specific to this implementation:
#
#  *  Element-wise residual.  R must have units of q/time for the ratio
#     R/‖q−⟨q⟩‖ to be a frequency and μ_res to come out as m²/s, and it
#     must be the ELEMENT's residual: rhs_el[K,i]/m_i^K, the element's weak
#     inviscid RHS over its lumped mass entry, against the assembled rate
#     ∂ₜq_i (see _dsgs_nodal_residual_1d! for why the assembled RHS is not
#     a residual at all with a lumped mass matrix).
#
#  *  Step-cadenced history.  params.qp.qnm1/qnm2 are advanced on every RK
#     *stage*, so they are stage snapshots, not states one Δt apart, and a
#     BDF2 stencil built on them does not approximate ∂q/∂t. DynSGS
#     therefore carries its own triple (params.dsgs_qn/qnm1/qnm2 = qⁿ, qⁿ⁻¹,
#     qⁿ⁻²), rolled once per time step by rhs!, and the stencil weights are
#     rebuilt at every stage from the stage time (rhs.jl, _dsgs_stencil):
#     a fixed BDF2 on (q_stage, qⁿ, qⁿ⁻¹) is right only at τ = Δt and reads
#     −∂ₜq/2 at the first stage, which fired the sensor on every smooth
#     moving structure (measured on the solitary wave of
#     ShallowWater/SoliWaveIslandDSGS: ν at the cap over the whole wave).
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
#     [1] ρ  : :μ[1]·μ                 — 0 (Marras) unless the case asks
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
# ⟨q⟩ and ‖q−⟨q⟩‖ are rank-local unless :dsgs_norms => "domain" — see
# _dsgs_norm_scope above. `comm` is what the global mode reduces over.
#
# Stratified atmospheres (problems/MHD/fluxEmergenceSon2025, eight decades
# of density between the photosphere and the corona) need two variants of
# the above, both off by default:
#
#  *  llocal_norms (:dsgs_norms => "element"). The residual of equation i is
#     normalized by the spread of q_i over the ELEMENT, ‖q_i − ⟨q_i⟩_e‖∞,e,
#     floored at local_rel (:dsgs_local_rel, default 1) times the
#     element-mean scales ρ_e, ρ_e c_e, ρ_e c_e², √ρ_e c_e, instead of the
#     domain spread. With the domain norm the dense bottom of the
#     atmosphere sets the scale of ρ, ρv and E, and a residual in the
#     corona — where those fields are 10⁻⁸ of it — is invisible: a
#     grid-scale sawtooth in the transition region grew unchecked with μ
#     at 10⁻¹¹ there. The element spread of a smooth stratified field is
#     O(q_i) (ρ changes by e⁻¹ across a 1 H₀ element), so the ratio stays
#     the relative under-resolution rate the model intends.
#
#  *  Cmin (:dsgs_Cmin). Background floor Cmin·Δ·(‖v‖+c_f) on μ, a fraction of the
#     Cmax cap, for the node-to-node modes the residual cannot sense (see the
#     kernel). 0 by default.
#
#  *  lconserved (:dsgs_conserved). Every slot receives the kinematic μ and
#     the case's user_primitives! hands the assembly the conserved variables
#     themselves, so the operator is a Laplacian on (ρ, ρv, E, B, ψ) — the
#     form in which a contact discontinuity (p continuous, ρ and T jumping)
#     diffuses consistently. rhs.jl drops the τ·u viscous-work term in this
#     mode (E already carries the dissipated kinetic energy) and the nodal-ρ
#     scaling below is not applied.
#
#  *  lnazarov_energy (:dsgs_nazarov_energy). The energy slot conducts heat
#     with Dao & Nazarov's κ = ρν/Pr (JSC 2022, §4.4) instead of the
#     Fourier-law κ = ρν γ/((γ−1)Pr) = c_p ρν/Pr of the default: in the
#     physical form the coefficient on ∇T (T = p/ρ) becomes ρν/Pr_t. In
#     the conserved form the energy flux is split (dsgs_split_energy):
#     the non-thermal part of E keeps ν, the thermal part p/(γ−1) gets
#     max(γ(γ−1)/Pr_t·ν_res, ν_floor), i.e. κ = ρν/Pr on the residual
#     viscosity with the Cmin floor intact. With :μ[1] = 1 (ν on ∇ρ)
#     the coefficients are then exactly their §4.4 set by equation: ν on
#     ρ, ρν on u, ρν/Pr on T, ν on B. (:dsgs_conserved_prandtl is accepted
#     as an alias.)
#
#  *  lnodal_rho (:dsgs_nodal_rho). Slots 2-5 receive the KINEMATIC μ (and
#     μγ/((γ−1)Pr_t) for E) and SGS_diffusion(::DSGS_MHD) multiplies by the
#     density OF THE QUADRATURE POINT. With the element mean ρ̄, the
#     effective diffusivity of u at the light side of an element is
#     (ρ̄/ρ)μ — up to 25μ across the chromosphere-corona transition — and
#     exceeds the explicit viscous stability limit as soon as the model
#     switches on there. The μ_dsgs output fields of slots 2-5 are then
#     kinematic too.
# ================================================================================
# TEMPORARY DIAGNOSTIC (JEXPRESSO_DSGS_DEBUG=1): per-equation maximum of the
# normalized residual, printed every 200 calls. Not for commit.
const _DSGS_DBG   = Ref(false)
const _DSGS_DBGN  = Ref(0)
const _DSGS_DBGV  = zeros(Float64, 8)
const _DSGS_DBGNU = Ref(0.0)
const _DSGS_DBGCAP = Ref(0.0)
const _DSGS_DBGLOC = zeros(Int, 4)     # (ie, i, j, ieq) of the largest ratio
const _DSGS_DBGTOP = Ref(0.0)

function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS_MHD, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 avg::AbstractVector{TT},
                                 denom::AbstractVector{TT},
                                 avg_e::AbstractVector{TT},
                                 den_e::AbstractVector{TT},
                                 wt::NTuple{3,TT},
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 γ::TT, Pr_t::TT, CR::TT, Cmax::TT,
                                 comm,
                                 nelem::Int, ngl::Int;
                                 lglobal_norms::Bool=false,
                                 llocal_norms::Bool=false,
                                 local_rel::TT=one(TT),
                                 lnodal_rho::Bool=false,
                                 lconserved::Bool=false,
                                 Cmin::TT=zero(TT),
                                 lnazarov_energy::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)          # residual max excludes the ψ slot
    rel  = TT(1.0e-3)            # floor fraction of the physical scales
    ldbg = get(ENV, "JEXPRESSO_DSGS_DEBUG", "") == "1"   # hoisted: no Ref read in the loops
    _DSGS_DBG[] = ldbg
    if ldbg
        fill!(_DSGS_DBGV, 0.0); _DSGS_DBGNU[] = 0.0; _DSGS_DBGCAP[] = 0.0
        _DSGS_DBGTOP[] = 0.0; fill!(_DSGS_DBGLOC, 0)
    end
    # avg_e / den_e: preallocated element mean / spread scratch (llocal_norms)
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
    @inbounds for ie = 1:nelem

        # Marras's element length scale: min edge / (N+1).
        Δ = Δelem[ie]/ngl

        ratio = zero(TT)      # max_i ‖R_i‖∞,e / denom_i
        wmax  = zero(TT)      # (‖v‖ + c_f)∞,e
        ρ_el  = zero(TT)      # element-mean density

        # Element-local normalization: ⟨q_i⟩_e and ‖q_i − ⟨q_i⟩_e‖∞,e, floored
        # at `rel` of the element-mean scales exactly as the domain norms are.
        if llocal_norms
            for ieq = 1:neqs
                avg_e[ieq] = zero(TT)
                den_e[ieq] = zero(TT)
            end
            for j = 1:ngl, i = 1:ngl
                ip = connijk[ie,i,j,1]
                for ieq = 1:neqs
                    avg_e[ieq] += q[ip,ieq]
                end
            end
            inv_ne = one(TT)/TT(ngl*ngl)
            for ieq = 1:neqs
                avg_e[ieq] *= inv_ne
            end
            for j = 1:ngl, i = 1:ngl
                ip = connijk[ie,i,j,1]
                for ieq = 1:neqs
                    den_e[ieq] = max(den_e[ieq], abs(q[ip,ieq] - avg_e[ieq]))
                end
            end
            # Floors at local_rel × the element's natural scales (ρ_e, ρ_e c_e,
            # ρ_e c_e², √ρ_e c_e). Unlike the domain norms, whose 10⁻³ floors
            # only guard against a degenerate spread, these floors ARE the
            # normalization of a quiescent element: the spread of ρv in an
            # atmosphere at rest is zero, and with a 10⁻³ floor a residual
            # of 2.5·10⁻³ ρc per unit time already drove μ to the wave-speed
            # cap over the whole quiet chromosphere of the flux-emergence
            # case (measured), whose sheet then eroded by resistive diffusion
            # (12% of its peak field in 2 τ₀) and sank. With local_rel = 1 the
            # ratio is the residual relative to the local physical rate ρc/τ,
            # which leaves a smooth settling flow at μ ~ 10⁻³ and still fires
            # on a grid-scale sawtooth (ratio ~ v_saw/Δ) or a shock (~ c/Δ).
            ρ_e = max(abs(avg_e[1]), eps)
            p_e = γm1*max(avg_e[4] - TT(0.5)*(avg_e[2]*avg_e[2] + avg_e[3]*avg_e[3] + avg_e[5]*avg_e[5])/ρ_e
                          - TT(0.5)*(avg_e[6]*avg_e[6] + avg_e[7]*avg_e[7] + avg_e[8]*avg_e[8]), zero(TT))
            c_e = sqrt(max(γ*p_e/ρ_e, eps))
            den_e[1] = max(den_e[1], local_rel*ρ_e)
            mom_e    = local_rel*ρ_e*c_e
            den_e[2] = max(den_e[2], mom_e)
            den_e[3] = max(den_e[3], mom_e)
            den_e[4] = max(den_e[4], local_rel*ρ_e*c_e*c_e)
            if neqs >= 5; den_e[5] = max(den_e[5], mom_e); end
            b_e = local_rel*sqrt(ρ_e)*c_e
            for ieq = 6:min(neqs,8)
                den_e[ieq] = max(den_e[ieq], b_e)
            end
            for ieq = 1:neqs
                den_e[ieq] += eps
            end
        end
        den = llocal_norms ? den_e : denom

        for j = 1:ngl
            for i = 1:ngl
                ip = connijk[ie,i,j,1]
                imK = one(TT)/(ω[i]*ω[j]*Je[ie,i,j])   # element lumped mass at the node

                for ieq = 1:NRES
                    R = abs((wt[1]*q[ip,ieq] + wt[2]*q1[ip,ieq] + wt[3]*q2[ip,ieq]) - imK*rhs_el[ie,i,j,ieq])
                    r = R/den[ieq]
                    ratio = max(ratio, r)
                    if ldbg
                        _DSGS_DBGV[ieq] = max(_DSGS_DBGV[ieq], Float64(r))
                        if Float64(r) > _DSGS_DBGTOP[]
                            _DSGS_DBGTOP[] = Float64(r)
                            _DSGS_DBGLOC[1] = ie; _DSGS_DBGLOC[2] = i
                            _DSGS_DBGLOC[3] = j;  _DSGS_DBGLOC[4] = ieq
                        end
                    end
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

        μ_res = CR*Δ*Δ*ratio
        μ_max = Cmax*Δ*wmax
        μ_c   = max(zero(TT), min(μ_max, μ_res))    # kinematic, m²/s (residual, capped)
        μ     = μ_c
        if ldbg
            _DSGS_DBGNU[]  = max(_DSGS_DBGNU[],  Float64(μ_c))
            _DSGS_DBGCAP[] = max(_DSGS_DBGCAP[], Float64(μ_max))
        end

        # Background floor Cmin·Δ·(‖v‖+c_f), a fraction of the wave-speed cap
        # (Cmin = 0 by default: pure Marras). The residual sensor is blind to a
        # node-to-node (checkerboard) mode — the discrete operator returns
        # nearly nothing on it, which is exactly why the CG discretization
        # leaves it undamped — and in the low-density corona of the
        # flux-emergence case such a mode grew from 0.02 to 0.5 C_s in three
        # τ₀ with μ_res ≈ 10⁻³ there. A floor of a few percent of the cap damps
        # it at rate Cmin Δ c (π/Δ)² ≈ 7/τ₀ for Cmin = 0.03 while diffusing a
        # resolved structure by only √(Cmin Δ c t) ≈ 1 H₀ over the whole run.
        μ_floor = Cmin > zero(TT) ? Cmin*Δ*wmax : zero(TT)
        μ = max(μ, μ_floor)

        # dynamic coefficient for u/v/w/T: ρ̄·μ with the element mean, or the
        # kinematic μ that SGS_diffusion(::DSGS_MHD) scales by the nodal ρ
        μ_dyn = lnodal_rho ? μ : ρ_el*μ

        # ρ: the user's :μ[1] times the kinematic μ (its primitive is ρ
        # itself, so ∇·(μ∇ρ) is a conservative mass diffusion). The MHD
        # cases at rest keep :μ[1] = 0; a stratified atmosphere whose
        # density drops 25× across an under-resolved transition region
        # needs it — the LGL undershoot of that contact goes below the 1e-8
        # of its light side as soon as it is displaced (fluxEmergenceSon2025).
        μ_dsgs[ie,1] = visc_coeff[1]*μ                             # ρ
        if lconserved
            # Laplacian on the CONSERVED variables (the case's user_primitives!
            # returns ρ, ρu, ρv, E, ρw, B, ψ themselves, rhs.jl drops the τ·u
            # term): one kinematic coefficient for every slot, so that an
            # isobaric contact diffuses consistently — ρ spreads, E (constant
            # across it) does not, and p stays what it was. Diffusing ρ alone
            # under a T-based energy closure had driven p negative within a
            # few τ₀ at the 25× density drop of the solar transition region.
            #
            # Coefficients by equation (Dao & Nazarov 2022, JSC 92:77, §4.4):
            # ONE kinematic ν from the max of the normalized residuals (their
            # eq. 4.8, this μ), and per equation ν on ∇ρ, μ = ρν in the
            # stress, κ = μ/Pr on ∇T, η = ν on B. In conserved variables
            # ν∇(ρu) ≈ ρν∇u and ν∇B are those already; the energy slot is
            # not: ν∇E carries the internal energy ρT/(γ(γ−1)) at ν, i.e. a
            # heat conduction κ = ρν/(γ(γ−1)), 19× Nazarov's ρν/Pr at
            # γ = 1.05, Pr = 1.
            # Energy slot. With lnazarov_energy the case splits the E
            # primitive (dsgs_split_energy above): this coefficient then
            # acts on the THERMAL part only — Dao & Nazarov's κ = ρν/Pr on
            # the residual viscosity, with the full background floor — and
            # the non-thermal part is diffused with the ρv slot's ν by
            # _expansion_visc!. Without it, ν on the whole of E.
            fE = lnazarov_energy ? γ*γm1/Pr_t : one(TT)
            μ_dsgs[ie,2] = visc_coeff[2]*μ                         # ρu
            μ_dsgs[ie,3] = visc_coeff[3]*μ                         # ρv
            μ_dsgs[ie,4] = lnazarov_energy ? visc_coeff[4]*max(fE*μ_c, μ_floor) :
                                             visc_coeff[4]*μ                       # E (thermal part if split)
            if neqs >= 5
                μ_dsgs[ie,5] = visc_coeff[5]*μ                     # ρw
            end
        else
            μ_dsgs[ie,2] = visc_coeff[2]*μ_dyn                     # ρu
            μ_dsgs[ie,3] = visc_coeff[3]*μ_dyn                     # ρv
            # κ∇T on the energy: the Fourier law κ = c_p ρν/Pr with
            # c_p = γ/(γ−1) (T = p/ρ), or Dao & Nazarov's κ = ρν/Pr
            μ_dsgs[ie,4] = lnazarov_energy ? visc_coeff[4]*μ_dyn/Pr_t :
                                             visc_coeff[4]*μ_dyn*γ/(γm1*Pr_t)   # E
            if neqs >= 5
                μ_dsgs[ie,5] = visc_coeff[5]*μ_dyn                 # ρw
            end
        end
        for ieq = 6:min(neqs,8)
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μ                     # B (resistivity)
        end
        if neqs >= 9
            μ_dsgs[ie,9] = visc_coeff[9]*μ                         # ψ
        end
    end

    if ldbg
        _DSGS_DBGN[] += 1
        if _DSGS_DBGN[] <= 12 || _DSGS_DBGN[] % 200 == 0
            edge = (_DSGS_DBGLOC[2] == 1 || _DSGS_DBGLOC[2] == ngl ||
                    _DSGS_DBGLOC[3] == 1 || _DSGS_DBGLOC[3] == ngl) ? "EDGE" : "int "
            @printf(" # DSGS dbg call %6d  nu_max=%.4e cap=%.4e  argmax: eq %d node (%d,%d) of %d %s  ratio by eq: %s   denom: %s\n",
                    _DSGS_DBGN[], _DSGS_DBGNU[], _DSGS_DBGCAP[],
                    _DSGS_DBGLOC[4], _DSGS_DBGLOC[2], _DSGS_DBGLOC[3], ngl, edge,
                    join((@sprintf("%.2e", _DSGS_DBGV[k]) for k = 1:NRES), " "),
                    join((@sprintf("%.2e", Float64(denom[k])) for k = 1:NRES), " "))
        end
    end
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_1D)
#
# The 1D version of the MHD kernel above, for the 8-variable system
# (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz) of a 1D ideal-MHD shock tube (Bx constant;
# a 9th ψ slot, if present, is carried but excluded from the residual max).
# Same residual (BDF2 history, max over the equations of the normalized
# residual — Dao & Nazarov 2022, eq. 4.8), same domain/element normalization,
# same cap Cmax·Δ·(|u| + c_f) with the fast magnetosonic speed and the same
# Cmin floor; see the 2D header for the meaning of every option. The
# coefficients by slot follow the 2D assignment: conserved form — one ν on
# every slot (the case's user_primitives! returns the conserved variables);
# physical form — ν on ρ, ρ̄ν on the momenta (u, v, w primitives),
# ρ̄ν γ/((γ−1)Pr_t) or, with lnazarov_energy, ρ̄ν/Pr_t on T = p/ρ, ν on B.
# The 1D viscous loop applies a scalar Laplacian per slot (no deviatoric
# stress, no τ·u), so the conserved form is the exactly conservative one and
# the default of problems/MHD/brioWu1d.
# ================================================================================
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS_MHD, ::NSD_1D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 avg::AbstractVector{TT},
                                 denom::AbstractVector{TT},
                                 avg_e::AbstractVector{TT},
                                 den_e::AbstractVector{TT},
                                 wt::NTuple{3,TT},
                                 connijk::AbstractArray{TI,4},
                                 Δx::AbstractVector{TT},
                                 γ::TT, Pr_t::TT, CR::TT, Cmax::TT,
                                 comm,
                                 nelem::Int, ngl::Int;
                                 lglobal_norms::Bool=false,
                                 llocal_norms::Bool=false,
                                 local_rel::TT=one(TT),
                                 lnodal_rho::Bool=false,
                                 lconserved::Bool=false,
                                 Cmin::TT=zero(TT),
                                 lnazarov_energy::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)
    rel  = TT(1.0e-3)
    γm1  = γ - one(TT)
    eps  = TT(1.0e-16)

    @inline function pres(ρ, mu, mv, E, mw, bx, by, bz)
        ρp = max(ρ, eps)
        return γm1*max(E - TT(0.5)*(mu*mu + mv*mv + mw*mw)/ρp - TT(0.5)*(bx*bx + by*by + bz*bz), zero(TT))
    end

    # --- Pass 1: means and spreads over the (rank-local or global) domain
    @inbounds for ieq = 1:neqs
        avg[ieq] = zero(TT); denom[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, i = 1:ngl
        ip = connijk[ie,i,1,1]
        for ieq = 1:neqs
            avg[ieq] += q[ip,ieq]
        end
    end
    inv_npts = one(TT)/max(TT(nelem*ngl), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(nelem*ngl), MPI.SUM, comm)
        MPI.Allreduce!(avg, MPI.SUM, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end
    @inbounds for ie = 1:nelem, i = 1:ngl
        ip = connijk[ie,i,1,1]
        for ieq = 1:neqs
            denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
        end
    end
    if lglobal_norms
        MPI.Allreduce!(denom, MPI.MAX, comm)
    end
    ρ_avg = max(abs(avg[1]), eps)
    p_avg = pres(avg[1], avg[2], avg[3], avg[4], (neqs >= 5 ? avg[5] : zero(TT)),
                 (neqs >= 6 ? avg[6] : zero(TT)), (neqs >= 7 ? avg[7] : zero(TT)), (neqs >= 8 ? avg[8] : zero(TT)))
    c_avg = sqrt(max(γ*p_avg/ρ_avg, eps))
    @inbounds begin
        denom[1] = max(denom[1], rel*ρ_avg)
        mom_fl   = rel*ρ_avg*c_avg
        denom[2] = max(denom[2], mom_fl)
        denom[3] = max(denom[3], mom_fl)
        denom[4] = max(denom[4], rel*ρ_avg*c_avg*c_avg)
        if neqs >= 5; denom[5] = max(denom[5], mom_fl); end
        b_fl = rel*sqrt(ρ_avg)*c_avg
        for ieq = 6:min(neqs,8)
            denom[ieq] = max(denom[ieq], b_fl)
        end
        for ieq = 1:neqs
            denom[ieq] += eps
        end
    end

    # --- Pass 2: per element ------------------------------------------
    @inbounds for ie = 1:nelem
        Δ = Δx[ie]/ngl
        ratio = zero(TT)
        wmax  = zero(TT)
        ρ_el  = zero(TT)
        if llocal_norms
            for ieq = 1:neqs
                avg_e[ieq] = zero(TT); den_e[ieq] = zero(TT)
            end
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                for ieq = 1:neqs
                    avg_e[ieq] += q[ip,ieq]
                end
            end
            inv_ne = one(TT)/TT(ngl)
            for ieq = 1:neqs
                avg_e[ieq] *= inv_ne
            end
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                for ieq = 1:neqs
                    den_e[ieq] = max(den_e[ieq], abs(q[ip,ieq] - avg_e[ieq]))
                end
            end
            ρ_e = max(abs(avg_e[1]), eps)
            p_e = pres(avg_e[1], avg_e[2], avg_e[3], avg_e[4], (neqs >= 5 ? avg_e[5] : zero(TT)),
                       (neqs >= 6 ? avg_e[6] : zero(TT)), (neqs >= 7 ? avg_e[7] : zero(TT)), (neqs >= 8 ? avg_e[8] : zero(TT)))
            c_e = sqrt(max(γ*p_e/ρ_e, eps))
            den_e[1] = max(den_e[1], local_rel*ρ_e)
            mom_e    = local_rel*ρ_e*c_e
            den_e[2] = max(den_e[2], mom_e)
            den_e[3] = max(den_e[3], mom_e)
            den_e[4] = max(den_e[4], local_rel*ρ_e*c_e*c_e)
            if neqs >= 5; den_e[5] = max(den_e[5], mom_e); end
            b_e = local_rel*sqrt(ρ_e)*c_e
            for ieq = 6:min(neqs,8)
                den_e[ieq] = max(den_e[ieq], b_e)
            end
            for ieq = 1:neqs
                den_e[ieq] += eps
            end
        end
        den = llocal_norms ? den_e : denom

        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            imK = one(TT)/(ω[i]*Je[ie,i])   # element lumped mass at the node
            for ieq = 1:NRES
                R = abs((wt[1]*q[ip,ieq] + wt[2]*q1[ip,ieq] + wt[3]*q2[ip,ieq]) - imK*rhs_el[ie,i,ieq])
                ratio = max(ratio, R/den[ieq])
            end
            ρl = max(q[ip,1], eps)
            ul = q[ip,2]/ρl
            vl = q[ip,3]/ρl
            wl = (neqs >= 5) ? q[ip,5]/ρl : zero(TT)
            bx = (neqs >= 6) ? q[ip,6] : zero(TT)
            by = (neqs >= 7) ? q[ip,7] : zero(TT)
            bz = (neqs >= 8) ? q[ip,8] : zero(TT)
            B2 = bx*bx + by*by + bz*bz
            pl = pres(q[ip,1], q[ip,2], q[ip,3], q[ip,4], (neqs >= 5 ? q[ip,5] : zero(TT)), bx, by, bz)
            # fast magnetosonic speed along x (the only propagation direction)
            a2  = γ*pl/ρl
            b2  = B2/ρl
            bx2 = bx*bx/ρl
            cf  = sqrt(max(TT(0.5)*(a2 + b2 + sqrt(max((a2 + b2)*(a2 + b2) - 4*a2*bx2, zero(TT)))), zero(TT)))
            wmax  = max(wmax, sqrt(ul*ul + vl*vl + wl*wl) + cf)
            ρ_el += ρl
        end
        ρ_el /= TT(ngl)

        μ_res = CR*Δ*Δ*ratio
        μ_max = Cmax*Δ*wmax
        μ_c   = max(zero(TT), min(μ_max, μ_res))
        μ_fl  = Cmin > zero(TT) ? Cmin*Δ*wmax : zero(TT)
        μ     = max(μ_c, μ_fl)
        μ_dyn = lnodal_rho ? μ : ρ_el*μ

        μ_dsgs[ie,1] = visc_coeff[1]*μ
        if lconserved
            μ_dsgs[ie,2] = visc_coeff[2]*μ
            μ_dsgs[ie,3] = visc_coeff[3]*μ
            μ_dsgs[ie,4] = visc_coeff[4]*μ
            if neqs >= 5; μ_dsgs[ie,5] = visc_coeff[5]*μ; end
        else
            μ_dsgs[ie,2] = visc_coeff[2]*μ_dyn
            μ_dsgs[ie,3] = visc_coeff[3]*μ_dyn
            μ_dsgs[ie,4] = lnazarov_energy ? visc_coeff[4]*μ_dyn/Pr_t : visc_coeff[4]*μ_dyn*γ/(γm1*Pr_t)
            if neqs >= 5; μ_dsgs[ie,5] = visc_coeff[5]*μ_dyn; end
        end
        for ieq = 6:min(neqs,8)
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*μ
        end
        if neqs >= 9
            μ_dsgs[ie,9] = visc_coeff[9]*μ
        end
    end
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity_nodal!(::DSGS_MHD, ::NSD_1D)
#
# The NODAL form of the 1D MHD kernel, i.e. the residual viscosity of Dao &
# Nazarov (2022, JSC 92:77, §4.2–4.3) as they define it: the residual is the
# assembled (lumped-mass) nodal residual R_i = |BDF2(q)_i − M⁻¹_i rhs_i|, the
# normalization n(w)_i = S̄(w)·(1 − C_l·(local range of w over the support
# of node i)/(global range of w)) (their eq. 4.7, S̄ the global spread; C_l = 0
# is the classical S̄), the residual ratio R_i = max over the equations of
# R_i/n_i (eq. 4.8, with the n²/(n²+ε) guard), and
#
#     ν_i = min(C_max h_i λ_max,i,  C_R h_i² R_i)      (eq. 4.10),   h_i = Δ_K/(k+1)
#
# at every node, floored at Cmin h_i λ_i as in the element kernel. ν is then a
# continuous (C⁰, "DSS'd") field: the element loop interpolates the nodal
# values, so the diffusive flux ∂x(ν ∂x q) has no jump at element interfaces.
# The element form above takes the maximum over each element and applies one
# ν per element; the resulting staircase in ν kinks the flux at every element
# boundary, which on the Brio–Wu tube appeared as one wiggle per element in
# the smooth plateau behind the compound wave (measured). The slot assignment
# is the element kernel's, with nodal ρ in the physical form. μ_dsgs[ie,:]
# receives the element mean of the nodal values for the output staircase.
# ================================================================================
function compute_dsgs_viscosity_nodal!(μ_dsgs::AbstractMatrix{TT},
                                       μ_pnode::AbstractMatrix{TT},
                                       ::DSGS_MHD, ::NSD_1D,
                                       q::AbstractMatrix{TT},
                                       q1::AbstractMatrix{TT},
                                       q2::AbstractMatrix{TT},
                                       rhs_el::AbstractArray{TT},
                                       ω::AbstractVector{TT},
                                       Je::AbstractArray{TT},
                                       visc_coeff::AbstractVector{TT},
                                       avg::AbstractVector{TT},
                                       denom::AbstractVector{TT},
                                       qmin::AbstractVector{TT},
                                       qmax::AbstractVector{TT},
                                       nmin::AbstractMatrix{TT},
                                       nmax::AbstractMatrix{TT},
                                       hnod::AbstractVector{TT},
                                       Rnod::AbstractMatrix{TT},
                                       mnod::AbstractVector{TT},
                                       wt::NTuple{3,TT},
                                       connijk::AbstractArray{TI,4},
                                       Δx::AbstractVector{TT},
                                       γ::TT, Pr_t::TT, CR::TT, Cmax::TT, Cl::TT,
                                       comm,
                                       nelem::Int, ngl::Int, npoin::Int;
                                       lglobal_norms::Bool=false,
                                       lconserved::Bool=false,
                                       Cmin::TT=zero(TT),
                                       lnazarov_energy::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)
    rel  = TT(1.0e-3)
    γm1  = γ - one(TT)
    eps  = TT(1.0e-16)
    k    = max(ngl - 1, 1)          # polynomial degree

    @inline function pres(ρ, mu, mv, E, mw, bx, by, bz)
        ρp = max(ρ, eps)
        return γm1*max(E - TT(0.5)*(mu*mu + mv*mv + mw*mw)/ρp - TT(0.5)*(bx*bx + by*by + bz*bz), zero(TT))
    end

    # --- global mean, spread S̄ and range per equation --------------------
    @inbounds for ieq = 1:neqs
        avg[ieq] = zero(TT); denom[ieq] = zero(TT)
        qmin[ieq] = typemax(TT); qmax[ieq] = typemin(TT)
    end
    @inbounds for ip = 1:npoin, ieq = 1:neqs
        avg[ieq] += q[ip,ieq]
        qmin[ieq] = min(qmin[ieq], q[ip,ieq]); qmax[ieq] = max(qmax[ieq], q[ip,ieq])
    end
    inv_npts = one(TT)/max(TT(npoin), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(npoin), MPI.SUM, comm)
        MPI.Allreduce!(avg, MPI.SUM, comm)
        MPI.Allreduce!(qmin, MPI.MIN, comm); MPI.Allreduce!(qmax, MPI.MAX, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end
    @inbounds for ip = 1:npoin, ieq = 1:neqs
        denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - avg[ieq]))
    end
    if lglobal_norms
        MPI.Allreduce!(denom, MPI.MAX, comm)
    end
    ρ_avg = max(abs(avg[1]), eps)
    p_avg = pres(avg[1], avg[2], avg[3], avg[4], (neqs >= 5 ? avg[5] : zero(TT)),
                 (neqs >= 6 ? avg[6] : zero(TT)), (neqs >= 7 ? avg[7] : zero(TT)), (neqs >= 8 ? avg[8] : zero(TT)))
    c_avg = sqrt(max(γ*p_avg/ρ_avg, eps))
    @inbounds begin
        denom[1] = max(denom[1], rel*ρ_avg)
        mom_fl   = rel*ρ_avg*c_avg
        denom[2] = max(denom[2], mom_fl)
        denom[3] = max(denom[3], mom_fl)
        denom[4] = max(denom[4], rel*ρ_avg*c_avg*c_avg)
        if neqs >= 5; denom[5] = max(denom[5], mom_fl); end
        b_fl = rel*sqrt(ρ_avg)*c_avg
        for ieq = 6:min(neqs,8)
            denom[ieq] = max(denom[ieq], b_fl)
        end
    end

    # --- local range over the support of each node (eq. 4.7) and h_i -------
    @inbounds for ip = 1:npoin
        hnod[ip] = zero(TT)
        for ieq = 1:neqs
            nmin[ip,ieq] = typemax(TT); nmax[ip,ieq] = typemin(TT)
        end
    end
    @inbounds for ie = 1:nelem
        # Mesh function of the node: the element form's Δ = Δ_K/(k+1)
        # (Marras), the same length the element kernel uses, so the two
        # forms differ only in where the coefficient lives. Dao & Nazarov's
        # h_K/k with h_K the circumradius is Δ_K/(√2 k) on a square, within
        # 12% of Δ_K/(k+1) at k = 4; Δ_K/k (25% larger) drove the nodal
        # coefficient of the rising-bubble case 1.56× above the element
        # form's at start-up and past the explicit viscous limit.
        h_e = Δx[ie]/TT(ngl)
        for ieq = 1:neqs
            emin = typemax(TT); emax = typemin(TT)
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                emin = min(emin, q[ip,ieq]); emax = max(emax, q[ip,ieq])
            end
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                nmin[ip,ieq] = min(nmin[ip,ieq], emin); nmax[ip,ieq] = max(nmax[ip,ieq], emax)
            end
        end
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            hnod[ip] = max(hnod[ip], h_e)
        end
    end

    # --- nodal viscosity ---------------------------------------------------
    _dsgs_nodal_residual_1d!(Rnod, mnod, q, q1, q2, wt, rhs_el, ω, Je, connijk, nelem, ngl, npoin, NRES)
    @inbounds for ip = 1:npoin
        ratio = _dsgs_nodal_ratio(Rnod, ip, NRES, denom, qmin, qmax, nmin, nmax, Cl, eps)
        ρl = max(q[ip,1], eps)
        ul = q[ip,2]/ρl
        vl = q[ip,3]/ρl
        wl = (neqs >= 5) ? q[ip,5]/ρl : zero(TT)
        bx = (neqs >= 6) ? q[ip,6] : zero(TT)
        by = (neqs >= 7) ? q[ip,7] : zero(TT)
        bz = (neqs >= 8) ? q[ip,8] : zero(TT)
        B2 = bx*bx + by*by + bz*bz
        pl = pres(q[ip,1], q[ip,2], q[ip,3], q[ip,4], (neqs >= 5 ? q[ip,5] : zero(TT)), bx, by, bz)
        a2  = γ*pl/ρl
        b2  = B2/ρl
        bx2 = bx*bx/ρl
        cf  = sqrt(max(TT(0.5)*(a2 + b2 + sqrt(max((a2 + b2)*(a2 + b2) - 4*a2*bx2, zero(TT)))), zero(TT)))
        λ   = sqrt(ul*ul + vl*vl + wl*wl) + cf
        h   = hnod[ip]

        ν_c = max(zero(TT), min(Cmax*h*λ, CR*h*h*ratio))
        ν   = Cmin > zero(TT) ? max(ν_c, Cmin*h*λ) : ν_c
        ν_d = ρl*ν

        μ_pnode[ip,1] = visc_coeff[1]*ν
        if lconserved
            μ_pnode[ip,2] = visc_coeff[2]*ν
            μ_pnode[ip,3] = visc_coeff[3]*ν
            μ_pnode[ip,4] = visc_coeff[4]*ν
            if neqs >= 5; μ_pnode[ip,5] = visc_coeff[5]*ν; end
        else
            μ_pnode[ip,2] = visc_coeff[2]*ν_d
            μ_pnode[ip,3] = visc_coeff[3]*ν_d
            μ_pnode[ip,4] = lnazarov_energy ? visc_coeff[4]*ν_d/Pr_t : visc_coeff[4]*ν_d*γ/(γm1*Pr_t)
            if neqs >= 5; μ_pnode[ip,5] = visc_coeff[5]*ν_d; end
        end
        for ieq = 6:min(neqs,8)
            μ_pnode[ip,ieq] = visc_coeff[ieq]*ν
        end
        if neqs >= 9
            μ_pnode[ip,9] = visc_coeff[9]*ν
        end
    end

    # --- element means for the staircase output ----------------------------
    inv_ngl = one(TT)/TT(ngl)
    @inbounds for ie = 1:nelem, ieq = 1:neqs
        m = zero(TT)
        for i = 1:ngl
            m += μ_pnode[connijk[ie,i,1,1], ieq]
        end
        μ_dsgs[ie,ieq] = m*inv_ngl
    end
    return nothing
end

# ================================================================================
# Nodal (Dao & Nazarov 2022) form in 2D — shared statistics
#
# Global mean, spread S̄ and range of every equation over the local nodes,
# the local range of every equation over the support of each node (the
# elements containing it, their eq. 4.7) and the nodal mesh function
# h_i = max Δ_K/(k+1) over those elements. All buffers are preallocated
# (params.dsgs_*), nothing is allocated here.
# ================================================================================
function _dsgs_nodal_stats_2d!(q::AbstractMatrix{TT},
                               connijk::AbstractArray{TI,4},
                               Δelem::AbstractVector{TT},
                               nelem::Int, ngl::Int, npoin::Int, neqs::Int, k::Int,
                               avg::AbstractVector{TT}, denom::AbstractVector{TT},
                               qmin::AbstractVector{TT}, qmax::AbstractVector{TT},
                               nmin::AbstractMatrix{TT}, nmax::AbstractMatrix{TT},
                               hnod::AbstractVector{TT},
                               lglobal_norms::Bool, comm;
                               qe::Union{Nothing,AbstractMatrix{TT}}=nothing) where {TT<:AbstractFloat, TI<:Integer}
    # qe given: every statistic is taken on the departure q − qe from the
    # reference state (the shallow-water kernel: the lake at rest)
    @inline dq(ip, ieq) = qe === nothing ? q[ip,ieq] : q[ip,ieq] - qe[ip,ieq]
    @inbounds for ieq = 1:neqs
        avg[ieq] = zero(TT); denom[ieq] = zero(TT)
        qmin[ieq] = typemax(TT); qmax[ieq] = typemin(TT)
    end
    @inbounds for ip = 1:npoin, ieq = 1:neqs
        v = dq(ip, ieq)
        avg[ieq] += v
        qmin[ieq] = min(qmin[ieq], v); qmax[ieq] = max(qmax[ieq], v)
    end
    inv_npts = one(TT)/max(TT(npoin), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(npoin), MPI.SUM, comm)
        MPI.Allreduce!(avg, MPI.SUM, comm)
        MPI.Allreduce!(qmin, MPI.MIN, comm); MPI.Allreduce!(qmax, MPI.MAX, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end
    @inbounds for ip = 1:npoin, ieq = 1:neqs
        denom[ieq] = max(denom[ieq], abs(dq(ip, ieq) - avg[ieq]))
    end
    if lglobal_norms
        MPI.Allreduce!(denom, MPI.MAX, comm)
    end
    @inbounds for ip = 1:npoin
        hnod[ip] = zero(TT)
        for ieq = 1:neqs
            nmin[ip,ieq] = typemax(TT); nmax[ip,ieq] = typemin(TT)
        end
    end
    @inbounds for ie = 1:nelem
        h_e = Δelem[ie]/TT(ngl)      # Δ_K/(k+1), the element form's Δ (see the 1D kernel)
        for ieq = 1:neqs
            emin = typemax(TT); emax = typemin(TT)
            for j = 1:ngl, i = 1:ngl
                ip = connijk[ie,i,j,1]
                v  = dq(ip, ieq)
                emin = min(emin, v); emax = max(emax, v)
            end
            for j = 1:ngl, i = 1:ngl
                ip = connijk[ie,i,j,1]
                nmin[ip,ieq] = min(nmin[ip,ieq], emin); nmax[ip,ieq] = max(nmax[ip,ieq], emax)
            end
        end
        for j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,1]
            hnod[ip] = max(hnod[ip], h_e)
        end
    end
    return nothing
end

# Nodal residual ratio max_i R_i/n_i (eq. 4.8) with the local-jump
# normalization (eq. 4.7) and the n²/(n²+ε) guard.
# ================================================================================
# The element-wise residual.
#
# With a lumped (LGL-collocated) mass matrix the assembled rate of change
# ∂ₜq_i = M⁻¹RHS_i IS what the integrator advances, so a residual built from
# the assembled RHS is the time-integration error and nothing else: it
# vanishes on a resolved AND on an under-resolved solution alike (measured on
# sod1d: ν at the shock 0.5 % of the cap). The residual of the method is the
# one-sided one of each element,
#
#     R_i^K = | ∂ₜq_i − rhs_el[K,i]/m_i^K |,    m_i^K = ω_i (ω_j) J_K,i(,j),
#
# rhs_el the element's own weak inviscid RHS (fluxes and sources) before the
# direct stiffness summation and m_i^K its lumped mass entry: at the interior
# nodes of K it equals the assembled one (the SBP property of the LGL
# operator makes rhs_el/m the strong nodal divergence), at the interface
# nodes it differs from the assembled rate by the mass-weighted JUMP of the
# flux divergence across the interface — O(h^k) where the solution is
# smooth, O(1/h) at a discontinuity, which is what the sensor is for. This
# is Dao & Nazarov's (1/m_i)∫|BDF(q) + ∇·f(q)|φ_i with the absolute value
# inside the integral, evaluated with the LGL rule; the element kernels take
# max_i R_i^K over the element, the nodal kernels the mass-weighted average
# Σ_K m_i^K R_i^K / Σ_K m_i^K over the elements that contain the node.
# ================================================================================
function _dsgs_nodal_residual_1d!(Rnod::AbstractMatrix{TT}, mnod::AbstractVector{TT},
                                  q::AbstractMatrix{TT}, q1::AbstractMatrix{TT}, q2::AbstractMatrix{TT},
                                  wt::NTuple{3,TT}, rhs_el::AbstractArray{TT},
                                  ω::AbstractVector{TT}, Je::AbstractArray{TT},
                                  connijk::AbstractArray{TI,4}, nelem::Int, ngl::Int, npoin::Int,
                                  NRES::Int) where {TT<:AbstractFloat, TI<:Integer}
    @inbounds for ip = 1:npoin
        mnod[ip] = zero(TT)
        for ieq = 1:NRES
            Rnod[ip,ieq] = zero(TT)
        end
    end
    @inbounds for ie = 1:nelem, i = 1:ngl
        ip = connijk[ie,i,1,1]
        m  = ω[i]*Je[ie,i]
        mnod[ip] += m
        for ieq = 1:NRES
            Rnod[ip,ieq] += abs((wt[1]*q[ip,ieq] + wt[2]*q1[ip,ieq] + wt[3]*q2[ip,ieq])*m - rhs_el[ie,i,ieq])
        end
    end
    @inbounds for ip = 1:npoin
        im = one(TT)/max(mnod[ip], TT(1.0e-300))
        for ieq = 1:NRES
            Rnod[ip,ieq] *= im
        end
    end
    return nothing
end

function _dsgs_nodal_residual_2d!(Rnod::AbstractMatrix{TT}, mnod::AbstractVector{TT},
                                  q::AbstractMatrix{TT}, q1::AbstractMatrix{TT}, q2::AbstractMatrix{TT},
                                  wt::NTuple{3,TT}, rhs_el::AbstractArray{TT},
                                  ω::AbstractVector{TT}, Je::AbstractArray{TT},
                                  connijk::AbstractArray{TI,4}, nelem::Int, ngl::Int, npoin::Int,
                                  NRES::Int) where {TT<:AbstractFloat, TI<:Integer}
    @inbounds for ip = 1:npoin
        mnod[ip] = zero(TT)
        for ieq = 1:NRES
            Rnod[ip,ieq] = zero(TT)
        end
    end
    @inbounds for ie = 1:nelem, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,1]
        m  = ω[i]*ω[j]*Je[ie,i,j]
        mnod[ip] += m
        for ieq = 1:NRES
            Rnod[ip,ieq] += abs((wt[1]*q[ip,ieq] + wt[2]*q1[ip,ieq] + wt[3]*q2[ip,ieq])*m - rhs_el[ie,i,j,ieq])
        end
    end
    @inbounds for ip = 1:npoin
        im = one(TT)/max(mnod[ip], TT(1.0e-300))
        for ieq = 1:NRES
            Rnod[ip,ieq] *= im
        end
    end
    return nothing
end

@inline function _dsgs_nodal_ratio(Rnod, ip, NRES, denom, qmin, qmax, nmin, nmax, Cl, eps)
    TT = eltype(denom)
    ratio = zero(TT)
    @inbounds for ieq = 1:NRES
        R = Rnod[ip,ieq]
        grange = qmax[ieq] - qmin[ieq]
        lfac   = grange > eps ? Cl*(nmax[ip,ieq] - nmin[ip,ieq])/grange : zero(TT)
        n      = denom[ieq]*(one(TT) - lfac)
        ratio  = max(ratio, R*n/(n*n + eps))
    end
    return ratio
end

# Element means of a nodal coefficient, for the output staircase.
function _dsgs_nodal_to_elements_2d!(μ_dsgs::AbstractMatrix{TT}, μ_pnode::AbstractMatrix{TT},
                                     connijk::AbstractArray{TI,4}, nelem::Int, ngl::Int) where {TT<:AbstractFloat, TI<:Integer}
    neqs = size(μ_dsgs, 2)
    inv_n = one(TT)/TT(ngl*ngl)
    @inbounds for ie = 1:nelem, ieq = 1:neqs
        m = zero(TT)
        for j = 1:ngl, i = 1:ngl
            m += μ_pnode[connijk[ie,i,j,1], ieq]
        end
        μ_dsgs[ie,ieq] = m*inv_n
    end
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity_nodal!(::DSGS_MHD, ::NSD_2D)
#
# The nodal form of the 2D MHD kernel (see the 1D header above): ν_i at every
# node from the assembled residual, the 2D fast speed √(γp/ρ + |B|²/ρ) in the
# cap, the Cmin floor, and the element kernel's slot assignment with the
# NODAL density in the physical form. μ_pnode receives the coefficients
# (no broadcast afterwards), μ_dsgs the element means.
# ================================================================================
function compute_dsgs_viscosity_nodal!(μ_dsgs::AbstractMatrix{TT},
                                       μ_pnode::AbstractMatrix{TT},
                                       ::DSGS_MHD, ::NSD_2D,
                                       q::AbstractMatrix{TT},
                                       q1::AbstractMatrix{TT},
                                       q2::AbstractMatrix{TT},
                                       rhs_el::AbstractArray{TT},
                                       ω::AbstractVector{TT},
                                       Je::AbstractArray{TT},
                                       visc_coeff::AbstractVector{TT},
                                       avg::AbstractVector{TT},
                                       denom::AbstractVector{TT},
                                       qmin::AbstractVector{TT},
                                       qmax::AbstractVector{TT},
                                       nmin::AbstractMatrix{TT},
                                       nmax::AbstractMatrix{TT},
                                       hnod::AbstractVector{TT},
                                       Rnod::AbstractMatrix{TT},
                                       mnod::AbstractVector{TT},
                                       wt::NTuple{3,TT},
                                       connijk::AbstractArray{TI,4},
                                       Δelem::AbstractVector{TT},
                                       γ::TT, Pr_t::TT, CR::TT, Cmax::TT, Cl::TT,
                                       comm,
                                       nelem::Int, ngl::Int, npoin::Int;
                                       lglobal_norms::Bool=false,
                                       lconserved::Bool=false,
                                       Cmin::TT=zero(TT),
                                       lnazarov_energy::Bool=false) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)
    rel  = TT(1.0e-3)
    γm1  = γ - one(TT)
    eps  = TT(1.0e-16)
    k    = max(ngl - 1, 1)

    _dsgs_nodal_stats_2d!(q, connijk, Δelem, nelem, ngl, npoin, neqs, k,
                          avg, denom, qmin, qmax, nmin, nmax, hnod, lglobal_norms, comm)

    ρ_avg = max(abs(avg[1]), eps)
    p_avg = γm1*max(avg[4] - TT(0.5)*(avg[2]*avg[2] + avg[3]*avg[3] + avg[5]*avg[5])/ρ_avg
                    - TT(0.5)*(avg[6]*avg[6] + avg[7]*avg[7] + avg[8]*avg[8]), zero(TT))
    c_avg = sqrt(max(γ*p_avg/ρ_avg, eps))
    @inbounds begin
        denom[1] = max(denom[1], rel*ρ_avg)
        mom_fl   = rel*ρ_avg*c_avg
        denom[2] = max(denom[2], mom_fl)
        denom[3] = max(denom[3], mom_fl)
        denom[4] = max(denom[4], rel*ρ_avg*c_avg*c_avg)
        if neqs >= 5; denom[5] = max(denom[5], mom_fl); end
        b_fl = rel*sqrt(ρ_avg)*c_avg
        for ieq = 6:min(neqs,8)
            denom[ieq] = max(denom[ieq], b_fl)
        end
    end

    fE = lnazarov_energy ? γ*γm1/Pr_t : one(TT)
    _dsgs_nodal_residual_2d!(Rnod, mnod, q, q1, q2, wt, rhs_el, ω, Je, connijk, nelem, ngl, npoin, NRES)
    @inbounds for ip = 1:npoin
        ratio = _dsgs_nodal_ratio(Rnod, ip, NRES, denom, qmin, qmax, nmin, nmax, Cl, eps)
        ρl = max(q[ip,1], eps)
        ul = q[ip,2]/ρl
        vl = q[ip,3]/ρl
        wl = (neqs >= 5) ? q[ip,5]/ρl : zero(TT)
        B2 = q[ip,6]*q[ip,6] + q[ip,7]*q[ip,7] + q[ip,8]*q[ip,8]
        pl = γm1*max(q[ip,4] - TT(0.5)*ρl*(ul*ul + vl*vl + wl*wl) - TT(0.5)*B2, zero(TT))
        cf = sqrt(max(γ*pl/ρl + B2/ρl, zero(TT)))
        λ  = sqrt(ul*ul + vl*vl + wl*wl) + cf
        h  = hnod[ip]

        ν_c = max(zero(TT), min(Cmax*h*λ, CR*h*h*ratio))
        ν_f = Cmin > zero(TT) ? Cmin*h*λ : zero(TT)
        ν   = max(ν_c, ν_f)
        ν_d = ρl*ν

        μ_pnode[ip,1] = visc_coeff[1]*ν
        if lconserved
            μ_pnode[ip,2] = visc_coeff[2]*ν
            μ_pnode[ip,3] = visc_coeff[3]*ν
            μ_pnode[ip,4] = lnazarov_energy ? visc_coeff[4]*max(fE*ν_c, ν_f) : visc_coeff[4]*ν
            if neqs >= 5; μ_pnode[ip,5] = visc_coeff[5]*ν; end
        else
            μ_pnode[ip,2] = visc_coeff[2]*ν_d
            μ_pnode[ip,3] = visc_coeff[3]*ν_d
            μ_pnode[ip,4] = lnazarov_energy ? visc_coeff[4]*ν_d/Pr_t : visc_coeff[4]*ν_d*γ/(γm1*Pr_t)
            if neqs >= 5; μ_pnode[ip,5] = visc_coeff[5]*ν_d; end
        end
        for ieq = 6:min(neqs,8)
            μ_pnode[ip,ieq] = visc_coeff[ieq]*ν
        end
        if neqs >= 9
            μ_pnode[ip,9] = visc_coeff[9]*ν
        end
    end
    _dsgs_nodal_to_elements_2d!(μ_dsgs, μ_pnode, connijk, nelem, ngl)
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity_nodal!(::DSGS, ::NSD_2D)
#
# The nodal form of the 2D Euler kernel, θ form (ρ, ρu, ρv, ρθ) or total
# energy form (ρ, ρu, ρv, ρE): per-equation denominators with the element
# kernels' floors, ν_i as above with |u| + c at the node, and the element
# kernels' slot assignment with the NODAL density — θ form: 0 on ρ, ρν on
# the momenta, (Pr/(γ−1))ρν on ρθ; energy form: (ν on ∇ρ), ρν on u, v,
# (Pr/(γ−1))ρν on T. (The element energy kernel normalizes the momentum
# residual as a vector; here every component is normalized on its own.)
# ================================================================================
function compute_dsgs_viscosity_nodal!(μ_dsgs::AbstractMatrix{TT},
                                       μ_pnode::AbstractMatrix{TT},
                                       ::DSGS, ::NSD_2D,
                                       q::AbstractMatrix{TT},
                                       q1::AbstractMatrix{TT},
                                       q2::AbstractMatrix{TT},
                                       rhs_el::AbstractArray{TT},
                                       ω::AbstractVector{TT},
                                       Je::AbstractArray{TT},
                                       visc_coeff::AbstractVector{TT},
                                       avg::AbstractVector{TT},
                                       denom::AbstractVector{TT},
                                       qmin::AbstractVector{TT},
                                       qmax::AbstractVector{TT},
                                       nmin::AbstractMatrix{TT},
                                       nmax::AbstractMatrix{TT},
                                       hnod::AbstractVector{TT},
                                       Rnod::AbstractMatrix{TT},
                                       mnod::AbstractVector{TT},
                                       wt::NTuple{3,TT},
                                       connijk::AbstractArray{TI,4},
                                       Δelem::AbstractVector{TT},
                                       PhysConst::PhysicalConst{TT},
                                       Pr::TT, CR::TT, Cmax::TT, Cl::TT,
                                       comm,
                                       nelem::Int, ngl::Int, npoin::Int;
                                       ltheta::Bool=true,
                                       lglobal_norms::Bool=false,
                                       Cmin::TT=zero(TT)) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 4)
    rel  = TT(1.0e-3)
    γ    = PhysConst.γ
    γm1  = PhysConst.γm1
    Cθ   = PhysConst.C0
    eps  = TT(1.0e-16)
    k    = max(ngl - 1, 1)

    _dsgs_nodal_stats_2d!(q, connijk, Δelem, nelem, ngl, npoin, neqs, k,
                          avg, denom, qmin, qmax, nmin, nmax, hnod, lglobal_norms, comm)

    ρ_avg = max(abs(avg[1]), eps)
    if ltheta
        θ_avg = avg[4]/ρ_avg
        p_avg = Cθ*(max(ρ_avg*θ_avg, zero(TT)))^γ
    else
        p_avg = γm1*max(avg[4] - TT(0.5)*(avg[2]*avg[2] + avg[3]*avg[3])/ρ_avg, zero(TT))
    end
    c_avg = sqrt(max(γ*p_avg/ρ_avg, eps))
    @inbounds begin
        denom[1] = max(denom[1], rel*ρ_avg)
        denom[2] = max(denom[2], rel*ρ_avg*c_avg)
        denom[3] = max(denom[3], rel*ρ_avg*c_avg)
        denom[4] = max(denom[4], ltheta ? rel*abs(avg[4]) : rel*ρ_avg*c_avg*c_avg)
    end

    _dsgs_nodal_residual_2d!(Rnod, mnod, q, q1, q2, wt, rhs_el, ω, Je, connijk, nelem, ngl, npoin, NRES)
    @inbounds for ip = 1:npoin
        ratio = _dsgs_nodal_ratio(Rnod, ip, NRES, denom, qmin, qmax, nmin, nmax, Cl, eps)
        ρl = max(q[ip,1], eps)
        ul = q[ip,2]/ρl
        vl = q[ip,3]/ρl
        if ltheta
            θl = q[ip,4]/ρl
            pl = Cθ*(max(ρl*θl, zero(TT)))^γ
            c  = sqrt(max(γ*pl/ρl, zero(TT)))
        else
            Tl = max(q[ip,4]/ρl - TT(0.5)*(ul*ul + vl*vl), zero(TT))
            c  = sqrt(γ*Tl)
        end
        λ  = sqrt(ul*ul + vl*vl) + c
        h  = hnod[ip]
        ν_c = max(zero(TT), min(Cmax*h*λ, CR*h*h*ratio))
        ν   = Cmin > zero(TT) ? max(ν_c, Cmin*h*λ) : ν_c
        μd  = ρl*ν
        if ltheta
            μ_pnode[ip,1] = zero(TT)
            μ_pnode[ip,2] = visc_coeff[2]*μd
            μ_pnode[ip,3] = visc_coeff[3]*μd
            μ_pnode[ip,4] = visc_coeff[4]*(Pr/γm1)*μd
        else
            μ_pnode[ip,1] = visc_coeff[1]*ν
            μ_pnode[ip,2] = visc_coeff[2]*μd
            μ_pnode[ip,3] = visc_coeff[3]*μd
            μ_pnode[ip,4] = visc_coeff[4]*(Pr/γm1)*μd
        end
        for ieq = 5:neqs                       # passive tracers: kinematic ν (as the element form)
            μ_pnode[ip,ieq] = visc_coeff[ieq]*ν
        end
    end
    _dsgs_nodal_to_elements_2d!(μ_dsgs, μ_pnode, connijk, nelem, ngl)
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D)
#
# Marras-Nazarov DynSGS for the 2D non-linear shallow-water system
#
#     q = (H, Hu, Hv),     H the water depth above the bathymetry,
#
# after Marras, Kopera, Constantinescu, Suckale, Giraldo, Adv. Water
# Resour. 114 (2018) 45-63 (residual-based shock capturing for the SWE with
# CG/DG spectral elements), in the form of the Euler and MHD kernels of
# this file:
#
#     ν_res|e = C_R · Δ² · max_i ‖R_i‖∞,e / ‖δq_i − ⟨δq_i⟩‖∞,Ω
#     ν_max|e = C_max · Δ · (‖v‖ + √(gH))∞,e
#     ν|e     = max(ν_floor, max(0, min(ν_max, ν_res))),  ν_floor = C_min Δ (‖v‖ + √(gH))∞,e
#
# with the BDF2 lumped-mass residual R_i of every equation (see the MHD
# header) and the statistics taken on the DEPARTURE δq = q − qe from the
# reference state qe = (He, 0, 0), the lake at rest of the case: the depth
# itself is cone-shaped over an island (‖H − ⟨H⟩‖ is the island, not the
# wave), while the solitary wave is a small perturbation of it. The
# denominators are floored at 10⁻³ of the still-water scales H̄, H̄√(gH̄)
# (H̄ the mean depth) so that a field that is uniform, or at rest,
# contributes 0/floor = 0.
#
# The velocity of the wave speed is Hu/max(H, h_min) with h_min the
# case's wet/dry threshold (:dsgs_swe_hmin), so a thin film does not
# produce |v| = |Hu|/ε; g is :dsgs_swe_g. The three slots receive the SAME
# kinematic ν (times the deck's :μ multipliers): the case's
# user_primitives! hands (H − He, Hu, Hv), so the operator is
# ∇·(ν∇(H − He)) on the continuity equation, which vanishes at the lake at
# rest, and the divergence of ν(∇(Hv) + ∇(Hv)ᵀ − ⅔∇·(Hv) I) on the momenta
# (the stress form of _expansion_visc! on the conserved momenta).
# ================================================================================
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS_SW, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs_el::AbstractArray{TT},
                                 ω::AbstractVector{TT},
                                 Je::AbstractArray{TT},
                                 visc_coeff::AbstractVector{TT},
                                 avg::AbstractVector{TT},
                                 denom::AbstractVector{TT},
                                 wt::NTuple{3,TT},
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 g::TT, hmin::TT, CR::TT, Cmax::TT,
                                 comm,
                                 nelem::Int, ngl::Int;
                                 lglobal_norms::Bool=false,
                                 Cmin::TT=zero(TT)) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 3)
    rel  = TT(1.0e-3)
    eps  = TT(1.0e-16)

    # --- Pass 1: means of δq = q − qe and of the depth --------------------
    @inbounds for ieq = 1:neqs
        avg[ieq] = zero(TT)
    end
    Hsum = zero(TT)
    @inbounds for ie = 1:nelem, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,1]
        for ieq = 1:neqs
            avg[ieq] += q[ip,ieq] - qe[ip,ieq]
        end
        Hsum += q[ip,1]
    end
    inv_npts = one(TT)/max(TT(nelem*ngl*ngl), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(nelem*ngl*ngl), MPI.SUM, comm)
        MPI.Allreduce!(avg, MPI.SUM, comm)
        Hsum      = MPI.Allreduce(Hsum, MPI.SUM, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end
    H_avg = max(Hsum*inv_npts, hmin)

    # --- Pass 2: L∞ spreads of δq, floored at the still-water scales -------
    @inbounds for ieq = 1:neqs
        denom[ieq] = zero(TT)
    end
    @inbounds for ie = 1:nelem, j = 1:ngl, i = 1:ngl
        ip = connijk[ie,i,j,1]
        for ieq = 1:neqs
            denom[ieq] = max(denom[ieq], abs(q[ip,ieq] - qe[ip,ieq] - avg[ieq]))
        end
    end
    if lglobal_norms
        MPI.Allreduce!(denom, MPI.MAX, comm)
    end
    c_avg = sqrt(g*H_avg)
    @inbounds begin
        denom[1] = max(denom[1], rel*H_avg) + eps
        for ieq = 2:NRES
            denom[ieq] = max(denom[ieq], rel*H_avg*c_avg) + eps
        end
    end

    # --- Pass 3: per element ------------------------------------------------
    @inbounds for ie = 1:nelem
        Δ     = Δelem[ie]/ngl
        ratio = zero(TT)
        wmax  = zero(TT)
        for j = 1:ngl, i = 1:ngl
            ip = connijk[ie,i,j,1]
            imK = one(TT)/(ω[i]*ω[j]*Je[ie,i,j])   # element lumped mass at the node
            for ieq = 1:NRES
                R = abs((wt[1]*q[ip,ieq] + wt[2]*q1[ip,ieq] + wt[3]*q2[ip,ieq]) - imK*rhs_el[ie,i,j,ieq])
                ratio = max(ratio, R/denom[ieq])
            end
            Hc = max(q[ip,1], zero(TT))
            Hd = max(q[ip,1], hmin)
            ul = q[ip,2]/Hd
            vl = q[ip,3]/Hd
            wmax = max(wmax, sqrt(ul*ul + vl*vl) + sqrt(g*Hc))
        end
        ν_res = CR*Δ*Δ*ratio
        ν_max = Cmax*Δ*wmax
        ν     = max(zero(TT), min(ν_max, ν_res))
        ν     = Cmin > zero(TT) ? max(ν, Cmin*Δ*wmax) : ν
        for ieq = 1:neqs
            μ_dsgs[ie,ieq] = visc_coeff[ieq]*ν
        end
    end
    return nothing
end

# ================================================================================
# compute_dsgs_viscosity_nodal!(::DSGS_SW, ::NSD_2D)
#
# The nodal (Dao & Nazarov) form of the shallow-water kernel: ν_i at every
# node from the assembled residual, the eq. 4.7 normalization with C_l on
# the departure δq = q − qe, h_i = Δ_K/(k+1), the wave speed |v_i| + √(gH_i),
# the C_min floor; one kinematic ν_i on the three slots. μ_dsgs receives the
# element means for the output staircase.
# ================================================================================
function compute_dsgs_viscosity_nodal!(μ_dsgs::AbstractMatrix{TT},
                                       μ_pnode::AbstractMatrix{TT},
                                       ::DSGS_SW, ::NSD_2D,
                                       q::AbstractMatrix{TT},
                                       q1::AbstractMatrix{TT},
                                       q2::AbstractMatrix{TT},
                                       qe::AbstractMatrix{TT},
                                       rhs_el::AbstractArray{TT},
                                       ω::AbstractVector{TT},
                                       Je::AbstractArray{TT},
                                       visc_coeff::AbstractVector{TT},
                                       avg::AbstractVector{TT},
                                       denom::AbstractVector{TT},
                                       qmin::AbstractVector{TT},
                                       qmax::AbstractVector{TT},
                                       nmin::AbstractMatrix{TT},
                                       nmax::AbstractMatrix{TT},
                                       hnod::AbstractVector{TT},
                                       Rnod::AbstractMatrix{TT},
                                       mnod::AbstractVector{TT},
                                       wt::NTuple{3,TT},
                                       connijk::AbstractArray{TI,4},
                                       Δelem::AbstractVector{TT},
                                       g::TT, hmin::TT, CR::TT, Cmax::TT, Cl::TT,
                                       comm,
                                       nelem::Int, ngl::Int, npoin::Int;
                                       lglobal_norms::Bool=false,
                                       Cmin::TT=zero(TT)) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 3)
    rel  = TT(1.0e-3)
    eps  = TT(1.0e-16)
    k    = max(ngl - 1, 1)

    _dsgs_nodal_stats_2d!(q, connijk, Δelem, nelem, ngl, npoin, neqs, k,
                          avg, denom, qmin, qmax, nmin, nmax, hnod, lglobal_norms, comm; qe=qe)

    # mean depth for the floors of the spreads
    Hsum = zero(TT)
    @inbounds for ip = 1:npoin
        Hsum += q[ip,1]
    end
    inv_npts = one(TT)/max(TT(npoin), one(TT))
    if lglobal_norms
        npts_glob = MPI.Allreduce(TT(npoin), MPI.SUM, comm)
        Hsum      = MPI.Allreduce(Hsum, MPI.SUM, comm)
        inv_npts  = one(TT)/max(npts_glob, one(TT))
    end
    H_avg = max(Hsum*inv_npts, hmin)
    c_avg = sqrt(g*H_avg)
    @inbounds begin
        denom[1] = max(denom[1], rel*H_avg)
        for ieq = 2:NRES
            denom[ieq] = max(denom[ieq], rel*H_avg*c_avg)
        end
    end

    _dsgs_nodal_residual_2d!(Rnod, mnod, q, q1, q2, wt, rhs_el, ω, Je, connijk, nelem, ngl, npoin, NRES)
    @inbounds for ip = 1:npoin
        ratio = _dsgs_nodal_ratio(Rnod, ip, NRES, denom, qmin, qmax, nmin, nmax, Cl, eps)
        Hc = max(q[ip,1], zero(TT))
        Hd = max(q[ip,1], hmin)
        ul = q[ip,2]/Hd
        vl = q[ip,3]/Hd
        λ  = sqrt(ul*ul + vl*vl) + sqrt(g*Hc)
        h  = hnod[ip]
        ν_c = max(zero(TT), min(Cmax*h*λ, CR*h*h*ratio))
        ν   = Cmin > zero(TT) ? max(ν_c, Cmin*h*λ) : ν_c
        for ieq = 1:neqs
            μ_pnode[ip,ieq] = visc_coeff[ieq]*ν
        end
    end
    _dsgs_nodal_to_elements_2d!(μ_dsgs, μ_pnode, connijk, nelem, ngl)
    return nothing
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
        κ_turb = μ_turb / (ρ * Pr_t)
        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
    else                                   # other scalars (moisture, species)
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
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
    Ri_crit = sgs.Ri_crit
    C_s2    = sgs.C_s2

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
            f_Ri_val = if Ri >= Ri_crit
                0.0
            elseif Ri >= 0.0
                ratio = Ri / Ri_crit
                (1.0 - ratio) * (1.0 - ratio)
            else
                min(sqrt(1.0 - 16.0*Ri), 3.0)
            end
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,m,1]
        sgs.μ_turb[ip] = ρ * C_s2 * Δ2 * Sij_val * f_Ri_val
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
    Ri_crit = sgs.Ri_crit
    C_vrem  = sgs.C_vrem
    eps_v   = eps(1.0)

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

        f_Ri_val = 1.0
        if lrichardson
            S11 = dudx;  S22 = dvdy;  S33 = dwdz
            S12 = 0.5*(dudy + dvdx)
            S13 = 0.5*(dudz + dwdx)
            S23 = 0.5*(dvdz + dwdy)
            S_ij_S_ij = S11*S11 + S22*S22 + S33*S33 + 2.0*(S12*S12 + S13*S13 + S23*S23)
            Sij2_val  = 2.0 * S_ij_S_ij
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = if Ri >= Ri_crit
                0.0
            elseif Ri >= 0.0
                ratio = Ri / Ri_crit
                (1.0 - ratio) * (1.0 - ratio)
            else
                min(sqrt(1.0 - 16.0*Ri), 3.0)
            end
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,m,1]
        μ_base = (u_ij_u_ij > eps_v && B_β > 0.0) ?
                 ρ * C_vrem * sqrt(B_β / u_ij_u_ij) : 0.0
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
    Ri_crit = sgs.Ri_crit
    C_s2    = sgs.C_s2

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
            f_Ri_val = if Ri >= Ri_crit
                0.0
            elseif Ri >= 0.0
                ratio = Ri / Ri_crit
                (1.0 - ratio) * (1.0 - ratio)
            else
                min(sqrt(1.0 - 16.0*Ri), 3.0)
            end
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,1]
        sgs.μ_turb[ip] = ρ * C_s2 * Δ2 * Sij_val * f_Ri_val
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
    Ri_crit = sgs.Ri_crit
    C_vrem  = sgs.C_vrem
    eps_v   = eps(1.0)

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

        f_Ri_val = 1.0
        if lrichardson
            Sij2_val = 2.0*(dudx*dudx + dvdy*dvdy + 2.0*(0.5*(dudy + dvdx))^2)
            Ri = Sij2_val > 1e-12 ? N2_val / Sij2_val : 0.0
            f_Ri_val = if Ri >= Ri_crit
                0.0
            elseif Ri >= 0.0
                ratio = Ri / Ri_crit
                (1.0 - ratio) * (1.0 - ratio)
            else
                min(sqrt(1.0 - 16.0*Ri), 3.0)
            end
        end
        sgs.f_Ri[ip] = f_Ri_val

        ρ = uprimitive[k,l,1]
        μ_base = (u_ij_u_ij > eps_v && B_β > 0.0) ?
                 ρ * C_vrem * sqrt(B_β / u_ij_u_ij) : 0.0
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
        κ_turb = μ_turb / (ρ * Pr_t)
        if ltheta_eqn
            return κ_turb * visc_coeffieq[ieq]
        else
            return (κ_mol + κ_turb) * visc_coeffieq[ieq]
        end
    else                                      # other scalars (moisture, species)
        κ_turb_scalar = μ_turb / (ρ * Sc_t)
        return (κ_mol + κ_turb_scalar) * visc_coeffieq[ieq]
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
                                          PhysConst, Δ2, ::SMAG)
    C_s2 = PhysConst.C_s * PhysConst.C_s
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
                                          PhysConst, Δ2, ::VREM)
    C_s2   = PhysConst.C_s * PhysConst.C_s
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
                                          PhysConst, Δ2, ::Any)
    return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

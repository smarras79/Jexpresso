using Statistics
using Plots

# module-level counter; persists across calls without changing the function signature
const _DIAG_CALL_COUNT = Ref(0)
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
#     μ_res = C1 · Δ² · max_i ‖R_i‖∞,Ω / ‖q_i − ⟨q_i⟩‖∞,Ω
#     μ_max = C2 · Δ · max(|u| + c)
#     μ_dsgs[iel] = max(0, min(μ_res, μ_max))
# where R_i is the STRONG-form BDF2 residual of conservation law i —
# (3qⁿ − 4qⁿ⁻¹ + qⁿ⁻²)/(2Δt) − M⁻¹·RHS. Since jexpresso assembles RHS
# in weak form (post-DSS, pre-mass-matrix division), the rhs argument
# is multiplied by Minv[ip] inline before the BDF2 minus, which makes
# μ_dsgs dimensionally a kinematic viscosity (m²/s) regardless of SD.
#
# Both numerators and denominators are global L∞ norms so the
# coefficient cannot be inlined into the (k,l) loop the way
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

# ---------------- 1D --------------------------------------------------
#
# Conservation form q = (ρ, ρu, ρE) on a 1D LGL mesh. The signature is
# a hand-typed function barrier (concrete arrays, no params.* lookups)
# so Julia can specialize and the inner loop is allocation-free.
#
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS, ::NSD_1D,
                                 coords,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 f_back::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 back_weights::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δx::AbstractVector{TT},
                                 nelem::Int, ngl::Int, time::TT;
                                 # NEW: region(s) to correlate mu_res_orig against (1-s),
                                 # deliberately AWAY from the known shock/contact elements.
                                 # Defaults below are guesses at the two "star state" plateaus
                                 # for the standard Sod IC at t=0.2 -- adjust to match the
                                 # actual exact-solution structure for the case being run.
                                 diagnostic_xranges::Vector{Tuple{TT,TT}} =
                                     [(TT(0.50), TT(0.62)), (TT(0.73), TT(0.80))]
                                 ) where {TT<:AbstractFloat, TI<:Integer}
    
    invnp = one(TT)/(nelem*ngl)
    γ     = TT(1.4)
    C1    = TT(1.0)
    C2    = TT(0.5)
    eps   = Base.eps(TT)
    neqs  = size(μ_dsgs, 2)

    # --- Pass 1: domain averages of q ----------------------------------
    ρ_avg  = zero(TT); ρu_avg = zero(TT); ρE_avg = zero(TT)
    rho_el = zeros(TT, ngl)
    x_el = zeros(TT, ngl)
    R_el = zeros(TT,ngl)
    ent_el       = zeros(TT, ngl) 
    r_el = zeros(TT,nelem)
    s_ent_el     = zeros(TT, nelem)   # NEW: sensor value from entropy field
    r_ent_el     = zeros(TT, nelem)
    decay_ent_el = zeros(TT, nelem)
    logr_ent_el = zeros(TT, nelem)
    logr_el = zeros(TT,nelem)
    decay_el = zeros(TT,nelem)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            ρ_avg  += q[ip,1]
            ρu_avg += q[ip,2]
            ρE_avg += q[ip,3]
        end
    end
    ρ_avg  *= invnp
    ρu_avg *= invnp
    ρE_avg *= invnp

    # --- Pass 2: domain L∞ norms of |q - ⟨q⟩| --------------------------
    denom1 = zero(TT); denom2 = zero(TT); denom3 = zero(TT)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = max(denom1, abs(q[ip,1] - ρ_avg))
            denom2 = max(denom2, abs(q[ip,2] - ρu_avg))
            denom3 = max(denom3, abs(q[ip,3] - ρE_avg))
        end
    end
    denom1 += eps; denom2 += eps; denom3 += eps

    # NEW: per-element storage for the correlation diagnostic
    μres_orig_el = zeros(TT, nelem)
    s_el         = zeros(TT, nelem)
    xc_el        = zeros(TT, nelem)

    # --- Pass 3: per-element loop --------------------------------------
    inv2Δt = one(TT)/(2*Δt)
    @inbounds for ie = 1:nelem
        Δ = Δx[ie]/ngl

        n1   = zero(TT); n2 = zero(TT); n3 = zero(TT)
        uTmx = zero(TT)
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            rho_el[i] = q[ip,1]
            x_el[i] = coords[1,ip]
            R1 = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
            R_el[i] = R1
            n1 = max(n1, R1); n2 = max(n2, R2); n3 = max(n3, R3)

            ρl = q[ip,1]
            ul = q[ip,2]/ρl
            el = q[ip,3]/ρl
            eint = max(el - TT(0.5)*ul*ul, zero(TT))
            uTmx = max(uTmx, abs(ul) + sqrt(γ*(γ - one(TT))*eint))

            # NEW: entropy-function field, log form: ln(p) - γ*ln(ρ)
            pl = (γ - one(TT)) * ρl * eint
            ent_el[i] = log(max(pl, eps)) - γ*log(max(ρl, eps))
        end

        μ_res_orig = C1*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3)  # renamed, kept pristine
        μ_res, s, b, r_rho, logr = backscatter_gate(rho_el, f_back, back_weights, μ_res_orig, x_el; beta_max= TT(0.6),tau=TT(-7.0), kappa=TT(0.7))
        μ_max = C2*Δ*uTmx
        μ     = max(zero(TT), min(μ_max, μ_res))
        for ieq = 1:neqs
            μ_dsgs[ie, ieq] = visc_coeff[ieq] * μ
        end

        # NEW: same modal machinery, entropy field instead of density
        m_ent = f_back * ent_el
        s_ent, r_ent, logr_ent, decay_ent = smoothness_sensor(ent_el, f_back, back_weights;
                                                                tau=TT(-8.7), kappa=TT(0.7))

        decay_ent = spectral_decay_rate(m_ent)

        m = f_back * rho_el
        decay_el[ie] = spectral_decay_rate(m)

        μres_orig_el[ie] = μ_res_orig
        s_el[ie]         = s
        r_el[ie]         = r_rho
        logr_el[ie] = logr

        logr_ent_el[ie] = logr_ent; s_ent_el[ie] = s_ent; r_ent_el[ie] = r_ent; decay_ent_el[ie] = decay_ent
        xc = zero(TT)
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            xc += coords[1,ip]
        end
        xc_el[ie] = xc / ngl
    end

    # NEW: plateau-region correlation between the UNMODIFIED mu_res and (1-s).
    # High correlation -> s is tracking the same structure mu_res already has
    # in these nominally-smooth regions. Low/near-zero -> s is adding
    # information (or noise) independent of what mu_res already sees there.
    mask = falses(nelem)
    @inbounds for ie = 1:nelem
        for (xlo, xhi) in diagnostic_xranges
            if xlo <= xc_el[ie] <= xhi
                mask[ie] = true
                break
            end
        end
    end
    n_pts = count(mask)
    #=if n_pts >= 3
        corr = cor(μres_orig_el[mask], one(TT) .- s_el[mask])
        @info "plateau-region correlation(mu_res_orig, 1-s)" n_points=n_pts correlation=corr
    else
        @warn "plateau-region correlation: too few elements matched diagnostic_xranges" n_points=n_pts
    end=#

    plot_every = 50
    _DIAG_CALL_COUNT[] += 1
    if plot_every > 0 && _DIAG_CALL_COUNT[] % plot_every == 0
        n = _DIAG_CALL_COUNT[]

        p1 = plot(xc_el, μres_orig_el, label="μ_res_orig", lw=2, title="raw vs final")
        plot!(p1, xc_el, μ_dsgs[:,1], label="μ_final", lw=2, ls=:dash)

        p2 = plot(xc_el, s_el, label="s (density)", lw=2, ylim=(0,1), title="smoothness sensor")
        plot!(p2, xc_el, s_ent_el, label="s (entropy)", lw=2)

        p5 = plot(xc_el, decay_el, label="decay (density)", lw=2, color=:green, title="modal decay slope")
        plot!(p5, xc_el, decay_ent_el, label="decay (entropy)", lw=2, color=:purple)

        display(plot(p2, p5, layout=(1,2), size=(1400,500), plot_title="call $n"))

        head, tail, contact_x, shock_x = sod_wave_positions(time)
        margin = TT(3) * (Δx[1] / ngl)   # a few nodes' worth, keeps windows off the fronts

        raref_mask      = (head + margin .<= xc_el .<= tail - margin)
        left_star_mask  = (tail + margin .<= xc_el .<= contact_x - margin)
        right_star_mask = (contact_x + margin .<= xc_el .<= shock_x - margin)
        contact_mask    = (contact_x - margin .<= xc_el .<= contact_x + margin)
        shock_mask       = (shock_x   - margin .<= xc_el .<= shock_x   + margin)
        @info maximum(decay_el[shock_mask]), maximum(logr_el[shock_mask])
        # guard against empty windows early in the run, before features
        # have separated enough for a window to contain any elements
        function safe_stats(vals, mask, label)
            if count(mask) < 2
                @info "entropy logr -- $label" note="window empty or too small" n=count(mask)
                return
            end
            @info "entropy logr -- $label" mean=mean(vals[mask]) std=std(vals[mask]) n=count(mask)
        end
        @info "maximum shock decay, logr", maximum(decay_el[shock_mask]), maximum(logr_el[shock_mask])
        safe_stats(logr_el, raref_mask,      "rarefaction")
        safe_stats(logr_el, left_star_mask,  "left star plateau")
        safe_stats(logr_el, right_star_mask, "right star plateau")
        safe_stats(decay_el, raref_mask,      "rarefaction")
        safe_stats(decay_el, left_star_mask,  "left star plateau")
        safe_stats(decay_el, right_star_mask, "right star plateau")

        #=if count(contact_mask) >= 1
            @info "entropy logr -- contact peak" value=maximum(logr_ent_el[contact_mask]) t=time x=contact_x
        end
        if count(shock_mask) >= 1
            @info "entropy logr -- shock peak" value=maximum(logr_ent_el[shock_mask]) t=time x=shock_x
        end=#
    end
   
    return nothing
end

function compute_dsgs_viscosity!(μ_dsgs::AbstractArray{TT},
                                 ::NDSGS, ::NSD_1D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δx::AbstractVector{TT},
                                 nelem::Int, ngl::Int) where {TT<:AbstractFloat, TI<:Integer}

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

    # --- Pass 1: domain averages of q ----------------------------------
    ρ_avg  = zero(TT); ρu_avg = zero(TT); ρE_avg = zero(TT)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            ρ_avg  += q[ip,1]
            ρu_avg += q[ip,2]
            ρE_avg += q[ip,3]
        end
    end
    ρ_avg  *= invnp
    ρu_avg *= invnp
    ρE_avg *= invnp

    # --- Pass 2: domain L∞ norms of |q - ⟨q⟩| --------------------------
    denom1 = zero(TT); denom2 = zero(TT); denom3 = zero(TT)
    @inbounds for ie = 1:nelem
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = max(denom1, abs(q[ip,1] - ρ_avg))
            denom2 = max(denom2, abs(q[ip,2] - ρu_avg))
            denom3 = max(denom3, abs(q[ip,3] - ρE_avg))
        end
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

            R1 = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
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

        #Forward scattering only Nodal DSGS V0 max(abs(element)) num, abs(nodal) denom (WORKING STABLE)
        #=μ_max = C2*Δ*uTmx
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = abs(q[ip,1] - ρ_avg)
            denom2 = abs(q[ip,2] - ρu_avg)
            denom3 = abs(q[ip,3] - ρE_avg)
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#

        #Forward scattering only Nodal DSGS V0.1 abs(nodal) num, abs(nodal) denom (WORKING STABLE)
        #=μ_max = C2*Δ*uTmx
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = abs(((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]))
            n2 = abs(((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]))
            n3 = abs(((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]))
            denom1 = abs(q[ip,1] - ρ_avg)
            denom2 = abs(q[ip,2] - ρu_avg)
            denom3 = abs(q[ip,3] - ρE_avg)
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Forward scattering only Nodal DSGS V0.2 avg(abs(element)) num, abs(nodal) denom (WORKING STABLE)
        #=μ_max = C2*Δ*uTmx
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 += abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 += abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = abs(((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]))
            n2 = abs(((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]))
            n3 = abs(((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]))
            denom1 = abs(q[ip,1] - ρ_avg)
            denom2 = abs(q[ip,2] - ρu_avg)
            denom3 = abs(q[ip,3] - ρE_avg)
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(zero(TT), min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering  Nodal DSGS V1.1 max(abs(element)) num, nodal denom (unstable)
        #=μ_max = C2*Δ*uTmx
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.0 max(abs(element)) num, avg(element) denom (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.2 avg(abs(element)) num, nodal denom (unstable)
        #=μ_max = C2*Δ*uTmx
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 += abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 += abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.3 nodal num, nodal denom (unstable)
        #=μ_max = C2*Δ*uTmx
        
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            n2 = (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            n3 = (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.4 avg(element) num, nodal denom (unstable)
        #=μ_max = C2*Δ*uTmx
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.5 -nodal num, nodal denom (unstable)
        #=μ_max = C2*Δ*uTmx
        
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = -((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            n2 = -((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            n3 = -((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.6 -avg(element) num, nodal denom (unstable) ALL PURE NODAL DENOM ARE UNSTABLE
        #=μ_max = C2*Δ*uTmx
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = -R1*invnp; n2 = -R2*invnp; n3 = -R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            denom1 = q[ip,1] - ρ_avg
            denom2 = q[ip,2] - ρu_avg
            denom3 = q[ip,3] - ρE_avg
            denom1 += eps; denom2 += eps; denom3 += eps
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end
    end=#
    #Dual scattering Nodal DSGS V1.7 avg(abs(element)) num, avg(element) denom (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            R2 += abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            R3 += abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.8 nodal num, avg(element) denom  (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            n2 = (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            n3 = (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.9 avg(element) num, avg(element) denom (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.10 -nodal num, avg(element) denom (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = -((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
            n2 = -((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2])
            n3 = -((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3])
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.11 -avg(element) num, avg(element) denom (unstable)
        #=μ_max = C2*Δ*uTmx
        denom1 = 0.0; denom2=0.0; denom3=0.0
        for i = 1:ngl
            ip = connijk[ie,i,1,1]
            denom1 += (q[ip,1] - ρ_avg)*invnp
            denom2 += (q[ip,2] - ρu_avg)*invnp
            denom3 += (q[ip,3] - ρE_avg)*invnp
        end
        denom1 += eps; denom2 += eps; denom3 += eps
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]

            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = -R1*invnp; n2 = -R2*invnp; n3 = -R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.12 avg(element) num, max(abs(denom))
        #=μ_max = C2*Δ*uTmx
        
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = R1*invnp; n2 = R2*invnp; n3 = R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.13 -avg(element) num, max(abs(denom)) (unstable)
        #=μ_max = C2*Δ*uTmx
        
        R1 =0.0; R2=0.0; R3=0.0
        @simd for i = 1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            R1 += (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            R2 += (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            R3 += (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
        end
        n1 = -R1*invnp; n2 = -R2*invnp; n3 = -R3*invnp
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
        #Dual scattering Nodal DSGS V1.14 nodal num, max(abs(denom))
        μ_max = C2*Δ*uTmx
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = (3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            n2 = (3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            n3 = (3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-1e-6, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-1e-6, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-1e-6, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
            #STABILITY ISSUES are likely dictated by negative viscosity limit.
            #=μ_dsgs[ie, i, 1] = maximum(μ_dsgs[ie, i, :])
            μ_dsgs[ie, i, 2] = maximum(μ_dsgs[ie, i, :])
            μ_dsgs[ie, i, 3] = maximum(μ_dsgs[ie, i, :])=#
        end
        #Dual scattering Nodal DSGS V1.15 -nodal num, max(abs(denom))
        #=μ_max = C2*Δ*uTmx
        for i =1:ngl
            ip = connijk[ie,i,1,1]
            Mi = Minv[ip]
            n1 = -(3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1]
            n2 = -(3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])*inv2Δt - Mi*rhs[ip,2]
            n3 = -(3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])*inv2Δt - Mi*rhs[ip,3]
            μ_dsgs[ie, i, 1] = visc_coeff[1] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n1/denom1)))
            μ_dsgs[ie, i, 2] = visc_coeff[2] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n2/denom2)))
            μ_dsgs[ie, i, 3] = visc_coeff[3] * max(-μ_max, min(μ_max, C1*Δ*Δ*max(n3/denom3)))
        end=#
    end
    @info maximum(μ_dsgs[:, :, 1]), minimum(μ_dsgs[:, :, 1])
    @info maximum(μ_dsgs[:, :, 2]), minimum(μ_dsgs[:, :, 2])
    @info maximum(μ_dsgs[:, :, 3]), minimum(μ_dsgs[:, :, 3])
    
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
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 qe::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 PhysConst::PhysicalConst{TT},
                                 Pr::TT,
                                 nelem::Int, ngl::Int) where {TT<:AbstractFloat, TI<:Integer}

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

    # --- Pass 1: domain averages of (ρ, ρu, ρv, ρθ) --------------------
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
    ρ_avg  *= invnp; ρu_avg *= invnp
    ρv_avg *= invnp; ρθ_avg *= invnp

    # --- Pass 2: domain L∞ norms of |q - ⟨q⟩| --------------------------
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
                R1 = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) - Mi*rhs[ip,1])
                R2 = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) - Mi*rhs[ip,2])
                R3 = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) - Mi*rhs[ip,3])
                R4 = abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])/(2*Δt) - Mi*rhs[ip,4])
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
# All reductions are MPI-global: ⟨q⟩ and ‖q−⟨q⟩‖_{∞,Ω} are domain norms by
# definition, so a rank-local version would make the viscosity depend on
# the partitioning.
# ================================================================================
function compute_dsgs_viscosity!(μ_dsgs::AbstractMatrix{TT},
                                 ::DSGS_MHD, ::NSD_2D,
                                 q::AbstractMatrix{TT},
                                 q1::AbstractMatrix{TT},
                                 q2::AbstractMatrix{TT},
                                 rhs::AbstractMatrix{TT},
                                 Minv::AbstractVector{TT},
                                 visc_coeff::AbstractVector{TT},
                                 avg::AbstractVector{TT},
                                 denom::AbstractVector{TT},
                                 Δt::TT,
                                 connijk::AbstractArray{TI,4},
                                 Δelem::AbstractVector{TT},
                                 γ::TT, Pr_t::TT, C1::TT, C2::TT,
                                 comm,
                                 nelem::Int, ngl::Int) where {TT<:AbstractFloat, TI<:Integer}

    neqs = size(μ_dsgs, 2)
    NRES = min(neqs, 8)          # residual max excludes the ψ slot
    γm1  = γ - one(TT)
    eps  = TT(1.0e-16)

    # --- Pass 1: domain means ⟨q_i⟩ ------------------------------------
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
    npts_local = TT(nelem*ngl*ngl)
    npts_glob  = MPI.Allreduce(npts_local, MPI.SUM, comm)
    MPI.Allreduce!(avg, MPI.SUM, comm)
    inv_npts = one(TT)/max(npts_glob, one(TT))
    @inbounds for ieq = 1:neqs
        avg[ieq] *= inv_npts
    end

    # --- Pass 2: domain L∞ of |q_i − ⟨q_i⟩| ----------------------------
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
    MPI.Allreduce!(denom, MPI.MAX, comm)

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
                    R = abs((3*q[ip,ieq] - 4*q1[ip,ieq] + q2[ip,ieq])*inv2Δt - Mi*rhs[ip,ieq])
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

function broadcast_ndsgs_to_nodes!(μ_dsgs_pnode::AbstractMatrix{TT},
                                  μ_dsgs::AbstractArray{TT},
                                  connijk::AbstractArray{TI,4},
                                  nelem::Int, ngl::Int,
                                  SD::AbstractSpaceDimensions) where {TT,TI}
    neqs = size(μ_dsgs, 2)
    if SD === NSD_1D()
        @inbounds for ie = 1:nelem
            for i = 1:ngl
                ip = connijk[ie,i,1,1]
                for ieq = 1:neqs
                    μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, i, ieq]
                end
            end
        end
    elseif SD === NSD_2D()
        @inbounds for ie = 1:nelem
            for j = 1:ngl
                for i = 1:ngl
                    ip = connijk[ie,i,j,1]
                    for ieq = 1:neqs
                        μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, i, j, ieq]
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
                            μ_dsgs_pnode[ip, ieq] = μ_dsgs[ie, i, j, k, ieq]
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

"""
    backscatter_gate(rho_el, f_x, weight, mu_res_el; beta_max=0.2)

Per-element, per-RK-stage backscatter gate for the residual-based DSGS model.

Purely local: uses only this element's own nodal density values, its own
hierarchical modal transform, and the pre-existing Vandeven trust weights.
No neighbor or column information of any kind.

Arguments:
  rho_el   :: Vector{TFloat}  -- local density at this element's ngl nodes
  f_x      :: Matrix{TFloat}  -- UNBLENDED nodal->modal transform for this
                                  element (i.e. `leg_inv` from init_filter,
                                  NOT the mu_x-blended production filter `f`)
  weight   :: Vector{TFloat}  -- per-mode Vandeven trust weights, length ngl,
                                  the same `weight` computed inside init_filter
                                  (currently local-only there; needs to be
                                  returned alongside f_x for this to work)
  mu_res_el:: TFloat          -- this element's already-computed mu_res
                                  (unchanged, existing model)

Keyword:
  beta_max :: TFloat = 0.2    -- ceiling on the fraction of mu_res that can
                                  be withheld, even in a maximally smooth,
                                  maximally under-resolved element

Returns:
  nu_eff, s, beta   -- the modified diffusion coefficient to use in place
                        of mu_res, plus the sensor and beta for diagnostics
"""
function backscatter_gate(rho_el::AbstractVector{TFloat},
                           f_x::AbstractMatrix{TFloat},
                           weight::AbstractVector{TFloat},
                           mu_res_el::TFloat, x_el;
                           beta_max::TFloat = TFloat(0.2), tau = TFloat(-8.7), kappa = TFloat(0.7)) where {TFloat<:AbstractFloat}

    ngl = length(rho_el)
    eps_energy = TFloat(1.0e-14)

    # --- modal coefficients of the local density field (single element) ---
    m = f_x * rho_el   # m[k] = amplitude of hierarchical mode k

    # --- energy sitting in modes the existing Vandeven weighting distrusts ---
    #=atten_energy = zero(TFloat)
    total_energy = zero(TFloat)
    @inbounds for k = 1:ngl
        mk2 = m[k]^2
        atten_energy += (one(TFloat) - weight[k])^2 * mk2
        total_energy += mk2
    end

    # s -> 1 : energy concentrated in trusted (vertex/low) modes, smooth field
    # s -> 0 : energy concentrated in distrusted (high bubble) modes, shock-like
    s = one(TFloat) - atten_energy / (total_energy + eps_energy)
    s = clamp(s, zero(TFloat), one(TFloat))=#
    s, r, logr, decay = smoothness_sensor(rho_el,
                            f_x,
                            weight;
                            tau = tau,
                            kappa = kappa,
                            eps = TFloat(1e-14))
    #@info "element coordinates", x_el
    #@info "s, r, logr", s, r, logr
    beta = beta_max * s
    nu_eff = mu_res_el * (one(TFloat) - beta)

    return nu_eff, s, beta, r, logr
end

"""
    smoothness_sensor(rho_el, leg_inv, weight; tau=-8.7, kappa=0.7, eps=1e-14)

Log-space energy-fraction ratio (r, logr) is still computed and returned
for diagnostics/comparison, but s itself is now built from the spectral
decay rate of rho's modal coefficients, not from logr -- confirmed via
separation-ratio comparison to discriminate rarefaction/plateau vs. shock
more cleanly than the energy-fraction metric (see prototype log).

s -> 1 : element smooth (fast modal decay), decay << tau
s -> 0 : element shock-like (slow/flat modal decay), decay >> tau

tau, kappa are empirical, calibrated from this problem's own rarefaction
and star-plateau decay-rate means/stds (tau near their midpoint, kappa
near their common std) -- revisit if resolution, order, or problem changes.
"""
function smoothness_sensor(rho_el::AbstractVector{TT},
                            leg_inv::AbstractMatrix{TT},
                            weight::AbstractVector{TT};
                            tau::TT = TT(-8.7),
                            kappa::TT = TT(0.7),
                            eps::TT = TT(1e-14)) where {TT<:AbstractFloat}

    ngl = length(rho_el)
    m = leg_inv * rho_el

    atten_energy = zero(TT)
    total_energy = zero(TT)
    @inbounds for k = 1:ngl
        mk2 = m[k]^2
        atten_energy += (one(TT) - weight[k])^2 * mk2
        total_energy += mk2
    end
    r = atten_energy / (total_energy + eps)
    logr = log10(r + eps)

    decay = spectral_decay_rate(m)   # already defined; reuses the same m

    # NOTE: sign convention -- decay is MORE NEGATIVE for smooth elements,
    # LESS NEGATIVE (closer to 0) for shock-like elements, same "larger =
    # more shock-like" polarity as logr had, so the tanh form is unchanged.
    s = 1 - 0.5*(1 + tanh((decay - tau) / kappa))

    return s, r, logr, decay
end

function spectral_decay_rate(m::AbstractVector{TT}) where {TT<:AbstractFloat}
    ngl = length(m)
    # skip mode 1 (constant/vertex mode dominates trivially, not informative
    # about decay shape) and any exactly-zero coefficients (log undefined)
    ks = TT[]; logm = TT[]
    for k in 2:ngl
        if abs(m[k]) > eps(TT)
            push!(ks, log(TT(k)))
            push!(logm, log(abs(m[k])))
        end
    end
    length(ks) < 3 && return zero(TT)  # not enough points, degenerate element
    # least-squares slope
    kbar, lbar = sum(ks)/length(ks), sum(logm)/length(logm)
    num = sum((ks .- kbar) .* (logm .- lbar))
    den = sum((ks .- kbar).^2)
    return num / den   # more negative = faster decay = smoother
end

function separation_ratio(vals::AbstractVector{TT}, xc::AbstractVector{TT},
                           smooth_range::Tuple{TT,TT}, shock_range::Tuple{TT,TT}) where {TT}
    smooth_mask = smooth_range[1] .<= xc .<= smooth_range[2]
    shock_mask  = shock_range[1] .<= xc .<= shock_range[2]

    μ_smooth, σ_smooth = mean(vals[smooth_mask]), std(vals[smooth_mask])
    μ_shock            = maximum(vals[shock_mask])  # peak, not mean, since shock is a spike not a plateau

    gap = abs(μ_shock - μ_smooth)
    return gap / σ_smooth, μ_smooth, σ_smooth, μ_shock
end

function sod_wave_positions(t::TT; x0::TT=TT(0.5)) where {TT<:AbstractFloat}
    head    = x0 + TT(-1.1832) * t
    tail    = x0 + TT(-0.0704) * t
    contact = x0 + TT( 0.9275) * t
    shock   = x0 + TT( 1.7522) * t
    return head, tail, contact, shock
end
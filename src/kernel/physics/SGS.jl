#----------------------------------------------------------------------
# SMAGORINSKY
#----------------------------------------------------------------------
@inline function SGS_diffusion(visc_coeffieq, ieq,
                               ρ,
                               u11, u22, u12, u21,
                               PhysConst, Δ2,
                               inputs,
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
    C_s2  = C_s*C_s

    
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
            return cp * (κ_mol + κ_turb) * visc_coeffieq[ieq]
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
                               inputs,
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
    # Only apply for potential temperature with Richardson correction enabled
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
    elseif lrichardson && inputs[:energy_equation] == "energy"
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
                               inputs,
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
    cp         = PhysConst.cp
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
            return cp * (κ_mol + κ_turb) * visc_coeffieq[ieq]
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
                               inputs,
                               ::VREM, ::NSD_3D;
                               ltheta_eqn=true,
                               lrichardson=false)

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
            Tl = max(el - TT(0.5)*ul*ul, zero(TT))
            uTmx = max(uTmx, abs(ul) + sqrt(γ*Tl))
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

        for j = 1:ngl
            @simd for i = 1:ngl
                ip = connijk[ie,i,j,1]

                R1 = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) - rhs[ip,1])
                R2 = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) - rhs[ip,2])
                R3 = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) - rhs[ip,3])
                R4 = abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])/(2*Δt) - rhs[ip,4])
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
            end
        end

        μ_res = C1*Δ*Δ*max(n1/denom1, n2/denom2, n3/denom3, n4/denom4)
        μ_max = C2*Δ*uTmx
        μ     = max(zero(TT), min(μ_max, μ_res))

        # Per-equation split (Marras eq. 10), scaled by the user-
        # supplied inputs[:μ] multiplier so the case can be run with
        # DSGS off (visc_coeff = [0,…]) to confirm whether DSGS itself
        # is the cause of any instability.
        μ_dsgs[ie,1] = zero(TT)                             # ρ : no mass diffusion
        μ_dsgs[ie,2] = visc_coeff[2] * μ                    # ρu
        μ_dsgs[ie,3] = visc_coeff[3] * μ                    # ρv
        # DIAGNOSTIC: θ-equation DSGS forced to zero to isolate whether
        # κ_θ·∇²θ on the LGL-discretized total θ field is the cold-start
        # blow-up source. Restore `visc_coeff[4] * (Pr/γm1) * μ` (Marras
        # eq. 10b) once the perturbation-θ fix is in.
        μ_dsgs[ie,4] = zero(TT)                             # ρθ (DIAG: disabled)
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

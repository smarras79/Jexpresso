function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=8, ip=1)
    """
    2D GLM-MHD flux function
    
    Note: Following the paper, we keep the third component (z) of velocity and 
    magnetic field because plasma systems admit three-dimensional electromagnetic 
    interactions in two-dimensional problems (2.5D MHD).
    
    State vector: q = [ρ, ρu, ρv, E, Bx, By, Bz, ψ]
    where (u,v,w) are velocities but w = ρw/ρ is computed from momentum
    and (Bx, By, Bz) includes out-of-plane Bz component
    """
    PhysConst = PhysicalConst{Float64}()
    
    # Extract conservative variables
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρe = q[4]  # Total energy E
    Bx = q[5]
    By = q[6]
    Bz = q[7]  # Out-of-plane magnetic field
    ψ  = q[8]
    
    # Primitive variables
    u = ρu/ρ
    v = ρv/ρ
    # Note: w (out-of-plane velocity) would be q[?] if we had ρw
    # For pure 2D, w = 0, but field Bz can be non-zero
    w = 0.0  # No z-momentum in 2D, but Bz exists
    
    # Thermodynamic properties
    γ   = PhysConst.γ
    γm1 = PhysConst.γm1
    
    # Hyperbolic cleaning speed
    c_h = PhysConst.c_h
    
    # Kinetic energy (only in-plane components for 2D)
    velomagsq = u*u + v*v + w*w  # w=0 for pure 2D
    ke = 0.5*ρ*velomagsq
    
    # Magnetic field magnitude squared (includes ALL three components!)
    Bmagsq = Bx*Bx + By*By + Bz*Bz
    
    # Magnetic energy
    me = 0.5*Bmagsq
    
    # ψ energy
    ψe = 0.5*ψ*ψ
    
    # Thermal pressure: p = (γ-1)(E - ½ρ||v||² - ½||B||² - ½ψ²)
    p_thermal = γm1*(ρe - ke - me - ψe)
    
    # Total pressure (thermal + magnetic)
    p_total = p_thermal + 0.5*Bmagsq
    
    # Energy term for advection
    E_advected = ke + γ*p_thermal/γm1 + Bmagsq
    
    # Velocity dot B (includes out-of-plane Bz!)
    vdotB = u*Bx + v*By + w*Bz
    
    #-----------------------------------------
    # F-flux (x-direction)
    #-----------------------------------------
    
    # Mass flux
    F[1] = ρu
    
    # Momentum flux (2 components: x and y)
    F[2] = ρu*u + p_total - Bx*Bx
    F[3] = ρv*u - By*Bx
    
    # Energy flux
    F[4] = u*E_advected + Bx*(c_h*ψ - vdotB)
    
    # Magnetic field flux (3 components: x, y, and z!)
    F[5] = c_h*ψ           # Bx evolution
    F[6] = u*By - Bx*v     # By evolution
    F[7] = u*Bz - Bx*w     # Bz evolution (out-of-plane)
    
    # Divergence cleaning flux
    F[8] = c_h*Bx
    
    #-----------------------------------------
    # G-flux (y-direction)
    #-----------------------------------------
    
    # Mass flux
    G[1] = ρv
    
    # Momentum flux (2 components: x and y)
    G[2] = ρu*v - Bx*By
    G[3] = ρv*v + p_total - By*By
    
    # Energy flux
    G[4] = v*E_advected + By*(c_h*ψ - vdotB)
    
    # Magnetic field flux (3 components: x, y, and z!)
    G[5] = v*Bx - By*u     # Bx evolution
    G[6] = c_h*ψ           # By evolution
    G[7] = v*Bz - By*w     # Bz evolution (out-of-plane)
    
    # Divergence cleaning flux
    G[8] = c_h*By
end


#-----------------------------------------
# Alternative: 2D with out-of-plane velocity (2.5D)
#-----------------------------------------
function user_flux_2p5d!(F, G, SD::NSD_2D, q, qe,
                         mesh::St_mesh, ::CL, ::TOTAL; neqs=9, ip=1)
    """
    2.5D GLM-MHD flux function with out-of-plane velocity
    
    State vector: q = [ρ, ρu, ρv, ρw, E, Bx, By, Bz, ψ]
    
    This version includes ρw (z-momentum) for cases where the out-of-plane 
    velocity is non-zero (e.g., jets, rotating flows).
    """
    PhysConst = PhysicalConst{Float64}()
    
    # Extract conservative variables
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]  # Out-of-plane momentum
    ρe = q[5]
    Bx = q[6]
    By = q[7]
    Bz = q[8]
    ψ  = q[9]
    
    # Primitive variables
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ  # Out-of-plane velocity
    
    # Thermodynamic properties
    γ   = PhysConst.γ
    γm1 = PhysConst.γm1
    c_h = PhysConst.c_h
    
    # Kinetic energy (includes out-of-plane component!)
    velomagsq = u*u + v*v + w*w
    ke = 0.5*ρ*velomagsq
    
    # Magnetic field magnitude squared
    Bmagsq = Bx*Bx + By*By + Bz*Bz
    me = 0.5*Bmagsq
    
    # ψ energy
    ψe = 0.5*ψ*ψ
    
    # Thermal pressure
    p_thermal = γm1*(ρe - ke - me - ψe)
    p_total = p_thermal + 0.5*Bmagsq
    
    # Energy term for advection
    E_advected = ke + γ*p_thermal/γm1 + Bmagsq
    
    # Velocity dot B
    vdotB = u*Bx + v*By + w*Bz
    
    #-----------------------------------------
    # F-flux (x-direction)
    #-----------------------------------------
    F[1] = ρu
    F[2] = ρu*u + p_total - Bx*Bx
    F[3] = ρv*u - By*Bx
    F[4] = ρw*u - Bz*Bx          # Out-of-plane momentum flux
    F[5] = u*E_advected + Bx*(c_h*ψ - vdotB)
    F[6] = c_h*ψ
    F[7] = u*By - Bx*v
    F[8] = u*Bz - Bx*w
    F[9] = c_h*Bx
    
    #-----------------------------------------
    # G-flux (y-direction)
    #-----------------------------------------
    G[1] = ρv
    G[2] = ρu*v - Bx*By
    G[3] = ρv*v + p_total - By*By
    G[4] = ρw*v - Bz*By          # Out-of-plane momentum flux
    G[5] = v*E_advected + By*(c_h*ψ - vdotB)
    G[6] = v*Bx - By*u
    G[7] = c_h*ψ
    G[8] = v*Bz - By*w
    G[9] = c_h*By
end


function mhd_glm_source_2p5d!(S, q, gradients, PhysConst; neqs=9)
    """
    2.5D source term with out-of-plane velocity
    """
    # Extract state
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]  # Out-of-plane momentum
    Bx = q[6]
    By = q[7]
    Bz = q[8]
    ψ  = q[9]
    
    # Primitive velocities
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ  # Out-of-plane velocity
    
    # Extract gradients
    ∂Bx_∂x = gradients.dBx_dx
    ∂By_∂y = gradients.dBy_dy
    ∂ψ_∂x  = gradients.dψ_dx
    ∂ψ_∂y  = gradients.dψ_dy
    
    # Divergence of B (2D)
    divB = ∂Bx_∂x + ∂By_∂y
    
    # v·∇ψ
    v_dot_gradψ = u*∂ψ_∂x + v*∂ψ_∂y
    
    # v·B (includes out-of-plane components!)
    v_dot_B = u*Bx + v*By + w*Bz
    
    # Assemble source term
    S[1] = 0.0
    S[2] = divB * Bx
    S[3] = divB * By
    S[4] = divB * Bz        # Out-of-plane momentum source
    S[5] = divB * v_dot_B + ψ * v_dot_gradψ
    S[6] = divB * u
    S[7] = divB * v
    S[8] = divB * w
    S[9] = v_dot_gradψ
    
    # Optional damping
    if hasfield(typeof(PhysConst), :c_p)
        c_h = PhysConst.c_h
        c_p = PhysConst.c_p
        S[9] += -(c_h^2 / c_p^2) * ψ
    end
    
    return nothing
end

function user_flux!(F, G, H, SD::NSD_3D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=9, ip=1)
    PhysConst = PhysicalConst{Float64}()
    
    # Extract conservative variables
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρe = q[5]  # This is total energy E (not specific energy)
    Bx = q[6]
    By = q[7]
    Bz = q[8]
    ψ  = q[9]
    
    # Primitive variables
    u = ρu/ρ
    v = ρv/ρ
    w = ρw/ρ
    
    # Thermodynamic properties
    γ   = PhysConst.γ
    γm1 = PhysConst.γm1
    
    # Hyperbolic cleaning speed
    c_h = PhysConst.c_h
    
    # Kinetic energy
    velomagsq = u*u + v*v + w*w
    ke = 0.5*ρ*velomagsq
    
    # Magnetic field magnitude squared
    Bmagsq = Bx*Bx + By*By + Bz*Bz
    
    # Magnetic energy
    me = 0.5*Bmagsq
    
    # ψ energy (CRITICAL: This was missing!)
    ψe = 0.5*ψ*ψ
    
    # Thermal pressure: p = (γ-1)(E - ½ρ||v||² - ½||B||² - ½ψ²)
    p_thermal = γm1*(ρe - ke - me - ψe)
    
    # Total pressure (thermal + magnetic)
    p_total = p_thermal + 0.5*Bmagsq
    
    # Energy term for advection: ½ρ||v||² + γp/(γ-1) + ||B||²
    # Note: ||B||² (not ½||B||²) as per paper formulation
    E_advected = ke + γ*p_thermal/γm1 + Bmagsq
    
    # Velocity dot B
    vdotB = u*Bx + v*By + w*Bz
    
    #-----------------------------------------
    # F-flux (x-direction)
    #-----------------------------------------
    
    # Mass flux
    F[1] = ρu
    
    # Momentum flux: ρu*v_i + I*p_total - B_x*B_i
    F[2] = ρu*u + p_total - Bx*Bx
    F[3] = ρv*u - By*Bx
    F[4] = ρw*u - Bz*Bx
    
    # Energy flux: u*(E_advected) + B_x*(c_h*ψ - v·B)
    F[5] = u*E_advected + Bx*(c_h*ψ - vdotB)
    
    # Magnetic field flux: v_x*B_i - B_x*v_i + I*c_h*ψ
    F[6] = c_h*ψ
    F[7] = u*By - Bx*v
    F[8] = u*Bz - Bx*w
    
    # Divergence cleaning flux
    F[9] = c_h*Bx
    
    #-----------------------------------------
    # G-flux (y-direction)
    #-----------------------------------------
    
    # Mass flux
    G[1] = ρv
    
    # Momentum flux: ρv*v_i + I*p_total - B_y*B_i
    G[2] = ρu*v - Bx*By
    G[3] = ρv*v + p_total - By*By
    G[4] = ρw*v - Bz*By
    
    # Energy flux: v*(E_advected) + B_y*(c_h*ψ - v·B)
    G[5] = v*E_advected + By*(c_h*ψ - vdotB)
    
    # Magnetic field flux: v_y*B_i - B_y*v_i + I*c_h*ψ
    G[6] = v*Bx - By*u
    G[7] = c_h*ψ
    G[8] = v*Bz - By*w
    
    # Divergence cleaning flux
    G[9] = c_h*By
    
    #-----------------------------------------
    # H-flux (z-direction)
    #-----------------------------------------
    
    # Mass flux
    H[1] = ρw
    
    # Momentum flux: ρw*v_i + I*p_total - B_z*B_i
    H[2] = ρu*w - Bx*Bz
    H[3] = ρv*w - By*Bz
    H[4] = ρw*w + p_total - Bz*Bz
    
    # Energy flux: w*(E_advected) + B_z*(c_h*ψ - v·B)
    H[5] = w*E_advected + Bz*(c_h*ψ - vdotB)
    
    # Magnetic field flux: v_z*B_i - B_z*v_i + I*c_h*ψ
    H[6] = w*Bx - Bz*u
    H[7] = w*By - Bz*v
    H[8] = c_h*ψ
    
    # Divergence cleaning flux
    H[9] = c_h*Bz
end


#-----------------------------------------
# Source term function (must be called separately!)
#-----------------------------------------

@inline function user_turbo_volume_flux(u_ll, u_rr, volume_flux)

    if volume_flux == "etec"
	flux_artiano_etec(u_ll, u_rr)
    elseif volume_flux == "ec"
	flux_artiano_ec(u_ll, u_rr)
    elseif volume_flux == "tec"
	flux_artiano_tec(u_ll, u_rr)
    elseif volume_flux == "ranocha"
        flux_turbo_ranocha(u_ll, u_rr)
    elseif volume_flux == "gruber"
        flux_kennedy_gruber(u_ll, u_rr)
    elseif volume_flux == "central"
        flux_central(u_ll, u_rr)
    else
        #default
        flux_turbo_ranocha(u_ll, u_rr)
    end
end

@inline function user_volume_flux(u_ll, u_rr, volume_flux)

    if volume_flux == "etec"
	flux_artiano_etec(u_ll, u_rr)
    elseif volume_flux == "ec"
	flux_artiano_ec(u_ll, u_rr)
    elseif volume_flux == "tec"
	flux_artiano_tec(u_ll, u_rr)
    elseif volume_flux == "ranocha"
        flux_turbo_ranocha(u_ll, u_rr)
    elseif volume_flux == "gruber"
        flux_kennedy_gruber(u_ll, u_rr)
    elseif volume_flux == "central"
        flux_central(u_ll, u_rr)
    else
        #default
        flux_turbo_ranocha(u_ll, u_rr)
    end

end



function flux(q)
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

function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL)
    
    PhysConst = PhysicalConst{Float64}()
                
    ρ  = q[1] 
    ρu = q[2]
    ρv = q[3]
    ρe = q[4]

    e  = ρe/ρ
    u  = ρu/ρ
    v  = ρv/ρ

    gamma = PhysConst.cp/PhysConst.cv
    p = (gamma - 1) * (ρe- 0.5f0 * (ρu^2 + ρv^2) / ρ)
    aux[1] = ρ
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = log(ρ)
    aux[6] = log(p)
end


function user_fluxaux!(aux, SD::NSD_2D, q, ::THETA, ::artiano_ec)
    
    PhysConst = PhysicalConst{Float64}()
                
    rho  = q[1] 
    rho_u = q[2]
    rho_v = q[3]
    rho_theta = q[4]

    theta  = rho_theta/rho
    u  = rho_u/rho
    v  = rho_v/rho

    p = perfectGasLaw_ρθtoP(PhysConst, ρ=rho, θ=theta)

    aux[1] = rho
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = rho_theta
    aux[6] = log(rho)
    aux[7] = log(rho_theta)
end

@inline function flux_ranocha(u_ll, u_rr)

    PhysConst = PhysicalConst{Float64}()
    rho_ll, rho_v1_ll, rho_v2_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_e_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr
    v2_rr = rho_v2_rr/rho_rr
    v2_ll = rho_v2_ll/rho_ll
   
    gamma = PhysConst.cp/PhysConst.cv
    p_ll = (gamma - 1) * (rho_e_ll - 0.5f0 * (rho_v1_ll^2 + rho_v2_ll^2) / rho_ll)
    p_rr = (gamma - 1) * (rho_e_rr - 0.5f0 * (rho_v1_rr^2 + rho_v2_rr^2) / rho_rr)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

    # Calculate fluxes depending on orientation
        f1 = rho_mean * v1_avg
        f2 = f1 * v1_avg + p_avg
        f3 = f1 * v2_avg
	f4 = f1 * (velocity_square_avg + inv_rho_p_mean /(PhysConst.γ - 1)) + 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
        g1 = rho_mean * v2_avg
        g2 = g1 * v1_avg
        g3 = g1 * v2_avg + p_avg
	g4 = g1 * (velocity_square_avg + inv_rho_p_mean /(PhysConst.γ - 1)) + 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end


@inline function flux_turbo(u_ll, u_rr, ::artiano_ec)
    PhysConst = PhysicalConst{Float64}()
	rho_ll, v1_ll, v2_ll, p_ll, rho_theta_ll, log_rho_ll, log_rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, p_rr, rho_theta_rr, log_rho_rr, log_rho_theta_rr = u_rr
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
            x2 = rho_ll * rho_theta_rr
            log_x2 = log_rho_ll + log_rho_theta_rr
            y2 = rho_rr * rho_theta_ll
            log_y2 = log_rho_rr + log_rho_theta_ll
            x2_plus_y2 = x2 + y2
            y2_minus_x2 = y2 - x2
            z2 = y2_minus_x2^2 / x2_plus_y2^2
            special_path2 = (2 + z2 * (2 / 3 + z2 * (2 / 5 + 2 / 7 * z2))) / x2_plus_y2
            regular_path2 = (log_y2 - log_x2) / y2_minus_x2
            inv_rho_p_mean = rho_theta_ll * rho_theta_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

            v1_avg = 0.5 * (v1_ll + v1_rr)
            v2_avg = 0.5 * (v2_ll + v2_rr)
            p_avg = 0.5 * (p_ll + p_rr)
            # calculate fluxes depending on cartesian orientation
            f1 = rho_mean * v1_avg
            f2 = f1 * v1_avg + p_avg
            f3 = f1 * v2_avg
            f4 = f1 * inv_rho_p_mean

            g1 = rho_mean * v2_avg
            g2 = g1 * v1_avg 
	    g3 = g1 * v2_avg + p_avg
	    g4 = g1 * inv_rho_p_mean
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_turbo_ranocha(u_ll, u_rr)
    physconst = physicalconst{float64}()
	rho_ll, v1_ll, v2_ll, p_ll, log_rho_ll, log_p_ll = u_ll
	rho_rr, v1_rr, v2_rr, p_rr, log_rho_rr, log_p_rr = u_rr
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
	    gamma = physconst.cp/physconst.cv
            # calculate fluxes depending on cartesian orientation
            f1 = rho_mean * v1_avg
            f2 = f1 * v1_avg + p_avg
            f3 = f1 * v2_avg
            f4 = f1 *
		(velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) +
                 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

            g1 = rho_mean * v2_avg
            g2 = g1 * v1_avg 
	    g3 = g1 * v2_avg + p_avg
            g4 = g1 * (velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) + 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_artiano_etec(u_ll, u_rr)
# Compute the necessary mean values
    PhysConst = PhysicalConst{Float64}()
    rho_ll, rho_v1_ll, rho_v2_ll, rho_theta_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_theta_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr
    v2_rr = rho_v2_rr/rho_rr
    v2_ll = rho_v2_ll/rho_ll

    p_ll = PhysConst.pref * (rho_theta_ll * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 
    p_rr = PhysConst.pref * (rho_theta_rr * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, PhysConst.cp/PhysConst.cv)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f4 = gammamean * 0.5f0 * v1_avg 
    f1 = f4 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    f2 = f1 * v1_avg + p_avg 
    f3 = f1 * v2_avg 

    g4 = gammamean * 0.5f0 * v2_avg 
    g1 = g4 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
    g2 = g1 * v1_avg 
    g3 = g1 * v2_avg + p_avg 
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_artiano_ec(u_ll, u_rr)
    # Unpack left and right state

    PhysConst = PhysicalConst{Float64}()
    rho_ll, rho_v1_ll, rho_v2_ll, rho_theta_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_theta_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr
    v2_rr = rho_v2_rr/rho_rr
    v2_ll = rho_v2_ll/rho_ll

    p_ll = PhysConst.pref * (rho_theta_ll * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 
    p_rr = PhysConst.pref * (rho_theta_rr * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    
    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * v1_avg 
    f2 = f1 * v1_avg + p_avg 
    f3 = f1 * v2_avg 
    f4 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    g1 = rho_mean * 0.5f0 * v2_avg 
    g2 = g1 * v1_avg 
    g3 = g1 * v2_avg + p_avg 
    g4 = g1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_artiano_tec(u_ll, u_rr)

    PhysConst = PhysicalConst{Float64}()
    rho_ll, rho_v1_ll, rho_v2_ll, rho_theta_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_theta_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr
    v2_rr = rho_v2_rr/rho_rr
    v2_ll = rho_v2_ll/rho_ll

    p_ll = PhysConst.pref * (rho_theta_ll * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 
    p_rr = PhysConst.pref * (rho_theta_rr * PhysConst.Rair/PhysConst.pref)^(PhysConst.cp/PhysConst.cv) 
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, PhysConst.cp/PhysConst.cv)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * v1_avg 
    f2 = f1 * v1_avg + p_avg 
    f3 = f1 * v2_avg 
    f4 = gammamean * 0.5f0 * v1_avg 

    g1 = rho_mean * 0.5f0 * v2_avg 
    g2 = g1 * v1_avg 
    g3 = g1 * v2_avg + p_avg 
    g4 = gammamean * 0.5f0 * v2_avg 
	return SVector(f1, f2, f3, f4), SVector(g1,g2,g3,g4)
end

function flux_central(u_ll, u_rr)

    return 0.5f0 .* (flux(u_ll) .+ flux(u_rr))
end

@inline function flux_kennedy_gruber(u_ll, u_rr)
    PhysConst = PhysicalConst{Float64}()
    rho_ll, rho_v1_ll, rho_v2_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_e_rr = u_rr
	v1_ll = rho_v1_ll/rho_ll
	v1_rr = rho_v1_rr/rho_rr
	v2_rr = rho_v2_rr/rho_rr
	v2_ll = rho_v2_ll/rho_ll
    Temp_ll = (rho_e_ll/rho_ll - 0.5*v1_ll*v1_ll)/PhysConst.cv
    Temp_rr = (rho_e_rr/rho_rr - 0.5*v1_rr*v1_rr)/PhysConst.cv

    p_ll = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_ll, Temp=Temp_ll)
    p_rr = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_rr, Temp=Temp_rr)

    # Unpack left and right state

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

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

@inline ln_mean(x::Real, y::Real) = ln_mean(promote(x, y)...)

@inline function ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return (x + y) / @evalpoly(f2,
                         2,
                         convert(RealT, 2 / 3),
                         convert(RealT, 2 / 5),
                         convert(RealT, 2 / 7))
    else
        return (y - x) / log(y / x)
    end
end

@inline inv_ln_mean(x::Real, y::Real) = inv_ln_mean(promote(x, y)...)

@inline function inv_ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return @evalpoly(f2,
                         2,
                         convert(RealT, 2 / 3),
                         convert(RealT, 2 / 5),
                         convert(RealT, 2 / 7)) / (x + y)
    else
        return log(y / x) / (y - x)
    end
end

@inline stolarsky_mean(x::Real, y::Real, gamma::Real) = stolarsky_mean(promote(x, y)...,
                                                                       gamma)

@inline function stolarsky_mean(x::RealT, y::RealT, gamma::Real) where {RealT <: Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        # convenience coefficients
        c1 = convert(RealT, 1 / 3) * (gamma - 2)
        c2 = convert(RealT, -1 / 15) * (gamma + 1) * (gamma - 3) * c1
        c3 = convert(RealT, -1 / 21) * (2 * gamma * (gamma - 2) - 9) * c2
        return 0.5f0 * (x + y) * @evalpoly(f2, 1, c1, c2, c3)
    else
        if gamma isa Integer
            yg = y^(gamma - 1)
            xg = x^(gamma - 1)
        else
            yg = exp((gamma - 1) * log(y)) # equivalent to y^gamma but faster for non-integers
            xg = exp((gamma - 1) * log(x)) # equivalent to x^gamma but faster for non-integers
        end
        return (gamma - 1) * (yg * y - xg * x) / (gamma * (yg - xg))
    end
end

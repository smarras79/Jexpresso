function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρe = q[4]

    e  = ρe/ρ
    u  = ρu/ρ
    v  = ρv/ρ

    γ   = PhysConst.γ
    γm1 = PhysConst.γm1
    
    velomagsq = (u*u + v*v)
    ke        = 0.5*ρ*velomagsq
    Pressure  = γm1*(ρe - ke)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρv*u
    F[4] = u*(ke + γ*Pressure/γm1)

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure
    G[4] = v*(ke + γ*Pressure/γm1)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1] + qe[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4] + qe[4]
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    Press = Press - qe[end]
    
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
end

function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::THETA; neqs=4, ip=1)
    
    PhysConst = PhysicalConst{Float64}()
                
    ρ  = q[1] 
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4] 
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
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

    γ   = PhysConst.γ
    γm1 = γ - 1.0
    
    velomagsq = (u*u + v*v)
    ke        = 0.5*ρ*velomagsq
    p  = γm1*(ρe - velomagsq)
    aux[1] = ρ
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = log(ρ)
    aux[6] = log(p)
end


function user_flux!(F, G, SD::NSD_2D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::NCL, ::AbstractPert; neqs=4, ip=1)
    
    PhysConst = PhysicalConst{Float64}()
                
    ρ = q[1]
    u = q[2]
    v = q[3]
    θ = q[4]
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρ*u
    F[2] = u
    F[3] = v
    F[4] = θ
    
    G[1] = ρ*v
    G[2] = u
    G[3] = v
    G[4] = θ
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    if (lpert)
        ρ  = q[1]+qe[1]
        ρu = q[2]+qe[2]
        ρv = q[3]+qe[3]
        ρθ = q[4]+qe[4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ) - qe[5]
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρθ*u), T(ρv),T(ρu*v),T(ρv*v + Pressure),T(ρθ*v)
    else
        ρ  = q[1]
        ρu = q[2]
        ρv = q[3]
        ρθ = q[4]
        θ  = ρθ/ρ
        u  = ρu/ρ
        v  = ρv/ρ
        Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
        return T(ρu), T(ρu*u + Pressure), T(ρv*u), T(ρθ*u), T(ρv),T(ρu*v),T(ρv*v + Pressure),T(ρθ*v)
    end
end


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

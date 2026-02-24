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

@inline function flux(q, ::central_theta) 
    PhysConst = PhysicalConst{Float64}()
                
    ρ  = q[1] 
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4] 
    
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    f1 = ρu
    f2 = ρu*u + Press
    f3 = ρv*u
    f4 = ρθ*u
    
    g1 = ρv
    g2 = ρu*v
    g3 = ρv*v + Press
    g4 = ρθ*v
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)

end

function user_fluxaux!(aux, SD::NSD_2D, q, ::THETA, ::central_theta)
    
    PhysConst = PhysicalConst{Float64}()

    aux[1] = q[1]
    aux[2] = q[2]
    aux[3] = q[3]
    aux[4] = q[4]
end


function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::central_euler)
    
    PhysConst = PhysicalConst{Float64}()

    aux[1] = q[1]
    aux[2] = q[2]
    aux[3] = q[3]
    aux[4] = q[4]
end

function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::kennedy_gruber)
    
    PhysConst = PhysicalConst{Float64}()
                
    rho  = q[1] 
    rho_u = q[2]
    rho_v = q[3]
    rho_e = q[4]

    u  = rho_u/rho
    v  = rho_v/rho
	
    γ   = PhysConst.γ
    gammam1 = γ - 1.0
    p = gammam1 * (rho_e - (0.5 * rho_u * u + 0.5 * v * rho_v)) 

    aux[1] = rho
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = rho_e/rho
end

function user_fluxaux!(aux, SD::NSD_2D, q, ::TOTAL, ::ranocha)
    
    PhysConst = PhysicalConst{Float64}()
                
    rho  = q[1] 
    rho_u = q[2]
    rho_v = q[3]
    rho_e = q[4]

    u  = rho_u/rho
    v  = rho_v/rho
	
    γ   = PhysConst.γ
    gammam1 = γ - 1.0
    p = gammam1 * (rho_e - (0.5 * rho_u * u + 0.5 * v * rho_v)) 

    aux[1] = rho
    aux[2] = u
    aux[3] = v
    aux[4] = p
    aux[5] = rho_e
    aux[6] = log(rho)
    aux[7] = log(p)
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

function user_fluxaux!(aux, SD::NSD_2D, q, ::THETA, ::artiano_tec)
    
    PhysConst = PhysicalConst{Float64}()
                
    rho  = q[1] 
    rho_u = q[2]
    rho_v = q[3]
    rho_theta = q[4]

    γ   = PhysConst.γ
    gammam1 = γ - 1.0
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
    aux[7] = rho_theta^gammam1
end
function user_fluxaux!(aux, SD::NSD_2D, q, ::THETA, ::artiano_etec)
    
    PhysConst = PhysicalConst{Float64}()
                
    rho  = q[1] 
    rho_u = q[2]
    rho_v = q[3]
    rho_theta = q[4]

    γ   = PhysConst.γ
    gammam1 = γ - 1.0
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
    aux[8] = rho_theta^gammam1
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
            f4 = f1 *
		(velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) +
                 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

            g1 = rho_mean * v2_avg
            g2 = g1 * v1_avg 
	    g3 = g1 * v2_avg + p_avg
            g4 = g1 * (velocity_square_avg + inv_rho_p_mean * 1/(gamma - 1)) + 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
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

            # algebraically equivalent to `inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)`
            # in exact arithmetic since
            #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/rho_thetaₗ - ϱᵣ/rho_thetaᵣ)
            #   = rho_thetaₗ rho_thetaᵣ log((ϱₗ rho_thetaᵣ) / (ϱᵣ rho_thetaₗ)) / (ϱₗ rho_thetaᵣ - ϱᵣ rho_thetaₗ)
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
            inv_rho_rho_theta_mean = rho_theta_ll * rho_theta_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

            v1_avg = 0.5 * (v1_ll + v1_rr)
            v2_avg = 0.5 * (v2_ll + v2_rr)
            p_avg = 0.5 * (p_ll + p_rr)
            # calculate fluxes depending on cartesian orientation
            f1 = rho_mean * v1_avg
            f2 = f1 * v1_avg + p_avg
            f3 = f1 * v2_avg
            f4 = f1 * inv_rho_rho_theta_mean

            g1 = rho_mean * v2_avg
            g2 = g1 * v1_avg 
	    g3 = g1 * v2_avg + p_avg
	    g4 = g1 * inv_rho_rho_theta_mean
    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

@inline function flux_turbo(u_ll, u_rr, ::artiano_tec)


	rho_ll, v1_ll, v2_ll, p_ll, rho_theta_ll, log_rho_ll, rho_theta_pow_gamma_m1_ll = u_ll
	rho_rr, v1_rr, v2_rr, p_rr, rho_theta_rr, log_rho_rr, rho_theta_pow_gamma_m1_rr = u_rr
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

x3 = rho_theta_ll
y3 = rho_theta_rr
x3_plus_y3 = x3 + y3
y3_minus_x3 = y3 - x3

numerator_f2 = x3 * (y3_minus_x3 - y3) + y3 * y3
denominator_f2 = x3 * (x3_plus_y3 + x3) + y3 * y3
f2 = numerator_f2 / denominator_f2

    PhysConst = PhysicalConst{Float64}()
                

    gamma   = PhysConst.γ
c1 = (gamma - 2) / 3
c2 = -(gamma + 1) * (gamma - 3) * c1 / 15
c3 = -(2 * gamma * (gamma - 2) - 9) * c2 / 21

special_path3 = 0.5 * x3_plus_y3 * (1 + f2 * (c1 + f2 * (c2 + f2 * c3)))

yg_minus_xg = rho_theta_pow_gamma_m1_rr - rho_theta_pow_gamma_m1_ll
regular_path3 = (gamma - 1) * (rho_theta_pow_gamma_m1_rr * y3 - rho_theta_pow_gamma_m1_ll * x3) / (gamma * yg_minus_xg)

gammamean = ifelse(f2 < 1.0e-4, special_path3, regular_path3)

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * v1_avg 
    f2 = f1 * v1_avg + p_avg 
    f3 = f1 * v2_avg 
    f4 = gammamean * v1_avg 

    g1 = rho_mean * v2_avg 
    g2 = g1 * v1_avg 
    g3 = g1 * v2_avg + p_avg 
    g4 = gammamean * v2_avg 
	return SVector(f1, f2, f3, f4), SVector(g1,g2,g3,g4)
end

@inline function flux_turbo(u_ll, u_rr, ::artiano_etec)
    rho_ll, v1_ll, v2_ll, p_ll, rho_theta_ll, log_rho_ll, log_rho_theta_ll, rho_theta_pow_gamma_m1_ll = u_ll
    rho_rr, v1_rr, v2_rr, p_rr, rho_theta_rr, log_rho_rr, log_rho_theta_rr, rho_theta_pow_gamma_m1_rr = u_rr
    
    PhysConst = PhysicalConst{Float64}()
    gamma = PhysConst.γ
    
   
	x_ratio = rho_ll / rho_theta_ll
log_x_ratio = log_rho_ll - log_rho_theta_ll  # log(a/b) = log(a) - log(b)

y_ratio = rho_rr / rho_theta_rr  
log_y_ratio = log_rho_rr - log_rho_theta_rr

# Calcolo ln_mean standard
x_ratio_plus_y = x_ratio + y_ratio
y_ratio_minus_x = y_ratio - x_ratio

# f2 per decidere special/regular path
numerator_f2 = x_ratio * (x_ratio - 2*y_ratio) + y_ratio * y_ratio
denominator_f2 = x_ratio * (x_ratio + 2*y_ratio) + y_ratio * y_ratio
f2_ratio = numerator_f2 / denominator_f2

# Special path (Taylor expansion)
special_ratio = x_ratio_plus_y / (2 + f2_ratio * (2/3 + f2_ratio * (2/5 + 2/7 * f2_ratio)))

# Regular path
regular_ratio = y_ratio_minus_x / (log_y_ratio - log_x_ratio)

ln_rho_over_theta = ifelse(f2_ratio < 1.0e-4, special_ratio, regular_ratio)

    x_theta = rho_theta_ll
    y_theta = rho_theta_rr
    x_theta_plus_y = x_theta + y_theta
    y_theta_minus_x = y_theta - x_theta
    
    numerator_f2 = x_theta * (x_theta - 2*y_theta) + y_theta * y_theta
    denominator_f2 = x_theta * (x_theta + 2*y_theta) + y_theta * y_theta
    f2 = numerator_f2 / denominator_f2
    
    c1 = (gamma - 2) / 3
    c2 = -(gamma + 1) * (gamma - 3) * c1 / 15
    c3 = -(2 * gamma * (gamma - 2) - 9) * c2 / 21
    
    special_gamma = 0.5 * x_theta_plus_y * (1 + f2 * (c1 + f2 * (c2 + f2 * c3)))
    
    yg_minus_xg = rho_theta_pow_gamma_m1_rr - rho_theta_pow_gamma_m1_ll
    regular_gamma = (gamma - 1) * (rho_theta_pow_gamma_m1_rr * y_theta - rho_theta_pow_gamma_m1_ll * x_theta) / (gamma * yg_minus_xg)
    
    gammamean = ifelse(f2 < 1.0e-4, special_gamma, regular_gamma)


    v1_avg = 0.5 * (v1_ll + v1_rr)
    v2_avg = 0.5 * (v2_ll + v2_rr)
    p_avg = 0.5 * (p_ll + p_rr)
    
    f4 = gammamean * v1_avg 
    f1 = f4 * ln_rho_over_theta  
    f2 = f1 * v1_avg + p_avg 
    f3 = f1 * v2_avg 
    
    g4 = gammamean * v2_avg 
    g1 = g4 * ln_rho_over_theta
    g2 = g1 * v1_avg 
    g3 = g1 * v2_avg + p_avg 
    
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

@inline function flux_turbo(u_ll, u_rr, sol_type::central_theta)
	return 0.5f0 .* (flux(u_ll,sol_type) .+ flux(u_rr, sol_type))
end

@inline function flux_turbo(u_ll, u_rr, sol_type::central_euler)
	return 0.5f0 .* (flux(u_ll,sol_type) .+ flux(u_rr, sol_type))
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

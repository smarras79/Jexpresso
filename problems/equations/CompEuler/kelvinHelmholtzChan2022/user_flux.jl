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
    γm1 = γ - 1.0
    
    velomagsq = (u*u + v*v)
    ke        = 0.5*ρ*velomagsq
    Pressure  = γm1*(ρe - velomagsq)
    
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
    Press = Press 
    
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


@inline function user_volume_flux(u_ll, u_rr)
	flux_artiano_ec(u_ll, u_rr)
	# flux_ranocha(u_ll, u_rr)
#	flux_kennedy_gruber(u_ll, u_rr)
#        flux_central(u_ll, u_rr)
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
    Pressure  = γm1*(ρe - velomagsq)
    
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
    f3 = f1 * v2_avg + p_avg 
    f4 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    g1 = rho_mean * 0.5f0 * v1_avg 
    g2 = f1 * v1_avg + p_avg 
    g3 = f1 * v2_avg + p_avg 
    g4 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end

function flux_central(u_ll, u_rr)

    return 0.5f0 .* (flux(u_ll) .+ flux(u_rr))
end
@inline function flux_ranocha(u_ll, u_rr)

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
        g2 = f1 * v1_avg
        g3 = f1 * v2_avg + p_avg
	g4 = f1 * (velocity_square_avg + inv_rho_p_mean /(PhysConst.γ - 1)) + 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
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

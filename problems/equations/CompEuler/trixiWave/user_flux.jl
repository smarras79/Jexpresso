function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]
    E  = ρE/ρ
    u  = ρu/ρ

    Temp = (E - 0.5*u*u)/PhysConst.cv
    Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)

    #@info " FLUX USER: " Temp Press ρ
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρE*u + Press*u

end

function flux(q)
    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρE = q[3]
    E  = ρE/ρ
    u  = ρu/ρ

    Temp = (E - 0.5*u*u)/PhysConst.cv
    Press = perfectGasLaw_ρTtoP(PhysConst; ρ=ρ, Temp=Temp)

    #@info " FLUX USER: " Temp Press ρ
    return SVector(ρu,
            ρu*u + Press,
            ρE*u + Press*u)
end

function user_volume_flux(u_ll, u_rr)
         flux_ranocha(u_ll, u_rr)
        # flux_kennedy_gruber(u_ll, u_rr)
        # flux_chandrashekar(u_ll, u_rr)
        # flux_central(u_ll, u_rr)
end

function flux_ranocha(u_ll, u_rr)

    PhysConst = PhysicalConst{Float64}()

    rho_ll, rho_v1_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_e_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr

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
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr)

    # Calculate fluxes
    # Ignore orientation since it is always "1" in 1D
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = f1 * (velocity_square_avg + inv_rho_p_mean * (PhysConst.cv/(PhysConst.cp - PhysConst.cv))) +
         0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)

    return SVector(f1, f2, f3)
end

function flux_central(u_ll, u_rr)

    return 0.5f0 * (flux(u_ll) + flux(u_rr))
end

@inline function flux_kennedy_gruber(u_ll, u_rr)
    PhysConst = PhysicalConst{Float64}()

    rho_ll, rho_v1_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_e_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr

    Temp_ll = (rho_e_ll/rho_ll - 0.5*v1_ll*v1_ll)/PhysConst.cv
    Temp_rr = (rho_e_rr/rho_rr - 0.5*v1_rr*v1_rr)/PhysConst.cv

    p_ll = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_ll, Temp=Temp_ll)
    p_rr = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_rr, Temp=Temp_rr)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Ignore orientation since it is always "1" in 1D
    f1 = rho_avg * v1_avg
    f2 = rho_avg * v1_avg * v1_avg + p_avg
    f3 = (rho_avg * e_avg + p_avg) * v1_avg

    return SVector(f1, f2, f3)
end

@inline function flux_chandrashekar(u_ll, u_rr)
    PhysConst = PhysicalConst{Float64}()

    rho_ll, rho_v1_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_e_rr = u_rr
    v1_ll = rho_v1_ll/rho_ll
    v1_rr = rho_v1_rr/rho_rr

    Temp_ll = (rho_e_ll/rho_ll - 0.5*v1_ll*v1_ll)/PhysConst.cv
    Temp_rr = (rho_e_rr/rho_rr - 0.5*v1_rr*v1_rr)/PhysConst.cv

    p_ll = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_ll, Temp=Temp_ll)
    p_rr = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_rr, Temp=Temp_rr)
    beta_ll = 0.5f0 * rho_ll / p_ll
    beta_rr = 0.5f0 * rho_rr / p_rr
    specific_kin_ll = 0.5f0 * (v1_ll^2)
    specific_kin_rr = 0.5f0 * (v1_rr^2)

    # Compute the necessary mean values
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    rho_mean = ln_mean(rho_ll, rho_rr)
    beta_mean = ln_mean(beta_ll, beta_rr)
    beta_avg = 0.5f0 * (beta_ll + beta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_mean = 0.5f0 * rho_avg / beta_avg
    velocity_square_avg = specific_kin_ll + specific_kin_rr

    # Calculate fluxes
    # Ignore orientation since it is always "1" in 1D
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_mean
    f3 = f1 * 0.5f0 * (1 / (PhysConst.cp/PhysConst.cv - 1) / beta_mean - velocity_square_avg) +
         f2 * v1_avg

    return SVector(f1, f2, f3)
end

## Functions from Trixi.jl

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

@inline function cons2entropy(u)
    PhysConst = PhysicalConst{Float64}()
    rho, rho_v1, rho_e = u
    gamma = PhysConst.cp/PhysConst.cv
    v1 = rho_v1 / rho
    v_square = v1^2
    p = (gamma- 1) * (rho_e - 0.5f0 * rho * v_square)
    s = log(p) - gamma * log(rho)
    rho_p = rho / p

    w1 = (gamma - s) * 1/(gamma-1) -
         0.5f0 * rho_p * v_square
    w2 = rho_p * v1
    w3 = -rho_p

    return SVector(w1, w2, w3)
end

@inline function entropy_thermodynamic(cons)
    PhysConst = PhysicalConst{Float64}()
    gamma = PhysConst.cp/PhysConst.cv
    p = (gamma - 1) * (cons[3] - 0.5f0 * (cons[2]^2) / cons[1])

    # Thermodynamic entropy
    s = log(p) - gamma * log(cons[1])

    return -s * cons[1] * 1/(gamma -1)
end
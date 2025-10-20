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

## Ranocha flux
function user_volume_flux(u_ll, u_rr)
    # Unpack left and right state

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
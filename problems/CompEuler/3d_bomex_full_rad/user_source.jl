function user_source!(S,
                      q,
                      qe,
                      npoin,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]
    ρhl = q[5]
    ρqt = q[6]
    ρqp = q[7]

    ρhl_ref = qe[5]
    ρqt_ref = qe[6]

    # Coriolis (f-plane, BOMEX latitude ≈ 15°N)
    u_geo = 0.0
    f0 = 0.376e-4  # s⁻¹
    buc = -f0 * ρv
    bvc =  f0 * (ρu - u_geo*ρ)

    # Rayleigh sponge layer above 2400 m
    z_sponge = 2400.0
    z_max    = 3000.0
    α_max    = 0.75
    γ        = 2.0
    β_sponge = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / 2)^γ
    end

    S[1] = 0.0
    S[2] = buc - β_sponge * ρu
    S[3] = bvc - β_sponge * ρv
    S[4] = -ρ * PhysConst.g - β_sponge * ρw
    S[5] = 0.0#-β_sponge*(ρhl-ρhl_ref)
    S[6] = 0.0#-β_sponge*(ρqt-ρqt_ref)
    S[7] = 0.0#-β_sponge*(ρqp)
end

function user_source!(S,
                      q,
                      qe,
                      npoin,
                      ::CL, ::PERT;
                      neqs=1,
                      x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    u_geo = 0.0
    f0 = 0.376e-4
    buc = -f0 * ρv
    bvc =  f0 * (ρu - u_geo*ρ)

    z_sponge = 2400.0
    z_max    = 3000.0
    α_max    = 0.75
    γ        = 2.0
    β_sponge = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / 2)^γ
    end

    S[1] = 0.0
    S[2] = buc - β_sponge * ρu
    S[3] = bvc - β_sponge * ρv
    S[4] = -ρ * PhysConst.g - β_sponge * ρw
    S[5] = 0.0
    S[6] = 0.0
    S[7] = 0.0
end

function user_source_gpu(q, qe, x, y, z, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    T = eltype(q)
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    f0 = T(0.376e-4)
    buc = T(-f0 * ρv)
    bvc = T( f0 * ρu)

    z_sponge = T(2400.0)
    z_max    = T(3000.0)
    α_max    = T(0.75)
    β_sponge = T(0.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / T(2))^T(2)
    end

    return T(0.0), T(buc - β_sponge*ρu), T(bvc - β_sponge*ρv),
           T(-ρ*PhysConst.g - β_sponge*ρw), T(0.0), T(0.0), T(0.0)
end

function user_scattering_functions(θ, θ1, ϕ, ϕ1, g)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)
    return (1 - g^2) / ((4*π) * (1 + g^2 - 2*g*cos_Θ))^(3/2)
end


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

    # Top Rayleigh sponge (5.5–8.1 km)
    z_sponge = 5500.0
    z_top    = 8100.0
    α_top    = 0.75
    γ        = 2.0
    β_top    = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / 2)^γ
    end

    # Lateral sponge: absorb outgoing gravity waves near the periodic boundaries
    # before they re-enter and form standing-wave interference patterns.
    # α_lat = 0.25 so that β_top + β_lat ≤ 1.0 everywhere (no corner instability).
    # β_x and β_y are kept separate and combined with max() to avoid doubling
    # in the x-y corner regions.
    α_lat    = 0.25
    lat_frac = 0.15
    β_x      = 0.0
    β_y      = 0.0
    dx = lat_frac * (xmax - xmin)
    dy = lat_frac * (ymax - ymin)
    if x < xmin + dx
        r = (xmin + dx - x) / dx
        β_x = α_lat * sinpi(r / 2)^2
    elseif x > xmax - dx
        r = (x - (xmax - dx)) / dx
        β_x = α_lat * sinpi(r / 2)^2
    end
    if y < ymin + dy
        r = (ymin + dy - y) / dy
        β_y = α_lat * sinpi(r / 2)^2
    elseif y > ymax - dy
        r = (y - (ymax - dy)) / dy
        β_y = α_lat * sinpi(r / 2)^2
    end
    β_lat = max(β_x, β_y)

    β = β_top + β_lat

    S[1] = 0.0
    S[2] = -β * ρu
    S[3] = -β * ρv
    S[4] = -ρ * PhysConst.g - β * ρw
    S[5] = 0.0
    S[6] = 0.0
    S[7] = 0.0
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

    # Top Rayleigh sponge (5.5–8.1 km)
    z_sponge = 5500.0
    z_top    = 8100.0
    α_top    = 0.75
    γ        = 2.0
    β_top    = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / 2)^γ
    end

    # Lateral sponge
    α_lat    = 0.25
    lat_frac = 0.15
    β_x      = 0.0
    β_y      = 0.0
    dx = lat_frac * (xmax - xmin)
    dy = lat_frac * (ymax - ymin)
    if x < xmin + dx
        r = (xmin + dx - x) / dx
        β_x = α_lat * sinpi(r / 2)^2
    elseif x > xmax - dx
        r = (x - (xmax - dx)) / dx
        β_x = α_lat * sinpi(r / 2)^2
    end
    if y < ymin + dy
        r = (ymin + dy - y) / dy
        β_y = α_lat * sinpi(r / 2)^2
    elseif y > ymax - dy
        r = (y - (ymax - dy)) / dy
        β_y = α_lat * sinpi(r / 2)^2
    end
    β_lat = max(β_x, β_y)

    β = β_top + β_lat

    S[1] = 0.0
    S[2] = -β * ρu
    S[3] = -β * ρv
    S[4] = -ρ * PhysConst.g - β * ρw
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

    # Top sponge
    z_sponge = T(5500.0)
    z_top    = T(8100.0)
    α_top    = T(0.75)
    β_top    = T(0.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / T(2))^T(2)
    end

    # Lateral sponge
    α_lat    = T(0.25)
    lat_frac = T(0.15)
    β_x      = T(0.0)
    β_y      = T(0.0)
    dx = lat_frac * (xmax - xmin)
    dy = lat_frac * (ymax - ymin)
    if x < xmin + dx
        r = (xmin + dx - x) / dx
        β_x = α_lat * sinpi(r / T(2))^T(2)
    elseif x > xmax - dx
        r = (x - (xmax - dx)) / dx
        β_x = α_lat * sinpi(r / T(2))^T(2)
    end
    if y < ymin + dy
        r = (ymin + dy - y) / dy
        β_y = α_lat * sinpi(r / T(2))^T(2)
    elseif y > ymax - dy
        r = (y - (ymax - dy)) / dy
        β_y = α_lat * sinpi(r / T(2))^T(2)
    end
    β_lat = max(β_x, β_y)

    β = β_top + β_lat

    return T(0.0), T(-β*ρu), T(-β*ρv),
           T(-ρ*PhysConst.g - β*ρw), T(0.0), T(0.0), T(0.0)
end

function user_scattering_functions(θ, θ1, ϕ, ϕ1, g)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)
    return (1 - g^2) / ((4*π) * (1 + g^2 - 2*g*cos_Θ))^(3/2)
end

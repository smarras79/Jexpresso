function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    #
    # Beta-plane Coriolis + sponge-layer damping for TC3:
    #   f(y) = f₀ + β₀ * y
    #   σ(x,y) = σ_max * max(rx, ry)  with cubic polynomial near boundaries
    #
    #   S₁ = -σ η
    #   S₂ = +f(y) v - σ u
    #   S₃ = -f(y) u - σ v
    #
    f0    = 1.0e-4    # base Coriolis parameter [s⁻¹]
    beta0 = 2.0e-11   # beta parameter [m⁻¹ s⁻¹]

    # Sponge layer parameters
    Lx = xmax - xmin
    Ly = ymax - ymin
    sponge_frac     = 0.12
    sponge_sigma_max = 0.0 #5.0e-3  # [s⁻¹]
    sx = sponge_frac * Lx
    sy = sponge_frac * Ly

    η = q[1]
    u = q[2]
    v = q[3]

    # Beta-plane Coriolis
    f = f0 + beta0 * y

    # Sponge coefficient at (x, y)
    dx_left  = x - xmin
    dx_right = xmax - x
    dy_bot   = y - ymin
    dy_top   = ymax - y

    rx = 0.0
    if dx_left < sx
        rx = max(rx, ((sx - dx_left) / sx)^3)
    end
    if dx_right < sx
        rx = max(rx, ((sx - dx_right) / sx)^3)
    end

    ry = 0.0
    if dy_bot < sy
        ry = max(ry, ((sy - dy_bot) / sy)^3)
    end
    if dy_top < sy
        ry = max(ry, ((sy - dy_top) / sy)^3)
    end

    σ = sponge_sigma_max * max(rx, ry)

    S[1] = -σ * η
    S[2] =  f * v - σ * u
    S[3] = -f * u - σ * v
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    f0    = 1.0e-4
    beta0 = 2.0e-11
    Lx = xmax - xmin
    Ly = ymax - ymin
    sponge_frac      = 0.12
    sponge_sigma_max = 5.0e-3
    sx = sponge_frac * Lx
    sy = sponge_frac * Ly

    η = q[1]
    u = q[2]
    v = q[3]

    f = f0 + beta0 * y

    dx_left  = x - xmin
    dx_right = xmax - x
    dy_bot   = y - ymin
    dy_top   = ymax - y

    rx = 0.0
    if dx_left < sx;  rx = max(rx, ((sx - dx_left)  / sx)^3); end
    if dx_right < sx; rx = max(rx, ((sx - dx_right) / sx)^3); end
    ry = 0.0
    if dy_bot < sy; ry = max(ry, ((sy - dy_bot) / sy)^3); end
    if dy_top < sy; ry = max(ry, ((sy - dy_top) / sy)^3); end

    σ = sponge_sigma_max * max(rx, ry)

    S[1] = -σ * η
    S[2] =  f * v - σ * u
    S[3] = -f * u - σ * v
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    f0    = T(1.0e-4)
    beta0 = T(2.0e-11)
    sponge_frac      = T(0.12)
    sponge_sigma_max = T(5.0e-3)

    Lx = xmax - xmin
    Ly = ymax - ymin
    sx = sponge_frac * Lx
    sy = sponge_frac * Ly

    η = q[1]; u = q[2]; v = q[3]
    f = f0 + beta0 * y

    dx_left  = x - xmin; dx_right = xmax - x
    dy_bot   = y - ymin; dy_top   = ymax - y

    rx = T(0.0)
    if dx_left  < sx; rx = max(rx, ((sx - dx_left)  / sx)^3); end
    if dx_right < sx; rx = max(rx, ((sx - dx_right) / sx)^3); end
    ry = T(0.0)
    if dy_bot < sy; ry = max(ry, ((sy - dy_bot) / sy)^3); end
    if dy_top < sy; ry = max(ry, ((sy - dy_top) / sy)^3); end

    σ = sponge_sigma_max * max(rx, ry)

    return T(-σ*η), T(f*v - σ*u), T(-f*u - σ*v)
end

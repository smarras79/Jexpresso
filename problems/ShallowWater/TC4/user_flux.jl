function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)
    #
    # Linearized rotating shallow water equations with a sloping bottom
    # (TC4: Bishnu et al. 2024, Sec. 2.8): topography b(x,y) = -α₀y, so the
    # background depth is D(y) = H₀ + α₀y (depth increases northward).
    # Folding D(y) directly into the continuity flux reproduces the
    # linearized topographic term exactly: ∂/∂y[D(y)v] = D(y)v_y + α₀v,
    # matching the "-∇·(bu)" term (Eq. 6b/10 in the paper) without a
    # separate source contribution.
    #
    #   ∂η/∂t + ∂(D(y)u)/∂x + ∂(D(y)v)/∂y = 0
    #   ∂u/∂t + ∂(gη)/∂x                   = f₀v     (source)
    #   ∂v/∂t                 + ∂(gη)/∂y    = -f₀u    (source)
    #
    # State vector: q = [η, u, v]
    #
    H0 = 1.0e3     # mean fluid depth at y=0 [m]
    g  = 9.81       # gravitational acceleration [m/s²]
    f0    = 1.0e-4    # base Coriolis parameter [s⁻¹]
    beta0 = 2.0e-11   # beta parameter [m⁻¹ s⁻¹], used only to set α₀
    α0 = beta0 * H0 / f0   # bottom slope (dimensionless)

    y = mesh.y[ip]
    D  = H0 + α0 * y

    η = q[1]
    u = q[2]
    v = q[3]

    # x-flux
    F[1] = D * u
    F[2] = g * η
    F[3] = 0.0

    # y-flux
    G[1] = D * v
    G[2] = 0.0
    G[3] = g * η
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    H0 = 1.0e3
    g  = 9.81
    f0    = 1.0e-4
    beta0 = 2.0e-11
    α0 = beta0 * H0 / f0

    y = mesh.y[ip]
    D  = H0 + α0 * y

    η = q[1]
    u = q[2]
    v = q[3]

    F[1] = D * u
    F[2] = g * η
    F[3] = 0.0

    G[1] = D * v
    G[2] = 0.0
    G[3] = g * η
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    H0 = T(1.0e3)
    g  = T(9.81)

    η = q[1]
    u = q[2]
    v = q[3]

    # NOTE: y is not available in this signature; the GPU path does not
    # currently support the y-dependent topographic flux used above.
    return T(H0*u), T(g*η), T(0.0), T(H0*v), T(0.0), T(g*η)
end

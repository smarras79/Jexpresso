function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)
    #
    # Linearized rotating shallow water equations (conservation-law form):
    #   ∂η/∂t + ∂(H₀u)/∂x + ∂(H₀v)/∂y = -σ η           (sponge in source)
    #   ∂u/∂t + ∂(gη)/∂x                = f(y)v - σ u    (source)
    #   ∂v/∂t                + ∂(gη)/∂y  = -f(y)u - σ v   (source)
    #
    # State vector: q = [η, u, v]
    #
    # Physical parameters (TC3: Bishnu et al. 2024)
    H0 = 1.0e3   # mean fluid depth [m]
    g  = 9.81     # gravitational acceleration [m/s²]

    η = q[1]
    u = q[2]
    v = q[3]

    # x-flux
    F[1] = H0 * u
    F[2] = g * η
    F[3] = 0.0

    # y-flux
    G[1] = H0 * v
    G[2] = 0.0
    G[3] = g * η
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    H0 = 1.0e3
    g  = 9.81

    η = q[1]
    u = q[2]
    v = q[3]

    F[1] = H0 * u
    F[2] = g * η
    F[3] = 0.0

    G[1] = H0 * v
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

    return T(H0*u), T(g*η), T(0.0), T(H0*v), T(0.0), T(g*η)
end

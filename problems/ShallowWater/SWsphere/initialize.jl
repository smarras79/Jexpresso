#---------------------------------------------------------------------------------
# SWsphere — initial condition.
#
# GALEWSKY BAROTROPICALLY UNSTABLE JET on the sphere
#   Galewsky, Scott & Polvani (2004), "An initial-value problem for testing
#   numerical models of the global shallow-water equations", Tellus 56A, 429-440.
#
# A zonal jet confined to a mid-latitude band, in exact geostrophic balance with
# the height field, plus a small localized height bump that seeds the barotropic
# instability. The balanced state on its own is a steady solution; the perturbed
# one rolls up into vortices over ~6 days. That is what makes it a test: the
# balanced field must stay put, and the perturbed field must produce the
# published relative-vorticity pattern.
#
# STATE VECTOR (tangent basis of the shell)
#
#   q = [h, u, v]
#
#     h  fluid depth                              [m]
#     u  ZONAL      velocity (eastward, +λ)       [m/s]
#     v  MERIDIONAL velocity (northward, +φ)      [m/s]
#
# (u, v) are components in the LOCAL TANGENT BASIS at each node,
#
#     e_λ = (-sin λ,           cos λ,          0     )
#     e_φ = (-sin φ cos λ,    -sin φ sin λ,    cos φ )
#
# so the Cartesian velocity is u·e_λ + v·e_φ. The writer emits that Cartesian
# vector alongside the components, because ParaView cannot interpret tangent
# components on its own.
#
# THE FIELDS
#
#   u(φ) = (u_max/e_n) exp[ 1/((φ-φ₀)(φ-φ₁)) ]   for φ₀ < φ < φ₁,  0 otherwise
#   v    = 0
#
#     φ₀ = π/7,  φ₁ = π/2 - φ₀,  u_max = 80 m/s,  e_n = exp[-4/(φ₁-φ₀)²]
#
#   u is C^∞: every derivative vanishes at φ₀ and φ₁, so the jet joins the
#   motionless fluid smoothly and the balanced height below is smooth too.
#
#   h(φ) = h₀ - (1/g) ∫_{-π/2}^{φ} a u(φ') [ f(φ') + tan(φ')/a · u(φ') ] dφ'
#
#   with f = 2Ω sin φ. The integrand is zero outside the jet band, so the
#   integral only ever runs over [φ₀, min(φ, φ₁)].
#
#   h₀ is fixed by requiring the GLOBAL MEAN depth to be 10 km. Because h is
#   zonally symmetric that mean is a 1-D integral, and swapping the order of
#   integration collapses the double integral to a single one:
#
#     mean(h) = ½ ∫_{-π/2}^{π/2} h(φ) cos φ dφ
#             = h₀ - (1/2g) ∫_{φ₀}^{φ₁} a u [f + tan(φ)u/a] (1 - sin φ) dφ
#
#   so h₀ is obtained in one quadrature (galewsky_h0 below) instead of being
#   hard-coded. With the published parameters it comes out as 10158.1861704546,
#   which is the constant quoted in the literature — a good check that the
#   quadrature and the parameters agree with the paper.
#
#   PERTURBATION (switch :lgalewsky_perturbation, on by default):
#
#     h'(λ,φ) = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²]
#
#     ĥ = 120 m, α = 1/3, β = 1/15, φ₂ = π/4, λ ∈ (-π, π] centred on λ = 0.
#
# The UNPERTURBED balanced state is stored in q.qe: it is the steady solution
# the scheme has to preserve, so it is the natural reference to difference
# against later.
#
# NOTE ON THE RADIUS: the balance uses the radius of the grid you actually
# built, mesh.radius. Set :sphere_radius => 6.37122e6 in user_inputs.jl (as the
# shipped deck does) so the shell is the Galewsky Earth whatever radius the
# .msh file happens to carry.
#---------------------------------------------------------------------------------

#
# Zonal jet profile. Returns 0 outside the band, which also keeps the
# 1/((φ-φ₀)(φ-φ₁)) singular at the endpoints from ever being evaluated.
#
function galewsky_ujet(φ, umax, φ0, φ1, en)
    (φ <= φ0 || φ >= φ1) && return zero(φ)
    return (umax/en) * exp(1.0/((φ - φ0)*(φ - φ1)))
end

#
# The integrand of the geostrophic balance: a u [ f + tan(φ) u / a ].
#
function galewsky_balance_integrand(φ, p)
    uu = galewsky_ujet(φ, p.umax, p.φ0, p.φ1, p.en)
    uu == 0 && return zero(φ)
    return p.a * uu * (2.0*p.Ω*sin(φ) + tan(φ)*uu/p.a)
end

#
# Composite Simpson. n must be even; the integrand is C^∞ on [φ₀, φ₁] so the
# 4th-order rule converges fast: n = 512 puts the error at ~5e-8 m on a 10 km
# depth (5e-12 relative), and the answer is unchanged to 10 digits from n = 128.
#
function galewsky_simpson(f, lo, hi, n::Int)
    hi <= lo && return 0.0
    h = (hi - lo)/n
    s = f(lo) + f(hi)
    for k = 1:n-1
        s += (isodd(k) ? 4.0 : 2.0) * f(lo + k*h)
    end
    return s*h/3.0
end

#
# h(φ) - h₀ : the balanced height relative to the constant of integration.
# Zero south of the jet, constant north of it.
#
function galewsky_hbalance(φ, p)
    hi = min(φ, p.φ1)
    hi <= p.φ0 && return 0.0
    return -galewsky_simpson(x -> galewsky_balance_integrand(x, p), p.φ0, hi, p.nquad) / p.g
end

#
# The constant of integration that puts the global mean depth at p.hmean.
# Single quadrature, from the order swap documented in the header.
#
function galewsky_h0(p)
    I = galewsky_simpson(x -> galewsky_balance_integrand(x, p)*(1.0 - sin(x)),
                         p.φ0, p.φ1, p.nquad)
    return p.hmean + I/(2.0*p.g)
end

#
# The localized height bump that triggers the instability.
#
function galewsky_hpert(λ, φ, p)
    return p.hhat * cos(φ) * exp(-(λ/p.α)^2) * exp(-((p.φ2 - φ)/p.β)^2)
end

#
# Every parameter of the test, in one place, each overridable from
# user_inputs.jl so the deck can sweep them without touching this file.
#
function galewsky_params(mesh, inputs)

    φ0   = get(inputs, :galewsky_phi0, π/7)
    φ1   = get(inputs, :galewsky_phi1, π/2 - φ0)
    umax = get(inputs, :galewsky_umax, 80.0)

    return (a     = mesh.radius,
            Ω     = get(inputs, :galewsky_Omega, 7.292e-5),
            g     = get(inputs, :galewsky_gravity, 9.80616),
            umax  = umax,
            φ0    = φ0,
            φ1    = φ1,
            en    = exp(-4.0/(φ1 - φ0)^2),
            hmean = get(inputs, :galewsky_hmean, 10000.0),
            hhat  = get(inputs, :galewsky_hhat, 120.0),
            α     = get(inputs, :galewsky_alpha, 1.0/3.0),
            β     = get(inputs, :galewsky_beta,  1.0/15.0),
            φ2    = get(inputs, :galewsky_phi2,  π/4),
            nquad = Int(get(inputs, :galewsky_nquad, 512)))
end


function initialize(SD::NSD_2D, PT, mesh::St_mesh_sphere, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if rank == 0
        println(" # ")
        println(" # INITIAL CONDITION: Galewsky et al. (2004) barotropically unstable jet ...")
    end

    #---------------------------------------------------------------------------------
    # Solution variables: q = [h, u, v]
    #---------------------------------------------------------------------------------
    qvars    = ["h", "u", "v"]
    qoutvars = ["h", "u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    p  = galewsky_params(mesh, inputs)
    h0 = galewsky_h0(p)

    lpert = get(inputs, :lgalewsky_perturbation, true)

    if rank == 0
        @printf(" #   sphere radius a = %.6e m ; Ω = %.6e 1/s ; g = %.5f m/s²\n", p.a, p.Ω, p.g)
        @printf(" #   jet: u_max = %.2f m/s in φ ∈ [%.4f, %.4f] rad = [%.2f°, %.2f°]\n",
                p.umax, p.φ0, p.φ1, p.φ0*180/π, p.φ1*180/π)
        @printf(" #   balanced height: h₀ = %.10f m for a global mean depth of %.1f m\n", h0, p.hmean)
        println(" #   perturbation: ", lpert ? "ON" : "OFF (balanced steady state only)")
        if abs(p.a - 6.37122e6)/6.37122e6 > 0.01
            @warn string(" # initialize.jl: the shell radius is ", p.a,
                         " m, not the Galewsky Earth radius 6.37122e6 m. The balance is ",
                         "computed consistently on YOUR sphere, but this is no longer the ",
                         "published test. Set :sphere_radius => 6.37122e6 in user_inputs.jl.")
        end
    end

    #---------------------------------------------------------------------------------
    # Fill q. h depends on latitude only; the perturbation breaks that symmetry.
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())

        for ip = 1:mesh.npoin

            λ = mesh.lon[ip]
            φ = mesh.lat[ip]

            hb = h0 + galewsky_hbalance(φ, p)
            uu = galewsky_ujet(φ, p.umax, p.φ0, p.φ1, p.en)

            # the balanced, unperturbed state: a steady solution, kept as the
            # reference the scheme has to preserve
            q.qe[ip, 1] = hb
            q.qe[ip, 2] = uu
            q.qe[ip, 3] = 0.0

            q.qn[ip, 1] = lpert ? hb + galewsky_hpert(λ, φ, p) : hb
            q.qn[ip, 2] = uu
            q.qn[ip, 3] = 0.0
        end
    else
        error(" # ERROR problems/ShallowWater/SWsphere/initialize.jl: only the CPU backend is implemented.")
    end

    if rank == 0
        hmin, hmax = extrema(@view q.qn[1:mesh.npoin, 1])
        uminv, umaxv = extrema(@view q.qn[1:mesh.npoin, 2])
        @printf(" #   h ∈ [%.4f, %.4f] m ; u ∈ [%.4f, %.4f] m/s ; v ≡ 0\n", hmin, hmax, uminv, umaxv)
        println(" # INITIAL CONDITION: Galewsky et al. (2004) barotropically unstable jet ... DONE")
    end

    return q
end

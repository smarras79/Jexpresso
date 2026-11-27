#-------------------------------------------------------------------------------------------------
#   ∂u/∂t = ν ∂²u/∂x²  on [0, L],    a u + b ν ∂u/∂n = g  at each end (n the outward normal)
#
# Imposed in user_bc.jl as the inward flux F = ν ∂u/∂n = (g - a u)/b.
#
# Exact solution. Each end is a Newton-cooling condition, -ν ∂u/∂n = h (u - u∞), with h = a/b
# and u∞ = g/a. With β = h/ν at the left (βL) and right (βR) ends:
#   u(x, t) = u_s(x) + A Y(x) exp(-ν λ² t)
#   u_s = c0 + c1 x      steady profile:  -ν c1 = hL (u∞L - c0),   ν c1 = hR (u∞R - c0 - c1 L)
#   Y   = cos(λx) + (βL/λ) sin(λx),  λ the first root of (λ² - βL βR) sin(λL) - λ (βL + βR) cos(λL)
#-------------------------------------------------------------------------------------------------
const R1G_NU   = Ref(0.1)
const R1G_L    = Ref(1.0)
const R1G_LEFT  = Ref((a = 0.5, b = 1.0, g = 0.5))
const R1G_RIGHT = Ref((a = 1.0, b = 4.0, g = 0.0))
const R1G_A    = 0.5
const R1G_LAM  = Ref(0.0)

r1g_h(c)    = c.a / c.b                  # equivalent transfer coefficient
r1g_uinf(c) = c.g / c.a                  # equivalent ambient value

function r1g_steady()
    ν, L = R1G_NU[], R1G_L[]
    hL, uL = r1g_h(R1G_LEFT[]),  r1g_uinf(R1G_LEFT[])
    hR, uR = r1g_h(R1G_RIGHT[]), r1g_uinf(R1G_RIGHT[])
    # [hL  -ν      ] [c0]   [hL uL]
    # [hR  hR L + ν] [c1] = [hR uR]
    M = [hL -ν; hR (hR * L + ν)]
    return M \ [hL * uL, hR * uR]
end

function r1g_lambda()
    βL = r1g_h(R1G_LEFT[]) / R1G_NU[]; βR = r1g_h(R1G_RIGHT[]) / R1G_NU[]; L = R1G_L[]
    f(λ) = (λ^2 - βL * βR) * sin(λ * L) - λ * (βL + βR) * cos(λ * L)
    a, b = 1.0e-12, π / L                # f < 0 near 0, f(π/L) = π(βL + βR)/L > 0
    for _ = 1:200
        m = (a + b) / 2
        (f(a) < 0) == (f(m) < 0) ? (a = m) : (b = m)
    end
    return (a + b) / 2
end

function r1g_exact(x, t)
    c0, c1 = r1g_steady()
    λ  = R1G_LAM[]; βL = r1g_h(R1G_LEFT[]) / R1G_NU[]
    return c0 + c1 * x + R1G_A * (cos(λ * x) + βL / λ * sin(λ * x)) * exp(-R1G_NU[] * λ^2 * t)
end

function initialize(SD::NSD_1D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    R1G_NU[]    = Float64(inputs[:μ][1])
    R1G_L[]     = Float64(inputs[:xmax] - inputs[:xmin])
    R1G_LEFT[]  = inputs[:robin_left]
    R1G_RIGHT[] = inputs[:robin_right]
    for c in (R1G_LEFT[], R1G_RIGHT[])
        (c.b != 0 && c.a / c.b > 0) || error("1d_robin: need b ≠ 0 and a/b > 0 at both ends, got $c")
    end
    R1G_LAM[] = r1g_lambda()
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs = length(qvars))
    for ip = 1:mesh.npoin
        q.qn[ip, 1] = r1g_exact(mesh.coords[1, ip] - inputs[:xmin], 0.0)
        q.qe[ip, 1] = 0.0
    end
    return q
end

#-------------------------------------------------------------------------------------------------
#   ∂u/∂t = ν ∂²u/∂x²  on [0, L]
#   left  (x = 0): Dirichlet,  u = u_left
#   right (x = L): radiation,  -ν ∂u/∂n = εσ (u⁴ - u∞⁴)     (nonlinear in u)
#
# The flux condition need not be linear: user_bc.jl returns F = εσ (u∞⁴ - u⁴) evaluated with the
# current wall value. Starting from u = u_left everywhere, the wall radiates and the solution
# relaxes to the steady profile
#   u_s(x) = u_left + c1 x,   ν c1 = εσ (u∞⁴ - (u_left + c1 L)⁴),
# a quartic for c1, solved below by bisection. u_s is linear, so the elements represent it
# exactly: at the final time the only error left is the decayed transient.
#-------------------------------------------------------------------------------------------------
const RAD_NU    = Ref(0.1)
const RAD_L     = Ref(1.0)
const RAD_COEFF = Ref(0.1)
const RAD_UINF  = Ref(1.0)
const RAD_ULEFT = Ref(2.0)

function rad_steady_slope()
    ν, L, ϵσ, u∞, uL = RAD_NU[], RAD_L[], RAD_COEFF[], RAD_UINF[], RAD_ULEFT[]
    f(c1) = ν * c1 - ϵσ * (u∞^4 - (uL + c1 * L)^4)
    a, b = (u∞ - uL) / L, 0.0           # f(a) = ν a < 0 and f(0) = ϵσ (uL⁴ - u∞⁴) > 0 for uL > u∞
    for _ = 1:200
        m = (a + b) / 2
        (f(a) < 0) == (f(m) < 0) ? (a = m) : (b = m)
    end
    return (a + b) / 2
end

rad_steady(x) = RAD_ULEFT[] + rad_steady_slope() * x

function initialize(SD::NSD_1D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    RAD_NU[]    = Float64(inputs[:μ][1])
    RAD_L[]     = Float64(inputs[:xmax] - inputs[:xmin])
    RAD_COEFF[] = Float64(inputs[:rad_coeff])
    RAD_UINF[]  = Float64(inputs[:rad_uinf])
    RAD_ULEFT[] = Float64(inputs[:u_left])
    RAD_ULEFT[] > RAD_UINF[] || error("1d_radiation: the example assumes u_left > rad_uinf")
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs = length(qvars))
    for ip = 1:mesh.npoin
        q.qn[ip, 1] = RAD_ULEFT[]
        q.qe[ip, 1] = 0.0
    end
    return q
end

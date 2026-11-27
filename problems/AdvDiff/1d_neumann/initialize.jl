#-------------------------------------------------------------------------------------------------
#   ∂u/∂t = ν ∂²u/∂x²   on [0, 1]
#   left  (x = 0): Neumann,   ν ∂u/∂n = g(t)       (a prescribed inward flux, varying in time)
#   right (x = 1): Dirichlet, u = u_R(t)
#
# Both data are taken from the exact solution
#
#   u(x, t) = 1 + A exp(-ν k² t) cos(k x + φ),
#
# which solves the equation for any k and φ, so the boundary data are whatever this u gives at
# the walls: g(t) = ν ∂u/∂n = ν (du/dx) n_x with the outward normal n_x = -1 at x = 0, and
# u_R(t) = u(1, t). k and φ are chosen so that neither boundary value is trivial.
#-------------------------------------------------------------------------------------------------
const N1_NU  = Ref(0.1)
const N1_A   = 0.5
const N1_K   = 2.5
const N1_PHI = 0.3

n1_exact(x, t) = 1 + N1_A * exp(-N1_NU[] * N1_K^2 * t) * cos(N1_K * x + N1_PHI)
n1_dudx(x, t)  = -N1_A * N1_K * exp(-N1_NU[] * N1_K^2 * t) * sin(N1_K * x + N1_PHI)

function initialize(SD::NSD_1D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    N1_NU[] = Float64(inputs[:μ][1])
    qvars = ["q"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs = length(qvars))
    for ip = 1:mesh.npoin
        q.qn[ip, 1] = n1_exact(mesh.coords[1, ip], 0.0)
        q.qe[ip, 1] = 0.0
    end
    return q
end

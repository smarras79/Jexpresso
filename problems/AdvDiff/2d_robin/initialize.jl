#-------------------------------------------------------------------------------------------------
# 2D advection-diffusion with Robin boundary conditions: a template case.
#
#   ∂u/∂t + c ∂u/∂x = ν ∇²u              on [0, Lx] x [0, H], periodic in x
#   -ν ∂u/∂n = h (u - u∞)                on the walls y = 0 ("bottom") and y = H ("top")
#
# n is the OUTWARD normal, so the Robin condition says: the diffusive flux leaving the domain is
# proportional to the difference between the wall value and an ambient value u∞ (Newton cooling,
# a heat-transfer coefficient h). h → ∞ recovers Dirichlet u = u∞, h = 0 an insulated wall.
# See user_bc.jl for how this is imposed and why it is a flux, not a value.
#
# Exact solution, used here as the initial condition (t = 0) and for checking the run:
#
#   u(x, y, t) = u_s(y) + A cos(k (x - c t)) Y(y) exp(-ν (k² + λ²) t)
#
#   u_s(y) = a + b y      the steady profile between the two ambient values:
#                         ν b = h (a - u∞_bot)  and  ν b = h (u∞_top - a - b H)
#   Y(y)   = cos(λ y) + (β/λ) sin(λ y),  β = h/ν,  satisfies the Robin condition with u∞ = 0
#                         on both walls when λ is a root of
#                         (λ² - β²) sin(λH) - 2 β λ cos(λH) = 0   (first root, in (0, π/H))
#   k      = 2π/Lx        one wavelength across the periodic direction
#
# The steady part is linear in y, so the spectral elements represent it exactly: any error in it
# is an error in the boundary flux itself (sign, h, ν or the edge integral). The transient part
# checks the convergence of the boundary treatment with resolution.
#
# Parameters (user_inputs.jl): ν = :μ[1], :robin_h, :robin_uinf_bottom, :robin_uinf_top, :adv_c,
# :robin_amp.
#-------------------------------------------------------------------------------------------------

# Case parameters, set from user_inputs in initialize() and read by user_bc.jl and user_flux.jl.
const ROBIN_NU     = Ref(0.1)    # diffusivity ν (= inputs[:μ][1])
const ROBIN_H      = Ref(0.2)    # heat-transfer coefficient h
const ROBIN_UINF_B = Ref(1.0)    # ambient value at the bottom wall
const ROBIN_UINF_T = Ref(0.0)    # ambient value at the top wall
const ROBIN_C      = Ref(1.0)    # advection speed along x
const ROBIN_AMP    = Ref(0.5)    # amplitude A of the decaying mode
const ROBIN_LX     = Ref(2.0)    # domain length in x (from the mesh)
const ROBIN_YMIN   = Ref(0.0)    # bottom wall (from the mesh)
const ROBIN_HGT    = Ref(1.0)    # domain height H (from the mesh)
const ROBIN_LAMBDA = Ref(0.0)    # eigenvalue λ (computed)

# First root of (λ² - β²) sin(λH) - 2 β λ cos(λH) in (0, π/H), by bisection
# (f < 0 near λ = 0 and f(π/H) = 2βπ/H > 0 for β > 0).
function robin_eigenvalue(β, H)
    f(λ) = (λ^2 - β^2) * sin(λ * H) - 2β * λ * cos(λ * H)
    a, b = 1.0e-12, π / H
    for _ = 1:200
        m = (a + b) / 2
        (f(a) < 0) == (f(m) < 0) ? (a = m) : (b = m)
        b - a < 1.0e-15 * b && break
    end
    return (a + b) / 2
end

# Steady profile u_s(y) = a + b (y - ymin).
function robin_steady_coeffs()
    ν, h, H = ROBIN_NU[], ROBIN_H[], ROBIN_HGT[]
    ub, ut  = ROBIN_UINF_B[], ROBIN_UINF_T[]
    b = h * (ut - ub) / (2ν + h * H)
    a = ub + ν * b / h
    return a, b
end

# Exact solution at (x, y, t).
function robin_exact(x, y, t)
    ν, c    = ROBIN_NU[], ROBIN_C[]
    λ, β    = ROBIN_LAMBDA[], ROBIN_H[] / ROBIN_NU[]
    k       = 2π / ROBIN_LX[]
    η       = y - ROBIN_YMIN[]
    a, b    = robin_steady_coeffs()
    Y       = cos(λ * η) + (β / λ) * sin(λ * η)
    return a + b * η + ROBIN_AMP[] * cos(k * (x - c * t)) * Y * exp(-ν * (k^2 + λ^2) * t)
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for 2D advection-diffusion with Robin walls ........ ")

    qvars = ("q",)
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs = length(qvars))

    # Parameters
    ROBIN_NU[]     = Float64(inputs[:μ][1])
    ROBIN_H[]      = Float64(get(inputs, :robin_h, 0.2))
    ROBIN_UINF_B[] = Float64(get(inputs, :robin_uinf_bottom, 1.0))
    ROBIN_UINF_T[] = Float64(get(inputs, :robin_uinf_top, 0.0))
    ROBIN_C[]      = Float64(get(inputs, :adv_c, 1.0))
    ROBIN_AMP[]    = Float64(get(inputs, :robin_amp, 0.5))
    ROBIN_H[] > 0 || error("2d_robin: :robin_h must be > 0 (h = 0 is a Neumann wall, not a Robin one)")

    # Domain extent from the mesh (global over MPI ranks)
    comm  = MPI.COMM_WORLD
    xmin  = MPI.Allreduce(minimum(@view(mesh.coords[1, :])), MPI.MIN, comm)
    xmax  = MPI.Allreduce(maximum(@view(mesh.coords[1, :])), MPI.MAX, comm)
    ymin  = MPI.Allreduce(minimum(@view(mesh.coords[2, :])), MPI.MIN, comm)
    ymax  = MPI.Allreduce(maximum(@view(mesh.coords[2, :])), MPI.MAX, comm)
    ROBIN_LX[]   = xmax - xmin
    ROBIN_YMIN[] = ymin
    ROBIN_HGT[]  = ymax - ymin
    ROBIN_LAMBDA[] = robin_eigenvalue(ROBIN_H[] / ROBIN_NU[], ROBIN_HGT[])

    a, b = robin_steady_coeffs()
    if MPI.Comm_rank(comm) == 0
        @info "2d_robin: ν = $(ROBIN_NU[]), h = $(ROBIN_H[]), u∞ (bottom, top) = ($(ROBIN_UINF_B[]), $(ROBIN_UINF_T[])), " *
              "c = $(ROBIN_C[]); λ = $(ROBIN_LAMBDA[]), steady u_s = $(a) + $(b) (y - ymin), " *
              "decay rate ν(k² + λ²) = $(ROBIN_NU[] * ((2π / ROBIN_LX[])^2 + ROBIN_LAMBDA[]^2)) 1/s"
    end

    for ip = 1:mesh.npoin
        q.qn[ip, 1] = robin_exact(mesh.coords[1, ip], mesh.coords[2, ip], 0.0)
        q.qe[ip, 1] = 0.0
    end

    println(" Initialize fields for 2D advection-diffusion with Robin walls ........ DONE")
    return q
end

#---------------------------------------------------------------------------------
# TC6: Non-linear manufactured solution (Bishnu et al. 2024, Sec. 2.10)
#
# Non-linear rotating shallow water equations, conservative form (paper Eq. 41,
# flat bottom), doubly periodic, f-plane:
#
#   ∂h/∂t  + ∂(hu)/∂x              + ∂(hv)/∂y               = s_h
#   ∂hu/∂t + ∂(hu² + g h²/2)/∂x    + ∂(huv)/∂y      - f hv  = s_hu
#   ∂hv/∂t + ∂(huv)/∂x             + ∂(hv² + g h²/2)/∂y + f hu = s_hv
#
# State q = [h, hu, hv] with h = H + η. Exact solution, phase φ = kx x + ky y - ω t:
#
#   u = η̂ cos φ,   v = η̂ cos φ,   η = η̂ sin φ
#
# NOTE on v. The paper prints v = 0 (Eq. 38b), but its source terms (Eqs. 40 and
# 42) are those of v = η̂ cos φ = u: s_h contains (kx + ky), which needs
# ∂(hv)/∂y ≠ 0, and s_v contains ω sin φ, which needs ∂v/∂t ≠ 0. Eq. 38b is a
# typo; v = u is the solution that is consistent with all of the printed
# sources, and Eq. 42 was re-derived symbolically for it (exact match).
#
# Parameters (paper Sec. 2.10): η̂ = 0.01 m, H = 1000 m, g = 10 m/s² (c = √(gH)
# = 100 m/s), f = 1e-4 s⁻¹, Lx = Ly = 1e7 m, kx = 2π/Lx, ky = 2π/Ly,
# ω = c √(kx² + ky²).
#---------------------------------------------------------------------------------

const TC6_H      = 1.0e3
const TC6_G      = 10.0
const TC6_F      = 1.0e-4
const TC6_ETAHAT = 0.01
const TC6_L      = 1.0e7                       # Lx = Ly (mesh SWE_TC2_periodic_5x5)
const TC6_KX     = 2π / TC6_L
const TC6_KY     = 2π / TC6_L
const TC6_OMEGA  = sqrt(TC6_G * TC6_H) * sqrt(TC6_KX^2 + TC6_KY^2)

# Exact solution (h, hu, hv) at (x, y, t).
@inline function tc6_exact(x, y, t)
    φ  = TC6_KX * x + TC6_KY * y - TC6_OMEGA * t
    h  = TC6_H + TC6_ETAHAT * sin(φ)
    u  = TC6_ETAHAT * cos(φ)
    return h, h * u, h * u        # v = u
end

# Manufactured sources (paper Eq. 42, valid for u = v = η̂ cos φ).
@inline function tc6_sources(x, y, t)
    A, H, g, f, ω = TC6_ETAHAT, TC6_H, TC6_G, TC6_F, TC6_OMEGA
    kx, ky = TC6_KX, TC6_KY
    K  = kx + ky
    φ  = kx * x + ky * y - ω * t
    s, c   = sincos(φ)
    s2, c2 = sincos(2φ)
    cubic  = A^3 * K * (c^3 - s * s2)
    s_hu = A * H * ((g * kx - f) * c + ω * s) + A^2 * (s2 * (-H * K + (g * kx - f) / 2) - ω * c2) + cubic
    s_hv = A * H * ((g * ky + f) * c + ω * s) + A^2 * (s2 * (-H * K + (g * ky + f) / 2) - ω * c2) + cubic
    s_h  = A * (-H * K * s - ω * c + A * K * c2)
    return s_h, s_hu, s_hv
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        println(" Initialize fields for 2D Shallow Water TC6 (non-linear manufactured solution) ...")
    end

    qvars    = ["h", "hu", "hv"]
    qoutvars = ["η", "u", "v"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        Lx = MPI.Allreduce(maximum(@view(mesh.coords[1,:])), MPI.MAX, comm) - MPI.Allreduce(minimum(@view(mesh.coords[1,:])), MPI.MIN, comm)
        Ly = MPI.Allreduce(maximum(@view(mesh.coords[2,:])), MPI.MAX, comm) - MPI.Allreduce(minimum(@view(mesh.coords[2,:])), MPI.MIN, comm)
        if rank == 0
            println("  Domain size: Lx = $Lx m, Ly = $Ly m (TC6 assumes $(TC6_L) m)")
            println("  ω = $(TC6_OMEGA) s⁻¹, period = $(2π/TC6_OMEGA/3600) h")
        end
        (isapprox(Lx, TC6_L; rtol=1e-8) && isapprox(Ly, TC6_L; rtol=1e-8)) ||
            error(" # ERROR TC6: the mesh is $(Lx) x $(Ly) m but TC6_L = $(TC6_L) m; the manufactured solution would not be periodic.")

        t0 = inputs[:tinit]
        for ip = 1:mesh.npoin
            h, hu, hv = tc6_exact(mesh.coords[1, ip], mesh.coords[2, ip], t0)
            q.qn[ip, 1] = h
            q.qn[ip, 2] = hu
            q.qn[ip, 3] = hv
            q.qe[ip, 1] = h
            q.qe[ip, 2] = hu
            q.qe[ip, 3] = hv
        end
    end

    if rank == 0
        println(" Initialize fields for 2D Shallow Water TC6 (non-linear manufactured solution) ... DONE")
    end
    return q
end

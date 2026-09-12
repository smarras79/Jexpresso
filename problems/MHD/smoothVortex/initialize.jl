#---------------------------------------------------------------------------------
# The SMOOTH (isentropic) MHD VORTEX — the accuracy test of
# Dao & Nazarov, J. Sci. Comput. 92:77 (2022), §5.1.1, in the form of
# Balsara, ApJS 151:149 (2004): a steady vortex superimposed on a uniform
# advection, so the exact solution at any time is the initial condition
# translated by v₀t on the doubly periodic domain. It is the one problem of
# that paper with a closed-form solution, which is why it is the one the
# convergence history of their Fig. 1 is measured on.
#
# Domain [-5, 5]², doubly periodic, γ = 5/3, v₀ = (1, 1), so the vortex
# returns to its starting point at t = 10.
#
#   r² = (x - x_c)² + (y - y_c)²,        f(r) = exp((1 - r²)/2)
#
#   ρ  = 1
#   v  = v₀ + (κ/2π) f (-(y - y_c), (x - x_c))
#   B  =      (μ/2π) f (-(y - y_c), (x - x_c))
#   p  = 1 + (1/8π²) ( μ²(1 - r²) - κ² ) exp(1 - r²)
#
# with κ = μ = 1 the vortex and magnetic strengths (`SV_KAPPA`, `SV_MU`).
#
# The pressure is the one that makes the vortex an EXACT steady solution in
# the co-moving frame, which is what the accuracy test needs. With
# v_θ = κ g r and B_θ = μ g r, g = f/2π, radial momentum balance reads
#
#   ρ v_θ²/r - B_θ²/r = d/dr ( p + B_θ²/2 ),
#
# whose left side is (κ² - μ²) g² r and whose magnetic term contributes
# μ² (1/4π²) e^{1-r²} r (1 - r²); the p above closes it identically:
#
#   dp/dr = (1/4π²) e^{1-r²} r [ κ² - μ²(2 - r²) ].
#
# ψ and the out-of-plane components w, B_z are zero and stay zero.
#---------------------------------------------------------------------------------
const SV_KAPPA = 1.0     # vortex strength
const SV_MU    = 1.0     # magnetic strength
const SV_XC    = 0.0     # vortex centre at t = 0
const SV_YC    = 0.0
const SV_U0    = 1.0     # uniform advection
const SV_V0    = 1.0
const SV_RHO0  = 1.0
const SV_P0    = 1.0

# The state of the exact solution at one point, as the conserved vector of
# this case. `xr, yr` are the coordinates relative to the vortex centre —
# the caller wraps them into [-L/2, L/2) so the exact solution is periodic.
@inline function sv_state(xr, yr, γ)
    r2 = xr*xr + yr*yr
    f  = exp(0.5*(1.0 - r2))
    g  = f/(2.0*π)

    ρ  = SV_RHO0
    u  = SV_U0 - SV_KAPPA*g*yr
    v  = SV_V0 + SV_KAPPA*g*xr
    w  = 0.0
    Bx = -SV_MU*g*yr
    By =  SV_MU*g*xr
    Bz = 0.0
    ψ  = 0.0
    p  = SV_P0 + (SV_MU*SV_MU*(1.0 - r2) - SV_KAPPA*SV_KAPPA)*exp(1.0 - r2)/(8.0*π*π)

    ρE = p/(γ - 1.0) + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ
    return (ρ, ρ*u, ρ*v, ρE, ρ*w, Bx, By, Bz, ψ, p)
end

# Wrap a displacement into [-L/2, L/2): the exact solution is periodic, so
# the vortex that leaves one side comes back on the other.
@inline sv_wrap(d, L) = d - L*round(d/L)

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (smooth vortex) ........... "

    #---------------------------------------------------------------------------------
    # Slot 4 MUST carry the total energy ρE — see the header of user_flux.jl.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    inputs[:backend] == CPU() ||
        error(" problems/MHD/smoothVortex: only the CPU backend is supported for now.")
    inputs[:SOL_VARS_TYPE] == TOTAL() ||
        error(" problems/MHD/smoothVortex: only SOL_VARS_TYPE = TOTAL() is supported.")

    γ  = γ_mhd
    Lx = mesh.xmax - mesh.xmin
    Ly = mesh.ymax - mesh.ymin

    ch_local = 0.0
    for ip = 1:mesh.npoin
        xr = sv_wrap(mesh.x[ip] - SV_XC, Lx)
        yr = sv_wrap(mesh.y[ip] - SV_YC, Ly)
        s  = sv_state(xr, yr, γ)
        for ieq = 1:length(qvars)
            q.qn[ip,ieq] = s[ieq]
            q.qe[ip,ieq] = s[ieq]      # background state, for perturbation output only
        end
        q.qn[ip,end] = s[10]
        q.qe[ip,end] = s[10]

        ρ  = s[1]
        u  = s[2]/ρ; v = s[3]/ρ; w = s[5]/ρ
        B2 = s[6]^2 + s[7]^2 + s[8]^2
        cf = sqrt(γ*s[10]/ρ + B2/ρ)
        ch_local = max(ch_local, sqrt(u*u + v*v + w*w) + cf)
    end

    # GLM divergence-cleaning speed: the maximum wave speed of the initial
    # condition over the whole (global) domain, constant in time.
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    if rank == 0
        # What box this run is on, and what it costs. The vortex is a
        # Gaussian, so on [-L/2, L/2]² the exact solution is NOT periodic:
        # its velocity perturbation at the middle of an edge is
        # (SV_KAPPA/2π)(L/2)exp((1 − (L/2)²)/2), and it has the opposite sign on
        # the opposite edge, so the initial condition jumps across the
        # periodic seam by twice that. That jump is a discontinuity in the
        # DATA, so it floors the error of every scheme, at every order, at
        # about its own size — 4.9e-06 on the L = 10 box of the published
        # figure, where a 6th-order element reaches the floor sooner than a
        # 4th and the convergence history flattens. Say so, with the number,
        # rather than leave it to be rediscovered from a flat curve.
        tail = (SV_KAPPA/(2.0*π))*(0.5*Lx)*exp(0.5*(1.0 - (0.5*Lx)^2))
        hint = 2.0*tail > 1.0e-12 ?
               "\n   which is the floor of any error study on this mesh, at every order (widen the box with JEXPRESSO_SV_L=20)" :
               "\n   which is below round-off, so the error study is limited by the scheme and not by the box"
        @info string(@sprintf(" box [%.3f, %.3f] x [%.3f, %.3f]: the periodic-seam jump of the exact solution is %.2e",
                              mesh.xmin, mesh.xmax, mesh.ymin, mesh.ymax, 2.0*tail), hint)
        if abs(Lx - Ly) > 1.0e-8
            @warn " problems/MHD/smoothVortex: the box is not square; the vortex is written for a square one."
        end
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (smooth vortex) ........... DONE"
    end

    return q
end

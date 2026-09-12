#---------------------------------------------------------------------------------
# THE CLASSICAL ISENTROPIC (SHU) VORTEX for the 2D compressible Euler
# equations — the standard smooth accuracy test, and here the HYDRODYNAMIC
# CONTROL for problems/MHD/smoothVortex.
#
# Same box, same meshes, same machinery, same measured quantity; the only
# difference is that there is no magnetic field, so there is no ∇·B and no
# GLM cleaning field ψ. Whatever error the MHD case shows beyond this one is
# the price of the divergence constraint, not of the discretization.
#
# Nondimensional, as the test is always written: ρ∞ = p∞ = T∞ = 1, R = 1,
# γ = 1.4 (PhysConst.γ). Domain [-5,5]², doubly periodic, v₀ = (1,1), so the
# vortex returns to its starting point at t = 10 and the exact solution at
# any time is the initial condition translated by v₀t.
#
#   r² = (x - x_c)² + (y - y_c)²
#
#   δu = -(β/2π) (y - y_c) exp((1 - r²)/2)
#   δv =  (β/2π) (x - x_c) exp((1 - r²)/2)
#   δT = -(γ-1)β²/(8γπ²) exp(1 - r²)
#
#   T  = 1 + δT,    ρ = T^{1/(γ-1)},    p = ρ^γ = T^{γ/(γ-1)}
#
# The state is isentropic (p/ρ^γ ≡ 1) and steady in the co-moving frame: the
# centrifugal force of the swirl is balanced by the pressure gradient the
# temperature dip produces, exactly, for any β.
#
# β = 5 is the strength of the classical test (Shu, ICASE 97-65; the density
# dips to 0.66 at the core). β = 1 makes the perturbation comparable to the
# MHD vortex of Dao & Nazarov (κ = μ = 1), which is the setting to use for a
# like-for-like comparison with the MHD case — JEXPRESSO_EV_BETA sets it.
#---------------------------------------------------------------------------------
const EV_BETA_DEFAULT = 5.0     # vortex strength (JEXPRESSO_EV_BETA)
const EV_XC = 0.0               # vortex centre at t = 0
const EV_YC = 0.0
const EV_U0 = 1.0               # uniform advection
const EV_V0 = 1.0

# The state of the exact solution at one point, as the conserved vector of
# this case plus the pressure. `xr, yr` are the coordinates relative to the
# vortex centre — the caller wraps them into [-L/2, L/2) so that the exact
# solution is periodic.
@inline function ev_state(xr, yr, γ, β)
    r2 = xr*xr + yr*yr
    f  = exp(0.5*(1.0 - r2))

    u  = EV_U0 - (β/(2.0*π))*yr*f
    v  = EV_V0 + (β/(2.0*π))*xr*f

    δT = -(γ - 1.0)*β*β*exp(1.0 - r2)/(8.0*γ*π*π)
    T  = 1.0 + δT
    ρ  = T^(1.0/(γ - 1.0))
    p  = ρ^γ

    ρE = p/(γ - 1.0) + 0.5*ρ*(u*u + v*v)
    return (ρ, ρ*u, ρ*v, ρE, p)
end

# Wrap a displacement into [-L/2, L/2): the exact solution is periodic, so
# the vortex that leaves one side comes back on the other.
@inline ev_wrap(d, L) = d - L*round(d/L)

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D CompEuler (isentropic vortex) ........... "

    qvars    = ["ρ", "ρu", "ρv", "ρE"]
    qoutvars = ["ρ", "u", "v", "p", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    inputs[:backend] == CPU() ||
        error(" problems/CompEuler/smoothVortex: only the CPU backend is supported for now.")
    inputs[:SOL_VARS_TYPE] == TOTAL() ||
        error(" problems/CompEuler/smoothVortex: only SOL_VARS_TYPE = TOTAL() is supported.")

    γ  = PhysicalConst{Float64}().γ
    β  = _ev_beta()
    Lx = mesh.xmax - mesh.xmin
    Ly = mesh.ymax - mesh.ymin

    for ip = 1:mesh.npoin
        xr = ev_wrap(mesh.x[ip] - EV_XC, Lx)
        yr = ev_wrap(mesh.y[ip] - EV_YC, Ly)
        s  = ev_state(xr, yr, γ, β)
        for ieq = 1:length(qvars)
            q.qn[ip,ieq] = s[ieq]
            q.qe[ip,ieq] = s[ieq]      # background state, for perturbation output only
        end
        q.qn[ip,end] = s[5]            # pressure slot of uaux
        q.qe[ip,end] = s[5]
    end

    if rank == 0
        # What box this run is on, and what it costs. The vortex is a
        # Gaussian, so on [-L/2, L/2]² the exact solution is NOT periodic:
        # its velocity perturbation at the middle of an edge is
        # (β/2π)(L/2)exp((1 − (L/2)²)/2), and it has the opposite sign on
        # the opposite edge, so the initial condition jumps across the
        # periodic seam by twice that. That jump is a discontinuity in the
        # DATA, so it floors the error of every scheme, at every order, at
        # about its own size — 4.9e-06 on the L = 10 box of the published
        # figure, where a 6th-order element reaches the floor sooner than a
        # 4th and the convergence history flattens. Say so, with the number,
        # rather than leave it to be rediscovered from a flat curve.
        tail = (β/(2.0*π))*(0.5*Lx)*exp(0.5*(1.0 - (0.5*Lx)^2))
        hint = 2.0*tail > 1.0e-12 ?
               "\n   which is the floor of any error study on this mesh, at every order (widen the box with JEXPRESSO_EV_L=20)" :
               "\n   which is below round-off, so the error study is limited by the scheme and not by the box"
        @info string(@sprintf(" box [%.3f, %.3f] x [%.3f, %.3f]: the periodic-seam jump of the exact solution is %.2e",
                              mesh.xmin, mesh.xmax, mesh.ymin, mesh.ymax, 2.0*tail), hint)
        if abs(Lx - Ly) > 1.0e-8
            @warn " problems/CompEuler/smoothVortex: the box is not square; the vortex is written for a square one."
        end
        ρmin = minimum(@view q.qn[1:mesh.npoin,1])
        @info @sprintf(" isentropic vortex: β = %.3f, min ρ = %.6f, advection (%.1f, %.1f)",
                       β, ρmin, EV_U0, EV_V0)
        @info " Initialize fields for 2D CompEuler (isentropic vortex) ........... DONE"
    end

    return q
end

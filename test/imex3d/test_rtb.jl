#=============================================================================
 test_rtb.jl -- the fully 3D IMEX on a rising thermal bubble.

     julia --project=<env> test/imex3d/test_rtb.jl
     mpiexecjl -n 3 julia --project=<env> test/imex3d/test_rtb.jl

 The same benchmark, the same fixture and the SAME MESH as
 test/hevi/test_rtb.jl -- deliberately, because the headline claim of this
 scheme is relative to HEVI and a claim measured on a different mesh is not a
 comparison. Both schemes are built here and both are run, so the three step
 sizes come out of one process on one grid.

 WHAT THE MESH DOES TO THE ARGUMENT
 ----------------------------------
 8 x 1 x 24 elements over a 1000 m cube: Δz is three times finer than Δx, and
 the smallest LGL gaps are 21.6 m and 7.2 m. HEVI's gain is bounded by exactly
 that ratio, because the vertical acoustic term is the only one it removes.
 IMEX3D removes the horizontal ones as well, so its limit is set by advection
 -- which in a 0.5 K bubble is two orders of magnitude slower than sound. The
 measured step sizes below are the whole point of the file.
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "..", "hevi", "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl", "acoustic.jl",
          "ark.jl", "hevi.jl", "krylov.jl", "imex3d.jl")
    include(joinpath(SRC, f))
end
include(joinpath(@__DIR__, "..", "hevi", "rtb_fixture.jl"))

say(args...) = RANK == 0 && (println(args...); flush(stdout))

const NELX, NELY, NELZ, P = 8, 1, 24, 4
const L = 1000.0
const Θ0, ΘC, R0 = 300.0, 0.5, 250.0

say("\n=== IMEX3D rising thermal bubble: $NR rank(s), $(NELX)x$(NELY)x$(NELZ) elements, p=$P ===")

params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lx = L, Ly = L, Lz = L, comm = COMM,
                           θ0 = Θ0, base = :neutral)
npoin = Int(params.mesh.npoin)

#--- setup, in the order build_imex3d does it inside Jexpresso ---------------
topo = build_column_topology(params.mesh, COMM)

# Lateral free-slip walls. The fixture imposes free slip on every face of the
# box by projecting the normal momentum out of the STATE, so the exact
# linearisation of its RHS is (acoustic operator) ∘ (that projection) -- which
# is what the wallx/wally flags apply. Without them the implicit operator lets
# sound run out through the side walls and the difference lands in f_exp as a
# stiff residual at the boundary nodes.
wallx, wally, wsrc = imex_lateral_walls(params, :box)
@assert wsrc === :box

opfull = build_hevi_fast_operator(params, topo; wallx = wallx, wally = wally)
# The preconditioner drops the rows that are exactly the identity here.
PVARS  = hevi_choose_vars(params.metrics, COMM)
opvert = build_hevi_operator(params, topo, PVARS; full = false)
owner, own = assign_column_owners(topo, COMM)
cc   = build_column_comm(topo, owner, own, COMM, length(PVARS))

# HEVI's own operator, for the side-by-side comparison further down.
hvars = hevi_choose_vars(params.metrics, COMM)
oph   = build_hevi_operator(params, topo, hvars)
cch   = build_column_comm(topo, owner, own, COMM, length(hvars))

# Step sizes. The two small ones are the numbers test/hevi/test_rtb.jl measured
# on this mesh: the explicit limit sits between 0.02 and 0.04, HEVI's between
# 0.06 and 0.08. DT_IMEX is a further 5x beyond HEVI's, which is where this
# scheme should be -- limited by advection, not by sound in any direction.
const DT_SMALL = 0.02      # every scheme is stable here
const DT_HEVI  = 0.05      # HEVI stable, explicit not
const DT_IMEX  = 0.40      # 20x the explicit limit, 8x HEVI's
const GAMMA    = ark_tableau(:ARS343).γ

facv = build_column_factorization(params, opvert, cc, topo, GAMMA * DT_SMALL)
pc   = IMEX3DPrecond(opvert, cc, facv, topo, PVARS)
inner = build_distributed_inner(opfull.dss_cache, npoin, COMM)
# rtol 1e-10 rather than the production 1e-8: this file measures the SCHEME,
# and a stage residual that is not comfortably below the step's own truncation
# error would show up in the order and step-size comparisons as if it were the
# integrator's doing.
ws   = GMRESWorkspace(npoin, 5, inner; m = 30, maxiter = 300, rtol = 1e-10)
# warm_start = false: check 4 below measures the iteration-count LAW
# (iterations ~ 25 x CFL_h), which is a property of the preconditioner and the
# operator's spectrum. A warm start is a property of the SEQUENCE of
# right-hand sides and would confound it.
imex = IMEX3DCache(topo, opfull, pc, nothing, ws, zeros(npoin, 5), zeros(npoin, 5),
                   imex3d_fimp!, imex3d_solve!, :ARS343, GAMMA * DT_SMALL,
                   :RS, 5, false, COMM)

fach = build_column_factorization(params, oph, cch, topo, GAMMA * DT_HEVI)
hevi = HEVICache(topo, oph, cch, fach, hevi_fimp!, hevi_column_solve!, :ARS343,
                 GAMMA * DT_HEVI, :RS, 5)

ph = merge(params, (rtb = RTBCache(params), imex = imex, hevi = hevi))

const HX = 0.1727 * L/NELX
const HZ = 0.1727 * L/NELZ
say(@sprintf("  mesh: %d points/rank, %d levels, %d columns; smallest gaps dx = %.1f m, dz = %.1f m (anisotropy %.2f)",
             npoin, topo.nlev, topo.ncol, HX, HZ, HX/HZ))
say(@sprintf("  lateral free-slip wall nodes: %d in x, %d in y", count(wallx), count(wally)))

@testset "IMEX3D rising thermal bubble" begin

#--- 1. the split does not disturb hydrostatic balance -----------------------
# f_imp acts on (u - qe) and the bubble-free state IS qe, so it vanishes there
# identically -- exactly, not nearly. Whatever discrete imbalance the fixture's
# strong-form discretisation carries stays entirely inside f_exp, so IMEX3D,
# HEVI and the explicit scheme all inherit the same one and it cancels out of
# every comparison below.
ubase = rtb_initial_state(params; θc = 0.0, θ0 = Θ0)
f0 = similar(ubase)
imex3d_fimp!(f0, copy(ubase), ph, 0.0)
fmax = MPI.Allreduce(maximum(abs, f0), MPI.MAX, COMM)
@test fmax == 0.0
say(@sprintf("  f_imp on the bubble-free base state: %.1e  (must be exactly 0)", fmax))

#--- 2. the bubble ------------------------------------------------------------
u0 = rtb_initial_state(params; θc = ΘC, r0 = R0, θ0 = Θ0)
d0 = rtb_diagnostics(u0, params; θ0 = Θ0)
@test isapprox(d0.θmax, ΘC; rtol = 2e-3)
say(@sprintf("  t=0: theta' max = %.4f K, centre of buoyancy z = %.1f m", d0.θmax, d0.zbar))

#--- 3. IMEX3D agrees with an explicit reference at a shared step ------------
# Compare on w: it starts at exactly zero and is what buoyancy drives, so it
# carries no large background to hide a discrepancy in.
TCMP = 2.0
prob = ODEProblem(rtb_rhs!, copy(u0), (0.0, TCMP), ph)
ue = solve(prob, Tsit5(); dt = DT_SMALL, adaptive = false, save_everystep = false).u[end]
ui = solve(prob, IMEX_ARK(:ARS343); dt = DT_SMALL, adaptive = false,
           save_everystep = false).u[end]
it_small = ws.niter / max(ws.nsolve, 1)

wof(u, ip) = u[3*npoin + ip] / u[0*npoin + ip]
function wgap(ua, ub)
    d = 0.0; s = 0.0
    for ip = 1:npoin
        d = max(d, abs(wof(ua, ip) - wof(ub, ip))); s = max(s, abs(wof(ub, ip)))
    end
    MPI.Allreduce(d, MPI.MAX, COMM) / MPI.Allreduce(s, MPI.MAX, COMM)
end

gap_scheme = wgap(ui, ue)
@test gap_scheme < 5e-3
say(@sprintf("  after %.0f s at dt=%.3g: |w_IMEX - w_explicit|/|w| = %.2e  (%.1f Krylov iterations/solve)",
             TCMP, DT_SMALL, gap_scheme, it_small))

#--- 4. and with itself at a step 20x larger ---------------------------------
ws.nsolve = 0; ws.niter = 0
uib = solve(prob, IMEX_ARK(:ARS343); dt = DT_IMEX, adaptive = false,
            save_everystep = false).u[end]
it_big = ws.niter / max(ws.nsolve, 1)
gap_step = wgap(uib, ui)
@test gap_step < 1e-1
say(@sprintf("  IMEX3D dt=%.3g vs dt=%.3g: |Δw|/|w| = %.2e  (%.1f Krylov iterations/solve)",
             DT_IMEX, DT_SMALL, gap_step, it_big))
# THE ITERATION COUNT IS THE ENTIRE COST OF THIS SCHEME, and it is not a
# constant. Measured on this mesh, over Δt spanning 8x, it is very nearly
#
#     iterations  ~  25 x CFL_h,     CFL_h = γΔt·c/h_x
#
# the HORIZONTAL acoustic Courant number that survives the vertical
# preconditioner:
#
#     Δt     γΔt      CFL_h   iterations (rtol 1e-8)
#     0.05   0.0218   0.35    12
#     0.10   0.0436   0.70    19
#     0.20   0.0872   1.40    35
#     0.40   0.1743   2.80    71
#
# Linear, because the operator is skew: the eigenvalues of I − γΔt A lie on a
# segment of a vertical line whose LENGTH is CFL_h, and GMRES on such a
# spectrum needs iterations proportional to that length. The same law comes out
# on the rtb_imex production mesh, whose h_x is 8x larger.
#
# The consequence is the thing to remember: cost per step grows linearly with
# Δt, so cost per unit SIMULATED time approaches a floor rather than falling.
# Raising Δt past the point where the 4 rhs! evaluations stop dominating buys
# nothing. See README_IMEX3D.md.
#
# Asserted with margin around the measured value, so both a regression in the
# preconditioner and an unexpected improvement are visible.
@test 50 < it_big < 95
@test ws.nunconverged == 0

#--- 5. the bubble physics at the big step -----------------------------------
TRUN = 60.0
ui2 = solve(ODEProblem(rtb_rhs!, copy(u0), (0.0, TRUN), ph), IMEX_ARK(:ARS343);
            dt = DT_IMEX, adaptive = false, save_everystep = false).u[end]
d = rtb_diagnostics(ui2, params; θ0 = Θ0)

@test all(isfinite, ui2)
@test d.zbar > d0.zbar + 5.0                 # the bubble went UP
@test 0.05 < d.wmax < 10.0                   # and at a physical speed
@test d.θmax < 1.5 * ΘC                      # without the peak blowing up
@test d.θmax > 0.5 * ΘC                      # or being diffused away
@test d.vmax < 1e-8 * max(d.wmax, 1.0)       # and stayed y-invariant
say(@sprintf("  IMEX3D at dt=%.3g for %.0f s: theta' max = %.4f K, z = %.1f m (+%.1f), w max = %.3f m/s, |v| max = %.1e",
             DT_IMEX, TRUN, d.θmax, d.zbar, d.zbar - d0.zbar, d.wmax, d.vmax))

#--- 6. neither of the other two schemes survives that step ------------------
# HEVI is included on purpose. It is not a straw man -- it is the scheme this
# one has to beat, it is stable at 0.05 on this very mesh (test/hevi/test_rtb.jl
# measures that), and the vertical acoustic term is the only one it removes.
# At 0.40 the HORIZONTAL acoustic CFL is what kills it, which is precisely the
# term IMEX3D takes implicitly and HEVI does not.
NBLOW = 150
function blows_up(alg, dt, nstep)
    uend = with_logger(NullLogger()) do
        solve(ODEProblem(rtb_rhs!, copy(u0), (0.0, nstep*dt), ph), alg;
              dt = dt, adaptive = false, save_everystep = false,
              unstable_check = (dt, u, p, t) -> false).u[end]
    end
    # BOTH reductions run UNCONDITIONALLY and are combined only afterwards.
    # The obvious `!all(isfinite, u) || MPI.Allreduce(...) > 1e3` deadlocks:
    # `||` short-circuits, so the first rank whose solution has already gone
    # non-finite skips the Allreduce every other rank is sitting in. It passes
    # in serial and hangs on 3 ranks after the assertions have all gone green.
    g      = maximum(abs, uend .- u0)
    nonfin = MPI.Allreduce((all(isfinite, uend) && isfinite(g)) ? 0.0 : 1.0, MPI.MAX, COMM)
    growth = MPI.Allreduce(isfinite(g) ? g : 0.0, MPI.MAX, COMM)
    return nonfin > 0.0 || growth > 1e3
end

refactorize!(fach, params, oph, cch, topo, GAMMA * DT_IMEX)
@test blows_up(HEVI_ARK(:ARS343), DT_IMEX, NBLOW)
say(@sprintf("  HEVI at dt=%.3g: diverged within %d steps -- the HORIZONTAL acoustic term,",
             DT_IMEX, NBLOW))
say( "     which HEVI leaves explicit and this scheme does not")
@test blows_up(Tsit5(), DT_IMEX, NBLOW)
say(@sprintf("  explicit at dt=%.3g: diverged too, as expected", DT_IMEX))

say(@sprintf("  -> IMEX3D ran at %.0fx the explicit step and %.0fx HEVI's, on a mesh with",
             DT_IMEX/DT_SMALL, DT_IMEX/DT_HEVI))
say(@sprintf("     %.1fx vertical anisotropy -- which bounds what HEVI can buy and not what this can.",
             HX/HZ))

end # testset

MPI.Barrier(COMM)
say("\n=== IMEX3D rising thermal bubble done ===")

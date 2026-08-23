#=============================================================================
 test_rtb.jl -- HEVI on a rising thermal bubble.

     julia --project=<env> test/hevi/test_rtb.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_rtb.jl

 Needs only MPI, LinearAlgebra, Printf and OrdinaryDiffEq -- see
 rtb_fixture.jl for what the fixture provides and why this benchmark is the
 right one for a vertically-implicit split.

 The mesh is deliberately anisotropic: 8 x 1 x 24 elements over a 1000 m cube,
 so Δz is four times finer than Δx. That is the situation HEVI exists for, in
 miniature, and it makes the step-size claim testable rather than asserted.
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
for f in ("columns.jl", "operator.jl", "factorize.jl", "ark.jl", "hevi.jl")
    include(joinpath(SRC, f))
end
include(joinpath(@__DIR__, "rtb_fixture.jl"))

say(args...) = RANK == 0 && (println(args...); flush(stdout))

const NELX, NELY, NELZ, P = 8, 1, 24, 4
const L = 1000.0
const Θ0, ΘC, R0 = 300.0, 0.5, 250.0

say("\n=== HEVI rising thermal bubble: $NR rank(s), $(NELX)x$(NELY)x$(NELZ) elements, p=$P ===")

params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lx = L, Ly = L, Lz = L, comm = COMM,
                           θ0 = Θ0, base = :neutral)

# --- HEVI setup, the same sequence build_hevi does inside Jexpresso ---------
topo  = build_column_topology(params.mesh, COMM)
vars  = hevi_choose_vars(params.metrics, COMM)
op    = build_hevi_operator(params, topo, vars)
owner, own = assign_column_owners(topo, COMM)
cc    = build_column_comm(topo, owner, own, COMM, length(vars))

# Measured on this mesh, not estimated: a sweep of both schemes to t = 30 s
# puts the explicit limit between 0.02 and 0.04 and the HEVI limit between 0.06
# and 0.08. DT_BIG is set below the measured HEVI limit with margin, and above
# the measured explicit one.
#
# A first-principles estimate from the LGL gaps put both limits ~1.8x higher.
# That factor is the same for both schemes -- the spectral radius of the SEM
# derivative grows like N²/h, not 1/Δ_min -- so it cancels out of the ratio,
# which is the number this case is about.
const DT_SMALL = 0.02     # both schemes stable
const DT_BIG   = 0.05     # HEVI stable with margin, explicit not
const GAMMA    = ark_tableau(:ARS343).γ

fac  = build_column_factorization(params, op, cc, topo, GAMMA * DT_SMALL)
hevi = HEVICache(topo, op, cc, fac, hevi_fimp!, hevi_column_solve!, :ARS343, GAMMA * DT_SMALL,
                 :RS, 5)
ph   = merge(params, (rtb = RTBCache(params), hevi = hevi))

npoin = params.mesh.npoin
say(@sprintf("  mesh: %d points/rank, %d levels, %d columns; dx_elem = %.0f m, dz_elem = %.0f m",
             npoin, topo.nlev, topo.ncol, L/NELX, L/NELZ))

# LGL clustering: the smallest gap in an order-4 element is 0.1727 h.
const HX = 0.1727 * L/NELX
const HZ = 0.1727 * L/NELZ
say(@sprintf("  smallest node gaps: dx = %.1f m, dz = %.1f m  (anisotropy %.2f)",
             HX, HZ, HX/HZ))

@testset "HEVI rising thermal bubble" begin

#--- 1. the base state is (nearly) in discrete hydrostatic balance -----------
# A strong-form collocation scheme does not balance dp/dz against -ρg exactly:
# the pressure profile is not a polynomial, so its derivative carries a
# truncation error. What matters is that the leftover acceleration is small
# next to the buoyant one the bubble is supposed to feel, g·θc/θ0.
ubase = rtb_initial_state(params; θc = 0.0, θ0 = Θ0)
dubase = similar(ubase)
rtb_rhs!(dubase, ubase, ph, 0.0)
accel = 0.0
for ip = 1:npoin
    accel = max(accel, abs(dubase[3*npoin + ip]) / ubase[0*npoin + ip])
end
accel = MPI.Allreduce(accel, MPI.MAX, COMM)
PC    = PhysicalConst{Float64}()
buoy  = PC.g * ΘC / Θ0
# Expressed against dp/dz the same residual is a 1e-4 relative error in the
# pressure gradient -- an accurate derivative. It only looks large next to the
# buoyant acceleration of a 0.5 K perturbation, which is what "not well
# balanced" means for a strong-form scheme in TOTAL variables. (Jexpresso's
# answer to this is :SOL_VARS_TYPE => PERT(), which subtracts the base state
# and cancels the two large terms analytically; CompEuler/theta uses it.)
#
# It does not affect the checks that follow: the imbalance lives in f_full, so
# HEVI and the explicit scheme inherit exactly the same one and it cancels
# out of any comparison between them.
@test accel < 0.15 * buoy
say(@sprintf("  hydrostatic residual: %.3e m/s^2 = %.1f%% of buoyant, = %.1e relative in dp/dz",
             accel, 100*accel/buoy, accel/PC.g))

#--- 2. the split does not disturb that balance ------------------------------
# f_imp acts on (u - qe) and the bubble-free state IS qe, so f_imp vanishes on
# it identically. HEVI must therefore reproduce the explicit scheme's drift
# exactly, not merely closely: any difference here is the split leaking into a
# state it should not touch at all.
fimp0 = similar(ubase)
hevi_fimp!(fimp0, ubase, ph, 0.0)
fmax = MPI.Allreduce(maximum(abs, fimp0), MPI.MAX, COMM)
@test fmax == 0.0
say(@sprintf("  f_imp on the bubble-free base state: %.1e  (must be exactly 0)", fmax))

#--- 3. HEVI and the explicit scheme agree on the bubble ---------------------
u0 = rtb_initial_state(params; θc = ΘC, r0 = R0, θ0 = Θ0)
d0 = rtb_diagnostics(u0, params; θ0 = Θ0)
# Slightly under θc: the bubble centre does not sit on an LGL node, so the
# nodal maximum undershoots the analytic peak.
@test isapprox(d0.θmax, ΘC; rtol = 2e-3)
say(@sprintf("  t=0: theta' max = %.4f K, centre of buoyancy z = %.1f m", d0.θmax, d0.zbar))

TCMP = 4.0
prob = ODEProblem(rtb_rhs!, copy(u0), (0.0, TCMP), ph)
ue  = solve(prob, Tsit5(); dt = DT_SMALL, adaptive = false, save_everystep = false).u[end]
uh  = solve(prob, HEVI_ARK(:ARS343); dt = DT_SMALL, adaptive = false,
            save_everystep = false).u[end]
refactorize!(fac, params, op, cc, topo, GAMMA * DT_BIG)
uhb = solve(prob, HEVI_ARK(:ARS343); dt = DT_BIG, adaptive = false,
            save_everystep = false).u[end]

# Compare on w: it starts at exactly zero and is what buoyancy drives, so it
# has no large background to hide a discrepancy in.
wof(u, ip) = u[3*npoin + ip] / u[0*npoin + ip]
function wgap(ua, ub)
    d = 0.0; s = 0.0
    for ip = 1:npoin
        d = max(d, abs(wof(ua, ip) - wof(ub, ip))); s = max(s, abs(wof(ub, ip)))
    end
    MPI.Allreduce(d, MPI.MAX, COMM) / MPI.Allreduce(s, MPI.MAX, COMM)
end

gap_scheme = wgap(uh, ue)     # HEVI vs explicit, same step
gap_step   = wgap(uhb, uh)    # HEVI at the big step vs HEVI at the small one
@test gap_scheme < 5e-3
@test gap_step   < 5e-2
say(@sprintf("  after %.0f s: |w_HEVI - w_explicit|/|w| = %.2e at dt=%.3g;  HEVI dt=%.3g vs %.3g: %.2e",
             TCMP, gap_scheme, DT_SMALL, DT_BIG, DT_SMALL, gap_step))

#--- 4. HEVI runs the bubble at a step the explicit scheme cannot take -------
TRUN = 60.0
uh2 = solve(ODEProblem(rtb_rhs!, copy(u0), (0.0, TRUN), ph), HEVI_ARK(:ARS343);
            dt = DT_BIG, adaptive = false, save_everystep = false).u[end]
d = rtb_diagnostics(uh2, params; θ0 = Θ0)

@test all(isfinite, uh2)
@test d.zbar > d0.zbar + 5.0                 # the bubble went UP
@test 0.05 < d.wmax < 10.0                   # and at a physical speed
@test d.θmax < 1.5 * ΘC                      # without the peak blowing up
@test d.θmax > 0.5 * ΘC                      # or being diffused away
@test d.vmax < 1e-8 * max(d.wmax, 1.0)       # and stayed y-invariant
say(@sprintf("  HEVI at dt=%.3g for %.0f s: theta' max = %.4f K, z = %.1f m (+%.1f), w max = %.3f m/s, |v| max = %.1e",
             DT_BIG, TRUN, d.θmax, d.zbar, d.zbar - d0.zbar, d.wmax, d.vmax))

#--- 5. the explicit scheme at that step does not survive --------------------
NBLOW = 400
ue2 = with_logger(NullLogger()) do
    solve(ODEProblem(rtb_rhs!, copy(u0), (0.0, NBLOW*DT_BIG), ph), Tsit5();
          dt = DT_BIG, adaptive = false, save_everystep = false,
          unstable_check = (dt, u, p, t) -> false).u[end]
end

# Both reductions run UNCONDITIONALLY, and are combined only afterwards.
#
# The obvious form -- `!all(isfinite, ue2) || MPI.Allreduce(...) > 1e3` --
# deadlocks: `||` short-circuits, so the first rank whose solution has already
# gone non-finite skips the Allreduce that every other rank is sitting in. It
# passes in serial and hangs on 3 ranks after the assertions have all gone
# green, which is about the worst way for a test to be wrong.
g       = maximum(abs, ue2 .- u0)
nonfin  = MPI.Allreduce((all(isfinite, ue2) && isfinite(g)) ? 0.0 : 1.0, MPI.MAX, COMM)
growth  = MPI.Allreduce(isfinite(g) ? g : 0.0, MPI.MAX, COMM)
@test nonfin > 0.0 || growth > 1e3
say(@sprintf("  explicit at dt=%.3g: diverged within %d steps, as expected", DT_BIG, NBLOW))
say(@sprintf("  -> HEVI took a %.1fx larger step, on a mesh with %.1fx vertical anisotropy,",
             DT_BIG/DT_SMALL, HX/HZ))
say(@sprintf("     and the two step sizes agree to %.1e in w", gap_step))

end # testset

MPI.Barrier(COMM)
say("\n=== rising thermal bubble: all checks passed on $NR rank(s) ===\n")

#=============================================================================
 runtests_standalone.jl -- IMEX3D unit tests that do NOT need the Jexpresso
 dependency stack.

     julia --project=<env> test/imex3d/runtests_standalone.jl
     mpiexecjl -n 3 julia --project=<env> test/imex3d/runtests_standalone.jl

 Needs only MPI, LinearAlgebra, Printf and OrdinaryDiffEq, and reuses the
 spectral-element fixture the HEVI tests are built on (test/hevi/mock_sem.jl).

 RUN IT ON SEVERAL RANK COUNTS. Two of the things under test exist ONLY in
 parallel and are invisible on one rank:

   * the multiplicity-weighted inner product. On one rank every node is held
     once and the weights are all 1, so a wrong weighting is a no-op; on three
     it is the difference between a norm and a number.
   * the preconditioner's gather/scatter, which on 3 ranks over the default
     2x2x5 mesh has to move columns that the partition cuts through.
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "..", "hevi", "mock_sem.jl"))
include(joinpath(SRC, "columns.jl"))
include(joinpath(SRC, "vdiffusion.jl"))
include(joinpath(SRC, "operator.jl"))
include(joinpath(SRC, "factorize.jl"))
include(joinpath(SRC, "acoustic.jl"))
include(joinpath(SRC, "ark.jl"))
include(joinpath(SRC, "hevi.jl"))
include(joinpath(SRC, "krylov.jl"))
include(joinpath(SRC, "imex3d.jl"))

say(args...) = RANK == 0 && (println(args...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
const LX, LY, LZ = 400.0, 400.0, 1000.0
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lx = LX, Ly = LY, Lz = LZ, comm = COMM)

const NPOIN = Int(params.mesh.npoin)
const GIP   = params.mesh.ip2gip
const QE    = params.qp.qe

say("\n=== IMEX3D standalone tests: $NR rank(s), $(NELX)x$(NELY)x$(NELZ) elements, p=$P ===")

@testset "IMEX3D" begin

topo = build_column_topology(params.mesh, COMM)

#--- 1. the 3D operator is the vertical one plus the horizontal sweeps -------
# Both halves of this matter. Reducing to the vertical operator on a
# horizontally uniform field tests that the new ξ and η sweeps contribute
# EXACTLY zero where they should -- a wrong metric index or flux slot there
# does not cancel. Responding in the horizontal momenta to an x-varying field
# tests that they contribute at all: an operator that had silently collapsed
# back to HEVI's would pass the first check perfectly.
opfull = build_hevi_fast_operator(params, topo)
opvert = build_hevi_operator(params, topo, [1, 2, 3, 4, 5]; full = false)
vf = hevi_verify_fast(params, opfull, opvert, topo)
@test vf.rel_vertical < 1e-12
@test vf.horiz_response > 0.0
@test vf.ok
say(@sprintf("  full vs vertical on a horizontally uniform field: %.2e", vf.rel_vertical))
say(@sprintf("  horizontal momentum response to a field varying in x: %.3g", vf.horiz_response))

#--- 2. lateral free-slip walls ----------------------------------------------
# Two-sided, because a one-sided test passes for an operator that ignores the
# flags entirely OR for one that zeroes far too much:
#
#   (a) on a field whose horizontal momenta ALREADY vanish at the wall nodes,
#       the walled operator must agree with the unwalled one to round-off --
#       the projection has nothing left to remove;
#   (b) on a field whose horizontal momenta do NOT vanish there, it must
#       differ, and by an amount of the operator's own size.
wx, wy, src = imex_lateral_walls(params, :box)
@test src === :box
@test any(wx) && any(wy)
opwall = build_hevi_fast_operator(params, topo; wallx = wx, wally = wy)

Vok = zeros(NPOIN, 5)
for ip = 1:NPOIN
    a = sinpi(1e-3 * GIP[ip])
    Vok[ip, 1] = a; Vok[ip, 4] = 0.3a; Vok[ip, 5] = 0.7a
    # tangential at the walls by construction
    Vok[ip, 2] = wx[ip] ? 0.0 : a
    Vok[ip, 3] = wy[ip] ? 0.0 : a
end
Rw = hevi_apply_A!(zeros(NPOIN, 5), Vok, params, opwall)
Rn = hevi_apply_A!(zeros(NPOIN, 5), Vok, params, opfull)
dagree = MPI.Allreduce(maximum(abs, Rw .- Rn), MPI.MAX, COMM)
scale  = MPI.Allreduce(maximum(abs, Rn),       MPI.MAX, COMM)
@test dagree <= 1e-12 * scale
say(@sprintf("  wall projection is a no-op on an already-tangential field: %.2e", dagree))

Vbad = copy(Vok)
for ip = 1:NPOIN
    Vbad[ip, 2] = sinpi(1e-3 * GIP[ip])
    Vbad[ip, 3] = sinpi(1e-3 * GIP[ip])
end
Rw2 = hevi_apply_A!(zeros(NPOIN, 5), Vbad, params, opwall)
Rn2 = hevi_apply_A!(zeros(NPOIN, 5), Vbad, params, opfull)
ddiff = MPI.Allreduce(maximum(abs, Rw2 .- Rn2), MPI.MAX, COMM)
scal2 = MPI.Allreduce(maximum(abs, Rn2),        MPI.MAX, COMM)
@test ddiff > 1e-6 * scal2
say(@sprintf("  and does bite on a field with normal momentum at the wall: %.2e (rel %.2e)",
             ddiff, ddiff / scal2))

#--- 2b. the columnar-partition guard ----------------------------------------
# Shared by HEVI and IMEX3D, and the reason it is tested here rather than
# trusted: it used to read `:lxy_partition` alone, which is a REQUEST the mesh
# code ignores whenever refinement is on. A deck with `:linitial_refine =>
# true` then ran on p4est's space-filling-curve partition -- where the
# operator is not skew and the solution grows at a few 1/s -- with the guard
# and every self-check green. See check_columnar_partition in columns.jl.
let
    cases = ((:defaults,           Dict{Symbol,Any}(),                                          false),
             (:flag_false,         Dict{Symbol,Any}(:lxy_partition => false),                   true),
             (:flag_true_refine,   Dict{Symbol,Any}(:lxy_partition => true,
                                                    :linitial_refine => true),                  true),
             (:refine_flag_absent, Dict{Symbol,Any}(:linitial_refine => true),                  true))
    for (name, deck, bad) in cases, scheme in ("HEVI", "IMEX3D")
        threw = try
            check_columnar_partition(deck, COMM, scheme); false
        catch
            true
        end
        # Serial is a single partition with no rank-shared nodes, so every
        # configuration is legitimate there; the guard must not fire.
        @test threw == (bad && NR > 1)
    end
    say(NR > 1 ? "  partition guard: refuses :lxy_partition=>false AND :linitial_refine=>true" :
                 "  partition guard: serial, every configuration allowed (correct)")
end

#--- 3. the distributed inner product ----------------------------------------
# The one quantity here with a known exact answer: <1, 1> is the number of
# GLOBAL degrees of freedom, whatever the rank count. Counting each rank's
# local nodes and summing would give more than that as soon as a node is
# shared, which is exactly the error the weighting exists to prevent.
inner = build_distributed_inner(opfull.dss_cache, NPOIN, COMM)
ones5 = ones(NPOIN, 5)
gdofs = ((NELX*P + 1) * (NELY*P + 1) * (NELZ*P + 1)) * 5
@test isapprox(ddot(inner, ones5, ones5), gdofs; rtol = 1e-12)
say(@sprintf("  <1,1> = %.1f, global dofs = %d", ddot(inner, ones5, ones5), gdofs))

# and it must not depend on the partition: a field that is a function of the
# global point id has a norm every rank count has to agree on.
Vg = [sinpi(1e-3 * (GIP[ip] + 11q)) for ip = 1:NPOIN, q = 1:5]
nrm = dnorm(inner, Vg)
@test MPI.Allreduce(nrm, MPI.MAX, COMM) ≈ MPI.Allreduce(nrm, MPI.MIN, COMM)
say(@sprintf("  ‖V‖ = %.12f  (identical on every rank)", nrm))

#--- 4. the preconditioner and the Krylov solve ------------------------------
GDT = 0.02
# The preconditioner drops the rows that are exactly the identity on a mesh
# whose ζ lines are vertical -- see IMEX3DPrecond. `hevi_choose_vars` reads the
# metrics, so on this Cartesian box it returns (1, 4, 5).
PVARS = hevi_choose_vars(params.metrics, COMM)
@test PVARS == [1, 4, 5]
opprec = build_hevi_operator(params, topo, PVARS; full = false)
owner, own = assign_column_owners(topo, COMM)
cc  = build_column_comm(topo, owner, own, COMM, length(PVARS))
fac = build_column_factorization(params, opprec, cc, topo, GDT)
pc  = IMEX3DPrecond(opprec, cc, fac, topo, PVARS)

matvec! = (W, V) -> begin
    hevi_apply_A!(W, V, params, opwall)
    @. W = V - GDT * W
    W
end

B = [sinpi(1e-3 * (GIP[ip] + 53q)) for ip = 1:NPOIN, q = 1:5]

function solve_and_check(ws, precon!)
    X = zeros(NPOIN, 5)
    it, rel, conv = gmres_solve!(X, B, ws, matvec!, precon!)
    R = matvec!(zeros(NPOIN, 5), X)
    err = MPI.Allreduce(maximum(abs, R .- B), MPI.MAX, COMM)
    scl = MPI.Allreduce(maximum(abs, B),      MPI.MAX, COMM)
    return (X, it, rel, conv, err / scl)
end

ws_pc = GMRESWorkspace(NPOIN, 5, inner; m = 30, maxiter = 300, rtol = 1e-10)
Xp, itp, relp, convp, resp = solve_and_check(ws_pc, Z -> imex3d_precond!(Z, pc, params, GDT))
@test convp
# The APPLIED residual, not the one the Arnoldi recursion reported. They differ
# if the recursion has drifted, and it is the applied one the integrator feels.
@test resp < 1e-8
say(@sprintf("  GMRES + column preconditioner: %d iterations, applied residual %.2e",
             itp, resp))

ws_no = GMRESWorkspace(NPOIN, 5, inner; m = 30, maxiter = 300, rtol = 1e-10)
Xn, itn, reln, convn, resn = solve_and_check(ws_no, Z -> Z)
@test convn
say(@sprintf("  GMRES unpreconditioned:        %d iterations, applied residual %.2e",
             itn, resn))

# Both must reach the SAME answer; a preconditioner changes the path, not the
# solution. This is the check that would catch a preconditioner that is not an
# approximate inverse of anything -- an unconverged or wrong M just costs
# iterations, but a right-preconditioned solve whose final `X += M^-1(Vy)` step
# used a DIFFERENT M than the iteration did would land somewhere else.
agree = MPI.Allreduce(maximum(abs, Xp .- Xn), MPI.MAX, COMM)
scl   = MPI.Allreduce(maximum(abs, Xp),       MPI.MAX, COMM)
@test agree < 1e-6 * max(scl, 1.0)
say(@sprintf("  preconditioned and unpreconditioned answers agree to %.2e", agree))

#--- 4b. the preconditioner where it is supposed to matter -------------------
# On the mesh above, Δx = Δy = Δz, and a vertical-line preconditioner on an
# ISOTROPIC mesh has almost nothing to remove -- the vertical acoustic term it
# inverts is no stiffer than the two horizontal ones it leaves to the Krylov
# iteration. Measured there it saves one iteration out of eight, which is a
# true number and a useless test: it would stay green if the preconditioner
# were quietly doing nothing at all.
#
# The regime this scheme is FOR is an atmospheric mesh, where Δz is several
# times finer than Δx and the vertical acoustic term dominates the spectrum.
# So the claim is tested at 8:1, where the column solve is inverting most of
# what makes the operator stiff and the iteration count has room to separate.
let
    pa = build_mock_params(nelx = 2, nely = 2, nelz = 20, p = P,
                           Lx = 800.0, Ly = 800.0, Lz = 1000.0, comm = COMM)
    np = Int(pa.mesh.npoin)
    tp = build_column_topology(pa.mesh, COMM)
    of = build_hevi_fast_operator(pa, tp)
    pv = hevi_choose_vars(pa.metrics, COMM)
    ov = build_hevi_operator(pa, tp, pv; full = false)
    ow, on = assign_column_owners(tp, COMM)
    ca = build_column_comm(tp, ow, on, COMM, length(pv))
    fa = build_column_factorization(pa, ov, ca, tp, GDT)
    pa_pc = IMEX3DPrecond(ov, ca, fa, tp, pv)
    ia = build_distributed_inner(of.dss_cache, np, COMM)

    mv! = (W, V) -> begin
        hevi_apply_A!(W, V, pa, of)
        @. W = V - GDT * W
        W
    end
    gipa = pa.mesh.ip2gip
    Ba = [sinpi(1e-3 * (gipa[ip] + 53q)) for ip = 1:np, q = 1:5]

    wa = GMRESWorkspace(np, 5, ia; m = 40, maxiter = 400, rtol = 1e-10)
    Xa = zeros(np, 5)
    ita, _, cva = gmres_solve!(Xa, Ba, wa, mv!, Z -> imex3d_precond!(Z, pa_pc, pa, GDT))

    wb = GMRESWorkspace(np, 5, ia; m = 40, maxiter = 400, rtol = 1e-10)
    Xb = zeros(np, 5)
    itb, _, cvb = gmres_solve!(Xb, Ba, wb, mv!, Z -> Z)

    @test cva && cvb
    say(@sprintf("  8:1 anisotropic mesh: %d iterations preconditioned, %d without", ita, itb))
    # The claim the whole design rests on: the column solve is a real
    # preconditioner for the 3D operator, not decoration. If this ever
    # regresses, the scheme still gives the right answer and costs several
    # times more per step than it should, with nothing else in the output to
    # say so.
    @test ita < itb
    @test ita <= itb ÷ 2
    dab = MPI.Allreduce(maximum(abs, Xa .- Xb), MPI.MAX, COMM)
    sab = MPI.Allreduce(maximum(abs, Xa),       MPI.MAX, COMM)
    @test dab < 1e-6 * max(sab, 1.0)
end

#--- 5. the stage solve, end to end ------------------------------------------
ws   = GMRESWorkspace(NPOIN, 5, inner; m = 30, maxiter = 300, rtol = 1e-12)
imex = IMEX3DCache(topo, opwall, pc, ws, zeros(NPOIN, 5), zeros(NPOIN, 5),
                   imex3d_fimp!, imex3d_solve!, :ARS343, GDT, :RS, 5, COMM)
ph   = merge(params, (imex = imex,))

srcv = zeros(NPOIN * 5)
for ieq = 1:5, ip = 1:NPOIN
    srcv[(ieq-1)*NPOIN + ip] = QE[ip, ieq] + 1e-3 * sinpi(1e-3 * (GIP[ip] + 17ieq))
end
dstv = similar(srcv)
imex3d_solve!(dstv, srcv, ph, GDT)

# (dst - qe) - γΔt A (dst - qe) must equal (src - qe)
D  = [dstv[(ieq-1)*NPOIN + ip] - QE[ip, ieq] for ip = 1:NPOIN, ieq = 1:5]
AD = hevi_apply_A!(zeros(NPOIN, 5), D, params, opwall)
res = 0.0; scl = 0.0
for ieq = 1:5, ip = 1:NPOIN
    want = srcv[(ieq-1)*NPOIN + ip] - QE[ip, ieq]
    got  = D[ip, ieq] - GDT * AD[ip, ieq]
    res = max(res, abs(got - want)); scl = max(scl, abs(want))
end
res = MPI.Allreduce(res, MPI.MAX, COMM); scl = MPI.Allreduce(scl, MPI.MAX, COMM)
@test res / scl < 1e-9
say(@sprintf("  stage solve residual: %.2e relative", res / scl))

# Single-valued across ranks. Nothing else here can see this failure: every
# other check feeds the machinery a field that is single-valued by
# construction, and a CG state that has gone double-valued is differentiated
# across the jump by the next explicit RHS.
sv = imex_single_valued_error(D, opwall, COMM)
@test sv < 1e-12
say(@sprintf("  cross-rank single-valuedness of the answer: %.2e", sv))

#--- 6. well-balancedness: f_imp(qe) is EXACTLY zero -------------------------
# Exactly, not approximately: the operator acts on u - qe, so the reference
# state is a steady state of the implicit half by construction. Whatever
# discrete hydrostatic imbalance the discretisation carries stays entirely in
# the explicit part, and no rearrangement of the split can create or destroy
# balance.
qeflat = vec([QE[ip, ieq] for ip = 1:NPOIN, ieq = 1:5])
du0 = zeros(NPOIN * 5)
hevi_fast_rhs!(du0, copy(qeflat), params, opwall)
@test MPI.Allreduce(maximum(abs, du0), MPI.MAX, COMM) == 0.0
say("  f_imp(qe) = 0 exactly")

#--- 7. the integrator: order, and stability past the explicit limit ---------
# A synthetic "full RHS" whose fast part is exactly the implicit operator plus
# a slow relaxation. This is the only check that would catch fimp! and solve!
# disagreeing with EACH OTHER -- a sign, a scaling, a stale γΔt. Every other
# self-check validates one half against itself and would pass.
#
# Run from inside a function, deliberately: a closure defined at the top level
# of a @testset captures non-const globals, every `p.imex` access becomes a
# dynamic dispatch, and the run takes minutes instead of seconds -- which
# under MPI is indistinguishable from a hang.
function integration_checks(params, ph, qeflat, npoin)
    f_full! = let qeflat = qeflat
        function (du, u, p, t)
            imex3d_fimp!(du, u, p, t)
            @. du += -0.1 * (u - qeflat)
            du
        end
    end
    gip = params.mesh.ip2gip
    u0 = qeflat .+ [1e-3 * sinpi(1e-3 * (gip[mod1(i, npoin)] + 7 * div(i-1, npoin)))
                    for i = 1:npoin*5]
    T = 2.0

    # Fixed step for the reference, deliberately. An ADAPTIVE reference
    # deadlocks on more than one rank: each rank runs its own error controller,
    # they disagree about the step size, and the MPI collectives inside f_imp
    # are then reached a different number of times per rank.
    ref = solve(ODEProblem(f_full!, copy(u0), (0.0, T), ph), Tsit5();
                dt = T/400, adaptive = false, save_everystep = false).u[end]

    errs = Float64[]
    for n in (40, 80, 160)
        sol = solve(ODEProblem(f_full!, copy(u0), (0.0, T), ph), IMEX_ARK(:ARS343);
                    dt = T/n, adaptive = false, save_everystep = false).u[end]
        push!(errs, MPI.Allreduce(maximum(abs, sol .- ref), MPI.MAX, COMM))
    end

    # Stability far past where an explicit scheme is not. The 3D operator
    # carries the acoustics in ALL directions, so the explicit limit here is
    # set by the horizontal as well as the vertical -- which is precisely what
    # this scheme removes and HEVI does not.
    λ = 0.0
    let V = [sinpi(1e-3 * (gip[ip] + 3q)) for ip = 1:npoin, q = 1:5]
        # a crude power-iteration bound on ‖A‖, enough to pick a "big" dt
        R = hevi_apply_A!(zeros(npoin, 5), V, params, ph.imex.op)
        λ = MPI.Allreduce(maximum(abs, R), MPI.MAX, COMM) /
            MPI.Allreduce(maximum(abs, V), MPI.MAX, COMM)
    end
    big = 8.0 * 3.34 / max(λ, eps())
    si = solve(ODEProblem(f_full!, copy(u0), (0.0, 40*big), ph), IMEX_ARK(:ARS343);
               dt = big, adaptive = false, save_everystep = false).u[end]
    # unstable_check off: the built-in abort fires when u goes non-finite,
    # which need not happen on the same step on every rank.
    se = with_logger(NullLogger()) do
        solve(ODEProblem(f_full!, copy(u0), (0.0, 40*big), ph), Tsit5();
              dt = big, adaptive = false, save_everystep = false,
              unstable_check = (dt, u, p, t) -> false).u[end]
    end
    imax = MPI.Allreduce(maximum(abs, si .- qeflat), MPI.MAX, COMM)
    emax = MPI.Allreduce(maximum(abs, se .- qeflat), MPI.MAX, COMM)
    return errs, big, λ, imax, emax
end

errs, big, λ, imax, emax = integration_checks(params, ph, qeflat, NPOIN)
ordr = log2(errs[2] / errs[3])
@test errs[end] < 1e-6
@test ordr > 2.5
@test imax < 1e-2                # bounded
@test !(emax < 1e-2)             # explicit is not
say(@sprintf("  ARS343 vs Tsit5 reference: err %.2e, observed order %.2f", errs[end], ordr))
say(@sprintf("  at dt = %.3g s (~8x an explicit scheme's 3D acoustic limit, ‖A‖ ≈ %.3g 1/s):",
             big, λ))
say(@sprintf("      IMEX3D |q'| = %.2e, explicit |q'| = %.2e", imax, emax))

end # testset

MPI.Barrier(COMM)
say("\n=== IMEX3D standalone tests done ===")

#=============================================================================
 runtests_standalone.jl -- HEVI unit tests that do NOT need the Jexpresso
 dependency stack.

     julia --project=<env> test/hevi/runtests_standalone.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/runtests_standalone.jl

 Needs only MPI, LinearAlgebra, Printf and OrdinaryDiffEq. Run it on several
 rank counts: with 3 ranks over the default 2x2x5 element mesh the partition
 cuts through columns, which is the case the redistribution exists for and
 the one a single-rank run cannot reach.
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
include(joinpath(SRC, "columns.jl"))
include(joinpath(SRC, "operator.jl"))
include(joinpath(SRC, "factorize.jl"))
include(joinpath(SRC, "ark.jl"))
include(joinpath(SRC, "hevi.jl"))

say(args...) = RANK == 0 && (println(args...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
const LZ = 1000.0
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P, Lz = LZ, comm = COMM)

say("\n=== HEVI standalone tests: $NR rank(s), $(NELX)x$(NELY)x$(NELZ) elements, p=$P ===")

@testset "HEVI" begin

#--- 1. column topology -------------------------------------------------------
topo = build_column_topology(params.mesh, COMM)
@test topo.nlev  == NELZ*P + 1
@test topo.ncolg == (NELX*P + 1) * (NELY*P + 1)
# every local node landed in exactly one slot
@test count(!iszero, topo.node) == params.mesh.npoin
@test sum(topo.nheld) == params.mesh.npoin
nsplit_g = MPI.Allreduce(topo.nsplit, MPI.SUM, COMM)
say(@sprintf("  topology: %d levels, %d global columns, %d local, %d split (global)",
             topo.nlev, topo.ncolg, topo.ncol, nsplit_g))

vars  = hevi_choose_vars(params.metrics, COMM)
@test vars == [1, 4, 5]          # dζdx = dζdy = 0 on a Cartesian box
op    = build_hevi_operator(params, topo, vars)
owner, own = assign_column_owners(topo, COMM)
@test MPI.Allreduce(length(own), MPI.SUM, COMM) == topo.ncolg   # every column owned once
cc    = build_column_comm(topo, owner, own, COMM, length(vars))

#--- 2. gather / scatter is the identity on a consistent field ----------------
# The field has to take the same value on every rank holding a shared node,
# so it is built from the global point id rather than from rand().
gip = params.mesh.ip2gip
V0  = [sinpi(1e-3*(gip[ip] + 31*q)) for ip = 1:params.mesh.npoin, q = 1:op.nimp]
V   = copy(V0)
column_gather!(cc, V, 1:op.nimp)
fill!(V, NaN)
column_scatter!(cc, V, 1:op.nimp)
@test maximum(abs, V .- V0) == 0.0
say("  gather/scatter round trip: exact")

#--- 3. the operator couples nodes only within one column --------------------
# Light up one column and check nothing outside it responds.
lev = level_of_node(topo, params.mesh.npoin)
target = topo.gcol[1]
target = MPI.Allreduce(target, MPI.MIN, COMM)
colof = zeros(Int, params.mesh.npoin)
for ic = 1:topo.ncol, il = 1:topo.nlev
    ip = topo.node[il, ic]; ip != 0 && (colof[ip] = topo.gcol[ic])
end
Vc = zeros(params.mesh.npoin, op.nimp)
for ip = 1:params.mesh.npoin
    colof[ip] == target && (Vc[ip, 3] = 1.0 + 0.1*lev[ip])
end
Rc = similar(Vc)
hevi_apply_A!(Rc, Vc, params, op)
leak = 0.0
for ip = 1:params.mesh.npoin, q = 1:op.nimp
    colof[ip] != target && (leak = max(leak, abs(Rc[ip, q])))
end
inside = maximum(abs, Rc)
leak = MPI.Allreduce(leak, MPI.MAX, COMM)
inside = MPI.Allreduce(inside, MPI.MAX, COMM)
@test leak <= 1e-12 * max(inside, 1.0)
say(@sprintf("  column locality: leakage %.2e vs in-column %.2e", leak, inside))

#--- 4. the assembled band IS the operator ------------------------------------
γdt = 0.4358665215084590 * 0.05
AB, kl, ku, n = assemble_column_band(params, op, cc, topo, γdt)
verr = hevi_verify_band(params, op, cc, topo, AB, kl, ku, n, γdt)
@test verr < 1e-12
say(@sprintf("  band == operator to %.2e relative", verr))

#--- 5. spectrum of one column: an acoustic operator is (nearly) skew ---------
# Recover A from W = I - γΔt A and look at its eigenvalues. A hyperbolic
# operator has essentially imaginary spectrum; a real part that is not tiny
# would mean the discretisation is spuriously damping or amplifying sound.
if !isempty(AB)
    d = kl + ku + 1
    W = zeros(n, n)
    for j = 1:n, i = max(1, j-ku):min(n, j+kl)
        W[i, j] = AB[1][d + i - j, j]
    end
    Amat = (I - W) ./ γdt
    λ = eigvals(Amat)
    reim = maximum(abs, real.(λ)) / maximum(abs, imag.(λ))
    PC = PhysicalConst{Float64}()
    c  = sqrt(PC.γ * perfectGasLaw_ρθtoP(PC, 1.2, 300.0) / 1.2)
    dzmin = LZ/NELZ * 0.5 * (1.0 - (-1.0))/2   # crude: min LGL gap for p=2 is h/2
    say(@sprintf("  column spectrum: max|Re|/max|Im| = %.2e, max|Im| = %.3g 1/s (c/dz ~ %.3g)",
                 reim, maximum(abs, imag.(λ)), c/dzmin))
    @test reim < 1e-6
    @test 0.2 < maximum(abs, imag.(λ)) / (c/dzmin) < 20.0
end

#--- 6. end-to-end stage solve ------------------------------------------------
fac = build_column_factorization(params, op, cc, topo, γdt)
hevi = HEVICache(topo, op, cc, fac, hevi_fimp!, hevi_column_solve!, :ARS343, γdt, :RS, 5)
ph   = merge(params, (hevi = hevi,))

npoin = params.mesh.npoin
qe    = params.qp.qe
src   = zeros(npoin*5)
for ieq = 1:5, ip = 1:npoin
    src[(ieq-1)*npoin + ip] = qe[ip, ieq] + 1e-3*sinpi(1e-3*(gip[ip] + 17*ieq))
end
dst = similar(src)
hevi_column_solve!(dst, src, ph, γdt)

# (dst - qe) - γΔt A (dst - qe) must equal (src - qe)
D = [dst[(ieq-1)*npoin + ip] - qe[ip, ieq] for ip = 1:npoin, ieq in op.vars]
AD = similar(D)
hevi_apply_A!(AD, D, params, op)
res = 0.0; scl = 0.0
for (q, ieq) in enumerate(op.vars), ip = 1:npoin
    want = src[(ieq-1)*npoin + ip] - qe[ip, ieq]
    got  = D[ip, q] - γdt*AD[ip, q]
    res = max(res, abs(got - want)); scl = max(scl, abs(want))
end
res = MPI.Allreduce(res, MPI.MAX, COMM); scl = MPI.Allreduce(scl, MPI.MAX, COMM)
@test res / scl < 1e-9
say(@sprintf("  stage solve residual: %.2e relative", res/scl))

# equations outside the implicit set must come through untouched
untouched = 0.0
for ieq in setdiff(1:5, op.vars), ip = 1:npoin
    untouched = max(untouched, abs(dst[(ieq-1)*npoin+ip] - src[(ieq-1)*npoin+ip]))
end
@test MPI.Allreduce(untouched, MPI.MAX, COMM) == 0.0

# the answer must be single-valued across ranks (CG consistency)
if NR > 1
    Dm = copy(D)
    asm = setup_assembler(NSD_3D(), Dm, params.mesh.ip2gip, params.mesh.gip2owner)
    S = copy(Dm); assemble_mpi!(S, asm)     # S = sum over holders
    cnt = ones(npoin, op.nimp); assemble_mpi!(cnt, asm)
    spread = maximum(abs, Dm .- S ./ cnt)
    spread = MPI.Allreduce(spread, MPI.MAX, COMM)
    @test spread < 1e-12 * max(scl, 1.0)
    say(@sprintf("  cross-rank consistency of the solve: %.2e", spread))
end

#--- 7. the factorised solve inverts the matrix it was built from ------------
serr = hevi_verify_solve(params, op, cc, fac)
@test serr < 1e-9
say(@sprintf("  factorised solve residual (setup self-check): %.2e", serr))

#--- 8. end to end: the integrator and the column solver agree ---------------
# A synthetic "full RHS" whose fast part is exactly the implicit operator plus
# a slow relaxation. If solve! and fimp! ever disagreed -- a sign, a scaling, a
# transposed band -- both earlier self-checks would still pass (each validates
# one half against itself) and only this comparison against an independent
# integrator would catch it.
#
# Run from inside a function, deliberately. Defining the RHS closure at the top
# level of the @testset makes it capture non-const globals, every `p.hevi`
# access goes through dynamic dispatch, and the run takes minutes instead of
# seconds -- which under MPI is indistinguishable from a hang.
function integration_checks(params, ph, gip, npoin, qe)
    qeflat = vec([qe[ip, ieq] for ip = 1:npoin, ieq = 1:5])
    f_full! = let qeflat = qeflat
        function (du, u, p, t)
            hevi_fimp!(du, u, p, t)
            @. du += -0.1 * (u - qeflat)
            du
        end
    end
    u0 = qeflat .+ [1e-3*sinpi(1e-3*(gip[mod1(i, npoin)] + 7*div(i-1, npoin)))
                    for i = 1:npoin*5]
    T = 2.0

    # Fixed step, deliberately. An ADAPTIVE reference deadlocks on more than one
    # rank: each rank runs its own error controller, they disagree on the step
    # size, and the MPI collectives inside f_imp are then reached a different
    # number of times per rank. The same hazard is why HEVI_ARK refuses
    # :ode_adaptive_solver => true. Order 5 at dt = T/400 leaves a reference
    # error around 3e-12, far below the 1e-8 being measured.
    ref = solve(ODEProblem(f_full!, copy(u0), (0.0, T), ph), Tsit5();
                dt = T/400, adaptive = false, save_everystep = false).u[end]

    errs = Float64[]
    for n in (40, 80, 160)
        sol = solve(ODEProblem(f_full!, copy(u0), (0.0, T), ph), HEVI_ARK(:ARS343);
                    dt = T/n, adaptive = false, save_everystep = false).u[end]
        push!(errs, MPI.Allreduce(maximum(abs, sol .- ref), MPI.MAX, COMM))
    end

    # and it is stable far beyond where the explicit scheme is not
    big = 8.0 * (3.34 / 5.17)   # ~8x a 5-stage explicit scheme's acoustic limit
    sh = solve(ODEProblem(f_full!, copy(u0), (0.0, 100*big), ph), HEVI_ARK(:ARS343);
               dt = big, adaptive = false, save_everystep = false).u[end]
    # unstable_check off: the built-in abort fires when u goes non-finite, which
    # need not happen on the same step on every rank.
    se = with_logger(NullLogger()) do
        solve(ODEProblem(f_full!, copy(u0), (0.0, 100*big), ph), Tsit5();
              dt = big, adaptive = false, save_everystep = false,
              unstable_check = (dt, u, p, t) -> false).u[end]
    end
    hmax = MPI.Allreduce(maximum(abs, sh .- qeflat), MPI.MAX, COMM)
    emax = MPI.Allreduce(maximum(abs, se .- qeflat), MPI.MAX, COMM)
    return errs, big, hmax, emax
end

errs, big, hmax, emax = integration_checks(params, ph, gip, npoin, qe)
ordr = log2(errs[2] / errs[3])
@test errs[end] < 1e-6          # solution deviation is 1e-3, so this is 1e-3 relative
@test ordr > 2.5
@test hmax < 1e-2               # bounded
@test !(emax < 1e-2)            # explicit is not
say(@sprintf("  ARS343 vs Tsit5 reference: err %.2e, observed order %.2f", errs[end], ordr))
say(@sprintf("  at dt = %.2f s (~8x the explicit acoustic limit): HEVI |q'| = %.2e, explicit |q'| = %.2e",
             big, hmax, emax))

#--- 9. the typed function barrier actually does its job ---------------------
# Jexpresso's St_metrics and St_mesh have NO field type annotations, so
# metrics.dζdz and mesh.connijk are ::Any in a production run and every index
# inside an element loop is a boxed dynamic dispatch. The mocks above are
# concretely typed, so nothing else in this file can see that -- and it cost a
# real case a 60x slowdown in hevi_apply_A!, enough to make HEVI slower than
# the explicit scheme it is supposed to beat.
#
# hevi_apply_A! hands the arrays to _hevi_A_elements! positionally so the
# callee specialises on their runtime types. This checks that the untyped case
# is no longer dramatically slower, AND that it gives the same answer.
let
    pu   = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                             Lz = LZ, comm = COMM, untyped_metrics = true)
    topu = build_column_topology(pu.mesh, COMM)
    opu  = build_hevi_operator(pu, topu, vars)
    Vt   = [sinpi(1e-3*(gip[ip] + 11*q)) for ip = 1:params.mesh.npoin, q = 1:op.nimp]
    Rt = similar(Vt); Ru = similar(Vt)
    hevi_apply_A!(Rt, Vt, params, op)      # typed metrics
    hevi_apply_A!(Ru, Vt, pu,     opu)     # ::Any metrics
    @test maximum(abs, Rt .- Ru) == 0.0    # same answer, bit for bit

    bench(pp, oo, RR) = (hevi_apply_A!(RR, Vt, pp, oo);
                         minimum(@elapsed(hevi_apply_A!(RR, Vt, pp, oo)) for _ = 1:20))
    tt = bench(params, op,  Rt)
    tu = bench(pu,     opu, Ru)
    @test tu < 5 * tt
    say(@sprintf("  typed function barrier: ::Any metrics cost %.2fx the typed ones (%.1e vs %.1e s)",
                 tu/tt, tu, tt))
end

#--- 10. refactorisation for an off-nominal γΔt ------------------------------
# Not an exotic path: OrdinaryDiffEq shortens the step that lands on a tstop,
# and Jexpresso passes every diagnostic and statistics time as a tstop, so γΔt
# changes a few hundred times in a normal run.
# (the count is already non-zero: every HEVI_ARK stage above ran at a γΔt
# different from the one the factorisation was built with, which is exactly
# the production situation)
n0   = fac.nrefactor[]
gdt2 = 0.7 * γdt
refactorize!(fac, params, op, cc, topo, gdt2)
@test fac.gdt[] == gdt2
@test fac.nrefactor[] == n0 + 1
@test hevi_verify_solve(params, op, cc, fac) < 1e-9
refactorize!(fac, params, op, cc, topo, γdt)          # and back
@test hevi_verify_solve(params, op, cc, fac) < 1e-9
say("  refactorisation for a shortened step: residual holds, storage reused")

#--- 11. tableaux ------------------------------------------------------------
for (nm, ord, rad) in ((:ARS232, 2, 1.732), (:ARS343, 3, 2.828))
    tab = ark_tableau(nm)
    @test tab.order == ord
    @test isapprox(ark_imaginary_radius(tab), rad; atol = 2e-3)
    @test isapprox(sum(tab.bE), 1.0; atol = 1e-13)
    @test isapprox(sum(tab.bI), 1.0; atol = 1e-13)
    @test isapprox(dot(tab.bE, tab.cE), 0.5; atol = 1e-9)   # 2nd order, explicit
    @test isapprox(dot(tab.bI, tab.cI), 0.5; atol = 1e-9)   # 2nd order, implicit
    @test isapprox(dot(tab.bE, tab.cI), 0.5; atol = 1e-9)   # 2nd order, coupling
    @test isapprox(dot(tab.bI, tab.cE), 0.5; atol = 1e-9)
end
say("  tableaux: order conditions and stability radii check out")

end # testset

MPI.Barrier(COMM)
say("\n=== all HEVI standalone tests passed on $NR rank(s) ===\n")

#=============================================================================
 test_schur_stage.jl -- STEP 6: the PRODUCTION stage solve through the scalar
 Schur system, against the production five-field one.

 THE ONE THING THIS HAS TO PROVE
 -------------------------------
 `imex3d_solve_schur!` and `imex3d_solve!` solve the SAME linear system by
 different algebra. They must return the same five fields -- not "close", the
 same to the tolerance the two Krylov solves were asked for. If they do not,
 the elimination is wrong somewhere, and a wrong elimination gives a solver
 that converges beautifully to the wrong state.

 This runs the REAL entry points, not a replica of them: both `IMEX3DCache`
 objects are built here, `params.imex` is merged in exactly as `params_setup`
 does it, and the two production functions are called on the same `src`. A
 replica would prove the algebra and miss the wiring -- which is where the
 remaining risk actually is, since step 4 already proved the algebra densely.

 WHY THIS IS THE TEST THAT CATCHES A CACHE CLOBBER
 -------------------------------------------------
 `schur_setup_rhs!` leaves `q_b` and `m_b` in the SchurState, and
 `schur_recover!` reads them back after the whole Krylov iteration has run --
 hundreds of `schur_H!` calls in between, each of which writes the transient
 momentum slots of the same object. Comparing the two solves AFTER a full GMRES
 run is what makes a clobber visible; the dense step-4 test could not see it,
 because there the operator was applied to unit vectors and the right-hand side
 was rebuilt immediately before the recovery.

     julia --project=<env> test/imex3d/test_schur_stage.jl
     mpiexecjl -n 3 julia --project=<env> test/imex3d/test_schur_stage.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "..", "hevi", "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl", "acoustic.jl",
          "ark.jl", "hevi.jl", "krylov.jl", "precond_api.jl",
          "schur.jl", "schur_precond.jl", "schur_stage.jl", "imex3d.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 10, 2
const GDT = 0.30            # gamma*dt; well inside the stiff regime on this mock

params0 = build_mock_params(nelx=NELX, nely=NELY, nelz=NELZ, p=P,
                            Lx=800.0, Ly=800.0, Lz=1000.0, comm=COMM,
                            base=:stratified)
topo  = build_column_topology(params0.mesh, COMM)
npoin = Int(params0.mesh.npoin)

# THE ADVECTIVE ROW ON BOTH. The two solves must face the SAME operator, or the
# comparison measures the splitting difference (0.06% -- test_theta_advective)
# rather than the reduction, and would "fail" for a reason that is not a bug.
op = build_hevi_fast_operator(params0, topo; lwall_flux = true,
                              theta_advective = true)
inner = build_distributed_inner(op.dss_cache, npoin, COMM)

# Tight, so the comparison is about the algebra and not about where two loose
# iterations happened to stop.
const RTOL = 1.0e-12
inputs = Dict{Symbol,Any}(:imex_rtol => RTOL, :imex_restart => 40,
                          :imex_maxiter => 400, :imex_precond => :column)

#-- the five-field path, built exactly as build_imex3d builds it ---------------
pvars = hevi_choose_vars(params0.metrics, COMM)
opv5  = build_hevi_operator(params0, topo, pvars; lwall_flux = true, full = false)
own5  = assign_column_owners(topo, COMM)
cc5   = build_column_comm(topo, own5[1], own5[2], COMM, length(pvars))
fac5  = build_column_factorization(params0, opv5, cc5, topo, GDT)
pc5   = IMEX3DPrecond(opv5, cc5, fac5, topo, copy(pvars))
ws5   = GMRESWorkspace(npoin, op.nimp, inner; m = 40, maxiter = 400,
                       rtol = RTOL, atol = 1.0e-30)

#-- the scalar path -----------------------------------------------------------
sch = build_imex3d_schur(params0, topo, COMM, inputs, op, GDT, inner)

mkcache(pc, schur, ws, solve!) =
    IMEX3DCache(topo, op, pc, schur, ws,
                zeros(Float64, npoin, op.nimp), zeros(Float64, npoin, op.nimp),
                imex3d_fimp!, solve!, :ARS343, GDT, :RS, 1,
                false,            # warm_start OFF: each call must stand alone
                COMM)

params_full  = merge(params0, (imex = mkcache(pc5, nothing, ws5, imex3d_solve!),))
params_schur = merge(params0, (imex = mkcache(nothing, sch, ws5, imex3d_solve_schur!),))

# A state that is NOT the reference state, so the deviation is non-trivial in
# every equation. Built from the GLOBAL node id, so it is the same field however
# the mesh was cut and the comparison means the same thing on any rank count.
qe  = params0.qp.qe
gip = params0.mesh.ip2gip
src = zeros(Float64, 5 * npoin)
@inbounds for ieq = 1:5, ip = 1:npoin
    src[(ieq-1)*npoin + ip] = qe[ip, ieq] +
        (0.02 + 0.01*ieq) * sinpi(1e-2 * (gip[ip] + 5ieq)) * max(qe[ip, 1], 1.0)
end

say("\n=== IMEX3D step 6: the Schur stage solve against the five-field one, $NR rank(s) ===")

@testset "Schur stage solve ($NR ranks)" begin

    dst_full  = imex3d_solve!(zeros(5*npoin), src, params_full,  GDT)
    dst_schur = imex3d_solve_schur!(zeros(5*npoin), src, params_schur, GDT)

    @testset "the two stage solves agree" begin
        names = ("rho", "rho_u", "rho_v", "rho_w", "rho_theta")
        worst = 0.0
        for ieq = 1:5
            o = (ieq-1)*npoin
            a = view(dst_full,  o+1:o+npoin)
            b = view(dst_schur, o+1:o+npoin)
            # against the DEVIATION, not the total: rho is O(1) and its
            # deviation is O(1e-2), so a relative error measured against the
            # total would hide a 100% error in the thing actually being solved.
            dev = maximum(abs, a .- view(qe, :, ieq))
            dev = MPI.Allreduce(dev, MPI.MAX, COMM)
            e   = MPI.Allreduce(maximum(abs, a .- b), MPI.MAX, COMM)
            rel = e / max(dev, eps())
            worst = max(worst, rel)
            say(@sprintf("  %-10s |full - schur| = %.3e  against a deviation of %.3e  -> %.2e",
                         names[ieq], e, dev, rel))
        end
        # Two independent Krylov solves to 1e-12; agreement is limited by that,
        # not by the algebra, which step 4 pinned at 1e-14 densely.
        @test worst < 1e-7
    end

    @testset "the scalar solve really did the work" begin
        # Guards against the comparison passing for the wrong reason -- e.g. if
        # the Schur path silently fell through to the five-field one, or if the
        # right-hand side were zero and both returned qe.
        say(@sprintf("  scalar GMRES: %d iterations to %.2e;  five-field: %d to %.2e",
                     sch.ws.last_iters, sch.ws.last_relres,
                     ws5.last_iters, ws5.last_relres))
        @test sch.ws.nsolve == 1
        @test sch.ws.last_iters >= 1
        @test ws5.last_iters >= 1
    end

    @testset "schur_H! does not touch what the recovery reads" begin
        # THE INVARIANT the ordering in imex3d_solve_schur! depends on, tested
        # DIRECTLY: setup_rhs leaves q_b and m_b in the state, and the whole
        # Krylov iteration runs between there and the recovery.
        #
        # The first version of this testset re-ran schur_recover! after the
        # solve and compared it against `dst - qe`. That proved nothing twice
        # over: a clobber during the iteration would corrupt BOTH recoveries
        # equally and cancel, and `(X + qe) - qe` is not exact in floating
        # point -- with qe ~ 360 for rho_theta it loses ~5.7e-14, which is
        # precisely the 2.841e-14 that "failure" reported. It was measuring
        # its own cancellation.
        st2 = SchurState(npoin, op.nimp)
        B2  = zeros(npoin, op.nimp)
        @inbounds for (q, ieq) in enumerate(op.vars), ip = 1:npoin
            B2[ip, q] = src[(ieq-1)*npoin + ip] - qe[ip, ieq]
        end
        schur_setup_rhs!(zeros(npoin), st2, B2, params0, op, GDT)
        keep = (copy(st2.qb), copy(st2.mbx), copy(st2.mby), copy(st2.mbz),
                copy(st2.cx), copy(st2.cy), copy(st2.cz))
        # as many applications as the solve above ran, then some
        P = [sinpi(1e-2*gip[i]) for i = 1:npoin]
        for _ = 1:20
            schur_H!(zeros(npoin), P, params0, op, GDT, st2)
        end
        now = (st2.qb, st2.mbx, st2.mby, st2.mbz, st2.cx, st2.cy, st2.cz)
        e = maximum(maximum(abs, a .- b) for (a, b) in zip(keep, now))
        e = MPI.Allreduce(e, MPI.MAX, COMM)
        say(@sprintf("  q_b, m_b and c after 20 applications of H: max change %.3e", e))
        # BITWISE unchanged, not "small": schur_H! writes only the transient
        # m/r/d slots, so anything non-zero here is a genuine aliasing bug.
        @test e == 0.0
    end

    @testset "relinearisation refreshes what it derives" begin
        # ark_relinearize! rewrites thetabar from the state. grad(thetabar) is
        # DERIVED from thetabar and the Schur preconditioner's band is assembled
        # FROM both, so a refresh that stops at thetabar leaves the advective row
        # differentiating the setup-time state and the preconditioner probing an
        # operator that no longer exists. Neither is loud: the solve still
        # converges, to a slightly wrong splitting.
        tb0 = copy(op.thetabar); dz0 = copy(op.dtbdz)
        nref0 = sch.pc === nothing ? 0 : sch.pc.nrefactor[]
        ark_relinearize!(params_schur, params_schur.imex, src)
        dtb = MPI.Allreduce(maximum(abs, op.thetabar .- tb0), MPI.MAX, COMM)
        ddz = MPI.Allreduce(maximum(abs, op.dtbdz   .- dz0), MPI.MAX, COMM)
        say(@sprintf("  thetabar moved %.3e, grad(thetabar) moved %.3e, precond rebuilds %d",
                     dtb, ddz, sch.pc === nothing ? 0 : sch.pc.nrefactor[] - nref0))
        @test dtb > 0                       # the state really is off-reference
        @test ddz > 0                       # ... and the gradient followed it
        sch.pc === nothing || @test sch.pc.nrefactor[] == nref0 + 1

        # And the two paths still agree on the RELINEARISED operator, which is
        # the property that matters: the refresh must leave a consistent system,
        # not merely a changed one.
        a = imex3d_solve!(zeros(5*npoin), src, params_full,  GDT)
        b = imex3d_solve_schur!(zeros(5*npoin), src, params_schur, GDT)
        e = MPI.Allreduce(maximum(abs, a .- b), MPI.MAX, COMM)
        dev = MPI.Allreduce(maximum(abs, a .- vec(qe)), MPI.MAX, COMM)
        say(@sprintf("  after relinearisation: |full - schur| = %.3e (deviation %.3e)",
                     e, dev))
        @test e < 1e-7 * dev
    end
end

MPI.Barrier(COMM)
say("=== step 6 done on $NR rank(s) ===\n")

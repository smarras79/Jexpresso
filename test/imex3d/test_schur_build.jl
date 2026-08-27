#=============================================================================
 test_schur_build.jl -- does `build_imex3d` wire the two stage solves up the
 way it says it does?

 WHY THIS FILE EXISTS AT ALL
 ---------------------------
 `build_imex3d` had NO test. Every other IMEX3D test assembles the pieces by
 hand and constructs an `IMEX3DCache` directly, which is the right way to test
 the solve but leaves the 300 lines that decide WHICH solve runs, with WHICH
 operator and WHICH preconditioner, covered by nothing. `:imex_schur` is
 exactly a change to those 300 lines, so it needed them exercised.

 What it checks is the wiring, not the numerics -- the numerics are
 test_schur_stage.jl's job:

   flag off  ->  imex3d_solve!        five-field preconditioner built
                                      flux-form Theta row
   flag on   ->  imex3d_solve_schur!  NO five-field preconditioner
                                      advective Theta row, scalar precond

 The "no five-field preconditioner" half is the one worth asserting rather than
 eyeballing. Building it when the Schur path is on is not a correctness bug --
 nothing would apply it -- so it would never show up as a wrong answer. It would
 show up as setup time and resident memory on a 256-rank job, which is the
 place it is most expensive to notice.

     julia --project=<env> test/imex3d/test_schur_build.jl
     mpiexecjl -n 3 julia --project=<env> test/imex3d/test_schur_build.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "..", "hevi", "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl", "acoustic.jl",
          "ark.jl", "hevi.jl", "cfl_diagnostics.jl", "krylov.jl", "precond_api.jl",
          "schur.jl", "schur_precond.jl", "schur_stage.jl", "imex3d.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

params = build_mock_params(nelx=2, nely=2, nelz=10, p=2, Lx=800.0, Ly=800.0,
                           Lz=1000.0, comm=COMM, base=:stratified)
# The CFL diagnostics index params.inputs directly rather than `get` it.
params.inputs[:lvisc] = false
params.inputs[:lsgs]  = false

npoin = Int(params.mesh.npoin)
base  = Dict{Symbol,Any}(:ode_solver => IMEX_ARK(:ARS343), :Δt => 0.5,
                         :imex_verify => false, :imex_rtol => 1.0e-8,
                         :lvisc => false, :lsgs => false,
                         :SOL_VARS_TYPE => TOTAL())

say("\n=== IMEX3D step 6: build_imex3d wiring, $NR rank(s) ===")

@testset "build_imex3d wiring ($NR ranks)" begin

    @testset "the flag off leaves the five-field path exactly as it was" begin
        imex = build_imex3d(params, merge(base, Dict(:imex_schur => false)))
        @test imex.schur === nothing
        @test imex.pc !== nothing                 # five-field precond still built
        @test imex.solve! === imex3d_solve!
        @test imex.op.theta_advective == false    # flux form, the default
        say("  flag off: solve! = imex3d_solve!, pc built, flux-form Theta row")
    end

    @testset "the flag on switches the solve, the row and the preconditioner" begin
        imex = build_imex3d(params, merge(base, Dict(:imex_schur => true)))
        @test imex.schur !== nothing
        @test imex.solve! === imex3d_solve_schur!
        # :imex_schur IMPLIES the advective row -- the reduction does not close
        # on one scalar without it, so the flag sets it rather than trusting a
        # second deck key to agree.
        @test imex.op.theta_advective == true
        # ... and the five-field preconditioner is NOT built: nothing applies it
        # here, and it is nimp times wider and nimp times longer than the one
        # that replaces it.
        @test imex.pc === nothing
        @test imex.schur.pc !== nothing
        # The scalar band is 2(ngl-1) either side, one field, one row per level.
        @test imex.schur.pc.kl == 2 * (Int(params.mesh.ngl) - 1)
        @test imex.schur.pc.ku == imex.schur.pc.kl
        @test imex.schur.pc.n  == imex.topo.nlev
        say(@sprintf("  flag on: solve! = imex3d_solve_schur!, no five-field pc, band %d, n %d",
                     imex.schur.pc.kl, imex.schur.pc.n))
    end

    @testset "the DEFAULT self-checks pass with the Schur path on" begin
        # `:imex_verify` defaults to TRUE, so this is what a production deck
        # hits on its first run -- and it did not pass when this was written.
        # The self-check builds a throw-away VERTICAL operator and asserts the
        # 3D one reduces to it on a horizontally uniform field, at 1e-6. Built
        # with the flux Theta row while `op` carries the advective one, it was
        # measuring the CONSISTENCY error between the two rows (2.79e-04, the
        # same 0.06%-of-flux that test_theta_advective.jl reports) and calling
        # it "a metric or index error in the xi/eta sweeps". build_imex3d threw,
        # so a correct operator aborted the run at setup.
        #
        # build_imex3d throws on a failed self-check, so reaching the assertion
        # below IS the test.
        imex = build_imex3d(params, merge(base, Dict(:imex_schur => true,
                                                     :imex_verify => true)))
        @test imex.schur !== nothing
        say("  self-checks pass with :imex_verify => true (the default)")
    end

    @testset "the cache it built actually solves" begin
        # Structural assertions can all pass on a cache that throws the moment
        # it is used, so drive the hook it installed once, through params.imex
        # the way the integrator reaches it.
        imex = build_imex3d(params, merge(base, Dict(:imex_schur => true)))
        p    = merge(params, (imex = imex,))
        qe   = params.qp.qe
        gip  = params.mesh.ip2gip
        src  = zeros(Float64, 5 * npoin)
        @inbounds for ieq = 1:5, ip = 1:npoin
            src[(ieq-1)*npoin + ip] = qe[ip, ieq] +
                0.01 * sinpi(1e-2 * (gip[ip] + 3ieq)) * max(qe[ip, 1], 1.0)
        end
        dst = imex.solve!(zeros(5*npoin), src, p, imex.gdt)
        @test all(isfinite, dst)
        @test imex.schur.ws.nsolve == 1
        moved = MPI.Allreduce(maximum(abs, dst .- src), MPI.MAX, COMM)
        say(@sprintf("  one stage solve: %d iterations, state moved %.3e",
                     imex.schur.ws.last_iters, moved))
        # It must have DONE something -- an implicit solve that returns its
        # input is the failure mode a finiteness check sails straight past.
        @test moved > 0
    end

    @testset "the :imex_schur + :implicit_vdiff contradiction repairs itself" begin
        # THE PAIR IS NOT A TASTE QUESTION. The Schur reduction eliminates the
        # momentum and Theta rows assuming their coefficient is the identity;
        # implicit vertical diffusion makes it (I - lam*d/dz(mu d/dz)) and the
        # elimination drops the operator. The term is still in f_imp, so the
        # explicit half has it subtracted and the stage solve never puts it
        # back -- an error that is ~0 on a laminar initial state and grows with
        # the boundary layer. That is why this is repaired rather than left to
        # a deck author to notice tens of seconds into a run.
        #
        # Tested on the resolver directly rather than through build_imex3d:
        # the decision is a pure function of two Bools and the deck, and
        # reaching it through a real :implicit_vdiff build would drag in the
        # closure and the vertical-diffusion coefficients for no extra coverage.

        # 1. the contradiction: Schur loses, diffusion is kept.
        IMEX_SCHUR_DEMOTED[] = false        # so the assertion below means something
        d = Dict{Symbol,Any}(:imex_schur => true)
        @test imex_resolve_schur_vdiff!(d, true, true, RANK) == false
        # ... and the deck itself is corrected, so nothing downstream can read
        # the value the run is no longer using.
        @test d[:imex_schur] == false
        # The flag imex3d_report reads to repeat the adjustment further down a
        # log that will be days long by the time anyone reads it.
        @test IMEX_SCHUR_DEMOTED[]

        # 2. every legal combination is left exactly alone -- no banner, no
        #    edit to the deck. A guard that fires on a correct deck is worse
        #    than no guard.
        for (a, b) in ((true, false), (false, true), (false, false))
            e = Dict{Symbol,Any}()
            @test imex_resolve_schur_vdiff!(e, a, b, RANK) == a
            @test isempty(e)
        end

        # 3. the escape hatch runs the pair as asked. It exists so the failure
        #    can be MEASURED without editing source, which is the same reason
        #    :imex_schur_kernel => false exists.
        o = Dict{Symbol,Any}(:imex_allow_schur_vdiff => true)
        @test imex_resolve_schur_vdiff!(o, true, true, RANK) == true
        @test !haskey(o, :imex_schur)      # untouched: the deck asked for this

        say("  :imex_schur + :implicit_vdiff -> :imex_schur => false, diffusion kept")
    end
end

MPI.Barrier(COMM)
say("=== build_imex3d wiring done on $NR rank(s) ===\n")

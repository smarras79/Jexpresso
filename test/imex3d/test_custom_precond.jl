#=============================================================================
 test_custom_precond.jl -- can a deck supply its OWN stage-solve
                           preconditioner without editing src?

 WHAT THIS PINS
 --------------
 `:imex_precond => :custom` plus `:imex_precond_build => f` is the whole
 extension point (src/kernel/solvers/hevi/precond_api.jl). Three claims are
 made for it, and each of them is a claim about behaviour, not about wiring:

   1. IT IS REACHED. A custom preconditioner is built, and it is what the
      stage solve applies -- not the built-in one, not the identity. Asserted
      by counting applications inside the object itself, because a hook that is
      built and never called is exactly what "the deck seemed to have no
      effect" looks like from outside.

   2. IT IS EQUIVALENT TO THE BUILT-IN PATH. Handing back, through the hook,
      the very column preconditioner `:imex_precond => :column` builds
      internally must reproduce that arm's ITERATION COUNT exactly. Anything
      less means the hook applies M differently from the way build_imex3d does
      -- at a different point in the cycle, to a different field, or with a
      stale gdt -- and every one of those is invisible in the answer.

   3. IT CANNOT CHANGE THE ANSWER. A preconditioner is a change of basis for
      the iteration and nothing more, so a deliberately BAD one (scalar
      Jacobi on a stiff acoustic operator) must still converge to the same
      state as the good one, to the solve tolerance. If it does not, the hook
      is being applied somewhere it changes the system rather than its
      conditioning.

 Both arms are covered: the five-field solve and, with `:imex_schur`, the
 scalar one, where the field is npoin x 1 and the operator is H rather than
 I - gdt*A.

 The FAILURE modes get their own set. `imex_precond_selfcheck` exists to turn
 the out-of-place mistake -- computing M^-1 V into a fresh array, which acts as
 the identity and costs a whole run before anyone notices -- into a setup
 error, and that is only worth having if it actually fires.

     julia --project=<env> test/imex3d/test_custom_precond.jl
     mpiexecjl -n 3 julia --project=<env> test/imex3d/test_custom_precond.jl
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

#-----------------------------------------------------------------------------
# Three preconditioners written the way a user would write them: a type, one
# `imex_precond_apply!` method, and a builder that takes the context.
#-----------------------------------------------------------------------------

"""
    WrappedColumn

The built-in column preconditioner, reached through the hook instead of
through `build_imex3d`. Its whole purpose is claim 2 above: same M, same
application point, so it must produce the same iteration count.

`napply` is the evidence for claim 1.
"""
mutable struct WrappedColumn{P}
    inner::P
    napply::Int
end

imex_precond_apply!(pc::WrappedColumn, V::AbstractMatrix, params, gdt::Real) =
    (pc.napply += 1; imex_precond_apply!(pc.inner, V, params, gdt))
imex_precond_describe(pc::WrappedColumn) =
    string("the built-in column solve, wrapped (", imex_precond_describe(pc.inner), ")")

function build_wrapped(ctx)
    inner = ctx.schur ?
        build_schur_column_precond(ctx.params, ctx.topo, ctx.comm, ctx.gdt;
                                   lwall_flux = ctx.lwall_flux) :
        begin
            pvars = hevi_choose_vars(ctx.params.metrics, ctx.comm)
            opv = build_hevi_operator(ctx.params, ctx.topo, pvars;
                                      lwall_flux = ctx.lwall_flux, full = false)
            owner, own = assign_column_owners(ctx.topo, ctx.comm)
            cc  = build_column_comm(ctx.topo, owner, own, ctx.comm, length(pvars))
            fac = build_column_factorization(ctx.params, opv, cc, ctx.topo, ctx.gdt)
            IMEX3DPrecond(opv, cc, fac, ctx.topo, copy(pvars))
        end
    return WrappedColumn(inner, 0)
end

"""
    DiagJacobi

A deliberately WEAK preconditioner: the reciprocal of the operator's own
diagonal, recovered by probing `M` with a field of ones. Node-local, so it is
trivially single-valued across ranks and needs no communication -- the
cheapest legal thing anyone would write first.

It is here to be bad. Claim 3 is that a bad preconditioner still lands on the
same answer.
"""
mutable struct DiagJacobi
    invd::Matrix{Float64}
    napply::Int
end

imex_precond_apply!(pc::DiagJacobi, V::AbstractMatrix, params, gdt::Real) =
    (pc.napply += 1; V .*= pc.invd; V)
imex_precond_describe(::DiagJacobi) = "diagonal Jacobi (test)"

function build_diag_jacobi(ctx)
    npoin = Int(ctx.params.mesh.npoin)
    ones_ = ones(Float64, npoin, ctx.nimp)
    out   = zeros(Float64, npoin, ctx.nimp)
    if ctx.schur
        schur_H!(vec(out), vec(ones_), ctx.params, ctx.op, ctx.gdt,
                 SchurState(npoin, ctx.op.nimp))
    else
        hevi_apply_A!(out, ones_, ctx.params, ctx.op)
        @. out = ones_ - ctx.gdt * out
    end
    # A row sum is not a diagonal, so this is an approximation of an
    # approximation -- which is the point. Guard the reciprocal: a row that
    # sums to zero would hand the iteration an Inf.
    invd = similar(out)
    @inbounds for i in eachindex(out)
        invd[i] = abs(out[i]) > 1.0e-12 ? 1.0 / out[i] : 1.0
    end
    return DiagJacobi(invd, 0)
end

"""
    OutOfPlace

The mistake `imex_precond_selfcheck` exists for: correct arithmetic, returned
as a NEW array, so `V` is untouched and the iteration silently preconditions
with the identity.
"""
struct OutOfPlace end
imex_precond_apply!(::OutOfPlace, V::AbstractMatrix, params, gdt::Real) = V ./ 2

struct NoMethodPrecond end     # no imex_precond_apply! at all

#-----------------------------------------------------------------------------

params = build_mock_params(nelx=2, nely=2, nelz=10, p=2, Lx=800.0, Ly=800.0,
                           Lz=1000.0, comm=COMM, base=:stratified)
params.inputs[:lvisc] = false
params.inputs[:lsgs]  = false

npoin = Int(params.mesh.npoin)
base  = Dict{Symbol,Any}(:ode_solver => IMEX_ARK(:ARS343), :Δt => 0.5,
                         :imex_verify => false, :imex_rtol => 1.0e-10,
                         :imex_maxiter => 400, :imex_restart => 40,
                         :lvisc => false, :lsgs => false,
                         :SOL_VARS_TYPE => TOTAL())

deck(d...) = merge(base, Dict{Symbol,Any}(d...))

"""
    stage_solve(inputs) -> (dst, imex)

Build a cache from `inputs` and drive ONE stage solve through the hook the
integrator uses, on a deterministic perturbation of the reference state.

Deterministic and built from the GLOBAL point id, so every rank holding a
shared node feeds in the same value -- see imex_verify_solve for why `rand()`
would break this on more than one rank rather than testing anything.
"""
function stage_solve(inputs)
    imex = build_imex3d(params, inputs)
    p    = merge(params, (imex = imex,))
    qe   = params.qp.qe
    gip  = params.mesh.ip2gip
    src  = zeros(Float64, 5 * npoin)
    @inbounds for ieq = 1:5, ip = 1:npoin
        src[(ieq-1)*npoin + ip] = qe[ip, ieq] +
            0.01 * sinpi(1e-2 * (gip[ip] + 3ieq)) * max(qe[ip, 1], 1.0)
    end
    return imex.solve!(zeros(5*npoin), src, p, imex.gdt), imex
end

gmax(x) = NR > 1 ? MPI.Allreduce(x, MPI.MAX, COMM) : x

say("\n=== IMEX3D: the custom-preconditioner hook, $NR rank(s) ===")

@testset "custom stage-solve preconditioner ($NR ranks)" begin

  for schur in (false, true)
    arm = schur ? "Schur (scalar)" : "five-field"
    @testset "$arm" begin

        # -- the reference arms, built the ordinary way ----------------------
        d_col, imex_col = stage_solve(deck(:imex_schur => schur,
                                           :imex_precond => :column))
        it_col = (schur ? imex_col.schur.ws : imex_col.ws).last_iters
        d_non, imex_non = stage_solve(deck(:imex_schur => schur,
                                           :imex_precond => :none))
        it_non = (schur ? imex_non.schur.ws : imex_non.ws).last_iters

        # -- 1 & 2: the hook is reached, and it is the same M ----------------
        d_wrp, imex_wrp = stage_solve(deck(:imex_schur => schur,
                                           :imex_precond => :custom,
                                           :imex_precond_build => build_wrapped))
        pc_wrp = schur ? imex_wrp.schur.pc : imex_wrp.pc
        ws_wrp = schur ? imex_wrp.schur.ws : imex_wrp.ws
        @test pc_wrp isa WrappedColumn
        # Built AND applied. The self-check at build time is one of them; the
        # Arnoldi loop and the end-of-cycle update are the rest.
        @test pc_wrp.napply > ws_wrp.last_iters
        # THE assertion: same M, applied the same way, so the same descent.
        @test ws_wrp.last_iters == it_col
        @test gmax(maximum(abs, d_wrp .- d_col)) < 1.0e-11 * gmax(maximum(abs, d_col))
        say(@sprintf("  %-14s :column %3d iters | wrapped through the hook %3d iters (%d applications)",
                     arm, it_col, ws_wrp.last_iters, pc_wrp.napply))

        # -- 3: a WEAK preconditioner changes the cost, not the answer -------
        d_jac, imex_jac = stage_solve(deck(:imex_schur => schur,
                                           :imex_precond => :custom,
                                           :imex_precond_build => build_diag_jacobi))
        pc_jac = schur ? imex_jac.schur.pc : imex_jac.pc
        ws_jac = schur ? imex_jac.schur.ws : imex_jac.ws
        @test pc_jac isa DiagJacobi
        @test pc_jac.napply > 0
        @test ws_jac.last_relres <= 1.0e-10
        # Same state as the column arm, to the tolerance BOTH solves converged
        # to. This is the claim that a preconditioner cannot change the answer,
        # and it is the one that would fail if the hook were applied as a
        # multiplication on the operator rather than on the search direction.
        scl = gmax(maximum(abs, d_col))
        @test gmax(maximum(abs, d_jac .- d_col)) < 1.0e-8 * scl
        say(@sprintf("  %-14s :none   %3d iters | diagonal Jacobi          %3d iters, |Δstate| %.2e",
                     arm, it_non, ws_jac.last_iters, gmax(maximum(abs, d_jac .- d_col)) / scl))

        # -- the report names what actually ran ------------------------------
        @test occursin("wrapped", imex_precond_describe(pc_wrp))
    end
  end

  @testset "the mistakes are refused at setup" begin
    # Every one of these is silent at run time: the run converges to the right
    # answer and only the cost is wrong, which is why they are errors here.
    @test_throws ErrorException build_imex3d(params,
        deck(:imex_precond => :custom))                       # no builder
    @test_throws ErrorException build_imex3d(params,
        deck(:imex_precond => :custom, :imex_precond_build => 42))
    @test_throws ErrorException build_imex3d(params,
        deck(:imex_precond => :custom, :imex_precond_build => (_ -> nothing)))
    @test_throws ErrorException build_imex3d(params,
        deck(:imex_precond => :custom, :imex_precond_build => (_ -> OutOfPlace())))
    @test_throws ErrorException build_imex3d(params,
        deck(:imex_precond => :custom, :imex_precond_build => (_ -> NoMethodPrecond())))
    # ... and an unknown mode is still refused, rather than falling back.
    @test_throws ErrorException build_imex3d(params, deck(:imex_precond => :blockjacobi))
    say("  a missing builder, a non-callable one, `nothing`, an out-of-place")
    say("  apply, a type with no method, and an unknown mode: all refused at setup")
  end
end

MPI.Barrier(COMM)
say("=== custom-preconditioner hook done on $NR rank(s) ===\n")

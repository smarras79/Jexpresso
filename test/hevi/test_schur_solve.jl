#=============================================================================
 test_schur_solve.jl -- STEP 4 of the Schur complement work.

 The reduction is exact algebra on the discrete blocks, so there is exactly one
 thing to prove: solving the SCALAR Schur system and recovering the five fields
 must give the SAME answer as solving the full 5-field stage system.

 Not "close" -- the same, to the accuracy of the two solves. If it is not, the
 elimination is wrong somewhere, and a wrong elimination produces a solver that
 converges beautifully to the wrong state.

 Both systems are formed DENSELY here and solved directly. That removes the
 Krylov method from the comparison entirely: any difference is the algebra,
 not the iteration. Only tractable on the mock mesh, which is what it is for.

     julia --project=<env> test/hevi/test_schur_solve.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_schur_solve.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
for f in ("columns.jl","vdiffusion.jl","operator.jl","schur.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
params = build_mock_params(nelx=NELX, nely=NELY, nelz=NELZ, p=P, Lz=1000.0,
                           comm=COMM, base=:stratified)
topo = build_column_topology(params.mesh, COMM)
# ADVECTIVE: the reduction closes on one scalar only with this row.
op = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux=true, full=true,
                         theta_advective=true)
const N = Int(params.mesh.npoin)
const S = op.slot
const LAM = 0.05
st = SchurState(N, op.nimp)

say("\n=== Schur step 4: the reduction against the full solve, $NR rank(s) ===")

@testset "Schur reduction ($NR ranks)" begin

    @testset "H is linear" begin
        p1 = [sinpi(1e-2*(params.mesh.ip2gip[i]+3)) for i=1:N]
        p2 = [cospi(2e-2*params.mesh.ip2gip[i]) for i=1:N]
        h1 = schur_H!(zeros(N), p1, params, op, LAM, st)
        h2 = schur_H!(zeros(N), p2, params, op, LAM, st)
        hs = schur_H!(zeros(N), p1 .+ p2, params, op, LAM, st)
        @test maximum(abs, hs .- (h1 .+ h2)) < 1e-10 * maximum(abs, h1)
    end

if NR == 1
    @testset "THE TEST: same answer as the full 5-field solve" begin
        # full system (I - lam*A) U = B, formed densely
        m = N * op.nimp
        M = zeros(Float64, m, m)
        V = zeros(Float64, N, op.nimp); W = zeros(Float64, N, op.nimp)
        for j = 1:m
            fill!(V, 0.0); V[j] = 1.0
            hevi_apply_A!(W, V, params, op)
            @inbounds for i = 1:m; M[i,j] = (i==j ? 1.0 : 0.0) - LAM*W[i]; end
        end
        B = reshape([sinpi(1e-2*(params.mesh.ip2gip[mod1(k,N)] + 7.0*fld(k-1,N))) for k=1:m], N, op.nimp)
        Ufull = reshape(M \ vec(B), N, op.nimp)

        # Schur route: scalar H, solved densely, then recover the five fields
        rhs = schur_setup_rhs!(zeros(N), st, B, params, op, LAM)
        Hm = zeros(Float64, N, N); e = zeros(N)
        for j = 1:N
            fill!(e, 0.0); e[j] = 1.0
            Hm[:, j] = schur_H!(zeros(N), e, params, op, LAM, st)
        end
        Pv = Hm \ rhs
        # setup_rhs must run again: schur_H! above overwrote the cached m_b
        schur_setup_rhs!(zeros(N), st, B, params, op, LAM)
        Uschur = schur_recover!(zeros(N, op.nimp), Pv, st, B, params, op, LAM)

        for (c, nm) in ((S[1],"rho"), (S[2],"rho_u"), (S[3],"rho_v"),
                        (S[4],"rho_w"), (S[5],"rho_theta"))
            d = maximum(abs, Ufull[:,c] .- Uschur[:,c])
            r = d / max(maximum(abs, Ufull[:,c]), eps())
            say(@sprintf("  %-10s max|full - schur| = %.3e   (relative %.3e)", nm, d, r))
            @test r < 1e-10
        end
        say(@sprintf("  H is %d x %d against the full %d x %d", N, N, m, m))
    end
end
end

MPI.Barrier(COMM)
say("=== step 4 done on $NR rank(s) ===\n")

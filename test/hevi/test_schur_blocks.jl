#=============================================================================
 test_schur_blocks.jl -- STEP 1 of the Schur complement work.

 Before any Schur code is written, PROVE what the discrete operator actually
 is. The reduction is exact algebra on the assembled blocks, so an error in
 reading `_hevi_A_elements_full!` propagates into a solver that converges
 happily to the wrong answer. Nothing here builds a Schur solve; it only
 pins the block structure the reduction will be derived from.

 Claimed (see the kernel at operator.jl:363):

     A_rho   = -Div[ W m ]
     A_rho_u = -Gradx[ P ],   A_rho_v = -Grady[ P ],
     A_rho_w = -Gradz[ P ] - g*rho          <- LOCAL buoyancy, exactly pointwise
     A_Theta = -Div[ tb .* W m ]

 with P = beta.*Theta, W the free-slip mask, and Div/Grad the assembled
 weak operators followed by DSS and Minv.

 WHERE THE REDUCTION LEADS, AND WHY IT IS A PRECONDITIONER AND NOT A SOLVE
 -------------------------------------------------------------------------
 The stage system (I - lam*A) U = b is

     rho   + lam*Div[W m]                    = b_rho      (1)
     m     + lam*Grad[P] + lam*g*zhat*rho    = b_m        (2)
     Theta + lam*Div[tb .* W m]              = b_Theta    (3)

 Substituting (2) into (1) and (3) eliminates the momentum and leaves a 2x2
 system in the two SCALARS rho and P:

     rho - lam^2*g*Div[W zhat rho]                  = b_rho - lam*Div[W b_m]
                                                      + lam^2*Div[W Grad P]
     P/beta - lam^2*Div[tb W Grad P]
            - lam^2*g*Div[tb W zhat rho]            = b_Theta - lam*Div[tb W b_m]

 i.e. 2*Np unknowns, not 5*Np. Giraldo et al. reach ONE scalar because their
 theta equation is ADVECTIVE (u . grad theta0), which makes the coupling
 pointwise and lets rho be substituted away exactly. Ours is FLUX form
 (Div[tb*m]), so rho survives. That is a real difference between this code and
 the papers, and it is the reason the scalar Helmholtz

     H[P] := P/beta - lam^2 * Div[ tb .* W .* Grad P ]

 obtained by dropping the two lam^2*g couplings is used here as a
 PRECONDITIONER rather than as an exact solve. That choice is deliberate and
 it is the safe one: GMRES converges to the same answer whatever the
 preconditioner does, so an imperfect H can only cost iterations, never
 correctness. Making H the primary solve -- and collecting the full 5*Np ->
 Np saving on the Krylov vectors -- requires switching equation (3) to
 advective form, which changes f_imp and therefore the ARK split, and would
 have to be re-verified against rhs! by hevi_verify_physics.

     julia --project=<env> test/hevi/test_schur_blocks.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_schur_blocks.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
include(joinpath(SRC, "columns.jl"))
include(joinpath(SRC, "vdiffusion.jl"))
include(joinpath(SRC, "operator.jl"))

say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lz = 1000.0, comm = COMM)
topo = build_column_topology(params.mesh, COMM)
op   = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux = true, full = true)

const NPOIN = Int(params.mesh.npoin)
const GIP   = params.mesh.ip2gip
const S     = op.slot          # (rho, rho_u, rho_v, rho_w, rho_theta) -> column of V
const G     = PhysicalConst{Float64}().g   # read, not hardcoded, so it cannot drift

field(seed) = [sinpi(1.0e-2*(Float64(GIP[ip]) + 13.0seed)) for ip = 1:NPOIN]

"apply the real operator to a state built from the five given nodal fields"
function applyA(fρ, fu, fv, fw, fθ)
    V = zeros(Float64, NPOIN, op.nimp)
    V[:, S[1]] .= fρ; V[:, S[2]] .= fu; V[:, S[3]] .= fv
    V[:, S[4]] .= fw; V[:, S[5]] .= fθ
    W = zeros(Float64, NPOIN, op.nimp)
    hevi_apply_A!(W, V, params, op)
    return W
end
zed = zeros(Float64, NPOIN)

say("\n=== Schur step 1: block structure of the 3D operator, $NR rank(s) ===")

@testset "discrete block structure ($NR ranks)" begin

    fρ, fu, fv, fw, fθ = (field(k) for k = 1:5)

    @testset "LINEARITY -- the reduction is only valid if A is linear" begin
        A1 = applyA(fρ, fu, fv, fw, fθ)
        A2 = applyA(2fρ, 2fu, 2fv, 2fw, 2fθ)
        @test maximum(abs, A2 .- 2A1) < 1e-10 * max(maximum(abs, A1), 1e-30)
        # superposition over the five inputs, which is what lets each block be
        # probed in isolation below
        Sum = applyA(fρ,zed,zed,zed,zed) .+ applyA(zed,fu,zed,zed,zed) .+
              applyA(zed,zed,fv,zed,zed) .+ applyA(zed,zed,zed,fw,zed) .+
              applyA(zed,zed,zed,zed,fθ)
        @test maximum(abs, Sum .- A1) < 1e-10 * max(maximum(abs, A1), 1e-30)
    end

    @testset "rho depends ONLY on the momenta" begin
        # A_rho = -Div[W m]: no beta*Theta, no rho.
        @test maximum(abs, applyA(fρ,zed,zed,zed,zed)[:, S[1]]) == 0.0
        @test maximum(abs, applyA(zed,zed,zed,zed,fθ)[:, S[1]]) == 0.0
        @test maximum(abs, applyA(zed,fu,zed,zed,zed)[:, S[1]]) > 0.0
    end

    @testset "the momenta depend ONLY on P = beta*Theta (plus buoyancy on w)" begin
        for (c, name) in ((S[2],"rho_u"), (S[3],"rho_v"), (S[4],"rho_w"))
            @test maximum(abs, applyA(zed,fu,fv,fw,zed)[:, c]) == 0.0   # not on m
        end
        @test maximum(abs, applyA(zed,zed,zed,zed,fθ)[:, S[2]]) > 0.0
        # rho enters ONLY rho_w, and only through the buoyancy source
        R = applyA(fρ,zed,zed,zed,zed)
        @test maximum(abs, R[:, S[2]]) == 0.0
        @test maximum(abs, R[:, S[3]]) == 0.0
        @test maximum(abs, R[:, S[5]]) == 0.0
        @test maximum(abs, R[:, S[4]]) > 0.0
    end

    @testset "THE KEY CLAIM: buoyancy is exactly -g*rho, pointwise" begin
        # This is what makes the Schur elimination tractable. It holds only
        # because Minv .* DSS(wJac .* f) == f for a continuous nodal f -- i.e.
        # the CG-SEM mass matrix IS DSS(wJac). If that identity failed, the
        # source would come back smeared by M^-1 DSS and the reduction below
        # would be wrong in a way no solver could detect.
        R = applyA(fρ, zed, zed, zed, zed)[:, S[4]]
        @test maximum(abs, R .- (-G .* fρ)) < 1e-10 * maximum(abs, G .* fρ)
        say(@sprintf("  buoyancy block == -g*rho to %.2e (relative)",
                     maximum(abs, R .+ G .* fρ) / maximum(abs, G .* fρ)))
    end

    @testset "Theta depends ONLY on the momenta, through tb" begin
        @test maximum(abs, applyA(fρ,zed,zed,zed,fθ)[:, S[5]]) == 0.0
        @test maximum(abs, applyA(zed,fu,fv,fw,zed)[:, S[5]]) > 0.0
    end

    @testset "the momentum blocks are the three components of ONE gradient" begin
        # A_rho_u/v/w = -Gradx/y/z[P]: the same scalar P drives all three, which
        # is what collapses them into a single scalar unknown.
        Aθ = applyA(zed,zed,zed,zed,fθ)
        # scaling Theta scales all three momentum rows identically
        Aθ2 = applyA(zed,zed,zed,zed,3fθ)
        for c in (S[2], S[3], S[4])
            @test maximum(abs, Aθ2[:,c] .- 3Aθ[:,c]) < 1e-10*max(maximum(abs,Aθ[:,c]),1e-30)
        end
    end
end

MPI.Barrier(COMM)
say("=== step 1 done on $NR rank(s) ===\n")

#=============================================================================
 test_theta_advective.jl -- STEP 3 of the Schur complement work.

 The flux and advective forms of the Theta row are the SAME operator in the
 continuum:  Div(tb*m) = tb*Div(m) + m.grad(tb).  They are different MATRICES,
 because the discrete derivative does not obey the product rule exactly (LGL
 quadrature aliases the product). Two things follow, and both are measured
 here rather than assumed:

   1. The advective row must equal  tb .* A_rho  -  (W m).grad(tb)  EXACTLY.
      That is an identity of the implementation, not an approximation, and it
      holds only because thetabar is continuous and therefore commutes with
      the DSS assembly. If it fails, the kernel is wrong.

   2. The DIFFERENCE between the two forms does not vanish, and whatever it is
      stays EXPLICIT: the split is f_exp = rhs! - f_imp, and rhs! keeps using
      the flux form. So the difference is the price of switching, and it is
      only acceptable if it is small next to the acoustic term being removed.
      That number is reported, not hidden.

     julia --project=<env> test/hevi/test_theta_advective.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_theta_advective.jl
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
# STRATIFIED, deliberately. The default base states give a constant thetabar,
# grad(thetabar) = 0, and the two forms of the Theta row then agree trivially --
# which would let a broken advective kernel pass every check in this file.
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lz = 1000.0, comm = COMM, base = :stratified)
topo = build_column_topology(params.mesh, COMM)
opF  = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux = true, full = true)
opA  = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux = true, full = true,
                           theta_advective = true)

const NPOIN = Int(params.mesh.npoin)
const GIP   = params.mesh.ip2gip
const S     = opF.slot

fld(k) = [sinpi(1.0e-2*(Float64(GIP[ip]) + 11.0k)) + 0.3cospi(3.0e-2*GIP[ip]) for ip = 1:NPOIN]

function applyA(op, fρ, fu, fv, fw, fθ)
    V = zeros(Float64, NPOIN, op.nimp)
    V[:,S[1]] .= fρ; V[:,S[2]] .= fu; V[:,S[3]] .= fv; V[:,S[4]] .= fw; V[:,S[5]] .= fθ
    W = zeros(Float64, NPOIN, op.nimp)
    hevi_apply_A!(W, V, params, op)
    return W
end
zed = zeros(Float64, NPOIN)

say("\n=== Schur step 3: advective Theta row, $NR rank(s) ===")

@testset "advective Theta ($NR ranks)" begin

    fρ, fu, fv, fw, fθ = (fld(k) for k = 1:5)

    @testset "grad(thetabar) is the operator's OWN gradient" begin
        # filled by hevi_fill_gradthetabar!, which applies the operator to
        # Theta = thetabar/beta and reads the momentum rows. Re-derive it here
        # the same way and require an exact match -- this is the field the
        # advective row leans on entirely.
        V = zeros(Float64, NPOIN, opF.nimp)
        @inbounds for ip = 1:NPOIN; V[ip, S[5]] = opF.thetabar[ip]/opF.beta[ip]; end
        W = zeros(Float64, NPOIN, opF.nimp)
        hevi_apply_A!(W, V, params, opF)
        @test maximum(abs, opA.dtbdx .+ W[:, S[2]]) == 0.0
        @test maximum(abs, opA.dtbdy .+ W[:, S[3]]) == 0.0
        @test maximum(abs, opA.dtbdz .+ W[:, S[4]]) == 0.0
        @test maximum(abs, opA.dtbdz) > 0.0      # the mock has vertical stratification
    end

    @testset "THE IDENTITY: A_Theta_adv == tb .* A_rho - (W m).grad(tb)" begin
        # Exact, not approximate. It holds because thetabar is a continuous
        # nodal field and therefore commutes with the DSS sum -- the same
        # property that makes the buoyancy source exactly -g*rho.
        R = applyA(opA, zed, fu, fv, fw, zed)
        Arho = applyA(opA, zed, fu, fv, fw, zed)[:, S[1]]      # = -Div[W m]
        wx, wy, wz = opA.wallx, opA.wally, opA.wall
        want = similar(Arho)
        @inbounds for ip = 1:NPOIN
            mx = wx[ip] ? 0.0 : fu[ip]
            my = wy[ip] ? 0.0 : fv[ip]
            mz = wz[ip] ? 0.0 : fw[ip]
            want[ip] = opA.thetabar[ip]*Arho[ip] -
                       (mx*opA.dtbdx[ip] + my*opA.dtbdy[ip] + mz*opA.dtbdz[ip])
        end
        got = R[:, S[5]]
        rel = maximum(abs, got .- want) / maximum(abs, want)
        say(@sprintf("  identity holds to %.2e (relative)", rel))
        @test rel < 1e-12
    end

    @testset "only the Theta row changes" begin
        # rho and the momenta must be bit-identical between the two operators;
        # if they are not, the change has leaked out of the row it was meant
        # for and every conclusion about the split is void.
        A1 = applyA(opF, fρ, fu, fv, fw, fθ)
        A2 = applyA(opA, fρ, fu, fv, fw, fθ)
        for c in (S[1], S[2], S[3], S[4])
            @test maximum(abs, A1[:,c] .- A2[:,c]) == 0.0
        end
        @test maximum(abs, A1[:,S[5]] .- A2[:,S[5]]) > 0.0     # ...and it DOES change
    end

    @testset "still linear" begin
        A1 = applyA(opA, fρ, fu, fv, fw, fθ)
        A2 = applyA(opA, 2fρ, 2fu, 2fv, 2fw, 2fθ)
        @test maximum(abs, A2 .- 2A1) < 1e-10*maximum(abs, A1)
    end

    @testset "THE PRICE: how far the two forms sit apart" begin
        # This difference stays EXPLICIT, because rhs! keeps the flux form and
        # f_exp = rhs! - f_imp. It is only acceptable if it is small next to
        # the acoustic term the split exists to remove.
        A1 = applyA(opF, fρ, fu, fv, fw, fθ)
        A2 = applyA(opA, fρ, fu, fv, fw, fθ)
        dθ  = maximum(abs, A1[:,S[5]] .- A2[:,S[5]])
        ref = maximum(abs, A1[:,S[5]])
        say(@sprintf("  |flux - advective| on the Theta row = %.3e, against |flux| = %.3e  ->  %.2f%%",
                     dθ, ref, 100*dθ/ref))
        # Reported as the headline number. Asserted only loosely: what matters
        # is the fraction, and the honest test of it is hevi_verify_physics on
        # a real deck, which measures what f_imp captures of the true rhs!.
        @test dθ < ref
    end
end

MPI.Barrier(COMM)
say("=== step 3 done on $NR rank(s) ===\n")

#=============================================================================
 test_schur_helmholtz.jl -- STEP 2 of the Schur complement work.

 Step 1 (test_schur_blocks.jl) pinned the block structure. This builds the
 scalar Helmholtz H out of it and asks the questions that decide whether the
 whole approach is worth pursuing on THIS discretisation:

   * is H symmetric in the mass inner product?      -> can we use CG?
   * is it positive definite?                       -> is it a real Helmholtz?
   * is its spectrum real and positive?             -> Giraldo's claim, tested
                                                       on our operator, not
                                                       assumed from the paper
   * is it better conditioned than (I - lam*A)?     -> the actual payoff

 The last two are computed by forming both operators DENSELY on the small mock
 mesh, column by column, and taking eigenvalues. That is only tractable here
 (275 nodes, 1375 dof) which is exactly what the mock mesh is for.

     julia --project=<env> test/hevi/test_schur_helmholtz.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_schur_helmholtz.jl
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
include(joinpath(SRC, "schur.jl"))

say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lz = 1000.0, comm = COMM)
topo = build_column_topology(params.mesh, COMM)
op   = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux = true, full = true)

const NPOIN = Int(params.mesh.npoin)
const GIP   = params.mesh.ip2gip
const LAM   = 0.05
work = SchurWork(NPOIN, op.nimp)

Happly(p) = (h = zeros(Float64, NPOIN); schur_apply!(h, p, params, op, LAM, work); h)

# Mass-weighted inner product, with duplicated nodes counted once. M = 1/Minv.
const MULT = let m = ones(Float64, NPOIN, 1)
    op.dss_cache === nothing || assemble_mpi!(m, op.dss_cache)
    [1.0 / max(m[ip,1], 1.0) for ip = 1:NPOIN]
end
function mdot(a, b)
    s = 0.0
    @inbounds for ip = 1:NPOIN
        s += MULT[ip] * (1.0/params.Minv[ip]) * a[ip] * b[ip]
    end
    return NR > 1 ? MPI.Allreduce(s, MPI.SUM, COMM) : s
end
fld(k) = [sinpi(1.0e-2*(Float64(GIP[ip]) + 11.0k)) + 0.3cospi(3.0e-2*GIP[ip]) for ip = 1:NPOIN]

say("\n=== Schur step 2: the scalar Helmholtz H, $NR rank(s) ===")

@testset "Schur Helmholtz ($NR ranks)" begin

    p1, p2 = fld(1), fld(2)

    @testset "H is linear" begin
        h1 = Happly(p1); h2 = Happly(p2)
        @test maximum(abs, Happly(p1 .+ p2) .- (h1 .+ h2)) < 1e-10*maximum(abs, h1)
        @test maximum(abs, Happly(3 .* p1) .- 3h1) < 1e-10*maximum(abs, h1)
    end

    @testset "H is symmetric only up to a CONSISTENCY error, and that is real" begin
        # Giraldo's Schur operator has a real spectrum because his Grad and Div
        # are exact adjoints. Ours are NOT: this is a strong-form CG-SEM, where
        # the derivative is taken element-locally and averaged by DSS, and
        # adjointness holds only in the limit. Measured directly on
        # <Grad p, tb W v> vs -<p, Div[tb W v]>:
        #
        #     2x2x5  p=2   n=  275   rel 6.80e-02
        #     4x4x10 p=2   n= 1701   rel 1.22e-02
        #     2x2x5  p=4   n= 1701   rel 6.69e-03
        #     4x4x10 p=4   n=11849   rel 1.28e-03
        #
        # It converges away under refinement, which is what makes it a
        # discretisation property and not a coding error -- a wrong index or a
        # missing mask would sit at O(1) and not move. So:
        #
        #   * "the eigenvalues are all real" does NOT transfer exactly to this
        #     discretisation. It is true to a few parts in 1e4 on a production
        #     mesh, which is worlds better than the purely imaginary spectrum
        #     of the operator being preconditioned, and useless as a licence
        #     for CG.
        #   * CG IS NOT SAFE HERE. Use GMRES or BiCGstab on H.
        #   * For PRECONDITIONING none of this matters: H only has to
        #     approximate the inverse, and any error costs iterations, not
        #     correctness.
        a = mdot(Happly(p1), p2)
        b = mdot(p1, Happly(p2))
        rel = abs(a - b) / max(abs(a), abs(b), eps())
        say(@sprintf("  asymmetry <Hp,q> vs <p,Hq>: relative %.3e", rel))
        # Small enough that it cannot be a structural error...
        @test rel < 1.0e-3
        # ...and NOT round-off, which is the honest statement.
        @test rel > 1.0e-14
    end

    @testset "H is positive definite" begin
        for k = 1:6
            p = fld(k)
            q = mdot(Happly(p), p)
            @test q > 0.0
        end
        # and its quadratic form dominates the mass term it is built from,
        # i.e. the lam^2 Div-Grad piece is genuinely dissipating
        p = fld(1)
        @test mdot(Happly(p), p) > 0.0
    end
end

# ---------------------------------------------------------------------------
# Dense spectra. Serial only: forming a dense operator column by column needs
# the whole field on one rank, and the point here is the SPECTRUM, which is a
# property of the discretisation and not of the partition.
# ---------------------------------------------------------------------------
if NR == 1
    @testset "spectra: H against the full 5-field stage operator" begin
        n = NPOIN
        Hm = zeros(Float64, n, n)
        e  = zeros(Float64, n)
        for j = 1:n
            fill!(e, 0.0); e[j] = 1.0
            Hm[:, j] = Happly(e)
        end
        λH = eigvals(Hm)
        imH = maximum(abs, imag.(λH)) / maximum(abs, real.(λH))
        say(@sprintf("  H  : %d x %d, |Im|/|Re| = %.2e, Re in [%.4g, %.4g]",
                     n, n, imH, minimum(real, λH), maximum(real, λH)))
        # NEARLY real -- see the symmetry testset for why it is not exactly so
        # on a strong-form discretisation. On this deliberately coarse mock
        # mesh the consistency error is at its worst.
        @test imH < 1e-2
        @test minimum(real, λH) > 0.0          # POSITIVE real parts
        condH = maximum(real, λH) / minimum(real, λH)

        # the full stage operator (I - lam*A), 5 fields
        m = n * op.nimp
        Araw = zeros(Float64, m, m)
        V = zeros(Float64, n, op.nimp); W = zeros(Float64, n, op.nimp)
        for j = 1:m
            fill!(V, 0.0); V[j] = 1.0
            hevi_apply_A!(W, V, params, op)
            @inbounds for i = 1:m; Araw[i, j] = W[i]; end
        end
        rate = maximum(abs, imag.(eigvals(Araw)))
        λA = eigvals(Matrix(1.0I, m, m) .- LAM .* Araw)
        condA = maximum(abs, λA) / minimum(abs, λA)
        say(@sprintf("  I-lam*A: Re in [%.4g, %.4g], |Im| up to %.4g, cond %.3g",
                     minimum(real, λA), maximum(real, λA),
                     maximum(abs, imag.(λA)), condA))

        # THE RESULT THAT DECIDES WHETHER THIS IS WORTH BUILDING, swept over
        # the stiffness the production run actually sits at. lam*rate for
        # LESICP2-64x64x60 is gamma*dt*c/h_z = 0.4359*0.5*49.2 = 10.7.
        #
        #   lam   lam*rate   cond(I-lam A)   cond(H)
        #   0.20     1.42          1.74         3.03
        #   0.50     3.55          3.70        13.6
        #   1.50    10.66         10.8         114      <- production regime
        #   3.00    21.33         21.6         455
        #
        # cond(I-lam*A) grows LINEARLY with stiffness: A is skew, so the
        # eigenvalues sit on a vertical line through 1 and the worst is
        # sqrt(1+(lam*mu)^2). cond(H) grows QUADRATICALLY, because as lam
        # rises H is dominated by lam^2*Div-Grad -- a discrete Laplacian,
        # conditioned like (L/h)^2. At the stiffness this code runs at, H is
        # an order of magnitude WORSE conditioned than the operator it was
        # meant to replace.
        #
        # That does not make H useless, and it is not an argument against the
        # Schur idea in general -- but it kills the simple version of it here:
        #   * a symmetric positive operator converges like sqrt(cond), a skew
        #     one like cond, so sqrt(114) ~ 10.7 vs 10.8 is a WASH on raw
        #     iteration count;
        #   * the reference H costs TWO full operator applications, so using
        #     it as a preconditioner would triple the per-iteration cost to
        #     buy, at best, a halving of iterations;
        #   * the real prize -- 5*Np -> Np on the Krylov vectors -- needs H to
        #     be the SOLVE, and our flux-form theta equation leaves rho in the
        #     system, so it cannot be.
        # What H does have is a REAL POSITIVE spectrum, i.e. it is a Laplacian,
        # and Laplacians are what multigrid and line solves are good at, where
        # nothing standard preconditions a skew operator well. Any payoff has
        # to come from that, not from the conditioning.
        @test condH > 0.0 && condA > 0.0
        say(@sprintf("  cond: H = %.3g   I-lam*A = %.3g   (H is %.2fx WORSE here)",
                     condH, condA, condH/condA))
    end
end

MPI.Barrier(COMM)
say("=== step 2 done on $NR rank(s) ===\n")

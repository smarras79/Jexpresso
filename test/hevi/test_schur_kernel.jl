#=============================================================================
 test_schur_kernel.jl -- the bespoke scalar sweeps against the reference form
 they replace.

 THIS IS THE CHECK schur.jl PROMISED. Its header says H is built out of two
 `hevi_apply_A!` calls on purpose, because a hand-written scalar kernel "would
 duplicate the metric handling, the free-slip masking, the DSS and the mass
 division -- every one of which is already verified, and every one of which is a
 place to introduce a bug that converges to a plausible wrong answer", and that
 an optimised kernel "can be added later and checked against THIS".

 So the reference is not dead code and must not be deleted: it is the only
 independent statement of what the kernel should compute. A fast kernel with no
 reference is a kernel nobody can debug.

 WHAT WOULD MAKE THIS TEST VACUOUS, and is therefore asserted rather than
 assumed:

   * `wallx`/`wally` DEFAULT TO falses(npoin). An operator built without them
     never exercises the lateral mask, so the divergence kernel's masking branch
     would be dead and a bug in it would pass silently. This file builds real
     lateral walls and checks that they are non-empty AND that removing them
     changes the answer.
   * A gradient of a field that happens to be constant, or a divergence of a
     field that happens to be solenoidal, compares zero against zero. The probe
     fields are built from the GLOBAL node id and their outputs are checked to
     be O(1) before the agreement is believed.

     julia --project=<env> test/hevi/test_schur_kernel.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_schur_kernel.jl
=============================================================================#

using Test, MPI, LinearAlgebra, Printf
MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "schur_kernel.jl", "schur.jl")
    include(joinpath(SRC, f))
end
say(a...) = RANK == 0 && (println(a...); flush(stdout))

# p = 4, as production runs. The kernel's whole advantage is in the derivative
# chains, whose length is ngl, so measuring it at p = 2 would understate it.
const NELX, NELY, NELZ, P = 3, 3, 6, 4
params = build_mock_params(nelx=NELX, nely=NELY, nelz=NELZ, p=P,
                           Lx=3000.0, Ly=3000.0, Lz=6000.0, comm=COMM,
                           base=:stratified)
topo  = build_column_topology(params.mesh, COMM)
const N = Int(params.mesh.npoin)
const GIP = params.mesh.ip2gip

# REAL lateral walls, or the mask in the divergence kernel is never touched.
xs, ys = params.mesh.coords[1, :], params.mesh.coords[2, :]
xmin = MPI.Allreduce(minimum(xs), MPI.MIN, COMM); xmax = MPI.Allreduce(maximum(xs), MPI.MAX, COMM)
ymin = MPI.Allreduce(minimum(ys), MPI.MIN, COMM); ymax = MPI.Allreduce(maximum(ys), MPI.MAX, COMM)
tol = 1.0e-8 * max(xmax - xmin, 1.0)
wx = [abs(xs[i] - xmin) < tol || abs(xs[i] - xmax) < tol for i = 1:N]
wy = [abs(ys[i] - ymin) < tol || abs(ys[i] - ymax) < tol for i = 1:N]

op = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux=true, full=true,
                         theta_advective=true, wallx=wx, wally=wy)

const LAM = 0.26
stref  = SchurState(N, op.nimp)                                # two hevi_apply_A!
stfast = SchurState(N, op.nimp; kern = SchurKernel(params))    # bespoke sweeps

probe(s) = [sinpi(1.0e-2 * (GIP[i] + 3s)) + 0.4cospi(7.0e-3 * (GIP[i] + s)) for i = 1:N]
relerr(a, b) = (d = maximum(abs, a .- b); s = maximum(abs, b);
                MPI.Allreduce(d, MPI.MAX, COMM) / max(MPI.Allreduce(s, MPI.MAX, COMM), eps()))

say("\n=== Schur kernel vs the reference form, $NR rank(s), p=$P ===")

@testset "schur kernel ($NR ranks)" begin

    @testset "the masks are real, or the comparison is vacuous" begin
        nwx = MPI.Allreduce(count(wx), MPI.SUM, COMM)
        nwy = MPI.Allreduce(count(wy), MPI.SUM, COMM)
        nwz = MPI.Allreduce(count(op.wall), MPI.SUM, COMM)
        say(@sprintf("  masked nodes: wallx %d, wally %d, wall(z) %d", nwx, nwy, nwz))
        @test nwx > 0 && nwy > 0 && nwz > 0
        # ... and they must CHANGE the divergence, or masking is untested even
        # though the flags are set.
        opn = build_hevi_operator(params, topo, [1,2,3,4,5]; lwall_flux=true,
                                  full=true, theta_advective=true)
        vx, vy, vz = probe(1), probe(2), probe(3)
        # schur_divW! writes into `out` and returns nothing -- it always has.
        dm = zeros(N); schur_divW!(dm, vx, vy, vz, params, op,  stref.w)
        dn = zeros(N); schur_divW!(dn, vx, vy, vz, params, opn, stref.w)
        say(@sprintf("  masked vs unmasked divergence differ by %.3e", relerr(dm, dn)))
        @test relerr(dm, dn) > 1.0e-3
    end

    @testset "gradient: kernel == reference" begin
        Pp = probe(11)
        gr = (zeros(N), zeros(N), zeros(N)); gf = (zeros(N), zeros(N), zeros(N))
        schur_grad!(gr..., Pp, params, op, stref.w)
        schur_grad!(gf..., Pp, params, op, stfast.w)
        for (c, nm) in zip(1:3, ("x", "y", "z"))
            # non-trivial output, then agreement -- in that order
            @test MPI.Allreduce(maximum(abs, gr[c]), MPI.MAX, COMM) > 1.0e-8
            e = relerr(gf[c], gr[c])
            say(@sprintf("  grad %s: rel |kernel - reference| = %.3e", nm, e))
            @test e < 1.0e-12
        end
    end

    @testset "divergence: kernel == reference" begin
        vx, vy, vz = probe(4), probe(5), probe(6)
        dr = zeros(N); schur_divW!(dr, vx, vy, vz, params, op, stref.w)
        df = zeros(N); schur_divW!(df, vx, vy, vz, params, op, stfast.w)
        @test MPI.Allreduce(maximum(abs, dr), MPI.MAX, COMM) > 1.0e-8
        e = relerr(df, dr)
        say(@sprintf("  div: rel |kernel - reference| = %.3e", e))
        @test e < 1.0e-12
    end

    @testset "THE OPERATOR: H via the kernel == H via the reference" begin
        Pp = probe(21)
        hr = schur_H!(zeros(N), Pp, params, op, LAM, stref)
        hf = schur_H!(zeros(N), Pp, params, op, LAM, stfast)
        @test MPI.Allreduce(maximum(abs, hr), MPI.MAX, COMM) > 1.0e-10
        e = relerr(hf, hr)
        say(@sprintf("  H: rel |kernel - reference| = %.3e", e))
        # The two paths do the same arithmetic in a different order, so the bar
        # is round-off and not bit-identity.
        @test e < 1.0e-12
    end

    @testset "H stays linear through the kernel" begin
        # Cheap, and it catches a whole class of indexing error that a single
        # probe field can miss: a kernel that reads the wrong node can still
        # agree on one input and stop being a linear operator.
        p1, p2 = probe(31), probe(32)
        h1 = schur_H!(zeros(N), p1, params, op, LAM, stfast)
        h2 = schur_H!(zeros(N), p2, params, op, LAM, stfast)
        hs = schur_H!(zeros(N), p1 .+ p2, params, op, LAM, stfast)
        e = relerr(hs, h1 .+ h2)
        say(@sprintf("  linearity: rel |H(a+b) - H(a) - H(b)| = %.3e", e))
        @test e < 1.0e-12
    end

if NR == 1
    @testset "and it is actually faster" begin
        # MINIMUM OF SEVERAL BATCHES, not one batch's mean. On a shared box the
        # mean is dragged by whatever else is running, and the kernel is the
        # fastest thing measured here so it suffers most: a single-batch reading
        # of it came back 3x its own minimum and reported 1.97x where careful
        # measurement gives 5.78x. The minimum is the robust estimator when the
        # noise is one-sided, which scheduler noise is.
        function bench(f, nb = 7, n = 60)
            f()
            m = Inf
            for _ = 1:nb
                t = time_ns(); for _ = 1:n; f(); end
                m = min(m, (time_ns() - t) * 1e-9 / n)
            end
            return m
        end
        Pp = probe(41); HP = zeros(N); W = zeros(N, op.nimp); V = zeros(N, op.nimp)
        t_ref  = bench(() -> schur_H!(HP, Pp, params, op, LAM, stref))
        t_fast = bench(() -> schur_H!(HP, Pp, params, op, LAM, stfast))
        t_A    = bench(() -> hevi_apply_A!(W, V, params, op))
        say(@sprintf("  H reference   %.5f s   = %.2f x one hevi_apply_A!", t_ref, t_ref/t_A))
        say(@sprintf("  H kernel      %.5f s   = %.2f x one hevi_apply_A!", t_fast, t_fast/t_A))
        say(@sprintf("  speedup       %.2f x", t_ref/t_fast))
        # The reference form is ~2x one full application; the point of the
        # kernel is to get BELOW one. Asserted loosely because this is a timing
        # on a shared machine -- the number is reported for the record, and only
        # a gross failure fails the test.
        @test t_fast < t_ref
    end
end
end

MPI.Barrier(COMM)
say("=== schur kernel checks done on $NR rank(s) ===\n")

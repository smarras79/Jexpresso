#=============================================================================
 tools/poisson3d_benchmark/mpi/verify_mpi.jl — correctness checks of the MPI
 solvers (run on any number of ranks):

   mpiexec -n 8 julia --project=tools/poisson3d_benchmark/mpi \
       tools/poisson3d_benchmark/mpi/verify_mpi.jl [--reference ref/results.csv]

 1. the three MPI solvers (MUMPS, BoomerAMG-CG, Jacobi-CG) give the same
    solution of the same discrete system: errors agree to 1e-10;
 2. h-convergence at N = 3, 4³ -> 16³ elements (two doublings; single
    doublings are still pre-asymptotic: 4.5, then 3.5): rate N+1 = 4, within 0.3;
 3. with --reference: the errors equal those of the serial, Jexpresso-based
    benchmark (rows of ../bench3d.jl, e.g. --solvers sem) at the same
    (ne, N), to 1e-10: same discretisation, same right-hand side, same
    exact solution, same pinning.
 Exits with status 1 if a check fails.
=============================================================================#
using MPI
MPI.Init()
using HYPRE
HYPRE.Init()
using Printf
include(joinpath(@__DIR__, "poisson3d_mpi.jl"))
using .P3DMPI

comm = MPI.COMM_WORLD; me = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
say(a...) = me == 0 && (println(a...); flush(stdout))
ok = true
ref = nothing
i = findfirst(==("--reference"), ARGS)
i === nothing || (ref = ARGS[i+1])
pmax = maximum(MPI.Dims_create(np, [0, 0, 0]))
say("verify_mpi: $np ranks, process grid ", join(MPI.Dims_create(np, [0, 0, 0]), "x"))

# 1. solver agreement
say("1. agreement of the MPI solvers")
for (ne, N) in ((4, 3), (6, 4))
    ne * N >= pmax || continue
    errs = [run_config(s, ne, N).linf for s in MPI_SOLVERS]
    spread = maximum(errs) - minimum(errs)
    say(@sprintf("   ne=%d N=%d  L∞: mumps %.12e  boomeramg %.12e  jacobi %.12e  spread %.1e", ne, N, errs..., spread))
    global ok &= spread < 1e-10
end

# 2. h-convergence
say("2. h-convergence, N = 3")
if 4 * 3 >= pmax
    e1 = run_config(:boomeramg, 4, 3).linf; e2 = run_config(:boomeramg, 16, 3).linf
    rate = log2(e1 / e2) / 2
    say(@sprintf("   L∞ %.3e (4³) -> %.3e (16³): rate %.2f (expected 4)", e1, e2, rate))
    global ok &= abs(rate - 4) < 0.3
end

# 3. against the serial benchmark
if ref !== nothing
    say("3. against the serial benchmark ($ref)")
    lines = readlines(ref); hdr = split(lines[1], ",")
    col(k) = findfirst(==(k), hdr)
    for l in lines[2:end]
        f = split(l, ",")
        f[col("status")] == "ok" || continue
        ne = parse(Int, f[col("ne")]); N = parse(Int, f[col("nop")]); es = parse(Float64, f[col("linf")])
        ne * N >= pmax || continue
        for s in MPI_SOLVERS
            em = run_config(s, ne, N).linf
            say(@sprintf("   ne=%d N=%d  serial %-6s %.12e   MPI %-9s %.12e   |Δ| = %.1e",
                         ne, N, f[col("solver")], es, s, em, abs(em - es)))
            global ok &= abs(em - es) < 1e-10
        end
    end
end

say(ok ? "PASS" : "FAIL")
MPI.Barrier(comm)
exit(ok ? 0 : 1)

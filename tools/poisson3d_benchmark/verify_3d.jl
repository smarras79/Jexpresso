#=============================================================================
 verify_3d.jl — correctness checks of the 3D solvers (small sizes, minutes)

   1. p-convergence: SEM (direct) on 4³ elements, N = 2..8: the error falls
      exponentially with N.
   2. h-convergence: SEM (direct) at N = 3, 4³ … 16³ elements: the L2 error
      falls like h^(N+1) (observed rates printed; the asymptotic rate is
      reached once the peak of u is resolved).
   3. agreement: the five SEM solvers give the same solution, and the
      pseudo-spectral and FFT solvers give the same solution.

     julia --project=. -t 4 tools/poisson3d_benchmark/verify_3d.jl
=============================================================================#
using Jexpresso, Printf
include(joinpath(@__DIR__, "poisson3d.jl"))
using .P3D

println("1. p-convergence, SEM direct, 4^3 elements")
for N in 2:8
    r = run_config(:sem, 3, 4, N)
    @printf("   N=%d  n=%8d  L∞=%.3e  L2rel=%.3e\n", N, r.n, r.linf, r.l2rel)
end

println("2. h-convergence, SEM direct, N = 3 (expected L2 rate N+1 = 4)")
prev = nothing; first6 = nothing
for ne in (4, 6, 8, 12, 16)
    r = run_config(:sem, 3, ne, 3)
    rate = prev === nothing ? "" : @sprintf("  rate %.2f", log(prev[2] / r.l2rel) / log(ne / prev[1]))
    @printf("   ne=%2d  n=%8d  L2rel=%.3e%s\n", ne, r.n, r.l2rel, rate)
    global prev = (ne, r.l2rel)
    ne == 6 && (global first6 = prev)
end
@printf("   overall rate 6^3 -> 16^3: %.2f (local rates oscillate while the peak of u is being resolved)\n",
        log(first6[2] / prev[2]) / log(prev[1] / first6[1]))

println("3. agreement at 6^3 elements, N = 4 (Fourier grid 24^3)")
rows = [run_config(s, 3, 6, 4) for s in P3D.SOLVERS]
for r in rows
    @printf("   %-11s L∞=%.15e  its=%d\n", r.solver, r.linf, r.iters)
end
sem = [r.linf for r in rows if r.solver in P3D.SEM_SOLVERS]
four = [r.linf for r in rows if !(r.solver in P3D.SEM_SOLVERS)]
dsem = maximum(sem) - minimum(sem); dfour = maximum(four) - minimum(four)
@printf("   spread of the SEM errors: %.1e   of the Fourier errors: %.1e\n", dsem, dfour)
println(dsem < 1e-9 && dfour < 1e-9 ? "PASS" : "FAIL", ": solvers of the same discretisation agree")

#=============================================================================
 verify_3d_jexpresso.jl — Jexpresso's own 3D solves (run_case on
 problems/Elliptic/poisson_periodic_sem_3d, jexpresso3d.jl) against the
 independent Kronecker assembly of poisson3d.jl, solver by solver:

   SEM direct, SEM AMG, SC direct, SC AMG   same SEM system: the L∞ errors
                                            must agree to round-off (1e-10)
   pseudo-spectral, FFT                     same Fourier collocation system:
                                            the same, to 1e-10

 on 4³ elements at N = 2, 3 and 6³ elements at N = 4 (the Fourier grids
 (ne·N)³). Exits with status 1 if any pair differs by more.

     julia --project=. -t 4 tools/poisson3d_benchmark/verify_3d_jexpresso.jl
=============================================================================#
using Jexpresso, Printf
include(joinpath(@__DIR__, "poisson3d.jl"))
include(joinpath(@__DIR__, "jexpresso3d.jl"))
using .P3D, .JX3D

worst = 0.0
println("  solver      ne  N   Jexpresso L∞             Kronecker L∞             |difference|")
for (ne, N) in ((4, 2), (4, 3), (6, 4)), s in JX3D.SOLVERS
    ej = JX3D.run_config(s, ne, N).linf
    ek = P3D.run_config(s, 3, ne, N; r = JX3D.R).linf
    dif = abs(ej - ek)
    global worst = max(worst, dif)
    @printf("  %-10s  %2d  %d   %.15e   %.15e   %.1e\n", s, ne, N, ej, ek, dif)
end
println(worst < 1e-10 ? "PASS" : "FAIL", ": largest absolute difference $worst (max |u| = 1)")
exit(worst < 1e-10 ? 0 : 1)

#=============================================================================
 verify_2d.jl — the Kronecker-assembled SEM system of poisson3d.jl, in 2D,
 against Jexpresso's own periodic SEM solve of the same problem, run live:
 Jexpresso.run_case("Elliptic", "poisson_periodic_sem") on 16×16 elements,
 N = 2..8, SEM direct (sem_setup → DSS_laplace_sparse → periodic_sem_system
 → periodic_sem_factorize / periodic_sem_direct_solve), its L∞ error read
 from Jexpresso.JX_LAST_SOLVE_ERR.

 Both solve the same discrete problem, so the errors must agree to
 round-off: an absolute difference of order 1e-12 (u peaks at 1), amplified
 by the condition number of K as N grows, not a relative one (at N = 8 the
 error itself is only 5e-5).

     julia --project=. tools/poisson3d_benchmark/verify_2d.jl [Nmax]
=============================================================================#
using Jexpresso, Printf
include(joinpath(@__DIR__, "poisson3d.jl"))
using .P3D

const NEL  = 16                       # the deck's Cartesian grid: 16×16 elements
const NMAX = isempty(ARGS) ? 8 : parse(Int, ARGS[1])

# Jexpresso's own periodic SEM direct solve of the deck, at order N
function jexpresso_linf(N)
    ov = Dict{Symbol, Any}(:nop => N, :linitial_refine => false,
                           :lfft => false, :lpseudospectral => false,
                           :linsolve_amg => false, :lstatic_condensation => false,
                           :luse_mesh_cache => false, :lbenchmark_solve => false,
                           :outformat => "none")
    Jexpresso.run_case("Elliptic", "poisson_periodic_sem"; inputs = ov)
    return Jexpresso.JX_LAST_SOLVE_ERR[].linf
end

worst = 0.0
println("  N    Jexpresso run_case L∞   Kronecker L∞            |difference|")
for N in 2:NMAX
    ej = jexpresso_linf(N)
    ek = run_config(:sem, 2, NEL, N).linf
    dif = abs(ek - ej)
    global worst = max(worst, dif)
    @printf("  %d    %.15e   %.15e   %.1e\n", N, ej, ek, dif)
end
println(worst < 1e-10 ? "PASS" : "FAIL", ": largest absolute difference $worst (max |u| = 1)")
exit(worst < 1e-10 ? 0 : 1)

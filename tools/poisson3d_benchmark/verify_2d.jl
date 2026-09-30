#=============================================================================
 verify_2d.jl — the Kronecker-assembled SEM system of poisson3d.jl, in 2D,
 against Jexpresso's own periodic SEM solve of the same problem
 (problems/Elliptic/poisson_periodic_sem, 16×16 elements, N = 2..8).

 Both solve the same discrete problem, so the errors must agree to
 round-off: an absolute difference of order 1e-12 (u peaks at 1), amplified
 by the condition number of K as N grows, not a relative one (at N = 8 the
 error itself is only 5e-5). The Jexpresso numbers are read from
 tools/periodic_poisson_benchmark/results.csv (the SEM direct rows, level 0).

     julia --project=. tools/poisson3d_benchmark/verify_2d.jl
=============================================================================#
using Jexpresso, Printf
include(joinpath(@__DIR__, "poisson3d.jl"))
using .P3D

csv = joinpath(@__DIR__, "..", "periodic_poisson_benchmark", "results.csv")
lines = readlines(csv); hdr = split(lines[1], ",")
col(r, k) = r[findfirst(==(k), hdr)]
ref = Dict{Int, Float64}()
for l in lines[2:end]
    r = split(l, ",")
    col(r, "solver") == "sem" && get(Dict(zip(hdr, r)), "level", "0") == "0" &&
        (ref[parse(Int, col(r, "nop"))] = parse(Float64, col(r, "linf")))
end

worst = 0.0
println("  N    Jexpresso L∞            Kronecker L∞            |difference|")
for N in sort(collect(keys(ref)))
    r = run_config(:sem, 2, 16, N)
    dif = abs(r.linf - ref[N])
    global worst = max(worst, dif)
    @printf("  %d    %.15e   %.15e   %.1e\n", N, ref[N], r.linf, dif)
end
println(worst < 1e-10 ? "PASS" : "FAIL", ": largest absolute difference $worst (max |u| = 1)")

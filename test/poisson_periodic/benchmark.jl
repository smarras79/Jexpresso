#=============================================================================
 test/poisson_periodic/benchmark.jl — solver-agnostic benchmark for the
 periodic Poisson / Laplace problem

     -∇²u = f   on a periodic box,   answer fixed to zero mean.

 The manufactured problems (PROBLEMS) are closed-form (u, f) pairs, so every
 solver — whatever grid or discretisation it uses — is scored against the SAME
 exact answers. Two scoring paths:

 1. UNIFORM-GRID solvers (the FFT) are called through the "solve function"

        solve(f::Array{Float64,ND}, Ls::NTuple{ND,Float64}) -> u

      f   right-hand side sampled on the uniform periodic grid
          x_d,i = x0_d + (i-1) L_d/n_d  (right-hand periodic image not stored)
      Ls  the periods
      u   any solution; the harness removes its mean before comparing, so a
          solver that pins a node instead of fixing the mean is not penalised.

    verify_periodic_poisson_solver(name, solve)
        correctness against the exact answers: Laplace equation, null space,
        spectral exactness and exponential convergence.
    compare_periodic_poisson_solvers(nameA, solveA, nameB, solveB)
        two uniform-grid solvers must return the same solution; timing table.

 2. NODAL solvers on Jexpresso's periodic SEM mesh (LGL nodes):
      • the standard direct SEM solve,
      • an AlgebraicMultigrid.jl solve of the SAME SEM system,
      • the element-learning static-condensation solve.
    They live on LGL nodes, not on the uniform grid, so they are scored with
    nodal_error_norms(p, X, u; w): node coordinates X (nsd × npoin, the layout
    of mesh.coords), the nodal solution u and the quadrature weights w (the
    lumped mass matrix), which fix the free constant by the w-weighted mean and
    weight the relative L² norm. The FFT result is scored the same way (on its
    own grid, uniform weights), which puts all four solvers in one table.
    The direct and AMG solves of the SEM system must additionally agree with
    each other to the AMG tolerance.

 Depends on Test, LinearAlgebra and Printf only.
=============================================================================#
module PeriodicPoissonBenchmark

using Test, LinearAlgebra, Printf

export PoissonProblem, PROBLEMS, problem, sample_problem, gauge, error_norms,
       grid_points, nodal_error_norms, reference_spectral_solve,
       verify_periodic_poisson_solver, compare_periodic_poisson_solvers

#-----------------------------------------------------------------------------
# Manufactured problems. `u` is the exact solution of -∇²u = f; `f` is its
# negative Laplacian, written out in closed form. `maxmode[d]` is the highest
# integer wavenumber (in units of 2π/L_d) present in u along axis d, or
# `nothing` when u is not band-limited.
#-----------------------------------------------------------------------------
struct PoissonProblem{ND}
    name    :: String
    Ls      :: NTuple{ND, Float64}
    x0s     :: NTuple{ND, Float64}
    u       :: Function
    f       :: Function
    maxmode :: Union{Nothing, NTuple{ND, Int}}
end

const PROBLEMS = Dict{String, PoissonProblem}()
_add!(p) = (PROBLEMS[p.name] = p)

# 2D, [0,2π]², a few low Fourier modes — the case deck's problem.
_add!(PoissonProblem("modes2d", (2π, 2π), (0.0, 0.0),
    (x, y) -> sin(2x)*cos(3y) + sin(x)*cos(y),
    (x, y) -> 13sin(2x)*cos(3y) + 2sin(x)*cos(y),
    (2, 3)))

# 2D, non-square, shifted box [-1,1) × [0.5,1.5): Lx = 2, Ly = 1. Anisotropic
# wavenumbers kx = π·m, ky = 2π·m; catches any Lx/Ly or x0 mix-up.
_add!(PoissonProblem("aniso2d", (2.0, 1.0), (-1.0, 0.5),
    (x, y) -> sin(π*x)*cos(4π*y) + 0.5cos(3π*x),
    (x, y) -> (π^2 + 16π^2)*sin(π*x)*cos(4π*y) + 0.5*9π^2*cos(3π*x),
    (1, 2)))

# 2D, smooth but NOT band-limited: u = exp(sin x + cos 2y). With
# g = sin x + cos 2y,  ∇²u = e^g (|∇g|² + ∇²g). Its mean is not zero (the
# harness gauges both sides), which also checks that the null space is handled.
_add!(PoissonProblem("smooth2d", (2π, 2π), (0.0, 0.0),
    (x, y) -> exp(sin(x) + cos(2y)),
    (x, y) -> -exp(sin(x) + cos(2y)) *
              (cos(x)^2 + 4sin(2y)^2 - sin(x) - 4cos(2y)),
    nothing))

# 3D, [0,2π]³.
_add!(PoissonProblem("modes3d", (2π, 2π, 2π), (0.0, 0.0, 0.0),
    (x, y, z) -> sin(x)*cos(2y)*sin(z) + cos(3z),
    (x, y, z) -> 6sin(x)*cos(2y)*sin(z) + 9cos(3z),
    (1, 2, 3)))

problem(name) = PROBLEMS[name]

grid_lines(p::PoissonProblem{ND}, dims::NTuple{ND, Int}) where {ND} =
    ntuple(d -> [p.x0s[d] + (i - 1) * p.Ls[d] / dims[d] for i in 1:dims[d]], ND)

function _sample(fun, lines::NTuple{ND}) where {ND}
    A = Array{Float64, ND}(undef, map(length, lines))
    @inbounds for I in CartesianIndices(A)
        A[I] = fun(ntuple(d -> lines[d][I[d]], ND)...)
    end
    return A
end

"(f, u_exact) sampled on the `dims` grid of problem `p`."
sample_problem(p::PoissonProblem{ND}, dims::NTuple{ND, Int}) where {ND} =
    (_sample(p.f, grid_lines(p, dims)), _sample(p.u, grid_lines(p, dims)))

"`u` with its discrete mean removed (the periodic solution's free constant)."
gauge(u) = u .- sum(u) / length(u)

"(L∞, relative discrete L²) error of `u` against `uex`, both gauged to zero mean."
function error_norms(u, uex)
    e, r = gauge(u) .- gauge(uex), gauge(uex)
    return (linf = maximum(abs, e), l2rel = norm(e) / max(norm(r), eps()))
end

"All points of the `dims` grid of `p` as an nsd × npoin matrix (mesh.coords layout)."
function grid_points(p::PoissonProblem{ND}, dims::NTuple{ND, Int}) where {ND}
    lines = grid_lines(p, dims)
    X = Matrix{Float64}(undef, ND, prod(dims))
    for (k, I) in enumerate(CartesianIndices(dims)), d = 1:ND
        X[d, k] = lines[d][I[d]]
    end
    return X
end

"""
    nodal_error_norms(p, X, u; w = nothing) -> (linf, l2rel)

Error of the nodal solution `u[ip]` at the nodes `X[:, ip]` (nsd × npoin)
against the exact solution of problem `p`, for solvers that do not work on the
uniform grid (SEM on LGL nodes). `w` are the quadrature weights of the nodes —
the diagonal of the lumped mass matrix — and default to uniform. The free
constant is fixed on both sides by the w-weighted mean; the relative L² norm is
w-weighted as well, so it approximates the continuous one.
"""
function nodal_error_norms(p::PoissonProblem{ND}, X::AbstractMatrix, u::AbstractVector;
                           w = nothing) where {ND}
    size(X, 1) == ND || throw(DimensionMismatch("X has $(size(X,1)) rows, problem $(p.name) is $(ND)D"))
    size(X, 2) == length(u) || throw(DimensionMismatch("X has $(size(X,2)) nodes, u has $(length(u))"))
    wt  = w === nothing ? ones(length(u)) : collect(Float64, w)
    uex = [p.u(ntuple(d -> X[d, ip], ND)...) for ip in axes(X, 2)]
    wm(v) = sum(wt .* v) / sum(wt)
    e   = (u .- wm(u)) .- (uex .- wm(uex))
    r   = uex .- wm(uex)
    return (linf = maximum(abs, e), l2rel = sqrt(sum(wt .* e .^ 2) / max(sum(wt .* r .^ 2), eps())))
end

#-----------------------------------------------------------------------------
# Independent reference: dense-DFT Fourier solve.
#-----------------------------------------------------------------------------
"""
    reference_spectral_solve(f, Ls) -> u

The Fourier spectral solution of -∇²u = f - mean(f), zero mean, computed with
EXPLICIT DFT matrices applied along each axis (O(n) work per point per axis),
no FFT library involved. Slow, but an independent implementation of the same
mathematics — the trusted baseline the FFT is compared against. Keep grids
modest (≲ 128 per axis).
"""
function reference_spectral_solve(f::Array{Float64, ND}, Ls::NTuple{ND, Float64}) where {ND}
    dims = size(f)
    Fwd  = [ComplexF64[cis(-2π * j * m / n) for m in 0:n-1, j in 0:n-1] for n in dims]
    Bwd  = [conj.(F) ./ size(F, 1) for F in Fwd]
    along(A, M, d) = mapslices(v -> M * v, A; dims = d)
    fh = ComplexF64.(f)
    for d = 1:ND
        fh = along(fh, Fwd[d], d)
    end
    for I in CartesianIndices(fh)
        k2 = 0.0
        for d = 1:ND
            m   = I[d] - 1
            m̃   = 2m <= dims[d] ? m : m - dims[d]
            k2 += (2π * m̃ / Ls[d])^2
        end
        fh[I] = k2 == 0 ? 0.0 : fh[I] / k2
    end
    for d = 1:ND
        fh = along(fh, Bwd[d], d)
    end
    return real.(fh)
end

#-----------------------------------------------------------------------------
# Correctness of ONE uniform-grid solver.
#-----------------------------------------------------------------------------
"Local algebraic order between successive (n, error) pairs."
_rates(ns, errs) = [log(errs[i-1] / errs[i]) / log(ns[i] / ns[i-1]) for i in 2:length(errs)]

"""
    verify_periodic_poisson_solver(name, solve; rtol_solver = 1e-10, verbose = true)

Run the benchmark's correctness checks on a uniform-grid `solve(f, Ls) -> u`
(see the file header). `rtol_solver` is the relative accuracy the solver
promises on its own discrete problem.
"""
function verify_periodic_poisson_solver(name::AbstractString, solve;
                                        rtol_solver::Real = 1e-10,
                                        verbose::Bool = true)
    @testset verbose = verbose "$name" begin

        @testset "Laplace equation (f ≡ 0) ⇒ u ≡ 0" begin
            for dims in ((16, 16), (15, 22), (8, 6, 10))
                Ls = ntuple(d -> 1.0 + d, length(dims))
                u  = solve(zeros(dims), Ls)
                @test size(u) == dims
                @test maximum(abs, u) < 1e-12
            end
        end

        @testset "linearity" begin
            p = problem("modes2d")
            f, _ = sample_problem(p, (32, 32))
            u  = solve(f, p.Ls)
            u2 = solve(2 .* f, p.Ls)
            @test maximum(abs, gauge(u2) .- 2 .* gauge(u)) < 10rtol_solver * maximum(abs, u)
        end

        @testset "incompatible RHS ⇒ solves the projected problem" begin
            # f + c has no periodic solution for c ≠ 0; the answer must be
            # that of f (the mean is projected out), not garbage.
            p = problem("aniso2d")
            f, _ = sample_problem(p, (40, 24))
            u  = solve(f,         p.Ls)
            uc = solve(f .+ 3.25, p.Ls)
            @test maximum(abs, gauge(uc) .- gauge(u)) < 10rtol_solver * maximum(abs, u)
        end

        @testset "band-limited u reproduced to round-off" begin
            for (pname, dimsets) in (("modes2d",  ((8, 8), (16, 16), (9, 11), (64, 32))),
                                     ("aniso2d",  ((8, 8), (7, 10), (24, 48))),
                                     ("modes3d",  ((8, 8, 8), (5, 7, 9), (16, 12, 10))))
                p = problem(pname)
                for dims in dimsets
                    # resolvable: every mode strictly below the Nyquist mode
                    all(d -> dims[d] > 2p.maxmode[d], eachindex(dims)) || continue
                    f, uex = sample_problem(p, dims)
                    e = error_norms(solve(f, p.Ls), uex)
                    @test e.linf < 1e-12 * max(1.0, maximum(abs, uex))
                end
            end
        end

        @testset "exponential convergence (non-band-limited u)" begin
            p  = problem("smooth2d")
            ns = [8, 12, 16, 24, 32]
            errs = map(ns) do n
                f, uex = sample_problem(p, (n, n))
                error_norms(solve(f, p.Ls), uex).linf
            end
            r = _rates(ns, errs)
            verbose && @info "  $name smooth2d: L∞ errors $(round.(errs; sigdigits=3)) at n = $ns, local orders $(round.(r; digits=1))"
            # Exponential, not algebraic: the local algebraic order
            # log(e_i-1/e_i)/log(n_i/n_i-1) of a fixed-order method tends to
            # a constant; for a spectral method it keeps GROWING with n.
            @test issorted(r)
            @test r[end] > 15
            f, uex = sample_problem(p, (48, 48))
            @test error_norms(solve(f, p.Ls), uex).linf < 1e-12
        end
    end
end

#-----------------------------------------------------------------------------
# Agreement of TWO uniform-grid solvers.
#-----------------------------------------------------------------------------
const DEFAULT_COMPARISON_CASES = (("modes2d", (64, 64)), ("aniso2d", (96, 48)),
                                  ("smooth2d", (64, 64)), ("smooth2d", (50, 45)),
                                  ("modes3d", (16, 16, 16)))

function _best_time(solve, f, Ls, reps)
    u = solve(f, Ls)                         # also warms up / compiles
    t = Inf
    for _ in 1:reps
        t = min(t, @elapsed solve(f, Ls))
    end
    return u, t
end

"""
    compare_periodic_poisson_solvers(nameA, solveA, nameB, solveB;
        cases = DEFAULT_COMPARISON_CASES, rtol = 1e-10, reps = 3, verbose = true) -> rows

For every (problem, grid) in `cases`, solve with both uniform-grid solvers,
check that the gauged solutions agree to `rtol` (relative to max|u|), and
report each solver's error against the exact solution and its best-of-`reps`
wall time. Returns the table rows (NamedTuples).
"""
function compare_periodic_poisson_solvers(nameA::AbstractString, solveA,
                                          nameB::AbstractString, solveB;
                                          cases = DEFAULT_COMPARISON_CASES,
                                          rtol::Real = 1e-10, reps::Int = 3,
                                          verbose::Bool = true)
    rows = NamedTuple[]
    @testset verbose = verbose "$nameA vs $nameB" begin
        for (pname, dims) in cases
            p = problem(pname)
            f, uex = sample_problem(p, dims)
            uA, tA = _best_time(solveA, f, p.Ls, reps)
            uB, tB = _best_time(solveB, f, p.Ls, reps)
            diff = maximum(abs, gauge(uA) .- gauge(uB)) / max(maximum(abs, gauge(uA)), eps())
            push!(rows, (problem = pname, dims = dims, npts = prod(dims), diff = diff,
                         errA = error_norms(uA, uex).linf, errB = error_norms(uB, uex).linf,
                         tA = tA, tB = tB))
            @test diff < rtol
        end
    end
    verbose && _print_comparison(nameA, nameB, rows)
    return rows
end

function _print_comparison(nameA, nameB, rows)
    println()
    println("  periodic Poisson, -∇²u = f")
    @printf("  %-9s %-13s %9s  %-10s %-10s %-10s  %-10s %-10s\n",
            "problem", "grid", "points", "‖uA-uB‖∞", "err A", "err B", "time A[s]", "time B[s]")
    println("  ", repeat("─", 94))
    for r in rows
        @printf("  %-9s %-13s %9d  %-10.2e %-10.2e %-10.2e  %-10.2e %-10.2e\n",
                r.problem, join(r.dims, "×"), r.npts, r.diff, r.errA, r.errB, r.tA, r.tB)
    end
    println("  A = $nameA,  B = $nameB;  errors are L∞ vs the exact solution, both gauged to zero mean")
    println()
end

end # module PeriodicPoissonBenchmark

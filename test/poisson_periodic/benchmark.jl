#=============================================================================
 test/poisson_periodic/benchmark.jl — solver-agnostic benchmark for the
 periodic Poisson / Laplace problem

     -∇²u = f   on a periodic box,   answer fixed to zero mean.

 Any solver is tested through ONE calling convention, the "solve function":

     solve(f::Array{Float64,ND}, Ls::NTuple{ND,Float64}, disc::Symbol) -> u

   f     right-hand side sampled on the uniform periodic grid
         x_d,i = x0_d + (i-1) L_d/n_d  (right-hand periodic image not stored)
   Ls    the periods
   disc  the discretisation of -∇² the solver is asked to invert:
           :spectral   Fourier (only a spectral solver can offer this)
           :fd2        2nd-order centred finite differences, i.e. the sparse
                       matrix returned by fd2_periodic_laplacian
   u     any solution; the harness removes its mean before comparing, so a
         solver that pins a node instead of fixing the mean is not penalised.

 A solver declares which discretisations it supports and is then run through
 the SAME checks as every other solver:

   verify_periodic_poisson_solver(name, solve; discretizations = (...))
       correctness of one solver against closed-form answers.

   compare_periodic_poisson_solvers(nameA, solveA, nameB, solveB; ...)
       agreement of two solvers of the SAME discrete system, plus a timing
       table. With :fd2 both must return the same discrete solution, so they
       have to agree to solver tolerance — far tighter than the truncation
       error they share with the exact solution.

 Step 1 (now): the FFTW solver, :spectral and :fd2, checked against the exact
 answers and against a sparse direct solve of the fd2 system.
 Step 2: an AlgebraicMultigrid.jl solve of fd2_periodic_laplacian plugs into
 verify_periodic_poisson_solver(...; discretizations = (:fd2,)) and into
 compare_periodic_poisson_solvers against the FFT with no change to this file.
 Step 3: Jexpresso's SEM solves — the standard direct SEM solve and the
 element-learning static-condensation solve — on a periodic mesh. They live on
 LGL nodes, not on this uniform grid, and discretise a different system, so
 they cannot be compared node-by-node with the FFT. They are compared through
 the SAME manufactured problems instead, each solver against the exact answer
 on its own nodes: nodal_error_norms(p, X, u; w) takes the node coordinates X
 (nsd × npoin, the layout of mesh.coords) and the quadrature weights w (the
 lumped mass matrix), fixes the free constant with the w-weighted mean and
 returns the same L∞ / relative-L² pair as error_norms. The problem data are
 closed-form functions, so a case deck evaluates them directly.

 Depends on Test, LinearAlgebra, SparseArrays and Printf only.
=============================================================================#
module PeriodicPoissonBenchmark

using Test, LinearAlgebra, SparseArrays, Printf

export PoissonProblem, PROBLEMS, problem, sample_problem, gauge, error_norms,
       grid_points, nodal_error_norms,
       fd2_periodic_laplacian, direct_fd2_solve,
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
# The fd2 discrete system.
#-----------------------------------------------------------------------------
"""
    fd2_periodic_laplacian(dims, Ls) -> SparseMatrixCSC

The periodic second-order centred finite-difference approximation of -∇² on
the uniform `dims` grid of periods `Ls`, acting on `vec(u)` (column-major).
Symmetric positive SEMI-definite; its null space is the constant vector.
Every axis needs at least 3 points (with 2, the ± neighbours coincide).
"""
function fd2_periodic_laplacian(dims::NTuple{ND, Int}, Ls::NTuple{ND, <:Real}) where {ND}
    all(>=(3), dims) || throw(ArgumentError("fd2 stencil needs ≥ 3 points per axis, got $dims"))
    N  = prod(dims)
    LI = LinearIndices(dims)
    I, J, V = Int[], Int[], Float64[]
    sizehint!(I, (2ND + 1) * N); sizehint!(J, (2ND + 1) * N); sizehint!(V, (2ND + 1) * N)
    for C in CartesianIndices(dims)
        row  = LI[C]
        diag = 0.0
        for d = 1:ND
            ih2 = (dims[d] / Ls[d])^2
            diag += 2ih2
            for s in (-1, 1)
                nb = ntuple(e -> e == d ? mod1(C[e] + s, dims[e]) : C[e], ND)
                push!(I, row); push!(J, LI[nb...]); push!(V, -ih2)
            end
        end
        push!(I, row); push!(J, row); push!(V, diag)
    end
    return sparse(I, J, V, N, N)
end

"""
    direct_fd2_solve(f, Ls) -> u

Reference solution of the fd2 system by sparse LU: the mean of f is removed
(compatibility), the first unknown is pinned to 0 to remove the null space,
and the result is gauged to zero mean. A trusted baseline for the comparison
harness, and a template for the AMG solve function.
"""
function direct_fd2_solve(f::Array{Float64, ND}, Ls::NTuple{ND, Float64}) where {ND}
    A  = fd2_periodic_laplacian(size(f), Ls)
    b  = vec(f) .- sum(f) / length(f)
    u  = zeros(length(b))
    u[2:end] = lu(A[2:end, 2:end]) \ b[2:end]
    return reshape(gauge(u), size(f))
end

#-----------------------------------------------------------------------------
# Correctness of ONE solver.
#-----------------------------------------------------------------------------
_fd2_residual(u, f, Ls) =
    maximum(abs, fd2_periodic_laplacian(size(f), Ls) * vec(u) .- (vec(f) .- sum(f) / length(f))) /
    max(maximum(abs, f), eps())

"Observed order of convergence between successive (h, error) pairs."
_rates(ns, errs) = [log(errs[i-1] / errs[i]) / log(ns[i] / ns[i-1]) for i in 2:length(errs)]

"""
    verify_periodic_poisson_solver(name, solve; discretizations = (:spectral, :fd2),
                                   rtol_solver = 1e-10, verbose = true)

Run the benchmark's correctness checks on `solve` (see the file header for
its calling convention) for every discretisation it supports. `rtol_solver`
is the relative accuracy the solver promises on its own discrete system — a
direct/FFT solver meets 1e-10 easily; an iterative one is run with a matching
stopping tolerance.
"""
function verify_periodic_poisson_solver(name::AbstractString, solve;
                                        discretizations = (:spectral, :fd2),
                                        rtol_solver::Real = 1e-10,
                                        verbose::Bool = true)
    @testset verbose = verbose "$name" begin

        for disc in discretizations
            @testset "$disc: Laplace equation (f ≡ 0) ⇒ u ≡ 0" begin
                for dims in ((16, 16), (15, 22), (8, 6, 10))
                    Ls = ntuple(d -> 1.0 + d, length(dims))
                    u  = solve(zeros(dims), Ls, disc)
                    @test size(u) == dims
                    @test maximum(abs, u) < 1e-12
                end
            end

            @testset "$disc: returned solution is gauge-consistent" begin
                p = problem("modes2d")
                f, _ = sample_problem(p, (32, 32))
                u  = solve(f, p.Ls, disc)
                u2 = solve(2 .* f, p.Ls, disc)        # linearity
                @test maximum(abs, gauge(u2) .- 2 .* gauge(u)) < 10rtol_solver * maximum(abs, u)
            end

            @testset "$disc: incompatible RHS ⇒ solves the projected problem" begin
                # f + c has no periodic solution for c ≠ 0; the answer must be
                # that of f (the mean is projected out), not garbage.
                p = problem("aniso2d")
                f, _ = sample_problem(p, (40, 24))
                u  = solve(f,         p.Ls, disc)
                uc = solve(f .+ 3.25, p.Ls, disc)
                @test maximum(abs, gauge(uc) .- gauge(u)) < 10rtol_solver * maximum(abs, u)
            end
        end

        if :fd2 in discretizations
            @testset "fd2: solves the sparse fd2 system (residual)" begin
                # Includes odd, non-power-of-2 and anisotropic grids.
                for (pname, dims) in (("modes2d", (32, 32)), ("aniso2d", (48, 30)),
                                      ("smooth2d", (27, 25)), ("modes3d", (12, 10, 9)))
                    p = problem(pname)
                    f, _ = sample_problem(p, dims)
                    u = solve(f, p.Ls, :fd2)
                    r = _fd2_residual(u, f, p.Ls)
                    verbose && @info @sprintf("  %-9s fd2 %-12s  relative residual %.2e", name, string(dims), r)
                    @test r < rtol_solver
                end
            end

            @testset "fd2: second-order convergence to the exact solution" begin
                for pname in ("smooth2d", "aniso2d")
                    p  = problem(pname)
                    ns = [32, 64, 128, 256]     # asymptotic range for both problems
                    errs = map(ns) do n
                        f, uex = sample_problem(p, (n, n))
                        error_norms(solve(f, p.Ls, :fd2), uex).linf
                    end
                    r = _rates(ns, errs)
                    verbose && @info "  $name fd2 $pname: L∞ errors $(round.(errs; sigdigits=3)), rates $(round.(r; digits=3))"
                    @test all(x -> 1.9 < x < 2.1, r)
                end
            end
        end

        if :spectral in discretizations
            @testset "spectral: band-limited u reproduced to round-off" begin
                for (pname, dimsets) in (("modes2d",  ((8, 8), (16, 16), (9, 11), (64, 32))),
                                         ("aniso2d",  ((8, 8), (7, 10), (24, 48))),
                                         ("modes3d",  ((8, 8, 8), (5, 7, 9), (16, 12, 10))))
                    p = problem(pname)
                    for dims in dimsets
                        # resolvable: every mode strictly below the Nyquist mode
                        all(d -> dims[d] > 2p.maxmode[d], eachindex(dims)) || continue
                        f, uex = sample_problem(p, dims)
                        e = error_norms(solve(f, p.Ls, :spectral), uex)
                        @test e.linf < 1e-12 * max(1.0, maximum(abs, uex))
                    end
                end
            end

            @testset "spectral: exponential convergence (non-band-limited u)" begin
                p  = problem("smooth2d")
                ns = [8, 12, 16, 24, 32]
                errs = map(ns) do n
                    f, uex = sample_problem(p, (n, n))
                    error_norms(solve(f, p.Ls, :spectral), uex).linf
                end
                r = _rates(ns, errs)
                verbose && @info "  $name spectral smooth2d: L∞ errors $(round.(errs; sigdigits=3)) at n = $ns, local orders $(round.(r; digits=1))"
                # Exponential, not algebraic: the local algebraic order
                # log(e_i-1/e_i)/log(n_i/n_i-1) of a fixed-order method tends to
                # a constant; for a spectral method it keeps GROWING with n.
                @test issorted(r)
                @test r[end] > 15
                f, uex = sample_problem(p, (48, 48))
                @test error_norms(solve(f, p.Ls, :spectral), uex).linf < 1e-12
            end
        end
    end
end

#-----------------------------------------------------------------------------
# Agreement of TWO solvers of the same discrete system.
#-----------------------------------------------------------------------------
"""
    compare_periodic_poisson_solvers(nameA, solveA, nameB, solveB;
        disc = :fd2, cases = DEFAULT_COMPARISON_CASES, rtol = 1e-8,
        reps = 3, verbose = true) -> rows

For every (problem, grid) in `cases`, solve with both solvers, check that the
gauged solutions agree to `rtol` (relative to max|u|), and report each
solver's error against the exact solution and its best-of-`reps` wall time.
Returns the table rows (NamedTuples) so a caller can post-process them.
"""
const DEFAULT_COMPARISON_CASES = (("modes2d", (64, 64)), ("aniso2d", (96, 48)),
                                  ("smooth2d", (128, 128)), ("smooth2d", (100, 90)),
                                  ("modes3d", (24, 24, 24)))

function _best_time(solve, f, Ls, disc, reps)
    u = solve(f, Ls, disc)                   # also warms up / compiles
    t = Inf
    for _ in 1:reps
        t = min(t, @elapsed solve(f, Ls, disc))
    end
    return u, t
end

function compare_periodic_poisson_solvers(nameA::AbstractString, solveA,
                                          nameB::AbstractString, solveB;
                                          disc::Symbol = :fd2,
                                          cases = DEFAULT_COMPARISON_CASES,
                                          rtol::Real = 1e-8, reps::Int = 3,
                                          verbose::Bool = true)
    rows = NamedTuple[]
    @testset verbose = verbose "$nameA vs $nameB ($disc)" begin
        for (pname, dims) in cases
            p = problem(pname)
            f, uex = sample_problem(p, dims)
            uA, tA = _best_time(solveA, f, p.Ls, disc, reps)
            uB, tB = _best_time(solveB, f, p.Ls, disc, reps)
            diff = maximum(abs, gauge(uA) .- gauge(uB)) / max(maximum(abs, gauge(uA)), eps())
            push!(rows, (problem = pname, dims = dims, npts = prod(dims), diff = diff,
                         errA = error_norms(uA, uex).linf, errB = error_norms(uB, uex).linf,
                         tA = tA, tB = tB))
            @test diff < rtol
        end
    end
    verbose && _print_comparison(nameA, nameB, disc, rows)
    return rows
end

function _print_comparison(nameA, nameB, disc, rows)
    println()
    println("  periodic Poisson, -∇²u = f, discretisation :$disc")
    @printf("  %-9s %-13s %9s  %-10s %-10s  %-10s %-10s  %-10s\n",
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

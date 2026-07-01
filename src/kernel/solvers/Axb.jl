# ---------------------------------------------------------------------------
# Uniform solver timing.
#
# Wrap the PDE-SOLVE kernel of each linear solver (element-learning inference,
# FFT/Fourier solve, Chebyshev solve) in
#
#     u = jx_time_solve("label", () -> <the solve call>)
#
# so all of them report wall-clock time the SAME way: one greppable line
# "# SOLVER TIMING [label]: <seconds> s". This isolates the comparable part of
# each method (the actual solve, not mesh read / setup / IO) so the spectral
# solvers can be timed head-to-head against element learning. The most recent
# elapsed time is also kept in JX_LAST_SOLVE_TIME[] for programmatic use.
#
# It is a plain function (not a macro) on purpose: it is called from
# elementLearningStructs.jl, which is `include`d before this file, and function
# calls resolve at run time (late binding), whereas a macro would need to exist
# at parse time.
# ---------------------------------------------------------------------------
const JX_LAST_SOLVE_TIME = Ref(0.0)

# Last element-learning run's SURROGATE-ONLY inference time (the elementLearning_infer!
# call), i.e. WITHOUT the one-time per-element block assembly that extracts the
# reference-element operators from the sparse matrix A. Set inside
# elementLearning_Axb! and read by the EL diagnostics so "inference" is compared
# against the direct solve on equal footing (matches the original @btime scope).
const JX_LAST_EL_INFER_TIME = Ref(0.0)

# Last solve's accuracy, recorded by the error reporters (fft_report_grid_error,
# the Chebyshev block, print_solution_L2_error). Read by compare_laplace_solvers.
const JX_LAST_SOLVE_ERR = Ref((linf = NaN, l2rel = NaN, npts = 0))

function jx_time_solve(label, f)
    t0  = time_ns()
    val = f()
    dt  = (time_ns() - t0) / 1e9
    JX_LAST_SOLVE_TIME[] = dt
    println(GREEN_FG(string(" # SOLVER TIMING [", label, "]: ", round(dt; sigdigits = 6), " s")))
    return val
end

# ---------------------------------------------------------------------------
# Robust timing via BenchmarkTools (@belapsed).
#
# A single @time/time_ns shot is noisy: it captures whatever GC pause, JIT, or OS
# scheduling hiccup happened on that one run, so repeated measurements scatter
# widely. BenchmarkTools instead runs the operation MANY times (auto-tuned to a
# time budget), excludes compilation, and reports the MINIMUM — which is what
# @btime prints and is far more repeatable.
#
# BenchmarkTools is loaded LAZILY, on first use, rather than with a top-level
# `using`: a package-wide import would add its load cost to the baseline of every
# run and every MPI rank (that baseline cost is exactly why `using BenchmarkTools`
# was removed from src/Jexpresso.jl). Here it is pulled in only when a diagnostics
# timer actually fires.
# ---------------------------------------------------------------------------
const _BENCHMARKTOOLS   = Ref{Union{Module,Nothing}}(nothing)
const _JX_BENCH_TARGET  = Ref{Any}(nothing)   # zero-arg fn the benchmarked expr calls

function _ensure_benchmarktools()
    if _BENCHMARKTOOLS[] === nothing
        m = Base.require(Base.PkgId(
                Base.UUID("6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"), "BenchmarkTools"))
        # Bind the name in THIS module so a later Core.eval can expand the
        # `BenchmarkTools.@belapsed` macro without a compile-time `using`.
        Core.eval(@__MODULE__, :(const BenchmarkTools = $m))
        _BENCHMARKTOOLS[] = m
    end
    return _BENCHMARKTOOLS[]
end

"""
    jx_belapsed(f; seconds = 2.0) -> Float64

Minimum wall-clock time (seconds) of the zero-arg function `f`, measured with
BenchmarkTools' `@belapsed` (many samples, compilation excluded). `f` is stashed
in a module global so the runtime-expanded macro can reach it. Falls back to a
warm best-of-N `time_ns` loop if BenchmarkTools cannot be loaded, so a timing
failure never breaks the run.
"""
function jx_belapsed(f; seconds::Real = 2.0)
    try
        _ensure_benchmarktools()
        _JX_BENCH_TARGET[] = f
        return Core.eval(@__MODULE__,
            :(BenchmarkTools.@belapsed (_JX_BENCH_TARGET[])() seconds = $(Float64(seconds))))::Float64
    catch e
        @warn "BenchmarkTools unavailable; using warm minimum-over-budget time_ns" exception = (e,)
        # Robust fallback: same idea as @belapsed — discard a warm-up run, then
        # sample repeatedly for `seconds` and return the MINIMUM. This is stable
        # run-to-run even without BenchmarkTools (a single shot is what scattered).
        f()                                   # warm-up (discarded)
        best   = Inf
        tstart = time_ns()
        while true
            t0   = time_ns(); f()
            best = min(best, (time_ns() - t0) / 1e9)
            (time_ns() - tstart) >= Float64(seconds) * 1e9 && break
        end
        return best
    end
end

"""
    jx_btime_solve(label, f; seconds = 2.0) -> Float64

Benchmark `f` with `jx_belapsed`, print a `# SOLVER TIMING [label]: … s` line
(the BenchmarkTools minimum), record it in `JX_LAST_SOLVE_TIME[]`, and return it.
The robust counterpart of `jx_time_solve`.
"""
function jx_btime_solve(label, f; seconds::Real = 2.0)
    t = jx_belapsed(f; seconds = seconds)
    JX_LAST_SOLVE_TIME[] = t
    println(GREEN_FG(string(" # SOLVER TIMING [", label, "]: ",
                            round(t; sigdigits = 6), " s  (@btime minimum, compilation excluded)")))
    return t
end

"""
    jx_robust_solve(label, f; robust = true, seconds = 2.0) -> value

Timing wrapper for the Laplace/Poisson comparison solvers (direct SEM, FFT,
Chebyshev, element learning) that, unlike `jx_time_solve`, reports a STABLE
number. When `robust`, it computes `f()` once (the real result, returned to the
caller) and then benchmarks `f` with BenchmarkTools (@belapsed minimum,
compilation excluded) — so re-running the same case reports essentially the same
time instead of the wide scatter a single `@time`/`time_ns` shot gives. When
`robust = false` it falls back to the single-shot `jx_time_solve`.

Callers gate `robust` on `get(inputs, :lbenchmark_solve, true)` and the budget on
`get(inputs, :EL_timing_seconds, 2.0)`.
"""
function jx_robust_solve(label, f; robust::Bool = true, seconds::Real = 2.0)
    robust || return jx_time_solve(label, f)
    val = f()                                   # real result for the caller
    t   = jx_belapsed(f; seconds = seconds)     # stable @btime minimum
    JX_LAST_SOLVE_TIME[] = t
    println(GREEN_FG(string(" # SOLVER TIMING [", label, "]: ",
                            round(t; sigdigits = 6), " s  (@btime minimum, compilation excluded)")))
    return val
end

# Stash the most recent solve's error norms for the side-by-side comparison.
function jx_record_solve_error(; linf = NaN, l2rel = NaN, npts = 0)
    JX_LAST_SOLVE_ERR[] = (linf = Float64(linf), l2rel = Float64(l2rel), npts = Int(npts))
    return nothing
end

function solveAx(L, RHS, linear_solver...)
    
    prob = LinearProblem(L, RHS);    
    sol = solve(prob, linear_solver)
    
    return sol
end


#function solveAx_sparse(A, b)
#     return A\b
#end

function standard_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    RHS   = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))
    Mdiag = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))
    
    #-----------------------------------------------------
    # L*q = M*RHS   See algo 12.18 of Giraldo's book
    #-----------------------------------------------------
    for ip =1:sem.mesh.npoin
        RHS[ip] = user_source!(RHS[ip],
                               params.qp.qn[ip],
                               params.qp.qe[ip],
                               sem.mesh.npoin,
                               inputs[:CL], inputs[:SOL_VARS_TYPE];
                               neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                               xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                               ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
    end
    RHS = sem.matrix.M.*RHS

    if inputs[:lsparse] ==  false
        for ip = 1:sem.mesh.npoin
            sem.matrix.L[ip,ip] += inputs[:rconst][1]
        end
    end

    apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                         0.0, params.qp.qe,
                                         params.mesh.coords,
                                         params.metrics.nx,
                                         params.metrics.ny,
                                         params.metrics.nz,
                                         sem.mesh.npoin,
                                         params.mesh.npoin_linear,
                                         params.mesh.poin_in_bdy_edge,
                                         params.mesh.poin_in_bdy_face,
                                         params.mesh.nedges_bdy,
                                         params.mesh.nfaces_bdy,
                                         params.mesh.ngl, params.mesh.ngr,
                                         params.mesh.nelem_semi_inf,
                                         params.basis.ψ, params.basis.dψ,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         RHS, 0.0, params.ubdy,
                                         params.mesh.connijk_lag,
                                         params.mesh.bdy_edge_in_elem,
                                         params.mesh.bdy_edge_type,
                                         params.ω, qp.neqs,
                                         params.inputs, params.AD, sem.mesh.SD)

    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage ..............")))

    #solAxb = @btime solveAx($sem.matrix.L, $RHS, inputs[:ode_solver])
    #sol = solAxb.u
    # Robust @btime timing (BenchmarkTools minimum) so re-running the same case
    # reports a stable time instead of the wide scatter of a single shot. Set
    # :lbenchmark_solve => false for one quick solve on very large problems.
    solAxb = jx_robust_solve("direct SEM (Ax=b)", () -> sem.matrix.L \ RHS;
                             robust  = get(inputs, :lbenchmark_solve, true),
                             seconds = Float64(get(inputs, :EL_timing_seconds, 2.0)))
    sol = solAxb

    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage .............. DONE")))
    args = (params.SD, sol, params.uaux, 1, 1,
            sem.mesh, nothing,
            nothing, nothing,
            0.0, 0.0, 0.0,
            OUTPUT_DIR, inputs,
            params.qp.qvars,
            params.qp.qoutvars,
            inputs[:outformat])

    # Automatic L2-error check against the exact field qe (e.g. a manufactured
    # solution). No-op when qe carries no exact field.
    print_solution_L2_error(sol, params.qp.qe, sem.matrix.M, sem.mesh.npoin;
                            label="direct solve")

    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe, metrics=params.metrics)

    return nothing
end

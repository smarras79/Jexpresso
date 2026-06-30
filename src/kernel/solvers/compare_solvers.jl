# =============================================================================
#  Side-by-side comparison of the Laplace/Poisson solvers
# =============================================================================
#
#  compare_laplace_solvers() runs the three comparison cases — all on the SAME
#  domain [-π,π]², the SAME manufactured solution u_ex = sin(x)cos(y), AND the
#  SAME number of degrees of freedom (64 points/dir = 4096 nodes) so the solve
#  timings are comparable — and prints ONE table of time and accuracy:
#
#      • Elliptic/laplace_periodic    — FFT / Fourier        (periodic)
#      • Elliptic/laplace_chebyshev   — Chebyshev collocation (Dirichlet)
#      • Elliptic/elementLearning_2pi — element learning / SEM (Dirichlet)
#
#  Each solver records its solve time (JX_LAST_SOLVE_TIME) and error norms
#  (JX_LAST_SOLVE_ERR) through jx_time_solve / jx_record_solve_error; this driver
#  reads them after each run. Each case is run `reps` times (default 2) and the
#  FASTEST solve is reported, so the one-time JIT compilation of the solve kernel
#  is not counted in the comparison.
#
#  Usage (from the Julia REPL, after `using Jexpresso`):
#      Jexpresso.compare_laplace_solvers()
#      Jexpresso.compare_laplace_solvers(reps = 3)
#
#  NOTE: as shipped, elementLearning_2pi runs the DIRECT SEM solve (its row is
#  labelled accordingly). To compare the element-learning INFERENCE instead,
#  enable :lelementLearning (+ a trained :NNfile) in that case — the table then
#  reports the EL inference time.
# =============================================================================

function compare_laplace_solvers(; reps::Int = 2,
        cases = [("FFT (Fourier)",          "Elliptic", "laplace_periodic"),
                 ("Chebyshev collocation",  "Elliptic", "laplace_chebyshev"),
                 ("element learning / SEM", "Elliptic", "elementLearning_2pi")])

    rows = Vector{NamedTuple}()
    for (label, eq, case) in cases
        times = Float64[]
        err   = (linf = NaN, l2rel = NaN, npts = 0)
        for r in 1:max(1, reps)
            try
                run_case(eq, case)
            catch e
                @warn "compare_laplace_solvers: case '$case' failed; skipping" exception = (e, catch_backtrace())
                break
            end
            push!(times, JX_LAST_SOLVE_TIME[])
            err = JX_LAST_SOLVE_ERR[]
        end
        push!(rows, (label = label, case = case,
                     seconds = isempty(times) ? NaN : minimum(times),
                     linf = err.linf, l2rel = err.l2rel, npts = err.npts))
    end

    _print_laplace_comparison(rows, reps)
    return rows
end

_fmt_num(x; sig = 4) = (x isa Real && isfinite(x)) ? string(round(x; sigdigits = sig)) : "n/a"

function _print_laplace_comparison(rows, reps)
    println()
    println(GREEN_FG(" ===================== Laplace solver comparison ====================="))
    println(GREEN_FG(string("   domain [-π,π]²,  u_ex = sin(x)·cos(y),  fastest of ", reps,
                            " run", reps == 1 ? "" : "s")))
    println("   ", rpad("method", 26), rpad("grid (nodes)", 14),
            rpad("solve time [s]", 16), rpad("‖e‖∞", 13), "rel ‖e‖₂")
    println("   ", repeat("─", 76))
    for r in rows
        println("   ", rpad(r.label, 26),
                rpad(r.npts == 0 ? "n/a" : string(r.npts), 14),
                rpad(_fmt_num(r.seconds), 16),
                rpad(_fmt_num(r.linf), 13),
                _fmt_num(r.l2rel))
    end
    println(GREEN_FG(" ====================================================================="))
    println()
    return nothing
end

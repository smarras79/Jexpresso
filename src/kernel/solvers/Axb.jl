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
# Warm-up + best-of-reps timing.
#
# The FIRST evaluation of `f` in a fresh Julia process pays one-time costs that
# do NOT reflect steady-state performance — Julia JIT-compiles the solve kernel,
# a sparse factorization builds its symbolic structure, an ONNX InferenceSession
# warms up. Timing that first call therefore over-reports the true solve cost.
#
# This variant runs `f` once as a throw-away warm-up (time discarded), then runs
# it `reps` more times and reports the FASTEST — the same "measure the second
# run" strategy compare_laplace_solvers uses across run_case repetitions, applied
# here at the individual-solve level. The best time is stored in
# JX_LAST_SOLVE_TIME[] and the last computed value is returned. Set warmup=false
# / reps=1 to recover a single-shot measurement.
# ---------------------------------------------------------------------------
function jx_time_solve_best(label, f; reps::Int = 2, warmup::Bool = true)
    if warmup
        f()                    # throw-away warm-up: compiles/initialises, not timed
    end
    n    = max(1, reps)
    best = Inf
    val  = nothing
    for _ in 1:n
        t0   = time_ns()
        val  = f()
        best = min(best, (time_ns() - t0) / 1e9)
    end
    JX_LAST_SOLVE_TIME[] = best
    note = warmup ? string(" (best of ", n, ", JIT/warm-up excluded)") :
                    string(" (best of ", n, ")")
    println(GREEN_FG(string(" # SOLVER TIMING [", label, "]: ",
                            round(best; sigdigits = 6), " s", note)))
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
    # Timed the same way as the FFT/Chebyshev/EL solvers (jx_time_solve) so the
    # direct SEM solve is comparable in compare_laplace_solvers.
    solAxb = jx_time_solve("direct SEM (Ax=b)", () -> sem.matrix.L \ RHS)
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

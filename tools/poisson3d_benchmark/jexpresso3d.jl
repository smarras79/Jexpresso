#=============================================================================
 tools/poisson3d_benchmark/jexpresso3d.jl — one configuration of the 3D
 periodic Poisson study run THROUGH JEXPRESSO:

     Jexpresso.run_case("Elliptic", "poisson_periodic_sem_3d"; inputs = overrides)

 the 3D deck (built-in Cartesian grid, ne³ hexahedra on [0,2π]³, periodic in
 x, y, z; the exact solution of poisson3d.jl with r = 0.5) and Jexpresso's
 own solvers, selected by the same deck flags as in 2D:

   :sem        sem_setup (DSS_laplace_sparse_3D) → periodic_sem_system →
               periodic_sem_factorize / periodic_sem_direct_solve
   :sem_amg    … → periodic_sem_amg_solve (jx_amg_setup / jx_amg_solve)
   :sc_direct  … → periodic_sem_sc_solve (elementLearning_Axb!, skeleton by
               el_skeleton_solve, CHOLMOD Cholesky)
   :sc_amg     the same, skeleton by AMG + CG
   :ps         pseudospectral_linsolve! (FourierCollocationPoissonSolver3D)
   :fft        fft_linsolve! (FFTPoissonSolver, FFTW.ESTIMATE)

 The row has the columns of poisson3d.jl's run_config, from Jexpresso's
 per-phase timers (JX_TIMINGS) of that run:
   assembly = :sem_setup (mesh, basis, metrics, mass and Laplacian), SEM only
   rhs, setup, solve as recorded by the solvers.
 The sparse direct solves (SEM direct, SC direct) factorise with METIS nested
 dissection (deck option :sparse_ordering => "metis"; ordering = :amd gives
 Jexpresso's default, CHOLMOD's own AMD ordering), as the Kronecker engine.
=============================================================================#
module JX3D

using Jexpresso

const SOLVERS = (:sem, :sem_amg, :pmg_amg, :pmg_gmg, :sc_direct, :sc_amg, :ps, :fft)
const R = 0.5                         # the deck's PPB3_R

function overrides(solver::Symbol, ne::Int, N::Int; rtol = 1e-12, amg_method = "sa", itmax = 100_000,
                   ordering::Symbol = :metis)
    return Dict{Symbol, Any}(
        :nop => N, :nelx => ne, :nely => ne, :nelz => ne,
        :fft_N => ne * N,                    # Fourier grids: (ne·N)³, the SEM's unknowns
        :lfft => solver === :fft, :lpseudospectral => solver === :ps,
        :linsolve_amg => solver === :sem_amg,
        :linsolve_pmg => solver === :pmg_amg ? "amg" : solver === :pmg_gmg ? "gmg" : "none",
        :lstatic_condensation => solver in (:sc_direct, :sc_amg),
        :EL_skeleton_solver => solver === :sc_amg ? "amg" : "direct",
        :amg_method => amg_method, :amg_rtol => rtol, :amg_itmax => itmax,
        :fft_plan => "estimate",
        :sparse_ordering => ordering === :metis ? "metis" : "cholmod",
        :luse_mesh_cache => false, :lbenchmark_solve => false, :outformat => "none")
end

"""
    run_config(solver, ne, N; rtol = 1e-12, amg_method = "sa") -> row

Run the 3D deck once with Jexpresso's `solver` on ne³ elements of order N
(the Fourier solvers on the (ne N)³ grid) and return the row.
"""
function run_config(solver::Symbol, ne::Int, N::Int; rtol = 1e-12, amg_method = "sa", ordering::Symbol = :metis)
    solver in SOLVERS || error("Jexpresso 3D solvers: $SOLVERS (got $solver)")
    Jexpresso.run_case("Elliptic", "poisson_periodic_sem_3d";
                       inputs = overrides(solver, ne, N; rtol = rtol, amg_method = amg_method,
                                          ordering = ordering))
    T = Jexpresso.JX_TIMINGS; err = Jexpresso.JX_LAST_SOLVE_ERR[]
    g(k) = Float64(get(T, k, 0.0))
    sem = solver in (:sem, :sem_amg, :pmg_amg, :pmg_gmg, :sc_direct, :sc_amg)
    Ng = ne * N; n = Ng^3
    solved = solver in (:sc_direct, :sc_amg) ? n - (ne * (N - 1))^3 : (sem ? n - 1 : n)
    asm = sem ? g(:sem_setup) : 0.0
    return (solver = solver, d = 3, r = R, ne = ne, nop = N, Ng = Ng, n = n, solved = solved,
            linf = err.linf, l2rel = err.l2rel,
            assembly = asm, rhs = g(:rhs), setup = g(:setup), solve = g(:solve),
            total = asm + g(:rhs) + g(:setup) + g(:solve),
            iters = solver in (:sem_amg, :sc_amg, :pmg_amg, :pmg_gmg) ? Jexpresso.JX_AMG_STATS[].iters : 0,
            nnz = 0, factor_nnz = 0, skeleton_nnz = 0,
            ordering = solver in (:sem, :sc_direct) ? ordering : :none)
end

end # module

# =============================================================================
#  Pseudo-spectral (Fourier collocation) solver for the periodic Laplace /
#  Poisson equation — the Jexpresso driver
# =============================================================================
#
#  Solves -∇²u = f on a PERIODIC rectangle with the Fourier collocation solver
#  of fourier_collocation.jl (Kopriva's FourierDerivativeMatrix, dense
#  physical-space operators, direct solve by matrix diagonalisation). Selected
#  in problems/drivers.jl by :lpseudospectral => true, next to :lfft (FFT) and
#  the SEM solve; it reads the SAME problem definition as the FFT solver:
#       user_fft_rhs(x,y)     the right-hand side f          (REQUIRED)
#       user_fft_exact(x,y)   the exact solution u_ex        (optional; enables
#                             the automatic error verification)
#
#  Case-deck inputs (user_inputs.jl); each falls back to the FFT's key:
#       :ps_N, :ps_M          collocation points per direction (EVEN; :ps_M
#                             defaults to :ps_N)            [default :fft_N, 64]
#       :ps_Lx, :ps_Ly        periods                       [:fft_Lx/:fft_Ly, 2π]
#       :ps_x0, :ps_y0        lower-left corner             [:fft_x0/:fft_y0, 0]
#
#  The reported SOLVER TIMING is the solve only (four dense N×N products); the
#  one-time eigen-decomposition of the 1-D operators is setup, like a
#  factorisation, and is excluded.
# =============================================================================

function pseudospectral_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    inputs[:backend] == CPU() ||
        error(" # pseudospectral_linsolve!: the pseudo-spectral solver is CPU-only.")
    isdefined(@__MODULE__, :user_fft_rhs) ||
        error(" # pseudospectral_linsolve!: define user_fft_rhs(x,y) in the case's user_source.jl")
    has_exact = isdefined(@__MODULE__, :user_fft_exact)

    N  = Int(get(inputs, :ps_N, get(inputs, :fft_N, 64)))
    M  = Int(get(inputs, :ps_M, get(inputs, :fft_M, N)))
    Lx = Float64(get(inputs, :ps_Lx, get(inputs, :fft_Lx, 2π)))
    Ly = Float64(get(inputs, :ps_Ly, get(inputs, :fft_Ly, 2π)))
    x0 = Float64(get(inputs, :ps_x0, get(inputs, :fft_x0, 0.0)))
    y0 = Float64(get(inputs, :ps_y0, get(inputs, :fft_y0, 0.0)))
    (N > 0 && M > 0 && iseven(N) && iseven(M)) ||
        error(" # pseudospectral_linsolve!: :ps_N=$N, :ps_M=$M must be positive and even.")
    x, y = periodic_grid_lines((N, M), (Lx, Ly), (x0, y0))

    println(YELLOW_FG(string(" # Solve -∇²u = f by pseudo-spectral Fourier collocation: ",
                             N, "×", M, " periodic grid ..............")))

    F = jx_phase(() -> _fft_sample_rhs((x, y)), :rhs)
    S = jx_phase(() -> FourierCollocationPoissonSolver((N, M), (Lx, Ly)), :setup)
    u = similar(F)
    jx_robust_solve("pseudo-spectral (Fourier collocation) solve",
                    () -> fourier_collocation_poisson_solve!(u, S, F);
                    robust  = get(inputs, :lbenchmark_solve, true),
                    seconds = Float64(get(inputs, :EL_timing_seconds, 2.0)))
    if abs(S.fmean[]) > 1e-10 * max(1.0, maximum(abs, F))
        println(string(" # pseudospectral_linsolve!: RHS mean = ", S.fmean[],
                       " ≠ 0; solved the projected problem -∇²u = f - mean(f) ",
                       "(periodic compatibility condition)."))
    end

    println(YELLOW_FG(string(" # Solve -∇²u = f by pseudo-spectral Fourier collocation ........... DONE")))

    uex = nothing; err = nothing
    has_exact && ((uex, err) = fft_report_grid_error(u, (x, y), (Lx, Ly);
                                                     label = "pseudo-spectral solve"))

    if !(inputs[:outformat] isa NONE)          # :outformat => "none" skips the file
        vtkpath = joinpath(OUTPUT_DIR, "pseudospectral_laplace.vtk")
        write_fft_vtk(vtkpath, x, y, u, uex, err)
        println(string(" # pseudo-spectral solution written to ", vtkpath))
    end
    return u
end

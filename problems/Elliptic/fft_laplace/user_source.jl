# =============================================================================
#  Manufactured periodic problem for the classical FFT Laplace/Poisson solver
# =============================================================================
#
#  Domain [0,2π]² with periodic BCs. Exact (manufactured) solution
#
#       u_ex(x,y) = sin(2x) cos(3y) + sin(x) cos(y)
#
#  is 2π-periodic and zero-mean, so it is the unique zero-mean solution of the
#  periodic Poisson problem  -∇²u = f  with
#
#       f(x,y) = -∇²u_ex = (2²+3²) sin(2x)cos(3y) + (1²+1²) sin(x)cos(y)
#              = 13 sin(2x)cos(3y) + 2 sin(x)cos(y).
#
#  Each term is a single Fourier mode, so the FFT solver reproduces u_ex to
#  machine precision (spectral exactness). These two functions are what
#  fft_linsolve! (src/kernel/solvers/fft_laplace.jl) samples on its grid:
#       user_fft_rhs   → the right-hand side f          (REQUIRED)
#       user_fft_exact → the exact solution u_ex        (enables the L2 check)
# =============================================================================

user_fft_exact(x, y) = sin(2x)*cos(3y) + sin(x)*cos(y)

user_fft_rhs(x, y)   = 13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y)


# ---------------------------------------------------------------------------
# Generic source hook required by the setup pipeline (unused by the FFT solve,
# which gets its RHS from user_fft_rhs above). Kept so the case is structurally
# identical to the other Elliptic problems.
# ---------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=1.0, y=1.0,
                      xmax=1.0, xmin=0.0,
                      ymax=1.0, ymin=0.0)
    return 0.0
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    return T(0.0)
end

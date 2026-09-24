# =============================================================================
#  Doubly periodic Poisson problem shared by the SEM and the FFT solvers
# =============================================================================
#
#  -∇²u = f  on [0,2π]², periodic in x and y, with the manufactured solution
#
#       u_ex(x,y) = sin(2x) cos(3y) + sin(x) cos(y)
#       f(x,y)    = -∇²u_ex = 13 sin(2x) cos(3y) + 2 sin(x) cos(y)
#
#  u_ex is zero-mean, so it is THE zero-mean periodic solution — the one both
#  solvers return (the periodic Laplacian determines u only up to a constant).
#  The same problem as Elliptic/fft_laplace, so the two runs are comparable.
#
#  One definition serves both solvers:
#       user_fft_rhs / user_fft_exact  → the FFT solver (:lfft => true)
#       user_source!                   → the SEM solver (:lfft => false)
#       initialize.jl sets qe = user_fft_exact → the SEM L2 error check
# =============================================================================

user_fft_exact(x, y) = sin(2x)*cos(3y) + sin(x)*cos(y)

user_fft_rhs(x, y)   = 13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y)

# SEM right-hand side: standard_linsolve! assembles  L u = M f  with L the
# (positive) stiffness matrix of -∇², so f enters with a + sign.
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=1.0, y=1.0,
                      xmax=1.0, xmin=0.0,
                      ymax=1.0, ymin=0.0)
    return user_fft_rhs(x, y)
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    return T(13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y))
end

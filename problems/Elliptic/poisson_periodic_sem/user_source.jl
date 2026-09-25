# =============================================================================
#  Doubly periodic Poisson problem shared by the SEM, pseudo-spectral and FFT
#  solvers
# =============================================================================
#
#  -∇²u = f  on [0,2π]², periodic in x and y, with the manufactured solution
#
#       u_ex(x,y) = A ( p(x) p(y) - 1/(c²-1) ),     p(s) = 1/(c - cos s),
#
#  the product of two periodic Poisson kernels. p has the Fourier series
#
#       p(s) = (1/√(c²-1)) (1 + 2 Σ_{k≥1} r^k cos ks),   r = c - √(c²-1) < 1,
#
#  so the Fourier coefficients of u decay like r^(|kx|+|ky|): geometrically,
#  NOT band-limited. Unlike a trigonometric polynomial, which a Fourier solver
#  reproduces to round-off as soon as the grid resolves its highest mode, this
#  u gives the FFT and pseudo-spectral solvers a grid-dependent error, about
#  r^(N_g/2) on an N_g × N_g grid (the aliased tail). PPB_R sets how hard the
#  problem is; with r = 0.8 the FFT error falls from O(1) on the 32² grid to the
#  round-off floor near N_g ≈ 300. u peaks sharply at the origin (width
#  ~√(2(c-1)) ≈ 0.22), which the SEM has to resolve too.
#
#  Normalisation: ∫u = 0 (the mean of p is 1/√(c²-1), so the mean of p(x)p(y)
#  is 1/(c²-1)), which makes u_ex THE zero-mean periodic solution — the one
#  every solver returns — and A = (c-1)²(c+1)/2 scales the peak to u_ex(0,0) = 1
#  (its minimum, at (π,π), is A(1/(c+1)² - 1/(c²-1)) ≈ -0.012).
#
#  f = -∇²u = -A ( p''(x) p(y) + p(x) p''(y) ),
#       p''(s) = -cos s/(c - cos s)² + 2 sin² s/(c - cos s)³ .
#
#  One definition serves every solver:
#       user_fft_rhs / user_fft_exact  → the FFT and pseudo-spectral solvers
#       user_source!                   → the SEM solvers
#       initialize.jl sets qe = user_fft_exact → the SEM error check
#  (Before this problem the deck used the trigonometric polynomial
#   u = sin 2x cos 3y + sin x cos y, still used by Elliptic/fft_laplace and
#   Elliptic/poisson_dirichlet_sc.)
# =============================================================================

const PPB_R = 0.8                                   # Fourier decay rate of p
const PPB_C = (PPB_R + 1 / PPB_R) / 2               # c = (r + 1/r)/2  (r = c - √(c²-1))
const PPB_A = (PPB_C - 1)^2 * (PPB_C + 1) / 2       # u_ex(0,0) = 1

_ppb_p(s)   = 1 / (PPB_C - cos(s))
_ppb_pss(s) = (d = PPB_C - cos(s); -cos(s) / d^2 + 2 * sin(s)^2 / d^3)

user_fft_exact(x, y) = PPB_A * (_ppb_p(x) * _ppb_p(y) - 1 / (PPB_C^2 - 1))

user_fft_rhs(x, y)   = -PPB_A * (_ppb_pss(x) * _ppb_p(y) + _ppb_p(x) * _ppb_pss(y))

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
    return T(user_fft_rhs(x, y))
end

# =============================================================================
#  Triply periodic Poisson problem: the 3D version of poisson_periodic_sem
# =============================================================================
#
#  -∇²u = f  on [0,2π]³, periodic in x, y and z, with the manufactured solution
#
#       u_ex = A ( p(x) p(y) p(z) - (c²-1)^(-3/2) ),   p(s) = 1/(c - cos s),
#
#  c = (r + 1/r)/2. Its Fourier coefficients decay like r^(|kx|+|ky|+|kz|):
#  geometrically, NOT band-limited, so the FFT and pseudo-spectral solvers have
#  a grid-dependent error too. (c²-1)^(-3/2) is the mean of p(x)p(y)p(z), so
#  u_ex has zero mean (the gauge every solver returns); A scales the peak to
#  u_ex(0,0,0) = 1. r = 0.5 here (the 2D deck uses 0.8): the product of three
#  sharper peaks is not resolved at feasible 3D sizes.
#
#  f = -∇²u = -A ( p''(x)p(y)p(z) + p(x)p''(y)p(z) + p(x)p(y)p''(z) ),
#       p''(s) = -cos s/(c - cos s)² + 2 sin² s/(c - cos s)³ .
#
#  user_fft_rhs / user_fft_exact → the FFT and pseudo-spectral solvers;
#  user_source! → the SEM solvers; initialize.jl sets qe = user_fft_exact.
# =============================================================================

const PPB3_R = 0.5
const PPB3_C = (PPB3_R + 1 / PPB3_R) / 2
const PPB3_M = (PPB3_C^2 - 1)^(-3 / 2)                     # mean of p(x)p(y)p(z)
const PPB3_A = 1 / ((PPB3_C - 1)^(-3) - PPB3_M)            # u_ex(0,0,0) = 1

_ppb3_p(s)   = 1 / (PPB3_C - cos(s))
_ppb3_pss(s) = (d = PPB3_C - cos(s); -cos(s) / d^2 + 2 * sin(s)^2 / d^3)

user_fft_exact(x, y, z) = PPB3_A * (_ppb3_p(x) * _ppb3_p(y) * _ppb3_p(z) - PPB3_M)

user_fft_rhs(x, y, z) = -PPB3_A * (_ppb3_pss(x) * _ppb3_p(y) * _ppb3_p(z) +
                                   _ppb3_p(x) * _ppb3_pss(y) * _ppb3_p(z) +
                                   _ppb3_p(x) * _ppb3_p(y) * _ppb3_pss(z))

# SEM right-hand side: the linear solve assembles L u = M f with L the
# (positive) stiffness matrix of -∇², so f enters with a + sign.
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=1.0, y=1.0, z=1.0,
                      xmax=1.0, xmin=0.0,
                      ymax=1.0, ymin=0.0,
                      zmax=1.0, zmin=0.0)
    return user_fft_rhs(x, y, z)
end

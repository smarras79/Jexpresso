# =============================================================================
#  Periodic Poisson problem  -∇²u = f  on [-π,π]²  (classical FFT solver)
# =============================================================================
#
#  Real (NON-manufactured) right-hand side: a smooth periodic dipole — a
#  positive and a negative Gaussian bump. This is a genuine source/sink Poisson
#  problem with no closed-form solution; the FFT solver still returns the exact
#  discrete periodic solution (up to the usual additive constant). The small
#  residual mean of f on the finite domain is removed automatically by
#  fft_linsolve! to satisfy the periodic compatibility condition ∫f = 0.
#
#  fft_linsolve! (src/kernel/solvers/fft_laplace.jl) samples user_fft_rhs on the
#  mesh grid. user_fft_exact is intentionally NOT defined here (no manufactured
#  solution ⇒ no L2-error print). To turn this into a verification case instead,
#  uncomment the manufactured block at the bottom.
# =============================================================================

function user_fft_rhs(x, y)
    A  = 5.0
    σ2 = 2*0.4^2
    xp, yp = -π/2, 0.0          # source (+)
    xm, ym =  π/2, 0.0          # sink   (−)
    return A*( exp(-((x-xp)^2 + (y-yp)^2)/σ2) - exp(-((x-xm)^2 + (y-ym)^2)/σ2) )
end


# ---------------------------------------------------------------------------
# Generic source hook required by the setup pipeline (unused by the FFT solve).
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


# ---------------------------------------------------------------------------
# OPTIONAL manufactured verification: a periodic exact solution on [-π,π]².
# Uncomment BOTH functions (and they override the dipole RHS above) to make
# fft_linsolve! print the L2/L∞ error against u_ex.
#
#   u_ex(x,y) = sin(2x)cos(3y) + sin(x)cos(y)          (2π-periodic, zero mean)
#   f(x,y)    = -∇²u_ex = 13 sin(2x)cos(3y) + 2 sin(x)cos(y)
# ---------------------------------------------------------------------------
# user_fft_rhs(x, y)   = 13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y)
# user_fft_exact(x, y) = sin(2x)*cos(3y) + sin(x)*cos(y)

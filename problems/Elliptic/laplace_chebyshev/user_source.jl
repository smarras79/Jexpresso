# =============================================================================
#  Non-periodic Laplace/Poisson problem for the Chebyshev spectral solver
# =============================================================================
#
#  As shipped this is a TRUE LAPLACE problem: u = exp(x)·cos(y) is harmonic
#  (∇²u = u_xx + u_yy = e^x cos y - e^x cos y = 0), so f = 0 and the solution is
#  driven entirely by the (non-homogeneous) Dirichlet boundary data g = u. The
#  Chebyshev collocation solver recovers it to spectral accuracy.
#
#  cheb_linsolve! (src/kernel/solvers/cheb_laplace.jl) samples these on the CGL
#  grid:
#       user_cheb_rhs(x,y)    the right-hand side f               (REQUIRED)
#       user_cheb_bc(x,y)     the Dirichlet boundary value g      (optional)
#       user_cheb_exact(x,y)  the exact solution u_ex             (optional; here
#                             it also supplies g, since user_cheb_bc is omitted)
# =============================================================================

user_cheb_exact(x, y) = exp(x)*cos(y)      # harmonic ⇒ also the Dirichlet data g

user_cheb_rhs(x, y)   = 0.0                # -∇²u = 0  (Laplace)


# ---------------------------------------------------------------------------
# Generic source hook required by the setup pipeline (unused by the Chebyshev
# solve, which gets its data from the user_cheb_* functions above).
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
# OPTIONAL Poisson (non-zero f) manufactured check. Uncomment to replace the
# Laplace case above with  u_ex = sin(πx)sin(πy)  (homogeneous Dirichlet) and
# f = -∇²u_ex = 2π² sin(πx)sin(πy):
# ---------------------------------------------------------------------------
# user_cheb_exact(x, y) = sin(π*x)*sin(π*y)
# user_cheb_rhs(x, y)   = 2*π^2*sin(π*x)*sin(π*y)

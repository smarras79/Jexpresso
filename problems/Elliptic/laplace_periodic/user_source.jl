# =============================================================================
#  Periodic Poisson problem  -∇²u = f  on [-π,π]ᵈ  (classical FFT solver)
# =============================================================================
#
#  Real (NON-manufactured) right-hand side: a smooth periodic dipole — a positive
#  and a negative Gaussian bump. A genuine source/sink Poisson problem with no
#  closed-form solution; the FFT solver still returns the exact discrete periodic
#  solution (up to the usual additive constant). The small residual mean of f on
#  the finite domain is removed automatically by fft_linsolve! to satisfy the
#  periodic compatibility condition ∫f = 0.
#
#  user_fft_rhs takes (x, y, z) with z defaulting to 0, so the SAME case runs on
#  a 2D mesh (the solver calls it with two args) and on a 3D mesh (three args).
#  fft_linsolve! (src/kernel/solvers/fft_laplace.jl) samples it on the mesh grid.
#  user_fft_exact is intentionally NOT defined here (no manufactured solution ⇒
#  no L2-error print). Uncomment the manufactured block at the bottom to verify.
# =============================================================================

function user_fft_rhs(x, y, z = 0.0)
    A  = 5.0
    σ2 = 2*0.4^2
    # dipole separated along x; on a 3D mesh the bumps are localized near z=0 too
    rp = (x + π/2)^2 + y^2 + z^2      # source (+) at (-π/2, 0, 0)
    rm = (x - π/2)^2 + y^2 + z^2      # sink   (−) at (+π/2, 0, 0)
    return A*( exp(-rp/σ2) - exp(-rm/σ2) )
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
# OPTIONAL manufactured verification: a periodic exact solution.
# Uncomment the pair for your mesh dimension to make fft_linsolve! print the
# L2/L∞ error against u_ex (each term is a single Fourier mode ⇒ exact to machine
# precision). -∇²u_ex for mode (kx,ky[,kz]) multiplies by (kx²+ky²[+kz²]).
#
# 2D (on [-π,π]²):
#   user_fft_exact(x, y)    = sin(2x)*cos(3y) + sin(x)*cos(y)
#   user_fft_rhs(x, y, z=0) = 13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y)
#
# 3D (on [-π,π]³):
#   user_fft_exact(x, y, z) = sin(x)*cos(y)*cos(z)
#   user_fft_rhs(x, y, z)   = 3.0*sin(x)*cos(y)*cos(z)
# ---------------------------------------------------------------------------

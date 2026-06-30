# =============================================================================
#  FFT (Fourier) member of the 3-method comparison trio.
#  Periodic Laplace problem  -∇²u = f  on [-π,π]²  (classical FFT solver).
# =============================================================================
#
#  Shared manufactured solution (same as laplace_chebyshev and elementLearning_2pi):
#
#       u_ex(x,y) = sin(x)·cos(y) ,   f = -∇²u_ex = 2·sin(x)·cos(y).
#
#  sin(x)cos(y) is 2π-periodic and zero-mean on [-π,π]², so it is exactly the
#  periodic solution recovered by the Fourier solver (it is a single Fourier mode
#  ⇒ machine-precision error). The other two cases impose it as Dirichlet data;
#  here it is carried by periodicity. Same domain, same PDE, same exact solution —
#  only the grid / BC treatment differs between the three methods.
#
#  fft_linsolve! (src/kernel/solvers/fft_laplace.jl) samples these on the grid.
#  user_fft_rhs takes (x,y,z) with z defaulting to 0 so the case also runs on a
#  3D mesh unchanged.
# =============================================================================

user_fft_exact(x, y, z = 0.0) = sin(x)*cos(y)

user_fft_rhs(x, y, z = 0.0)   = 2.0*sin(x)*cos(y)


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

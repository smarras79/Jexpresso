# =============================================================================
#  Chebyshev member of the 3-method comparison trio.
#  Non-periodic Laplace problem on [-π,π]² (Dirichlet), Chebyshev collocation.
# =============================================================================
#
#  Shared manufactured solution (same as laplace_periodic and elementLearning_2pi):
#
#       u_ex(x,y) = sin(x)·cos(y) ,   f = -∇²u_ex = 2·sin(x)·cos(y).
#
#  Here u_ex is imposed as the (non-homogeneous) Dirichlet boundary data g = u_ex
#  on ∂[-π,π]²; the FFT case carries the same u_ex by periodicity instead. Same
#  domain, same PDE, same exact solution — only the grid / BC treatment differs.
#
#  cheb_linsolve! (src/kernel/solvers/cheb_laplace.jl) samples these on the CGL
#  grid:
#       user_cheb_rhs(x,y)    the right-hand side f               (REQUIRED)
#       user_cheb_exact(x,y)  the exact solution u_ex (also the Dirichlet g,
#                             since user_cheb_bc is omitted)
# =============================================================================

user_cheb_exact(x, y) = sin(x)*cos(y)       # also supplies the Dirichlet data g

user_cheb_rhs(x, y)   = 2.0*sin(x)*cos(y)   # -∇²u_ex


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

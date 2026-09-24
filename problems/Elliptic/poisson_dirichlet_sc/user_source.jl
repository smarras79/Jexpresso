# =============================================================================
#  Dirichlet Poisson problem for the static-condensation / AMG comparison
# =============================================================================
#
#  -∇²u = f  on [0,2π]²,  u = u_ex on the boundary (Dirichlet), with the SAME
#  manufactured solution as Elliptic/poisson_periodic_sem:
#
#       u_ex(x,y) = sin(2x) cos(3y) + sin(x) cos(y)
#       f(x,y)    = -∇²u_ex = 13 sin(2x) cos(3y) + 2 sin(x) cos(y)
#
#  The static condensation of element learning (elementLearning_Axb!) is built
#  for Dirichlet problems: the skeleton splits into the boundary Γ (data g)
#  and the internal skeleton ∂O (unknowns). This case is the periodic
#  benchmark's problem on the same 16×16 mesh, with Dirichlet data instead of
#  periodicity, so the full and the condensed SEM solves can be compared.
#
#  user_fft_exact / user_fft_rhs keep the names of the periodic case so the
#  two decks share initialize.jl.
# =============================================================================

user_fft_exact(x, y) = sin(2x)*cos(3y) + sin(x)*cos(y)

user_fft_rhs(x, y)   = 13.0*sin(2x)*cos(3y) + 2.0*sin(x)*cos(y)

# SEM right-hand side of L u = M f (L = stiffness of -∇²).
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

# ============================================================================
#  NON-CONSTANT DIFFUSIVITY test  (derived from Elliptic/case1)
#
#  Variable-coefficient elliptic problem solved with the standard direct
#  solver (Ax=b, standard_linsolve!):
#
#        -∇·(a(x,y) ∇u) = f ,        u = g on ∂Ω.
#
#  Verified by the Method of Manufactured Solutions: choose a(x,y) and an
#  exact u_ex(x,y), set  f = -∇·(a ∇u_ex)  and  g = u_ex|_{∂Ω}.  The post-solve
#  L2 error (printed automatically by standard_linsolve!) should be small and
#  converge under refinement.
#
#    Diffusivity :  a(x,y) = A0 + A1·x        (strictly positive on the domain)
#    Exact soln  :  u_ex(x,y) = sin(kx·x) cos(ky·y)
#    Source      :  f = -∇·(a∇u)
#                     = -A1·kx·cos(kx x)cos(ky y)
#                       + (A0 + A1 x)(kx² + ky²) sin(kx x) cos(ky y)
#                   (a_y = 0 since a depends on x only).
#
#  The SAME a(x,y) is given to the operator assembly through
#      inputs[:diffusivity] => nonconstant_diffusivity      (see user_inputs.jl)
#  so the discrete operator -∇·(a∇·) and the manufactured source f are
#  consistent. This case is the verifiable reference for the element-learning
#  non-constant-diffusivity path.
# ============================================================================

const NCD_A0 = 1.0     # diffusivity offset   a = A0 + A1 x
const NCD_A1 = 0.2     # diffusivity slope    (A0 + A1·xmax > 0 required)
const NCD_KX = 1.0     # manufactured-solution wavenumbers
const NCD_KY = 1.0

"Non-constant scalar diffusivity a(x,y) = A0 + A1·x  (> 0 on the domain)."
nonconstant_diffusivity(x, y) = NCD_A0 + NCD_A1 * x

"Exact manufactured solution u_ex(x,y) = sin(kx x) cos(ky y)."
ncd_u(x, y) = sin(NCD_KX * x) * cos(NCD_KY * y)

"Manufactured source f = -∇·(a∇u_ex)."
ncd_f(x, y) = -NCD_A1 * NCD_KX * cos(NCD_KX * x) * cos(NCD_KY * y) +
              (NCD_A0 + NCD_A1 * x) * (NCD_KX^2 + NCD_KY^2) * sin(NCD_KX * x) * cos(NCD_KY * y)


function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=1.0, y=1.0,
                      xmax=1.0, xmin=0.0,
                      ymax=1.0, ymin=0.0)

    PhysConst = PhysicalConst{Float64}()

    # Right-hand side f of  -∇·(a∇u) = f  (manufactured).
    return ncd_f(x, y)
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    return T(-NCD_A1 * NCD_KX * cos(NCD_KX * x) * cos(NCD_KY * y) +
             (NCD_A0 + NCD_A1 * x) * (NCD_KX^2 + NCD_KY^2) * sin(NCD_KX * x) * cos(NCD_KY * y))
end

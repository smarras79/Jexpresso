# ============================================================================
#  Manufactured solution for verifying the element-learning Laplace solver.
#
#  Governing equation (constant diffusivity a = 1, the "standard Laplacian"
#  the element-learning network is trained on):
#
#        -∇·(∇u) = -Δu = f ,        u = g  on ∂Ω.
#
#  Method of Manufactured Solutions (MMS): pick a smooth exact solution
#  u_ex(x,y), insert it into the PDE to obtain the (non-constant) source
#  f = -Δu_ex, and impose the Dirichlet data g = u_ex|_{∂Ω}. The discrete
#  solution must then converge to u_ex. See e.g.
#    • P. J. Roache, "Code Verification by the Method of Manufactured
#      Solutions", J. Fluids Eng. 124 (2002) 4–10.
#    • K. Salari & P. Knupp, "Code Verification by the Method of Manufactured
#      Solutions", SAND2000-1444, Sandia (2000).
#
#  Chosen exact solution (a Laplacian eigenfunction — smooth, with a genuinely
#  NON-CONSTANT right-hand side and NON-ZERO, non-trivial boundary data on
#  every edge):
#
#        u_ex(x,y) = A · sin(kx·x) · cos(ky·y)
#        f(x,y)    = -Δu_ex = A·(kx² + ky²)·sin(kx·x)·cos(ky·y) = (kx²+ky²)·u_ex
#
#  These closed forms are defined once here and reused by user_bc.jl (for g)
#  and initialize.jl (for the exact field qe). user_source.jl is `include`d
#  before those files, so the helpers below are in scope there.
# ============================================================================

# ── Manufactured-solution parameters (smooth, well resolved at nop = 6) ──────
const MMS_A  = 1.0
const MMS_KX = 1.0
const MMS_KY = 1.0

"Exact manufactured solution  u_ex(x,y) = A sin(kx x) cos(ky y)."
manufactured_u(x, y) = MMS_A * sin(MMS_KX * x) * cos(MMS_KY * y)

"Manufactured source  f = -Δu_ex = A (kx²+ky²) sin(kx x) cos(ky y)."
manufactured_f(x, y) = MMS_A * (MMS_KX^2 + MMS_KY^2) * sin(MMS_KX * x) * cos(MMS_KY * y)

# ── Source / test mode selector ─────────────────────────────────────────────
# Edit el_source_mode() below to choose the source term used by this problem
# (self-contained here — no user_inputs.jl plumbing needed):
#   :mms    → manufactured solution above  (non-constant f, exact g, exact qe)
#   :const  → constant source  f = FCONST
#   :custom → arbitrary user-defined source f(x,y)  — see user_defined_source
#             below; use this to solve  A u = f  with ANY right-hand side
#             (not necessarily a manufactured one).
#   :zero   → homogeneous problem  f = 0    (original element-learning case;
#             reproduces the result obtained before the source term was added)
el_source_mode() = :mms
const FCONST = 1.0

# ── Arbitrary, user-defined source f(x,y) for the  :custom  mode ─────────────
# Edit this freely to solve  A u = f  with any source term. The domain bounds
# (xmin/xmax/ymin/ymax) are passed in so the expression can be made
# domain-relative if desired. The template below is the original
# (non-manufactured) source of this problem: with the coefficients at 0 it
# returns f = 0; set alpha/beta/gamma to activate the spatially-varying terms.
function user_defined_source(x, y, xmin, xmax, ymin, ymax)
    L     = abs(xmax - xmin)

    alpha = 0.0
    beta  = 0.0
    gamma = 0.0

    f   = -beta*( cos(x/L)*exp(-x/L)*cos(y)/L + sin(x/L)*exp(-x/L)*cos(y) )
    u_e = gamma*sin(x/L)*exp(-x/L)*cos(y)

    return f - alpha*u_e
end


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

    #
    # S(q(x)) — right-hand-side / source term f of  -Δu = f.
    #
    mode = el_source_mode()
    if mode == :mms
        return manufactured_f(x, y)          # non-constant manufactured source
    elseif mode == :const
        return FCONST                        # constant source
    elseif mode == :custom
        return user_defined_source(x, y, xmin, xmax, ymin, ymax)  # arbitrary f(x,y)
    else
        return 0.0                           # homogeneous (f = 0)
    end

end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    # Mirror of the CPU manufactured source (constant a = 1).
    return T(MMS_A * (MMS_KX^2 + MMS_KY^2) * sin(MMS_KX * x) * cos(MMS_KY * y))
end

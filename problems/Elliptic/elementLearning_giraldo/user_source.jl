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
    # S(q(x)) — right-hand-side / source term f of the governing equation
    #
    #            -∇·(a∇u) = f ,   u = g on ∂Ω.
    #
    # ────────────────────────────────────────────────────────────────────────
    #  CONSTANT source knob
    # ────────────────────────────────────────────────────────────────────────
    #   fconst = 0.0  ->  homogeneous problem (f = 0). Reproduces the original
    #                     element-learning result and is the case used to verify
    #                     that adding the non-zero source term f did NOT break it.
    #   fconst ≠ 0.0  ->  constant, non-zero right-hand side  f = fconst  (e.g.
    #                     set fconst = 1.0 to exercise the new source-term path).
    # ────────────────────────────────────────────────────────────────────────
    fconst = 0.0

    #L     = sqrt((xmax-xmin)^2 + (ymax-ymin)^2)
    L     = abs(xmax-xmin)

    alpha = 0.0
    beta  = 1.0
    gamma = 0.0

    # Optional manufactured (spatially-varying) source. It is fully gated by
    # `beta`, so with the defaults below f reduces exactly to the constant
    # `fconst` (and to 0 when fconst = 0).

    c = 1.0
    f = fconst - beta*(−2(c*π)^2 * sin(c*π*x)*sin(c*π*y))
    
    return f
    
end

function user_source_gpu(q, qe, x, y)
    T = eltype(q)
    L = 2
    alpha = 10
    # Constant source knob (mirror of the CPU path above).
    fconst = T(0.0)
    f   = fconst #T(- (cos(x/L) * exp(-x/L)*cos(y))/L - sin(x/L)*exp(-x/L)*cos(y))
    u_e = T(0.0) #T(sin(x/L)*exp(-x/L)*cos(y))

    return T(f - alpha*u_e)
end

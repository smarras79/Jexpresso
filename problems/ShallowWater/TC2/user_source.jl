function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    #
    # Coriolis source terms for f-plane:
    #   S₁ = 0
    #   S₂ = +f₀ v
    #   S₃ = -f₀ u
    #
    f0 = 1.0e-4  # Coriolis parameter [s⁻¹]

    u = q[2]
    v = q[3]

    S[1] = 0.0
    S[2] =  f0 * v
    S[3] = -f0 * u
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    f0 = 1.0e-4

    u = q[2]
    v = q[3]

    S[1] = 0.0
    S[2] =  f0 * v
    S[3] = -f0 * u
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    f0 = T(1.0e-4)

    u = q[2]
    v = q[3]

    return T(0.0), T(f0*v), T(-f0*u)
end

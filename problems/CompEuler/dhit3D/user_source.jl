function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=5,
                      x=0.0,
                      y=0.0,
                      z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)
    # Taylor-Green vortex: no body force, no gravity.
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0
    S[5] = 0.0
end

function user_source_gpu(q, qe, x, y, z, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
end

#---------------------------------------------------------------------------------
# No source terms: the homogeneous compressible Navier-Stokes system.
# There is no gravity in a shock-tunnel compression ramp problem (the
# domain is 60 mm tall and the flow crosses it in 35 microseconds), and the
# paper's equations (2.1)-(2.4) carry none either.  user_inputs.jl sets
# :lsource => false, so these are never called; they exist because the
# solver requires the methods to be defined.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0)
end

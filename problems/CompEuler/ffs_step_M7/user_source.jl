#---------------------------------------------------------------------------------
# No source terms: this is the homogeneous compressible Euler system. There
# is no gravity in the forward-facing-step problem (the tunnel is 1 m tall;
# the hydrostatic variation over that height is ~10 Pa against a 101325 Pa
# free stream, and neither the reference solutions nor the tutorial include
# it), so every slot is zero and user_inputs.jl sets :lsource => false.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=1.0, xmin=0.0, xmax=3.0)

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
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=1.0, xmin=0.0, xmax=3.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0)
end

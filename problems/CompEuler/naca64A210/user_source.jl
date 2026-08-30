#---------------------------------------------------------------------------------
# No source terms: this is the homogeneous compressible Euler system. Gravity
# over a 2 m tunnel is a ~24 Pa hydrostatic variation against a 101325 Pa free
# stream and plays no part in a Mach-3 aerofoil calculation, so every slot is
# zero and user_inputs.jl sets :lsource => false.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=-1.0, ymax=1.0, xmin=0.0, xmax=3.0)

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
                      neqs=4, x=0.0, y=0.0, ymin=-1.0, ymax=1.0, xmin=0.0, xmax=3.0)

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0)
end

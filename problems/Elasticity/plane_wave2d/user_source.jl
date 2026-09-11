#---------------------------------------------------------------------------------
# 2D elastodynamics has no source: every term of the system is a divergence.
# The deck sets :lsource => false, so this is never called; it exists so that
# this case's `include` overwrites whatever a previously-run case left behind
# at these signatures (see the note in problems/AdvDiff/kopriva/user_source.jl).
#
# A body force (gravity, a distributed load) would go in S[1] and S[2].
#---------------------------------------------------------------------------------
function user_source!(S, q, qe, npoin::Int64, ::CL, ::TOTAL;
                      neqs=5, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    for ieq = 1:neqs
        S[ieq] = 0.0
    end
end

function user_source!(S, q, qe, npoin::Int64, ::CL, ::PERT;
                      neqs=5, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    for ieq = 1:neqs
        S[ieq] = 0.0
    end
end

function user_source!(S, q, qe, npoin::Int64, ::NCL, ::AbstractPert;
                      neqs=5, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    for ieq = 1:neqs
        S[ieq] = 0.0
    end
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
end

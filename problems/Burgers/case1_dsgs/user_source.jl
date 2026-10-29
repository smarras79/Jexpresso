#
# NO SOURCE: inviscid Burgers, one equation.
#
# The PERT / NCL spellings exist only so that this file's `include` really does
# overwrite whatever a previously-run case left behind at those signatures (see
# the note in problems/AdvDiff/kopriva/user_source.jl). This case runs with
# CL() and TOTAL().
#
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::NCL, ::AbstractPert;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
end

function user_source_gpu(q, qe, x, PhysConst, xmax, xmin, lpert)

    T = eltype(q)
    return T(0.0)
end

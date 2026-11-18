# No source term.
function user_source!(S, q, qe, npoin::Int64, ::CL, ::TOTAL; neqs = 1, kwargs...)
    S[1] = 0.0
end
function user_source!(S, q, qe, npoin::Int64, ::CL, ::PERT; neqs = 1, kwargs...)
    S[1] = 0.0
end

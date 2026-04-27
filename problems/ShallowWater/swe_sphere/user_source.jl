function user_source!(S, q, qe, npoin, ::CL, ::TOTAL;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    fill!(S, 0.0)
end

function user_source!(S, q, qe, npoin, ::CL, ::PERT;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    fill!(S, 0.0)
end

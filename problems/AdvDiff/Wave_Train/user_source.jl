function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::AbstractPert; #for this wave train there is no difference between PERT() and TOTAL() in the source
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    S[1] = 0.0
    S[2] = 0.0

end

function user_source_gpu(q,qe,x,y,PhysConst, xmax, xmin, ymax, ymin,lpert)

    T = eltype(q)

    return T(0.0), T(0.0)
end

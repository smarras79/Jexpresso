function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    nothing
end

function user_source_gpu(q,qe,x,PhysConst, xmax, xmin,lpert)

    T = eltype(q)

    return T(0.0), T(0.0)
end

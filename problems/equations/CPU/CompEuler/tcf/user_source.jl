function user_source!(S::SubArray{TFloat},
                      q::SubArray{TFloat}, 
                      qe::SubArray{TFloat},
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax =0.0)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = 0.0
    #
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0
    S[5] = 0.0
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax =0.0, )

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = 0.0
    #
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = 0.0
    S[5] = 0.0
   
end

function user_source_gpu(q,qe,x,y,PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)

    T = eltype(q)
    
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
    
end

function user_source!(S::SubArray{TFloat},
                      q::SubArray{TFloat}, 
                      qe::SubArray{TFloat},
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()

    dπdz = x #the naming is unfortunate, but I want to avoid adding extra params here
    
    #
    # S(q(x)) = -ρg
    #
    ρ = q[1]
    θ = q[4]/ρ

    
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.cp*θ*dπdz -ρ*PhysConst.g
    S[4] = 0.0
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ = q[1] #- qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::NCL,
                      ::AbstractPert; #for NCL() there is no differece between PERT() and TOTAL() in the source
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -g
    #    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0
    
end

function user_source_gpu(q,qe,x,y,PhysConst, xmax, xmin, ymax, ymin,lpert)

    T = eltype(q)
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    return T(0.0), T(0.0), T(-ρ*PhysConst.g), T(0.0)
end

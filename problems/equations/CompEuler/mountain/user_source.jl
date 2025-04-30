function user_source!(S, q, ρref::Float64, npoin::Int64, ::CL; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    dρ  = q[1] - ρref
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -dρ*PhysConst.g
    S[4] = 0.0
    
end

function user_source!(S, q, qref::Float64, npoin::Int64, ::NCL; neqs=1)

    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -g
    #
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -PhysConst.g
    S[4] = 0.0
    
end

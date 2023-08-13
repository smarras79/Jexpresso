
function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    
end

function user_source(q::Array, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
    S = zeros(Float64, neqs)
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    
    return  S
    
end


function olduser_source(T, q::Array, npoin::Int64; neqs=1)

    PhysConst = PhysicalConst{Float64}()
    S = zeros(T, npoin, neqs)
    
    #
    # S(q(x)) = -ρg
    #
    for ip=1:npoin
        ρ  = q[ip,1]
        
        S[ip,1] = 0.0
        S[ip,2] = 0.0
        S[ip,3] = -ρ*PhysConst.g
        S[ip,4] = 0.0
    end
    
    return  S
    
end

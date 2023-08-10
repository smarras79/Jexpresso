
function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin; neqs=1, y=0.0, ymin=0.0, ymax=10000.0, ngl=5, nely=10)
   
    PhysConst = PhysicalConst{Float64}()
    
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0

    #### SPONGE

    #
    # clateral
    nsponge_points = 12

    # distance from the boundary. xs in Restelli's thesis
    dsy = (ymax - ymin)/(nely*(ngl - 1))# equivalent grid spacing
    dbl = ymax - y

    if (abs(dbl) <= nsponge_points * dsy) #&& dbl >= 0.0)
        beta_coe =  1.0 - tanh(dbl/(nsponge_points * dsy))
    else
        beta_coe = 0.0
    end
    cside= 0.75*beta_coe
    
    θref=300.0
    pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³
    
    #@info "β x: " beta_coe, y, ymin, ymax, dsy, dbl
    S[1] -= (cside)*(q[1]-ρref)
    S[2] -= (cside)*q[2]
    S[3] -= (cside)*q[3]
    S[4] -= (cside)*(q[4]-ρref*θref)
    
    
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

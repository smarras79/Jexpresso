
function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin; neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, ngl=5, nely=10,xmin = -100000, xmax =100000)
   
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
    nsponge_points = 8

    # distance from the boundary. xs in Restelli's thesis
    dsy = (ymax - ymin)/(nely*(ngl - 1))# equivalent grid spacing
    dbl = ymax - y
   
    dsx = (xmax - xmin)/(nely*(ngl - 1))# equivalent grid spacing
    dbx = min(xmax - x,x-xmin) 

    if (abs(dbl) <= 5000.0)#nsponge_points * dsy) #&& dbl >= 0.0)
        betay_coe =  1.0 - tanh(dbl/5000.0)#(nsponge_points * dsy))
    else
        betay_coe = 0.0
    end
    ctop= 0.0075*betay_coe
   
    #=if (abs(dbx) <= 20000) #&& dbl >= 0.0)
        betax_coe =  1.0 - tanh(dbx/20000.0)#(nsponge_points * dsy))
    else
        betax_coe = 0.0
    end
    cside= 0.0075*beta_coe=#
 
    θref=250.0
    pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³
    
    #@info nsponge_points * dsy
     
    #@info "β x: " beta_coe, y, ymin, ymax, dsy, dbl
    S[1] -= (ctop)*(q[1]-ρref)
    S[2] -= (ctop)*(q[2]-ρref*20.0)
    S[3] -= (ctop)*q[3]
    S[4] -= (ctop)*(q[4]-ρref*θref)
    
    
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

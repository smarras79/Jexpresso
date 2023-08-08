function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin; neqs=1, x=0.0, xmin=0.0, xmax=0.0, ngl=2, nelx=1)
    
    #
    # clateral
    nsponge_points = 5
    
    # distance from the boundary. xs in Restelli's thesis
    dsx = (xmax - xmin)/(nelx*(ngl - 1))# equivalent grid spacing
    dbl = min(x - xmin, xmax - x)
    
    if (dbl <= nsponge_points * dsx)    
        beta_coe =  1.0 - tanh(dbl/(nsponge_points * dsx))
    else
        beta_coe = 0.0
    end
    cside= beta_coe
    
    #@info "β x: " beta_coe, x, xmin, xmax
    S .-= 20 .*(cside).*q
    
end

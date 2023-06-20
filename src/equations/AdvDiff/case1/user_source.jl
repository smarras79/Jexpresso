function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin; neqs=1, x=0.0, xmin=0.0, xmax=0.0, ngl=2, nelx=1)
    
    #
    # clateral
    nsponge_points = 12
    
    # distance from the boundary. xs in Restelli's thesis
    dsx = (xmax - xmin)/(nelx*(ngl - 1)) # equivalent grid spacing
    dbl = min(x - xmin, xmax - x)
    
    beta_coe =  1.0 - tanh(dbl/(nsponge_points * dsx))
    cside= beta_coe
    
    #@info "β x: " beta_coe x
    S .= -(1.0 - cside).*q
    
end

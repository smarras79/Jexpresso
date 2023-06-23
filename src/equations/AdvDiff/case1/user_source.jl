function user_source!(S::SubArray{Float64}, q::SubArray{Float64}, npoin; neqs=1, x=0.0, xmin=0.0, xmax=0.0, ngl=2, nelx=1)
    
    #
    # clateral
    nsponge_points = 12
    α = 25.0
    
    # distance from the boundary. xs in Restelli's thesis
    dsx = (xmax - xmin)/(nelx*(ngl - 1)) # equivalent grid spacing
    dbl = min(x - xmin, xmax - x)
    
    β = α*tanh(dbl/(nsponge_points * dsx))
    
    #@info "β x: " beta_coe x
    S[:] .= β.*q[:]
    
end

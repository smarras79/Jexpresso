function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe::SubArray{Float64},::TOTAL)
    
    qbdy[1] = 0.0
    
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)

    qbdy[1] = 0.0
    
end


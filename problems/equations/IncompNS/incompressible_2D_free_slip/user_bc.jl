function user_bc_dirichlet!(q::SubArray{TFloat}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe::SubArray{TFloat},::TOTAL)
    
    qbdy[1] = 0.0
    
end

function user_bc_dirichlet!(q::SubArray{TFloat}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{TFloat},::PERT)

    qbdy[1] = 0.0
    
end


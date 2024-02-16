function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe::SubArray{Float64}, ::TOTAL)

    
    qbdy[2] = 0.0
    qbdy[3] = 0.0
    qbdy[4] = 0.0

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

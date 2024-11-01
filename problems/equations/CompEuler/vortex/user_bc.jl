function user_bc_dirichlet!(q::SubArray{TFloat}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe::SubArray{TFloat},::TOTAL)
    nothing   
end

function user_bc_dirichlet!(q::SubArray{Float64}, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe::SubArray{Float64},::PERT)
#    if (tag == "free_slip")
    nothing    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs::Dict)
    nothing
end

function user_bc_dirichlet_gpu(q,qe,x,y,t,nx,ny,qbdy,lpert)
    nothing
end

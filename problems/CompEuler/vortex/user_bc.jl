function user_bc_dirichlet!(q, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe,::TOTAL)
    nothing   
end

function user_bc_dirichlet!(q, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat,qe,::PERT)
#    if (tag == "free_slip")
    nothing    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, t::AbstractFloat, inputs)
    nothing
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    x, y = coords[1], coords[2]
    nothing
end

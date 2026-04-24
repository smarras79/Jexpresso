"""
    For the doubly-periodic 2D viscous Burgers case no Dirichlet/Neumann data
    are needed, so the boundary-condition callbacks are empty. Periodicity is
    handled at the mesh level (connijk identifies matching boundary nodes on
    opposite sides), so the IMEX Laplacian assembly also becomes periodic
    automatically.
"""
function user_bc_dirichlet!(q, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    nothing
end

function user_bc_dirichlet!(q, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, x, y, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    return T(0.0)
end

"""
    For the periodic 1D viscous Burgers case no Dirichlet/Neumann data are
    needed, so the boundary-condition callbacks are empty. Periodicity is
    handled at the mesh level (connijk identifies the first and last global
    node), so the IMEX Laplacian assembly also becomes periodic automatically.
"""
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    nothing
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs::Dict)
    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    return T(q[1])
end

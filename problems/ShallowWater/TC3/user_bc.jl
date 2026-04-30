function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    # TC3 uses sponge-layer absorbing BCs.
    # Free-slip wall: zero normal velocity, free tangential velocity.
    qnl = nx * q[2] + ny * q[3]
    qbdy[2] = q[2] - qnl * nx
    qbdy[3] = q[3] - qnl * ny
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    qnl = nx * q[2] + ny * q[3]
    qbdy[2] = q[2] - qnl * nx
    qbdy[3] = q[3] - qnl * ny
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs::Dict)
    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    qnl = nx * q[2] + ny * q[3]
    u = q[2] - qnl * nx
    v = q[3] - qnl * ny
    return T(qbdy[1]), T(u), T(v)
end

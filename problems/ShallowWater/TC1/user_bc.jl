function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    # TC1: solid-wall coast at x=0 ("left") and open boundary at x=Lx ("right").
    # The exact solution has u≡0 everywhere, so a free-slip wall (zero normal
    # velocity) is satisfied exactly on both x-boundaries. Top/bottom are
    # periodic (enforced at the mesh level), so this is never invoked there.
    qnl = nx * q[2] + ny * q[3]
    qbdy[2] = q[2] - qnl * nx
    qbdy[3] = q[3] - qnl * ny
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    qnl = nx * q[2] + ny * q[3]
    qbdy[2] = q[2] - qnl * nx
    qbdy[3] = q[3] - qnl * ny
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2), 1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
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

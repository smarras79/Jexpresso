#
# Free-slip reflective walls on all four sides of the closed basin.
#
# In conservative variables (H, Hu, Hv) the no-penetration condition is
# expressed by zeroing the normal momentum component:
#
#   (Hu, Hv) <- (Hu, Hv) - ((Hu) nx + (Hv) ny) (nx, ny)
#
# H is left untouched.
#
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    Hu = q[2]
    Hv = q[3]
    qn = nx * Hu + ny * Hv
    qbdy[2] = Hu - qn * nx
    qbdy[3] = Hv - qn * ny
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    user_bc_dirichlet!(q, coords, t, tag, qbdy, nx, ny, qe, TOTAL())
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
    T  = eltype(q)
    Hu = q[2]
    Hv = q[3]
    qn = nx * Hu + ny * Hv
    u  = Hu - qn * nx
    v  = Hv - qn * ny
    return T(qbdy[1]), T(u), T(v)
end

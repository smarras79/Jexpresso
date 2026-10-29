#---------------------------------------------------------------------------------
# Periodic domain (:lperiodic_1d => true): no boundary data.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    nothing
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    return T(q[1]), T(q[2])
end

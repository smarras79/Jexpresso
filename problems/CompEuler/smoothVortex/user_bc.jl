#---------------------------------------------------------------------------------
# Boundary conditions.
#
# The vortex is doubly periodic on [-5,5]²; periodicity is handled at the
# mesh level (the "periodicx"/"periodicy" physical-curve tags of
# vortex_periodic.geo), so these functions only fire if the mesh carries
# non-periodic boundary tags. In that case we fall back to free-slip walls:
# project the normal momentum out and leave ρ and ρE alone.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    qnl     = nx*q[2] + ny*q[3]
    qbdy[2] = q[2] - qnl*nx
    qbdy[3] = q[3] - qnl*ny

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

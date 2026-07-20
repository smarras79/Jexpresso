#---------------------------------------------------------------------------------
# Boundary conditions.
#
# The magnetized KHI of Chan et al. (2022), Section 3.2.1, is doubly periodic
# on [-1,1]²; periodicity is handled at the mesh level (periodic edge tags),
# so these functions only fire if the mesh carries non-periodic boundary
# tags. In that case we fall back to free-slip walls: project the normal
# in-plane momentum out and leave everything else (ρ, ρE, ρw, B, ψ) alone —
# the same treatment used by the companion CompEuler KHI case.
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

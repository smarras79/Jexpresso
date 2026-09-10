#---------------------------------------------------------------------------------
# Boundary conditions (Son, Jang & Magara 2025, Section 2.1):
#
#   "we apply periodic conditions horizontally, symmetric conditions at the
#    bottom (z = 0), and free conditions with an absorbing layer at the top
#    boundary (z = Z_max = 35 H₀)"
#
# - Horizontal periodicity is handled at the mesh level (the "periodicx"
#   physical-curve tag of FE_80x35.geo), so no function here sees the left
#   or right edge.
#
# - "bottom": symmetric boundary. For the in-plane fields this is the
#   reflecting wall of a horizontal flux sheet: the normal velocity and the
#   normal magnetic field vanish (B_z is antisymmetric about z = 0, B_x, ρ,
#   p, V_x are symmetric). Implemented as the free-slip projection of the
#   momentum and of B onto the wall.
#
# - "top": the paper's boundary is open, with the absorbing layer of
#   user_source.jl doing the real work. The strong-form CG discretization
#   has no boundary flux to prescribe, so a genuinely "free" edge would be
#   the do-nothing condition; a free-slip wall behind a 5 H₀ sponge is the
#   robust choice used throughout Jexpresso's stratified-atmosphere cases and
#   is what is applied here (normal momentum removed, everything else free).
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    # Zero normal momentum (free-slip) at both walls
    qnl     = nx*q[2] + ny*q[3]
    qbdy[2] = q[2] - qnl*nx
    qbdy[3] = q[3] - qnl*ny

    if tag == "bottom"
        # Symmetric wall: zero normal magnetic field
        bnl     = nx*q[6] + ny*q[7]
        qbdy[6] = q[6] - bnl*nx
        qbdy[7] = q[7] - bnl*ny
    end

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

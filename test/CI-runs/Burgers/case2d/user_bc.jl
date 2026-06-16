"""
    For the doubly-periodic 2D viscous Burgers case no Dirichlet/Neumann data
    are needed, so the boundary-condition callbacks are empty. Periodicity is
    handled at the mesh level (connijk identifies matching boundary nodes on
    opposite sides), so the IMEX Laplacian assembly also becomes periodic
    automatically.
"""
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny,qe,::TOTAL)
    nothing
end

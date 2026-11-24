#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani_DynSGS — boundary conditions.
#
# There are none: a spherical shell is a CLOSED 2D manifold. It has no
# boundary, so there is no boundary condition to impose anywhere — the
# "periodicity" is the topology itself, enforced by the node numbering built
# when the grid is read (every edge shared by exactly two elements, every seam
# node created once and re-used by both neighbours).
#
# These stubs exist only because problems/drivers.jl loads a user_bc.jl for
# every case.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    # No boundary on a closed shell.
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    # No boundary on a closed shell.
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
    return T(qbdy[1]), T(qbdy[2]), T(qbdy[3])
end

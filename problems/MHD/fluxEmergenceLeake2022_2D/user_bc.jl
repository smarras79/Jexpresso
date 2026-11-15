#---------------------------------------------------------------------------------
# Boundary conditions (as problems/MHD/fluxEmergenceLeake2022; Leake et al. 2013, 2014,
# 2022: "all components of the velocity, and the gradients of magnetic field, gas
# density, and specific energy density are set to zero"). All four sides:
#
#   ρu = ρv = ρw = 0   imposed (no flow through or along the boundary: line-tied)
#   ψ = 0              the one incoming condition of the GLM pair (Bₙ, ψ)
#   ρ, E, B            left free: the strong-form CG boundary nodes evolve with
#                      their own element derivatives, and the natural condition of
#                      the dissipative terms is a zero normal gradient
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray,
                            nx, ny, qe, ::TOTAL)
    qbdy[2] = 0.0
    qbdy[3] = 0.0
    qbdy[5] = 0.0
    qbdy[9] = 0.0
    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray,
                            nx, ny, qe, ::PERT)
    error(" problems/MHD/fluxEmergenceLeake2022_2D: PERT() solution variables are not supported.")
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

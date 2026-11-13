#---------------------------------------------------------------------------------
# Boundary conditions. The paper (Leake, Linton & Antiochos 2010) does not state
# them; its authors' Lare runs of the same problem in 3D (Leake, Linton & Török
# 2013, Sec. 2.3) use "all components of the velocity, and the gradients of
# magnetic field, gas density, and specific energy density are set to zero" on
# every boundary, which is what is applied here, on the four walls of LK10.geo
# ("bottom", "top", "left", "right"):
#
#   ρu = ρv = ρw = 0   imposed (closed, line-tied walls)
#   ψ = 0              the incoming condition of the GLM pair (Bₙ, ψ)
#   ρ, E, B            free: the CG boundary nodes evolve with their own element
#                      derivatives, and the natural condition of the dissipative
#                      terms is a zero normal gradient
#
# Waves are absorbed before they reach the side and top walls by the velocity
# damping layers of user_source.jl.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    qbdy[2] = 0.0
    qbdy[3] = 0.0
    qbdy[5] = 0.0
    qbdy[9] = 0.0
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

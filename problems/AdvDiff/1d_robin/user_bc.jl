#-------------------------------------------------------------------------------------------------
# General Robin condition  a u + b ν ∂u/∂n = g  (n outward) at the "left" and "right" ends.
#
# Solved for the inward diffusive flux, ν ∂u/∂n = (g - a u)/b, which is what user_bc_neumann!
# returns (needs :bdy_fluxes => true). Nothing is set in user_bc_dirichlet!: the wall value is
# part of the solution, not prescribed.
#-------------------------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::TOTAL)
    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::PERT)
    return nothing
end

function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst;
                          t = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, kwargs...)
    u = SOL_VARS_TYPE == PERT() ? q[1] + qe[1] : q[1]
    c = tag == "left" ? R1G_LEFT[] : R1G_RIGHT[]
    F[1] = (c.g - c.a * u) / c.b
    return nothing
end

#-------------------------------------------------------------------------------------------------
# Boundary conditions: Robin walls on "bottom" and "top" (tags from robin_2d.geo); "periodicx" is
# handled by the mesh and needs nothing here.
#
# HOW A ROBIN CONDITION IS IMPOSED IN JEXPRESSO (CG)
#
#   -ν ∂u/∂n = h (u - u∞),   n the outward normal.
#
# The diffusion term is assembled in weak form,
#
#   ∫ ψ ∂u/∂t = -∫ ∇ψ · ν∇u dΩ  +  ∮ ψ ν ∂u/∂n dS ,
#
# and the kernel adds the boundary integral ∮ ψ F dS for whatever F user_bc_neumann! returns
# (requires :bdy_fluxes => true). F is therefore the diffusive flux ν ∂u/∂n along the OUTWARD
# normal, i.e. what ENTERS the domain through the wall. A Robin condition is a natural condition
# whose flux depends on the solution itself:
#
#   F = ν ∂u/∂n = h (u∞ - u)        (positive when the ambient is warmer: heat enters)
#
# so it is imposed ONLY through user_bc_neumann!, with the current boundary value u = q[1].
# Nothing is prescribed in user_bc_dirichlet!: the wall value is free and settles where the
# flux balance puts it. (Writing a value there as well would turn the wall into a Dirichlet
# boundary and the flux would be ignored, since a constrained node's tendency is zeroed.)
#
# Adapting this to another case:
#   - several equations: set F[ieq] for each one that has a Robin wall, using q[ieq]
#     (q is the state in the solution variables; in PERT form add qe[ieq] for the full value);
#   - a space-dependent u∞ or h: use coords (x = coords[1], y = coords[2]);
#   - a pure Neumann wall: F = the prescribed inward flux, independent of q;
#   - a Dirichlet wall: set qbdy[ieq] in user_bc_dirichlet! for that tag instead.
#-------------------------------------------------------------------------------------------------

# Dirichlet hook: nothing is prescribed on any tag, so no boundary value is constrained.
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::PERT)
    return nothing
end

# Neumann hook: called at every boundary node of every non-periodic edge, at every RHS evaluation.
#   F          inward flux at this node, one entry per equation (zeroed by the kernel beforehand)
#   q          state at the boundary node; q1 the state at the next node inward (unused here)
#   qe         reference state at the node (zero in this case)
#   tag        the physical name of the boundary edge in the mesh
#   coords     (x, y) of the node
#   t          time (keyword)
#   nx, ny, nz outward unit normal at the node (keywords; nz = 0 in 2D)
# τ_f, wθ, PhysConst and the other keywords serve the atmospheric cases; ignore them here.
# Example of a time-dependent ambient value: u∞(t) = 1 + 0.1 sin(2π t), i.e.
#     F[1] = ROBIN_H[] * (1 + 0.1*sinpi(2t) - u)
function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst;
                          t = 0.0, nx = 0.0, ny = 0.0, nz = 0.0, kwargs...)
    u = SOL_VARS_TYPE == PERT() ? q[1] + qe[1] : q[1]
    if tag == "bottom"
        F[1] = ROBIN_H[] * (ROBIN_UINF_B[] - u)
    elseif tag == "top"
        F[1] = ROBIN_H[] * (ROBIN_UINF_T[] - u)
    end
    return nothing
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    return T(0.0)
end

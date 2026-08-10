#---------------------------------------------------------------------------------
# SWsphere — fluxes.
#
# PLACEHOLDER. The run stops after the grid is built (:lgrid_only => true), so
# nothing here is ever called. The signatures are kept so that the case loads
# through problems/drivers.jl exactly like every other case.
#
# What will go here: the shallow water fluxes written on the TANGENT PLANE of
# the shell. On a manifold they are not simply the Cartesian (F, G) pair of the
# flat cases — the divergence has to be taken with the surface metric, and the
# momentum equation carries the curvature term that keeps the velocity tangent
# to the sphere. Both need the metric terms of the manifold, which do not exist
# yet either.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)

    error(" # ERROR problems/ShallowWater/SWsphere/user_flux.jl: not implemented yet (grid-only case).")
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)

    error(" # ERROR problems/ShallowWater/SWsphere/user_flux.jl: not implemented yet (grid-only case).")
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
end

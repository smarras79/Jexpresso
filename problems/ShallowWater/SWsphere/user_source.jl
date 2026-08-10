#---------------------------------------------------------------------------------
# SWsphere — source terms.
#
# PLACEHOLDER, see user_flux.jl. What will go here: the full Coriolis term of
# the rotating sphere, f = 2Ω sin φ, with φ the latitude that the grid already
# stores per node (mesh.lat), plus the curvature/metric term that keeps the
# momentum tangent to the shell.
#---------------------------------------------------------------------------------
function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    error(" # ERROR problems/ShallowWater/SWsphere/user_source.jl: not implemented yet (grid-only case).")
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    error(" # ERROR problems/ShallowWater/SWsphere/user_source.jl: not implemented yet (grid-only case).")
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T = eltype(q)
    return T(0.0), T(0.0), T(0.0)
end

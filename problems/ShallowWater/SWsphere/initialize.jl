#---------------------------------------------------------------------------------
# SWsphere — initial condition.
#
# NOT WRITTEN YET, and deliberately so: with :lgrid_only => true in
# user_inputs.jl the driver returns as soon as the spherical shell grid has
# been built, verified and written to VTK, and never reaches this file.
#
# When the shallow-water kernels on the manifold are in place, this is where
# the classic Williamson et al. (1992) test cases go. The grid already carries
# what they need: mesh.lon / mesh.lat (radians) for every node, mesh.radius,
# and the Cartesian mesh.x/y/z of the shell. A Williamson-2 (steady zonal
# geostrophic flow) IC would read, in spherical components,
#
#   h(λ,φ) = h₀ - (R Ω u₀ + u₀²/2) sin²φ / g
#   u(λ,φ) = u₀ cos φ ,  v = 0
#
# and then be pushed to the tangent-plane Cartesian velocity the SEM kernels
# will use. Nothing here can be exercised until those kernels exist, so this
# function fails loudly rather than returning a plausible-looking zero state.
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    error(" # ERROR problems/ShallowWater/SWsphere/initialize.jl: " *
          "the shallow water equations on the sphere have no initial condition yet. " *
          "Keep :lgrid_only => true in user_inputs.jl — the run stops after the grid is built.")

end

#---------------------------------------------------------------------------------
# SWsphere_visc — the SAME equations as ShallowWater/SWsphere.
#
# This case differs from SWsphere in its user_inputs.jl ALONE: the artificial
# diffusion of Eq. (8b) replaces the modal filter as the stabilisation. Nothing
# about the equations, the initial condition or the boundary treatment changes,
# so this file is the SWsphere one — included rather than copied, so the two
# cases cannot drift apart.
#
# (src/run.jl includes exactly six user_*.jl per case directory and does not
# look for them anywhere else, which is why the file has to exist at all.)
#---------------------------------------------------------------------------------
include(joinpath(@__DIR__, "..", "SWsphere", "initialize.jl"))

#---------------------------------------------------------------------------------
# fluxEmergenceSon2025DSGS is the flux-emergence problem of Son, Jang &
# Magara (2025) with the SAME physics, mesh, initial and boundary conditions
# as problems/MHD/fluxEmergenceSon2025, differing only in how positivity is
# kept (see README.md). This file therefore just loads the sibling's.
#---------------------------------------------------------------------------------
include(joinpath(@__DIR__, "..", "fluxEmergenceSon2025", "initialize.jl"))

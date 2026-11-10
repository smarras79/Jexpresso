#---------------------------------------------------------------------------------
# MHD blast wave: the GLM-MHD flux of problems/MHD/orszagTangBormanis2024, with γ = 1.4.
#---------------------------------------------------------------------------------
if !@isdefined(γ_mhd)
    const γ_mhd = 1.4
end
if abs(γ_mhd - 1.4) > 1e-12
    error(" problems/MHD/blastBalsaraSpicer1999: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 1.4). Restart Julia before running this case.")
end
include(joinpath(@__DIR__, "..", "orszagTangBormanis2024", "user_flux.jl"))

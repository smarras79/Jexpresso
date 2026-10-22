#---------------------------------------------------------------------------------
# Shock / "MHD" clouds (Dai & Woodward 1998, Fig. 20): the GLM-MHD flux of
# problems/MHD/orszagTangBormanis2024 with its γ = 5/3, plus the two states of the shock and the letters.
#---------------------------------------------------------------------------------
if !@isdefined(γ_mhd)
    const γ_mhd = 5.0/3.0
end
if abs(γ_mhd - 5.0/3.0) > 1e-12
    error(" problems/MHD/shockCloudDaiWoodward1998: γ_mhd = $(γ_mhd) is already defined by another MHD case loaded in this Julia session (this case needs γ = 5/3). Restart Julia before running this case.")
end
include(joinpath(@__DIR__, "..", "orszagTangBormanis2024", "user_flux.jl"))
include(joinpath(@__DIR__, "mhd_letters.jl"))

# DW's states (ρ, p, u, v, w, Bx, By, Bz), Gaussian B as printed: magnetic pressure B²/8π.
# Heaviside-Lorentz B (½|B|²) is B/√(4π); ρ, p and v are the same in both.
if !@isdefined(SC_PRE)
    const SC_PRE  = (1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0)
    const SC_POST = (3.91109, 2.53724e2, 1.38388e1, -5.00186e-2, -5.00186e-2, 2.0, 7.84323, 7.84323)
    const SC_RHO_CLOUD = 50.0          # "50 times denser than that in the preshock state", in p and B equilibrium
end

# Conserved state (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ) of a DW primitive tuple, B converted to Heaviside-Lorentz.
@inline function sc_conserved(ρ, p, u, v, w, Bx, By, Bz)
    s = 1.0/sqrt(4.0*π)
    Bx, By, Bz = s*Bx, s*By, s*Bz
    ρE = p/(γ_mhd - 1.0) + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz)
    return (ρ, ρ*u, ρ*v, ρE, ρ*w, Bx, By, Bz, 0.0)
end
@inline sc_conserved(t::NTuple{8,Float64}) = sc_conserved(t...)

# |v| + c_f bound of a DW primitive tuple (c_f ≤ √((γp + |B|²)/ρ), HL |B|²)
@inline sc_wave_speed(ρ, p, u, v, w, Bx, By, Bz) =
    sqrt(u*u + v*v + w*w) + sqrt((γ_mhd*p + (Bx*Bx + By*By + Bz*Bz)/(4.0*π))/ρ)

#
# Non-linear shallow water equations with reduced gravity term-- conservative form.
# This model assumes a layer of water sitting on a deep bottom layer whic is assumed motionless. 
# The reduced gravity term is used to account for the density difference between the two layers.
# See Choi, Iskandarani, Levin & Haidvogel (2004), “A Spectral Finite-Volume Method 
#           for the Shallow Water Equations,” Monthly Weather Review, 132, 1777–1791.
#
# Adapted form the SoliWave test case in problems/ShallowWater/SoliWaveIsland/user_inputs.jl
# The flux term does not change except for teh reduced gravity parameter
# Wetting and drying is important to handle outcropping.
#
# Convention used by every file of this problem (it matches Marras et al.
# 2018, Sec. 2, where the total water surface is H_total = Hs + Hb):
#
#   z = 0           : flat bottom of the basin;
#   Hb(x,y) >= 0    : bathymetry (the conical island) measured UP from z=0;
#   q[1] = H        : water depth measured from the TOP of the bathymetry
#                     (Hs in the paper's notation);
#   z  = Hb + H     : position of the free surface;
#   qe[1] = He      : lake-at-rest depth max(h0 - Hb, H_dry) set in
#                     initialize.jl, so He + Hb = h0 wherever the cone is
#                     submerged and He = H_dry (thin dry film) where it
#                     pierces the surface.
#
# Equations:
#
#   ∂H/∂t  + ∂(Hu)/∂x  + ∂(Hv)/∂y                                 = 0
#   ∂Hu/∂t + ∂(Hu² + g H²/2)/∂x + ∂(Huv)/∂y                       = -g H ∂Hb/∂x
#   ∂Hv/∂t + ∂(Huv)/∂x          + ∂(Hv² + g H²/2)/∂y              = -g H ∂Hb/∂y
#
# Well-balanced (perturbation) split. The raw pair {flux gH²/2, source
# -gH∇Hb} balances the lake at rest only in exact arithmetic; the discrete
# residual at the wet/dry ring is amplified by the 1/H division and blows
# the run up at the first time step. Using the identity
#
#   ∇(g H²/2) + g H ∇Hb = ∇(g (H²-He²)/2) + g (H-He) ∇Hb + g He ∇(He+Hb)
#
# the last term vanishes identically on wet nodes (He + Hb = h0) and is
# O(g H_dry ∇Hb) on dry nodes, where it only pressurizes the artificial
# thin film, so we drop it and advance
#
#   pressure flux  :  g (H² - He²)/2        (here)
#   bathymetry     : -g (H - He) ∇Hb        (user_source.jl)
#
# which makes the rest state an exact equilibrium of the discrete operators.
#
const _G_SWE      = 0.015   # reduced gravity g'= g*(Δρ/ρ) [m/s^2] - value matches Choi et al. (2024)
const _H_WET_SWE  = 1.0e-3  # thin-film depth = wet/dry threshold (Marras et al. 2018, Sec. 5.4-5.5)

@inline function _swe_uvel(Hc, Hu)
    # Desingularized velocity (Kurganov & Petrova 2007): equals Hu/Hc for
    # Hc >= threshold and decays smoothly to zero on thinner layers instead
    # of switching off discontinuously -- the moving wet/dry front crosses
    # the threshold every step and a hard switch lets u jump to |Hu|/ε.
    # Combined with the dry-node momentum relaxation in user_source.jl this
    # realises the paper's "velocity forced to zero on the dry nodes".
    ε  = _H_WET_SWE
    H4 = max(Hc, ε)
    return sqrt(2.0) * Hc * Hu / sqrt(Hc^4 + H4^4)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)

    # The CG path has no positivity limiter (cf. Xing et al. 2010 used in
    # the paper), so H can undershoot below zero at the receding front:
    # clamp the depth used by the fluxes -- a negative layer must neither
    # carry mass nor exert pressure.
    Hc = max(q[1], 0.0)
    He = qe[1]

    u = _swe_uvel(Hc, q[2])
    v = _swe_uvel(Hc, q[3])

    # perturbation pressure: zero at lake-at-rest by construction
    p = 0.5 * _G_SWE * (Hc * Hc - He * He)

    # x-flux
    F[1] = Hc * u
    F[2] = Hc * u * u + p
    F[3] = Hc * u * v

    # y-flux
    G[1] = Hc * v
    G[2] = Hc * v * u
    G[3] = Hc * v * v + p
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    user_flux!(F, G, SD, q, qe, mesh, CL(), TOTAL(); neqs=neqs, ip=ip)
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T  = eltype(q)
    Hc = max(q[1], T(0.0))
    He = qe[1]
    ε  = T(_H_WET_SWE)
    H4 = max(Hc, ε)
    de = sqrt(T(2.0)) * Hc / sqrt(Hc^4 + H4^4)
    u  = de * q[2]
    v  = de * q[3]
    p  = T(0.5) * T(_G_SWE) * (Hc * Hc - He * He)

    return T(Hc*u),      T(Hc*u*u + p), T(Hc*u*v),
           T(Hc*v),      T(Hc*v*u),     T(Hc*v*v + p)
end

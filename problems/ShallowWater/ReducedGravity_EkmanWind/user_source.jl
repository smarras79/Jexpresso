#
# Sources of the reduced-gravity double gyre (Coriolis, wind stress, linear
# bottom drag, thin-film momentum relaxation), as in ReducedGravity/user_source.jl,
# with ONE change: the thickness scale over which the wind stress is handed to
# the active layer.
#
# Why. The wind stress τ is a force per unit AREA. Applied to a layer of
# thickness H it accelerates the layer at τ/(ρH), which is unbounded as the
# layer outcrops (H -> 0): the reduced-gravity momentum equation has a source
# that does NOT vanish with the layer, unlike every term of the SoliWaveIsland
# case (pressure g H ∇H, bathymetry g (H - He) ∇Hb), where the velocity stays
# bounded by the shallow-water characteristics however thin the layer gets.
# In the ReducedGravity deck the stress is ramped in over 1 mm < H < 0.1 m, so
# a 10 cm layer already carries the full Ekman transport τ0/(ρ f0) = 0.36 m²/s,
# i.e. u = 3.6 m/s (Fr ≈ 90); the upwelling that produced the thin layer keeps
# thinning it at that transport, and u = Hu/H crosses Δx_min/Δt = 23 m/s
# (CFL = 1 for the 25x25, N = 4 mesh and Δt = 600 s) near H ≈ 1.5 cm -- an
# order of magnitude ABOVE the 1 mm clip, which is why clipping does not help.
#
# Physically the stress is carried by an Ekman layer of depth H_E ~ 10-50 m;
# where the active layer is thinner than that only a fraction ~H/H_E of the
# stress goes into it (the rest into the deep layer, assumed motionless), so
# the acceleration is bounded by τ/(ρ H_E) and the layer velocity by
# τ0/(ρ f0 H_E) = 1.8 cm/s for H_E = 20 m (Fr ≈ 0.03). The smoothstep
# ramp of the deck is kept; only its upper end moves from 0.1 m to H_E.
# The thin-film relaxation keeps its own 0.1 m scale (_H_RELAX_SWE), so that
# this case differs from ReducedGravity in the wind forcing alone.
#

const _SIGMA_DRY_SWE = 0.25 / 1200   # [1/s] max relaxation rate of the thin-film momentum (σΔt = 0.125 at Δt = 600 s)
const _F0_SWE        = 8.37e-5 # Coriolis parameter [1/s]  (Choi et al. 2004, Sec. 3.1)
const _BETA_SWE      = 1.8e-11 # Coriolis beta parameter [1/(m s)] (Choi et al. 2004, Sec. 3.1)
const _TAU0_SWE      = 0.03    # maximum wind stress [N/m²]
const _RHO_SWE       = 1000.0  # reference water density [kg/m³]
const _KAPPA_SWE     = 1.0e-8  # linear bottom drag [1/s]

# Ekman-depth scale [m] over which the wind stress is handed to the layer
# (full stress for H >= _H_EKMAN_SWE). ReducedGravity uses 0.1 m here.
const _H_EKMAN_SWE   = 20.0
# Thickness [m] below which the thin-film momentum relaxation ramps in
# (unchanged from ReducedGravity, where it is also the wind scale).
const _H_RELAX_SWE   = 0.1

@inline function _swe_smoothstep(H, Hlo, Hhi)
    T = typeof(H)
    s = clamp((H - T(Hlo)) / (T(Hhi) - T(Hlo)), zero(T), one(T))
    return s*s*(T(3) - T(2)*s)
end

# fraction of the wind stress handed to the layer: 0 on the thin film, 1 for H >= H_E
@inline _swe_wind_fraction(H) = _swe_smoothstep(H, _H_WET_SWE, _H_EKMAN_SWE)

# fraction of the maximum relaxation rate applied: 1 on the thin film, 0 for H >= 0.1 m
@inline _swe_relax_fraction(H) = one(H) - _swe_smoothstep(H, _H_WET_SWE, _H_RELAX_SWE)

#No bathymetry gradient in the reduced gravity case, so we can just return 0.0, 0.0
#Keeping it here for later when we have multiple layers
@inline function _swe_bathy_grad(x, y)
    return 0.0, 0.0
end

#Calculate coriolis parameter with beta-plane approximation f(y) = f0 + β (y - yc) [1/s]
@inline function _swe_coriolis(y, ymin, ymax)
    y_center = 0.5 * (ymin + ymax)
    return _F0_SWE + _BETA_SWE * (y - y_center)
end

#Calculate wind stress with a zonal cosine profile in y-direction
@inline function _swe_wind_stress(y, ymin, ymax)
    Ly = ymax - ymin
    Ly > 0 || throw(ArgumentError(
        "Invalid y bounds: ymin=$ymin, ymax=$ymax"
    ))
    yn = (y - ymin) / Ly

    τx = - _TAU0_SWE * cospi(2.0 * yn)
    τy = 0.0

    return τx, τy
end

function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    H  = q[1]
    Hu = q[2]
    Hv = q[3]

    dH = max(H, 0.0) - qe[1]   # clamped depth perturbation w.r.t. lake at rest
    dHbdx, dHbdy = _swe_bathy_grad(x, y)

    τx, τy = _swe_wind_stress(y, ymin, ymax)

    f = _swe_coriolis(y, ymin, ymax)

    #Bathymetry (pressure) source term
    S[1] = 0.0
    S[2] = -_G_SWE * dH * dHbdx
    S[3] = -_G_SWE * dH * dHbdy

    #Coriolis term
    S[2] += f * Hv
    S[3] -= f * Hu

    #Wind shear stress term, handed to the layer over the Ekman-depth scale
    wet = _swe_wind_fraction(H)
    S[2] += wet * τx / _RHO_SWE
    S[3] += wet * τy / _RHO_SWE

    #bottom drag term
    S[2] -= _KAPPA_SWE * Hu
    S[3] -= _KAPPA_SWE * Hv

    # thin-film momentum relaxation (same scale as ReducedGravity: 1 mm - 0.1 m)
    σ = _swe_relax_fraction(H) * _SIGMA_DRY_SWE
    S[2] -= σ * Hu
    S[3] -= σ * Hv
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    user_source!(S, q, qe, npoin, CL(), TOTAL();
                 neqs=neqs, x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)
    T  = eltype(q)
    H  = q[1]
    Hu = q[2]
    Hv = q[3]

    #Wind stress
    Ly = ymax - ymin
    yn = (y - ymin) / Ly
    τx = -T(_TAU0_SWE) * cospi(T(2.0) * yn)
    τy = T(0.0)

    wet = _swe_wind_fraction(H)

    #Coriolis
    yc = T(0.5) * (ymin + ymax)
    f = T(_F0_SWE) + T(_BETA_SWE) * (y - yc)

    #Bottom drag
    ρ = T(_RHO_SWE)
    κ = T(_KAPPA_SWE)

    σ = _swe_relax_fraction(H) * T(_SIGMA_DRY_SWE)
    return T(0.0),
           f*Hv + wet*τx/ρ - (σ + κ) * Hu ,
          -f*Hu + wet*τy/ρ - (σ + κ) * Hv #pressure term removed - flat bathymetry assumed
end

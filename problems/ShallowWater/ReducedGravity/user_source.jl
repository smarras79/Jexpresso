#
# Bathymetry source term for the non-linear shallow water equations,
# in the well-balanced perturbation split documented in user_flux.jl:
#
#   S = (0, -g (H - He) ∂Hb/∂x, -g (H - He) ∂Hb/∂y)
#
# where H = q[1] is the water depth above the bathymetry and He = qe[1]
# is the lake-at-rest depth. At rest H = He, so S vanishes node-by-node
# exactly as the perturbation pressure flux g(H² - He²)/2 does: the
# discrete equilibrium leaves nothing for the wet/dry 1/H division to
# amplify.
#
# Here Hb is flat, but we keep it for the future multi-layer model
#

# Dry-node momentum relaxation rate [1/s]. On layers thinner than the
# wet/dry threshold the momentum is relaxed to zero (the paper forces the
# velocity to zero on dry nodes); σΔt ≈ 0.25 keeps the term well inside
# the explicit SSPRK54 stability region.

const _SIGMA_DRY_SWE = 25.0 / 1200   # NOTE - need to choose it appropriately with dt
const _F0_SWE        = 8.37e-5 # Coriolis parameter [1/s]  (Choi et al. 2004, Sec. 3.1)
const _BETA_SWE      = 1.8e-11 # Coriolis beta parameter [1/(m s)] (Choi et al. 2004, Sec. 3.1)
const _TAU0_SWE      = 0.03    # maximum wind stress [N/m²]
const _RHO_SWE       = 1000.0  # reference water density [kg/m³]
const _KAPPA_SWE     = 1.0e-8  # linear bottom drag [1/s]

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

    τx = _TAU0_SWE * cospi(2.0 * yn)
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

    τx, τy = _swe_wind_stress(y, ymin, ymax) #NOTE - make sure ymin and ymax are passed correctly to this function

    f = _swe_coriolis(y, ymin, ymax)

    #Bathymetry (pressure) source term
    S[1] = 0.0
    S[2] = -_G_SWE * dH * dHbdx 
    S[3] = -_G_SWE * dH * dHbdy 

    #Coriolis term
    S[2] += f * Hv
    S[3] -= f * Hu

    #Wind shear stress term
    wet = H >= _H_WET_SWE ? one(H) : zero(H) #apply wind stress only to wet nodes
    S[2] += wet * τx / _RHO_SWE
    S[3] += wet * τy / _RHO_SWE

    #bottom drag term
    S[2] -= _KAPPA_SWE * Hu
    S[3] -= _KAPPA_SWE * Hv
    #S[2] -= _C_D_SWE * Hu * sqrt(Hu^2 + Hv^2) / H
    #S[3] -= _C_D_SWE * Hv * sqrt(Hu^2 + Hv^2) / H


    # dry-node momentum relaxation (see _SIGMA_DRY_SWE)
    if H < _H_WET_SWE
        S[2] -= _SIGMA_DRY_SWE * q[2]
        S[3] -= _SIGMA_DRY_SWE * q[3]
    end
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
    τx = T(_TAU0_SWE) * cospi(T(2.0) * yn)
    τy = T(0.0)

    # Do not apply wind stress to the artificial dry film.
    wet = ifelse(H >= T(_H_WET_SWE), one(T), zero(T))

    #Coriolis
    yc = T(0.5) * (ymin + ymax)
    f = T(_F0_SWE) + T(_BETA_SWE) * (y - yc)

    #Bottom drag
    ρ = T(_RHO_SWE)
    κ = T(_KAPPA_SWE)

    g = T(_G_SWE)
    σ = H < T(_H_WET_SWE) ? T(_SIGMA_DRY_SWE) : T(0.0)
    #return T(0.0), T(-g * dH * dHbdx + f*q[3] + τx/ρ - κ*q[2] - σ * q[2]  ), T(-g * dH * dHbdy - f*q[2] + τy/ρ - κ*q[3] - σ * q[3])
    return T(0.0), 
           f*Hv + wet*τx/ρ - (σ + κ) * Hu , 
          -f*Hu + wet*τy/ρ - (σ + κ) * Hv #pressure term removed - flat bathymetry assumed
end

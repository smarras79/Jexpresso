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

const _SIGMA_DRY_SWE = 25.0

@inline function _swe_bathy_grad(x, y)
    dx = x - _XC_CONE_SWE
    dy = y - _YC_CONE_SWE
    r  = sqrt(dx*dx + dy*dy)
    if r < _RC_CONE_SWE && r > 1.0e-12
        slope = -_HC_CONE_SWE / (_RC_CONE_SWE * r)
        return slope * dx, slope * dy
    else
        return 0.0, 0.0
    end
end

function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=3, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    H  = q[1]
    dH = max(H, 0.0) - qe[1]   # clamped depth perturbation w.r.t. lake at rest
    dHbdx, dHbdy = _swe_bathy_grad(x, y)

    S[1] = 0.0
    S[2] = -_G_SWE * dH * dHbdx
    S[3] = -_G_SWE * dH * dHbdy

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
    dH = max(H, T(0.0)) - qe[1]

    dx = x - T(_XC_CONE_SWE)
    dy = y - T(_YC_CONE_SWE)
    r  = sqrt(dx*dx + dy*dy)

    dHbdx = T(0.0)
    dHbdy = T(0.0)
    if r < T(_RC_CONE_SWE) && r > T(1.0e-12)
        slope = -T(_HC_CONE_SWE) / (T(_RC_CONE_SWE) * r)
        dHbdx = slope * dx
        dHbdy = slope * dy
    end

    g = T(_G_SWE)
    σ = H < T(_H_WET_SWE) ? T(_SIGMA_DRY_SWE) : T(0.0)
    return T(0.0), T(-g * dH * dHbdx - σ * q[2]), T(-g * dH * dHbdy - σ * q[3])
end

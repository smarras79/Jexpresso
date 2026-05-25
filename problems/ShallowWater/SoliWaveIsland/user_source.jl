#
# Bathymetry source term for the non-linear shallow water equations.
#
#   S = (0, -g H ∂Hb/∂x, -g H ∂Hb/∂y)
#
# Hb(x,y) is the conical island of Marras et al. 2018, Sec. 5.5:
#
#   Hb(x,y) = 0.93 * (1 - r/rc)     if r ≤ rc
#           = 0                       otherwise
#
# with rc = 3.6 m and r = sqrt((x - xc_cone)² + (y - yc_cone)²).
# We place the cone in the middle of the basin at (xc_cone, yc_cone) = (12.5, 0).
#
const _XC_CONE_SWE =  12.5
const _YC_CONE_SWE =   0.0
const _RC_CONE_SWE =   3.6
const _HC_CONE_SWE =   0.93

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

    H = q[1]
    dHbdx, dHbdy = _swe_bathy_grad(x, y)

    S[1] = 0.0
    S[2] = -_G_SWE * H * dHbdx
    S[3] = -_G_SWE * H * dHbdy
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
    T = eltype(q)
    H  = q[1]

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
    return T(0.0), T(-g * H * dHbdx), T(-g * H * dHbdy)
end

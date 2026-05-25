#
# Non-linear shallow water equations -- conservative form.
#
#   ∂H/∂t  + ∂(Hu)/∂x  + ∂(Hv)/∂y                                = 0
#   ∂Hu/∂t + ∂(Hu² + g H²/2)/∂x + ∂(Huv)/∂y                       = -g H ∂Hb/∂x
#   ∂Hv/∂t + ∂(Huv)/∂x          + ∂(Hv² + g H²/2)/∂y              = -g H ∂Hb/∂y
#
# State: q = [H, Hu, Hv]. Bathymetry Hb(x,y) is time-independent and shows
# up only through the momentum source term (computed in user_source!).
#
# Marras et al. 2018, Sec. 2 and Sec. 5.5.
#
const _G_SWE      = 9.81
const _H_MIN_SWE  = 1.0e-4   # wet/dry safety floor used only when dividing by H

@inline function _swe_uvel(H, Hu)
    h = max(H, _H_MIN_SWE)
    return Hu / h
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)

    H  = q[1]
    Hu = q[2]
    Hv = q[3]

    u = _swe_uvel(H, Hu)
    v = _swe_uvel(H, Hv)

    gH2_half = 0.5 * _G_SWE * H * H

    # x-flux
    F[1] = Hu
    F[2] = Hu * u + gH2_half
    F[3] = Hu * v

    # y-flux
    G[1] = Hv
    G[2] = Hv * u
    G[3] = Hv * v + gH2_half
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    user_flux!(F, G, SD, q, qe, mesh, CL(), TOTAL(); neqs=neqs, ip=ip)
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T  = eltype(q)
    H  = q[1]
    Hu = q[2]
    Hv = q[3]
    h  = max(H, T(_H_MIN_SWE))
    u  = Hu / h
    v  = Hv / h
    g  = T(_G_SWE)
    p  = T(0.5) * g * H * H

    return T(Hu),          T(Hu*u + p),  T(Hu*v),
           T(Hv),          T(Hv*u),      T(Hv*v + p)
end

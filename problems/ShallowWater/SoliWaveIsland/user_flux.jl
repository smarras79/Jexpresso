#
# Non-linear shallow water equations -- conservative form.
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
const _G_SWE      = 9.81
const _H_WET_SWE  = 1.0e-3   # thin-film depth = wet/dry threshold (paper Sec. 5.4-5.5)

@inline function _swe_uvel(H, Hu)
    # Dry nodes carry no velocity (Marras et al. 2018: "the velocity in such
    # elements is forced to zero only on the dry nodes"); never divide
    # momentum noise by the thin-film depth.
    return H > _H_WET_SWE ? Hu / H : zero(H)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=3, ip=1)

    H  = q[1]
    He = qe[1]

    u = _swe_uvel(H, q[2])
    v = _swe_uvel(H, q[3])

    # perturbation pressure: zero at lake-at-rest by construction
    p = 0.5 * _G_SWE * (H * H - He * He)

    # x-flux
    F[1] = H * u
    F[2] = H * u * u + p
    F[3] = H * u * v

    # y-flux
    G[1] = H * v
    G[2] = H * v * u
    G[3] = H * v * v + p
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=3, ip=1)
    user_flux!(F, G, SD, q, qe, mesh, CL(), TOTAL(); neqs=neqs, ip=ip)
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T  = eltype(q)
    H  = q[1]
    He = qe[1]
    u  = H > T(_H_WET_SWE) ? q[2] / H : T(0.0)
    v  = H > T(_H_WET_SWE) ? q[3] / H : T(0.0)
    p  = T(0.5) * T(_G_SWE) * (H * H - He * He)

    return T(H*u),       T(H*u*u + p),  T(H*u*v),
           T(H*v),       T(H*v*u),      T(H*v*v + p)
end

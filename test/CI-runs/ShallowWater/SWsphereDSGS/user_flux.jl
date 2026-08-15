#---------------------------------------------------------------------------------
# SWsphereDSGS — fluxes of the shallow water equations on a spherical shell.
#
# Marras, Kopera & Giraldo (2015), QJRMS 141: 1727-1739,
#   "Simulation of shallow-water jets with a unified element-based
#    continuous/discontinuous Galerkin model with grid flexibility on the
#    sphere", Eq. (8):
#
#   ∂φ/∂t     + ∇·(φu)     = 0
#   ∂(φu)/∂t  + ∇·(φu⊗u)   = -φ∇φ - f(x × φu) + μx + δν∇²(φu)
#
# with φ = gh the geopotential, u = (u,v,w) the FULL CARTESIAN velocity,
# ∇ = (∂x, ∂y, ∂z)ᵀ, x the position vector, f = 2ωz/r² and μ the Lagrange
# multiplier that keeps the fluid on the shell (see user_source.jl).
#
# STATE VECTOR  q = [φ, φu, φv, φw]ᵀ   — 4 equations, conservative variables.
#
# Why Cartesian and not the tangent-plane pair (u,v): the whole point of the
# Lagrange multiplier is to remove the velocity component NORMAL to the shell.
# In a two-component tangent basis that component does not exist, so the
# constraint would be vacuous and there would be nothing for μ to do. The price
# is a fourth equation and a redundant degree of freedom that must be projected
# out — which is exactly what the paper does.
#
# CONSERVATIVE vs NON-CONSERVATIVE PRESSURE TERM
#
# The paper carries the pressure as the source -φ∇φ and leaves the flux as
# φu⊗u. Jexpresso's user_source! is a POINTWISE callback — it receives q at one
# node and no derivatives — so -φ∇φ cannot be evaluated there. The identity
#
#   ∇·( (φ²/2) I ) = ∇(φ²/2) = φ∇φ
#
# moves it into the flux instead, giving the standard conservative momentum
# flux tensor
#
#   T_ij = φ u_i u_j + (φ²/2) δ_ij ,
#
# which is algebraically the same equation for smooth solutions and is the form
# every element-based Galerkin SWE code on the sphere uses (e.g. Giraldo &
# Warburton, 2008). It also fits Jexpresso's split exactly: everything under a
# divergence goes in user_flux!, everything pointwise goes in user_source!.
#
# On the manifold this stays correct because the RHS assembly differentiates
# the Cartesian components with SURFACE derivatives: ∂ⱼ(δᵢⱼ φ²/2) = ∂ᵢ(φ²/2) is
# then the surface gradient, which is tangential by construction — no spurious
# normal pressure force appears.
#
# THREE flux components on a TWO-dimensional manifold: F, G, H are the fluxes
# in the x, y and z directions, and the surface divergence
#
#   ∇·F = (∂ξ/∂x ∂F/∂ξ + ∂η/∂x ∂F/∂η)
#       + (∂ξ/∂y ∂G/∂ξ + ∂η/∂y ∂G/∂η)
#       + (∂ξ/∂z ∂H/∂ξ + ∂η/∂z ∂H/∂η)
#
# uses the two reference directions (ξ,η) of the shell with the 3×2 metric of
# the manifold, built by src/kernel/mesh/sphere_metrics.jl and applied by
# src/kernel/operators/sphere_rhs.jl.
#---------------------------------------------------------------------------------

function user_flux!(F, G, H, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    u = φu/φ
    v = φv/φ
    w = φw/φ

    # the isotropic part of the momentum flux, φ²/2, whose divergence is φ∇φ
    p = 0.5*φ*φ

    # mass
    F[1] = φu
    G[1] = φv
    H[1] = φw

    # momentum: T_ij = φ u_i u_j + (φ²/2) δ_ij
    F[2] = φu*u + p
    G[2] = φu*v
    H[2] = φu*w

    F[3] = φv*u
    G[3] = φv*v + p
    H[3] = φv*w

    F[4] = φw*u
    G[4] = φw*v
    H[4] = φw*w + p
end


function user_flux!(F, G, H, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::PERT; neqs=4, ip=1)
    #
    # No perturbation form for this case: the Galewsky initial state is stored
    # as a full field (q.qe holds the unperturbed BALANCED state, which is a
    # reference to difference against, not a hydrostatic background to linearize
    # around). Route PERT to the same fluxes so a deck that sets
    # :SOL_VARS_TYPE => PERT() still gets the right equations.
    #
    user_flux!(F, G, H, SD, q, qe, mesh, CL(), TOTAL(); neqs=neqs, ip=ip)
end


function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)

    φ  = q[1]
    φu = q[2]
    φv = q[3]
    φw = q[4]

    u = φu/φ
    v = φv/φ
    w = φw/φ
    p = T(0.5)*φ*φ

    return T(φu), T(φv), T(φw),
           T(φu*u + p), T(φu*v),     T(φu*w),
           T(φv*u),     T(φv*v + p), T(φv*w),
           T(φw*u),     T(φw*v),     T(φw*w + p)
end

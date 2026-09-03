#---------------------------------------------------------------------------------
# SWsphere_ScottPolvani_DynSGS — fluxes of the shallow water equations on a spherical
# shell. IDENTICAL to problems/ShallowWater/SWsphere/user_flux.jl, which has the
# full derivation; in short, Marras, Kopera & Giraldo (2015) Eq. (8) with the
# pressure term moved into the flux,
#
#   ∂φ/∂t     + ∇·(φu)                   = 0
#   ∂(φu)/∂t  + ∇·(φu⊗u + (φ²/2) I)      = -f(x × φu) + μx + φu_F - ν_r φu + δν∇²(φu)
#
# STATE VECTOR  q = [φ, φu, φv, φw]ᵀ, φ = gh, (u,v,w) the FULL CARTESIAN velocity.
#
# The forcing φu_F and the Rayleigh friction are nodal terms added after
# assembly (src/kernel/operators/sphere_forcing.jl), not fluxes.
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
    # No perturbation form for this case: q.qe holds the REST state (uniform
    # depth, no motion), which is the reference the radiative relaxation pulls
    # the height back to, not a background to linearise around. Route PERT to
    # the same fluxes so a deck that sets :SOL_VARS_TYPE => PERT() still gets
    # the right equations.
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

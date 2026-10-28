#---------------------------------------------------------------------------------
# Korteweg-de Vries equation through its hyperbolic approximation (KdVH).
#
# KdV, in the scaling of Biswas, Ketcheson, Ranocha & Schütz, J. Sci. Comput.
# 103:90 (2025), arXiv:2412.17117, eq. (1):
#
#     ∂ₜη + η ∂ₓη + ∂ₓ³η = 0
#
# The third derivative has no place in a first-order flux, so the case solves
# their hyperbolic relaxation, eqs. (2a)-(2c), instead:
#
#     ∂ₜu + u ∂ₓu + ∂ₓw = 0
#     τ ∂ₜv = ∂ₓv - w
#     τ ∂ₜw = -(∂ₓu - v)
#
# As τ → 0 the two auxiliary variables relax onto the first and second
# derivatives of u,
#
#     v → ∂ₓu         w → ∂ₓv = ∂ₓ²u
#
# and the first row becomes KdV with u → η. In conservation form
# ∂ₜU + ∂ₓF(U) = S, with U = (u, v, w)ᵀ,
#
#     F(U) = (u²/2 + w, -v/τ, u/τ)ᵀ
#     S(U) = (0, -w/τ, v/τ)ᵀ
#
# The derivatives sit in the flux and the undifferentiated relaxation terms in
# the source, the same split as the Timoshenko beam (problems/Elasticity/
# timoshenko1d). The flux Jacobian at u = 0 has eigenvalues -1/τ and ±1/√τ
# (paper, Remark 2), so the system is hyperbolic, and the explicit time step
# shrinks like τ: that is the price of the relaxation.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Relaxation parameter and soliton, in one place. Every other user_*.jl of this
# case reads them from here. A method rather than `const`s for the same reason
# as timoshenko_properties(): run.jl re-`include`s these files when cases are
# switched within a session, and a const would collide.
#
#   τ   relaxation parameter. The KdVH soliton differs from the KdV one by
#       O(cτ), so τ = 1e-3 puts that model error near 0.1% of the amplitude.
#   c   soliton speed; the amplitude is 3c.
#   x0  soliton centre at t = 0.
#---------------------------------------------------------------------------------
@inline function kdv_parameters()
    return (τ  = 1.0e-3,
            c  = 1.0,
            x0 = -10.0)
end

#---------------------------------------------------------------------------------
# KdV soliton, eq. (3) of the paper,
#
#     η(x,t) = 3c sech²(√c (x - ct - x0)/2)
#
# together with its exact first and second x-derivatives, which are what v and
# w converge to. With ξ = √c (x - ct - x0)/2, s = sech ξ, θ = tanh ξ:
#
#     η    = 3c s²
#     ∂ₓη  = -3c√c s² θ
#     ∂ₓ²η = -(3c²/2) s² (s² - 2θ²)
#---------------------------------------------------------------------------------
@inline function kdv_soliton(x, t)
    p = kdv_parameters()
    c = p.c

    ξ = 0.5*sqrt(c)*(x - c*t - p.x0)
    s = sech(ξ)
    θ = tanh(ξ)

    η    = 3.0*c*s*s
    ηx   = -3.0*c*sqrt(c)*s*s*θ
    ηxx  = -1.5*c*c*s*s*(s*s - 2.0*θ*θ)

    return η, ηx, ηxx
end

#---------------------------------------------------------------------------------
# F(U) = (u²/2 + w, -v/τ, u/τ)ᵀ.   G is unused in 1D.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=3, ip=1)

    τ = kdv_parameters().τ

    u = q[1]
    v = q[2]
    w = q[3]

    F[1] = 0.5*u*u + w
    F[2] = -v/τ
    F[3] =  u/τ
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    τ = kdv_parameters().τ

    return T(0.5*q[1]*q[1] + q[3]), T(-q[2]/τ), T(q[1]/τ)
end

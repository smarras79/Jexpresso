#---------------------------------------------------------------------------------
# Kuramoto-Sivashinsky equation in conservation form.
#
# The KS equation for a front φ(x,t),
#
#     ∂ₜφ + ∂ₓ²φ + ∂ₓ⁴φ + ½(∂ₓφ)² = 0,
#
# becomes a conservation law for its slope u = ∂ₓφ (differentiate once in x):
#
#     ∂ₜu + ∂ₓ(u²/2) + ∂ₓ²u + ∂ₓ⁴u = 0
#
# which is what this case solves. ∂ₓ²u is anti-diffusion, the long-wave
# instability, and ∂ₓ⁴u the short-wave damping: a mode e^{ikx} grows at
# σ(k) = k² - k⁴, so 0 < k < 1 is unstable and the fastest mode is k = 1/√2.
#
# Jexpresso applies second derivatives through its viscous operator, which
# adds ∂ₓ(ν_ieq ∂ₓ uprimitive[ieq]) to equation ieq for whatever variable the
# case puts in uprimitive[ieq] (user_primitives.jl). The fourth derivative is
# reached by relaxing a second field w onto ∂ₓ²u with a short time τ,
#
#     ∂ₜu + ∂ₓ(u²/2) = -w - ∂ₓ²w
#     τ ∂ₜw          = ∂ₓ²u - w
#
# As τ → 0, w → ∂ₓ²u and the first row is KS. Split for the solver, with
# U = (u, w)ᵀ:
#
#     inviscid flux   F = (u²/2, 0)ᵀ           (user_flux!, this file)
#     source          S = (-w, -w/τ)ᵀ          (user_source.jl)
#     viscous         ν = (-1, 1/τ) on uprimitive = (w, u)
#                                              (user_inputs.jl :μ,
#                                               user_primitives.jl)
#
# Both KS terms of the u row are written with w, not -∂ₓ²(u + w): the
# linearised system then obeys
#
#     τ s² + s + k⁴ - k² = 0,
#
# whose slow root is s ≈ σ(1 - τσ) and whose fast roots have Re s = -1/(2τ)
# for every k, so it is stable. Writing ∂ₓ²u directly instead leaves the
# anti-diffusion unbalanced once τk² > 1 (Re s → k²/2): unstable at the mesh
# scale.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Parameters, in one place. A method rather than `const`s for the same reason
# as timoshenko_properties(): run.jl re-`include`s these files when cases are
# switched within a session, and a const would collide.
#
#   L   domain length, 32π: the classic run of Kassam & Trefethen, SIAM J. Sci.
#       Comput. 26:1214 (2005), Fig. 4.1 (their kursiv.m), with
#       u(x,0) = cos(x/16)(1 + sin(x/16)).
#   τ   relaxation time of w. It shifts each growth rate by ≈ -τσ², which is
#       below 1e-5 on the unstable band 0 < k < 1.
#---------------------------------------------------------------------------------
@inline function ks_parameters()
    return (L = 32.0π,
            τ = 1.0e-4)
end

#---------------------------------------------------------------------------------
# The initial condition and its exact second derivative, which initialize.jl
# gives to w so the run starts on the relaxed manifold:
#
#     u     = cos(x/16) (1 + sin(x/16)) = cos(x/16) + ½ sin(x/8)
#     ∂ₓ²u  = -cos(x/16)/256 - sin(x/8)/128
#---------------------------------------------------------------------------------
@inline function ks_initial(x)
    u   =  cos(x/16.0) + 0.5*sin(x/8.0)
    uxx = -cos(x/16.0)/256.0 - sin(x/8.0)/128.0
    return u, uxx
end

#---------------------------------------------------------------------------------
# F(U) = (u²/2, 0)ᵀ.   G is unused in 1D.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=2, ip=1)

    F[1] = 0.5*q[1]*q[1]
    F[2] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(0.5*q[1]*q[1]), T(0.0)
end

#---------------------------------------------------------------------------------
# Kuramoto-Sivashinsky equation in conservation form.
#
# The KS equation for a front φ(x,t),
#
#     ∂ₜφ + ∂ₓ²φ + ∂ₓ⁴φ + ½(∂ₓφ)² = 0,
#
# becomes a conservation law for its slope u = ∂ₓφ (differentiate once in x):
#
#     ∂ₜu + ∂ₓ(u²/2) + ∂ₓ²u + ∂ₓ⁴u = 0,
#
# which is what this case solves. ∂ₓ²u is anti-diffusion, the long-wave
# instability, and ∂ₓ⁴u the short-wave damping: a mode e^{ikx} grows at
# σ(k) = k² - k⁴, so 0 < k < 1 is unstable and the fastest mode is k = 1/√2.
#
# The equation is split as du/dt = L u + N(u):
#
#     L u  = -∂ₓ²u - ∂ₓ⁴u       stiff and linear: user_inputs.jl declares it
#                               (:lsplit_linear, :split_d2 = :split_d4 = -1)
#                               and the time integrator treats it exactly
#     N(u) = -∂ₓ(u²/2)          this file: the ordinary inviscid flux
#
# L is built on the spectral element mesh by kernel/solvers/split_linear.jl:
# ∂ₓ² is the SEM weak Laplacian D2 = -M⁻¹K (the operator the viscous term
# applies) and ∂ₓ⁴ = D2², the C⁰ mixed fourth derivative. ETDRK4 (Cox &
# Matthews 2002, as used by Kassam & Trefethen 2005) integrates L exactly,
# through its φ-functions, so the step is set by the nonlinear dynamics and
# not by the ∂ₓ⁴ eigenvalues: Δt = 1/4, the step of Kassam & Trefethen,
# where an explicit scheme needs Δt ≈ 2e-4.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Parameters, in one place. A method rather than `const`s for the same reason
# as timoshenko_properties(): run.jl re-`include`s these files when cases are
# switched within a session, and a const would collide.
#
#   L   domain length, 32π: the run of Kassam & Trefethen, SIAM J. Sci.
#       Comput. 26:1214 (2005), Fig. 4.1 (their kursiv.m).
#---------------------------------------------------------------------------------
@inline ks_parameters() = (L = 32.0π,)

# Kassam & Trefethen's initial condition: u(x,0) = cos(x/16)(1 + sin(x/16)).
@inline ks_initial(x) = cos(x/16.0)*(1.0 + sin(x/16.0))

#---------------------------------------------------------------------------------
# F(u) = u²/2.   G is unused in 1D.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    F[1] = 0.5*q[1]*q[1]
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(0.5*q[1]*q[1])
end

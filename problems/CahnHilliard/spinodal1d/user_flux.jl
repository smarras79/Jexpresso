#---------------------------------------------------------------------------------
# Cahn-Hilliard equation through a relaxation of the chemical potential.
#
#     ∂ₜc = ∂ₓ²μ,        μ = f'(c) - ε² ∂ₓ²c,        f'(c) = c³ - c
#
# Two second derivatives, of μ and of c, make a fourth-order equation. The
# solver applies the second derivatives through its viscous operator, which
# adds ∂ₓ(ν_ieq ∂ₓ uprimitive[ieq]) to equation ieq for ANY variable the case
# puts in uprimitive[ieq] (user_primitives.jl). μ is not an evolved quantity,
# so it is relaxed onto its definition with a rate 1/τ:
#
#     ∂ₜc                 = ∂ₓ²μ
#     τ ∂ₜμ               = c³ - c - ε² ∂ₓ²c - μ
#
# As τ → 0, μ → c³ - c - ε²∂ₓ²c and the first row is Cahn-Hilliard. Split for
# the solver, with U = (c, μ)ᵀ:
#
#     inviscid flux   F = (0, 0)ᵀ                  (user_flux!, this file)
#     source          S = (0, (c³ - c - μ)/τ)ᵀ     (user_source.jl)
#     viscous         ν = (1, -ε²/τ) on uprimitive = (μ, c)
#                                                  (user_inputs.jl :μ,
#                                                   user_primitives.jl)
#
# The negative coefficient is not anti-diffusion of μ: it multiplies ∂ₓ²c in
# the μ row. Linearised about c = 0, a mode e^{ikx+st} obeys
#
#     τ s² + s + k²(ε²k² - 1) = 0
#
# whose slow root is s ≈ σ(k)(1 - τσ(k)), σ(k) = k²(1 - ε²k²) the
# Cahn-Hilliard rate, and whose roots for large k have Re s = -1/(2τ) and
# |Im s| ≈ εk²/√τ. So the relaxation is stable, and it replaces the explicit
# fourth-order limit Δt ~ 1/(ε²k⁴) by an oscillation of frequency εk²/√τ.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Model and test parameters, in one place. Every other user_*.jl of this case
# reads them from here. A method rather than `const`s for the same reason as
# timoshenko_properties(): run.jl re-`include`s these files when cases are
# switched within a session, and a const would collide.
#
#   ε   interface width parameter. ε = 1/(4√2) puts the most unstable
#       wavenumber, k* = 1/(√2 ε), at exactly k = 4, with rate σ* = 1/(4ε²) = 8.
#   τ   relaxation time of μ. The growth rate is off by ≈ τσ, i.e. 0.08% at
#       τ = 1e-4.
#   k   wavenumber of the initial perturbation, the fastest growing one.
#   δ   its amplitude, small enough for the linear phase to last until t ≈ 0.5.
#---------------------------------------------------------------------------------
@inline function ch_parameters()
    return (ε = 1.0/(4.0*sqrt(2.0)),
            τ = 1.0e-4,
            k = 4.0,
            δ = 1.0e-3)
end

# Cahn-Hilliard growth rate of the mode cos(kx) about c = 0.
@inline function ch_growth_rate(k)
    ε = ch_parameters().ε
    return k*k*(1.0 - ε*ε*k*k)
end

#---------------------------------------------------------------------------------
# The perturbation and its exact first and second derivatives, together with
# the chemical potential they define:
#
#     c     = δ cos(kx)
#     ∂ₓc   = -δ k sin(kx)
#     ∂ₓ²c  = -δ k² cos(kx)
#     μ     = c³ - c - ε² ∂ₓ²c
#
# initialize.jl starts μ from this μ, so the run starts on the relaxed
# manifold and excites no O(τ) initial layer.
#---------------------------------------------------------------------------------
@inline function ch_initial(x)
    p = ch_parameters()

    c   =  p.δ*cos(p.k*x)
    cx  = -p.δ*p.k*sin(p.k*x)
    cxx = -p.δ*p.k*p.k*cos(p.k*x)
    μ   =  c*c*c - c - p.ε*p.ε*cxx

    return c, cx, cxx, μ
end

#---------------------------------------------------------------------------------
# Amplitude of the steady state the run saturates into.
#
# By symmetry the mode cos(kx) saturates into a steady profile of the same
# period 2π/k with μ = 0, i.e. a periodic orbit of
#
#     ε² c'' = c³ - c,       c(0) = A,  c'(0) = 0
#
# With F(c) = c⁴/4 - c²/2, energy conservation gives the period
#
#     P(A) = 4ε ∫₀ᴬ dc / √(2(F(c) - F(A)))
#          = 4ε ∫₀^{π/2} dθ / √(1 - A²(1 + sin²θ)/2)        (c = A sin θ)
#
# which grows from 2πε at A → 0 to ∞ at A → 1. Bisection on P(A) = 2π/k
# gives the amplitude: 0.80078 for the shipped parameters. Midpoint rule on
# the smooth θ-integrand; called only at output times.
#---------------------------------------------------------------------------------
function ch_equilibrium_amplitude()
    p = ch_parameters()

    n = 4000
    h = 0.5π/n
    period(A) = 4.0*p.ε*h*sum(1.0/sqrt(1.0 - 0.5*A*A*(1.0 + sin((i - 0.5)*h)^2)) for i = 1:n)

    target = 2π/p.k
    period(0.0) >= target && return 0.0     # k ≥ 1/ε: c = 0 is stable

    lo, hi = 0.0, 1.0 - 1.0e-12
    for _ = 1:60
        m = 0.5*(lo + hi)
        period(m) < target ? (lo = m) : (hi = m)
    end
    return 0.5*(lo + hi)
end

#---------------------------------------------------------------------------------
# F = 0: every derivative of this system goes through the viscous operator.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=2, ip=1)

    F[1] = 0.0
    F[2] = 0.0
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    return T(0.0), T(0.0)
end

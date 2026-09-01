#---------------------------------------------------------------------------------
# 2D linear elastodynamics, velocity-stress form.
#
# This is the genuine 2D counterpart of the Timoshenko beam: a first-order
# symmetric hyperbolic system in conservation form, with no reduced-model
# kinematic assumptions in it at all.
#
# State (plane strain), per unit volume:
#
#     U = (ρu, ρv, σxx, σyy, σxy)ᵀ
#
#     u, v            particle velocity
#     σxx, σyy, σxy   Cauchy stress
#
# Balance of linear momentum and the rate form of Hooke's law:
#
#     ∂ₜ(ρu)  = ∂ₓσxx + ∂_y σxy
#     ∂ₜ(ρv)  = ∂ₓσxy + ∂_y σyy
#     ∂ₜσxx   = (λ+2μ) ∂ₓu + λ ∂_y v
#     ∂ₜσyy   = λ ∂ₓu + (λ+2μ) ∂_y v
#     ∂ₜσxy   = μ (∂_y u + ∂ₓv)
#
# i.e. ∂ₜU + ∂ₓF(U) + ∂_y G(U) = 0 with
#
#     F(U) = -(σxx, σxy, (λ+2μ)u, λu,      μv)ᵀ
#     G(U) = -(σxy, σyy, λv,      (λ+2μ)v, μu)ᵀ
#
# NOTE THE SOURCE IS EMPTY. Unlike Timoshenko — whose shear coupling κₛGAγ and
# whose -ω are undifferentiated and therefore have to sit in S — every term of
# 2D elastodynamics is a divergence. That is not a detail: a Dirichlet-pinned
# boundary variable here has ∂ₜ(·) = 0 at the boundary for the exact solution,
# because the flux divergence that feeds it vanishes there too.
#
# The flux Jacobian in any direction n has eigenvalues
#
#     ±cp = ±√((λ+2μ)/ρ)   pressure (P) wave
#     ±cs = ±√(μ/ρ)        shear (S) wave
#      0                   the stress combination carrying no motion along n
#
# so the system is hyperbolic, and symmetrizable by the mechanical energy
#
#     ½ρ(u² + v²) + ½ σ : C⁻¹ : σ
#
# which is the convex entropy; its flux is -σ·v, the mechanical power. For
# ν = 1/4 the two speeds sit at cp/cs = √3; here ν = 3/10 gives 1.8708.
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Material, in one place. Every other user_*.jl of this case reads it from here.
#
# A METHOD and not a set of `const`s on purpose: src/run.jl loads a case by
# `include`ing its user_*.jl files into the Jexpresso module, so a `const`
# would collide with the sibling case (Elasticity/beam2d) the first time
# the two are run in the same session. Redefining a method at an identical
# signature is exactly what an `include` is supposed to do.
#
# Non-dimensional: unit density and Young's modulus, ν = 0.3 — the same
# material as the Timoshenko cases, so the two families are comparable.
#---------------------------------------------------------------------------------
@inline function elastic_properties()

    ρ = 1.0        # density
    E = 1.0        # Young's modulus
    ν = 0.3        # Poisson's ratio

    λ  = E*ν/((1.0 + ν)*(1.0 - 2.0*ν))   # Lamé λ   (plane strain)
    μs = E/(2.0*(1.0 + ν))               # Lamé μ = shear modulus.  NOT `μ`:
                                         # :μ is also an inputs key, and μ is
                                         # easy to confuse with the viscosity
                                         # the solver carries under that name.

    return (ρ = ρ, E = E, ν = ν, λ = λ, μ = μs, λ2μ = λ + 2.0*μs,
            cp = sqrt((λ + 2.0*μs)/ρ),   # P-wave speed
            cs = sqrt(μs/ρ))             # S-wave speed
end

#---------------------------------------------------------------------------------
# The plane P-wave this case is initialised with, and the one it must
# reproduce after every whole period.
#
# For a wave vector k = (2π p/L, 2π q/L) on a doubly periodic square of side L,
# with n = k/|k| and φ = k·x - ωt, ω = cp|k|:
#
#     displacement  U n sin(φ)
#     velocity      -U ω n cos(φ)
#     stress        σij = U|k| (λ δij + 2μ ni nj) cos(φ)
#
# Integer p, q make k·x periodic in both directions, so the wave simply
# re-enters and there are no boundaries anywhere in the problem. p = q = 1
# puts it on the diagonal, which is what makes it a real 2D test: all five
# components are active and both flux directions carry the wave.
#---------------------------------------------------------------------------------
@inline function elastic_plane_wave()

    p = 1                    # wave numbers, in units of 2π/L
    q = 1
    L = 1.0                  # domain side (must match the mesh)
    A = 1.0e-4               # displacement amplitude: strains stay ~1e-3

    P = elastic_properties()

    kx = 2.0*π*p/L
    ky = 2.0*π*q/L
    k  = sqrt(kx*kx + ky*ky)
    nx = kx/k
    ny = ky/k
    ω  = P.cp*k

    return (kx = kx, ky = ky, k = k, nx = nx, ny = ny, ω = ω, A = A, L = L)
end

#---------------------------------------------------------------------------------
# The five state components of that wave at (x, y, t).
#---------------------------------------------------------------------------------
@inline function elastic_plane_wave_state(x, y, t)

    P = elastic_properties()
    W = elastic_plane_wave()

    c = cos(W.kx*x + W.ky*y - W.ω*t)
    a = W.A*W.k*c

    return (ρu  = -P.ρ*W.A*W.ω*W.nx*c,
            ρv  = -P.ρ*W.A*W.ω*W.ny*c,
            σxx = a*(P.λ + 2.0*P.μ*W.nx*W.nx),
            σyy = a*(P.λ + 2.0*P.μ*W.ny*W.ny),
            σxy = a*2.0*P.μ*W.nx*W.ny)
end

function user_flux!(F, G, SD::NSD_2D, q, qe,
                    mesh::St_mesh, ::CL, ::TOTAL; neqs=5, ip=1)

    P = elastic_properties()

    u   = q[1]/P.ρ
    v   = q[2]/P.ρ
    σxx = q[3]
    σyy = q[4]
    σxy = q[5]

    F[1] = -σxx
    F[2] = -σxy
    F[3] = -P.λ2μ*u
    F[4] = -P.λ*u
    F[5] = -P.μ*v

    G[1] = -σxy
    G[2] = -σyy
    G[3] = -P.λ*v
    G[4] = -P.λ2μ*v
    G[5] = -P.μ*u
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)
    P = elastic_properties()

    u = q[1]/P.ρ
    v = q[2]/P.ρ

    return (T(-q[3]), T(-q[5]), T(-P.λ2μ*u), T(-P.λ*u), T(-P.μ*v),
            T(-q[5]), T(-q[4]), T(-P.λ*v),   T(-P.λ2μ*v), T(-P.μ*u))
end

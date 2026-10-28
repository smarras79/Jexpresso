# KdV

The Korteweg-de Vries equation

```
∂ₜη + η ∂ₓη + ∂ₓ³η = 0
```

has a third derivative, which a first-order flux cannot carry. These cases
solve the hyperbolic approximation KdVH of Biswas, Ketcheson, Ranocha & Schütz,
*J. Sci. Comput.* 103:90 (2025), [arXiv:2412.17117](https://arxiv.org/abs/2412.17117),
eqs. (2a)-(2c):

```
∂ₜu + u ∂ₓu + ∂ₓw = 0
τ ∂ₜv = ∂ₓv - w
τ ∂ₜw = -(∂ₓu - v)
```

As τ → 0, v → ∂ₓu and w → ∂ₓ²u, so u solves KdV. In conservation form
`∂ₜU + ∂ₓF(U) = S`, with `U = (u, v, w)ᵀ`:

```
F(U) = (u²/2 + w, -v/τ, u/τ)ᵀ
S(U) = (0, -w/τ, v/τ)ᵀ
```

At u = 0 the flux Jacobian has eigenvalues -1/τ and ±1/√τ, so the system is
hyperbolic. The explicit time step has to resolve the fastest one, which makes
it O(τ).

```bash
julia --project=. src/Jexpresso.jl KdV soliton1d
```

| case | dim | eqs | grid | boundaries |
|---|---|---|---|---|
| `soliton1d` | 1D | 3 | built-in, 120 elements, p = 4, on [-30, 30] | periodic |

## soliton1d

A KdV soliton, eq. (3) of the paper:

```
η(x,t) = 3c sech²(√c (x - ct - x0)/2)
```

with c = 1 and x0 = -10, run to t = 20 with τ = 1e-3. The two auxiliary fields
start from the soliton's exact first and second derivatives (`kdv_soliton` in
`user_flux.jl`). These are the well-prepared data the relaxation converges to,
so no initial layer is excited. Every output frame overlays the exact KdV
soliton and its derivatives on u, v and w (`user_analytic.jl`).

The overlay is the KdV solution, not the KdVH one. The gap between the curves
is therefore discretisation error plus the O(cτ) model error of the
relaxation. τ, c and x0 are set in `kdv_parameters()` in `user_flux.jl`. If
you lower τ, lower `:Δt` in proportion.

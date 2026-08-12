# SWsphere_visc — the Galewsky jet stabilised by artificial diffusion

The same test as [`../SWsphere`](../SWsphere/README.md) — same equations, same
cubed sphere, same Galewsky et al. (2004) barotropically unstable jet — with the
**artificial diffusion** `δν∇²(φu)` of Marras, Kopera & Giraldo (2015) Eq. (8b)
in place of the modal filter.

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_visc")
```

## What differs from SWsphere

Two lines of `user_inputs.jl`, and nothing else:

|  | `SWsphere` | `SWsphere_visc` |
|---|---|---|
| `:lfilter` | `true` | **`false`** |
| `:lvisc`, `:μ` | `false`, `0.0` | **`true`, `1.0e5`** |

The other five `user_*.jl` files in this directory are one-line `include`s of the
`SWsphere` ones, and `:gmsh_filename` points at that case's `cubed_sphere.msh`,
so the equations, the initial condition and the grid are the same **by
construction** rather than by having been copied and kept in sync. (`src/run.jl`
includes exactly six `user_*.jl` per case directory and looks nowhere else,
which is why the five stubs have to exist.)

## Why it is a separate case rather than a comment

The shell needs **at least one** of the two stabilisation mechanisms: the
inviscid, unfiltered high-order solution on the cubed sphere grows grid-scale
modes and blows up — which is what the paper reports in section 4.2, and what
you get if you switch the filter off without switching anything on in its place.
`SWsphere` keeps the filter, because that is the paper's own choice for the
published test. This case exercises the other branch, so both are covered by
something runnable.

Note that a run with **both** on is perfectly legitimate — the two mechanisms
compose, they are simply both dissipative. It just would not isolate either one.

## The term

The surface (Laplace-Beltrami) Laplacian, built from the same 3×2 manifold
metrics as the divergence and assembled in weak form — integration by parts
leaves no boundary term on a closed shell:

```
∫ ψ ν∇ₛ²q dΩ = −∫ ∇ₛψ·(ν ∇ₛq) dΩ ,     ∇ₛq = a¹ ∂q/∂ξ + a² ∂q/∂η
```

One derivative of the C⁰ solution instead of two, and symmetric negative
semi-definite by construction, so the term can only ever remove energy.
`:ivisc_equations => [2,3,4]` puts it on the momentum only, where Eq. (8b) has
it; the continuity equation carries no diffusion.

It is **not** the viscous path of `src/kernel/operators/rhs.jl`. That one builds
`∇q` from the flat 2×2 inverse Jacobian, which has no third metric direction —
dropping that component is not an approximation but a different operator, one
that differentiates along the projection of the shell onto `z = 0` and is
singular at the equator. See `src/kernel/operators/sphere_rhs.jl`.

## Cost

None, at this resolution. Diffusion is a second derivative, so it is
explicit-stable only for `Δt ~ Δ²/ν`, and `sphere_cfl_dt` takes the min of that
and the wave-speed limit — but `Δmin²/ν ≈ 1e5 s` against a gravity-wave step of
`~1e2 s`, so the waves still set `Δt` by three orders of magnitude. `ν = 1e5
m²/s` is nowhere near the stability boundary for this grid. Raising `ν` far
enough does start to cut the step, which is the point of taking the min: the
alternative is a silent blow-up.

## Output

`problems/ShallowWater/SWsphere_visc/output/`, same fields and same conventions
as the parent case — `sphere_grid_ho.vtu` plus `sphere_0001.vtu` … The field the
test is judged on is `vorticity`; `h` barely moves while the instability
develops. See the parent README for the full list and for how to view it.

## Unit test of the operator

`test/test_sphere_visc.jl` checks the assembled Laplacian against the
spherical-harmonic eigenvalues `Δₛ Y_l = −l(l+1)/R²·Y_l` for `l = 1…4`, plus
annihilation of constants, symmetry and negative-definiteness. A term that
merely damped would keep this case from blowing up and still fail that test.

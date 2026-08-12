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
Worth knowing: the paper itself runs **both**, filtering at every step (§4.2)
*and* carrying `ν = 1e5`. Viscosity alone, as this deck has it, is a deliberate
isolation of the term, not a reproduction of the published configuration.

## What a healthy run looks like

One simulated day at the shipped resolution, `:lfilter => false`, `ν = 1e5`,
against the same run with **neither** mechanism:

| step (t) | `δE/E` with ν | `δE/E` with neither | `max\|ζ\|` with ν | `max\|ζ\|` with neither |
|---|---|---|---|---|
| 200 (0.17 d) | −2.07e-05 | +3.98e-08 | 1.128e-04 | 1.148e-04 |
| 600 (0.52 d) | −6.09e-05 | +4.96e-07 | 1.099e-04 | 1.235e-04 |
| 1000 (0.87 d) | −9.96e-05 | +1.32e-05 | 1.122e-04 | 4.866e-04 |
| 1153 (1.00 d) | −1.14e-04 | +2.37e-04 | 1.125e-04 | 1.880e-03 |

The **sign of `δE/E` is the whole story**. With the diffusion, energy decreases
monotonically — that is what a correctly-signed dissipative term does, and it is
the negative semi-definiteness of the weak form showing up in a run. With
neither, energy is *created*, roughly ten-fold per 0.17 d by the end, and
`max|ζ|` has grown 16× while the viscous run's is flat. That is the grid-scale
instability, and it is what eventually blows up.

Mass is conserved to 3e-12 and the constraint `max|(φu)·x̂|` stays at 4e-10
throughout the viscous run.

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

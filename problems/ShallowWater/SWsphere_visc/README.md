# SWsphere_visc — the Galewsky jet with filter *and* artificial diffusion

The same test as [`../SWsphere`](../SWsphere/README.md) — same equations, same
cubed sphere, same Galewsky et al. (2004) barotropically unstable jet — with the
**artificial diffusion** `δν∇²(φu)` of Marras, Kopera & Giraldo (2015) Eq. (8b)
added on top of the modal filter. That pairing is the paper's own
configuration: it filters at every step (§4.2) *and* carries `ν = 1e5 m²/s`.

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_visc")
```

## What differs from SWsphere

One line of `user_inputs.jl`:

|  | `SWsphere` | `SWsphere_visc` |
|---|---|---|
| `:lfilter` | `true` | `true` |
| `:lvisc`, `:μ` | `false`, `0.0` | **`true`, `1.0e5`** |

The other five `user_*.jl` files in this directory are one-line `include`s of the
`SWsphere` ones, and `:gmsh_filename` points at that case's `cubed_sphere.msh`,
so the equations, the initial condition and the grid are the same **by
construction** rather than by having been copied and kept in sync. (`src/run.jl`
includes exactly six `user_*.jl` per case directory and looks nowhere else,
which is why the five stubs have to exist.)

## Why not viscosity alone

This case first shipped with `:lfilter => false`, on the reasoning that the
diffusion is a stabilisation mechanism in its own right and the case should
isolate it. **That deck reliably blew up at 2 days.** Measured, 3 simulated
days, everything else equal:

| stabilisation | outcome | `max\|ζ\|`, 0 → 3 d | `δE/E` at 3 d |
|---|---|---|---|
| `ν = 1e5`, filter **off**, momentum only | **NaN at 2.005 d** | — | — |
| `ν = 1e5`, filter **off**, all 4 equations | **NaN at 2.009 d** | — | — |
| `ν = 5e5`, filter **off** | completes | 1.12e-4 → **5.65e-5 (−50%)** | −1.09e-03 |
| `ν = 1e6`, filter **off** | completes | 1.12e-4 → **3.78e-5 (−66%)** | −1.57e-03 |
| **`ν = 1e5`, filter on** ← this deck | completes | 1.11e-4 → **9.58e-5 (−13%)** | −3.12e-04 |

Three things are worth taking from that table.

**It is not the undiffused φ.** The obvious suspicion about momentum-only
diffusion is that it leaves the continuity equation `∂φ/∂t + ∇·(φu) = 0` — pure
advection, and CG has no upwinding — with no dissipation at all. Adding φ to
`:ivisc_equations` tests that directly, and the run still dies, within 0.4% of
the same time. The instability is not living in the height field.

**ν = 1e5 is simply too weak for this grid.** 10 elements per panel at `nop=5`
puts the coarse elements at ~219 km effective resolution, so the grid Reynolds
number is `uΔ/ν = 80 × 219162 / 1e5 ≈ 175`. Viscosity alone holds the grid only
from about `ν = 5e5` upward.

**And raising ν buys stability by erasing the answer.** ν = 5e5 completes the
three days with half of `max|ζ|` gone; ν = 1e6 completes with two thirds gone.
Relative vorticity is the field the Galewsky test is judged on, so that is a
stable run of the wrong problem. The trend is monotone and there is no window
where viscosity alone is both stable and faithful. The filter, which damps only
the top modes rather than every scale, costs 13% over the same interval.

The lesson generalises: on this shell a blow-up is almost always **too little
dissipation**, not too large a step. `ν = 1e5` sits three orders below the
diffusive stability ceiling (`Δmin²/ν ≈ 1e5 s` against a 75 s step), so lowering
`:cfl` does nothing for it — which is exactly the symptom that identifies it.

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

No time step, at this resolution. Diffusion is a second derivative, so it is
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

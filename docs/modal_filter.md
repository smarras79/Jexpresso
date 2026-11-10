# Modal filter (`:lfilter`)

The modal filter damps the highest polynomial modes inside each spectral element. It's the main control for grid-scale noise in the non-dissipative CG discretisation, for example in a stable inversion.

## Inputs

| Key | Default | Meaning |
|---|---|---|
| `:lfilter` | `false` | turn the filter on |
| `:mu_x`, `:mu_y`, `:mu_z` | `0.0` | strength per direction: the filter is `mu·F + (1 − mu)·I` |
| `:filter_type` | `"erf"` | transfer function (weights per mode): `"erf"` (Boyd–Vandeven), `"quad"`, `"exp"` |
| `:filter_basis` | `"legendre"` | modal basis the weights act on: `"legendre"` or `"hierarchical"` |

With nop = 4 and `"erf"`, modes 0–2 have weight 1, mode 3 about 0.996 and mode 4 weight 0. So in practice the filter removes a fraction `mu` of the degree-4 content per application.

## When it is applied

The filter is applied **once per completed time step** to the solution (`integrator.u`). It runs as a `DiscreteCallback` registered in `time_loop!` (`src/kernel/solvers/TimeIntegrators.jl`), which calls `filter_step!` (`src/kernel/operators/filter.jl`).

Before this change, the filter ran at the start of every right-hand-side evaluation and overwrote the stage input in place, i.e. at every Runge–Kutta stage. That altered the stage states the scheme combines, including the step's base state, which later stages reuse. The step was then no longer the RK scheme's update. Filtering once per step keeps every step consistent, and it is also cheaper.

**Retuning:** the same `mu` now acts once per step instead of once per stage (five times per step for SSPRK54). For small `mu`, roughly `mu_new ≈ (number of stages) × mu_old` gives the old damping per unit time. A case that was stable only thanks to the per-stage filtering may need a larger `mu`.

## Filter basis

The weights are applied to modal coefficients. The choice of modes matters:

- **`"legendre"` (default)**: orthogonal Legendre polynomials.
  - With GLL quadrature, the filter never increases the solution's mass-weighted L2 norm, including after the DSS step.
  - It changes element-end values; the DSS averaging restores continuity.
  - Element integrals, and therefore mass, are preserved.
- **`"hierarchical"`**: the two linear vertex modes plus the bubbles `L_j − L_{j−2}`, which vanish at the element ends. This was the only option before.
  - The filter leaves element-end values unchanged.
  - The basis is not orthogonal, so one application can **increase** the solution's energy. With nop = 4 the maximum increase is about +0.17% at mu = 0.01, +0.9% at mu = 0.05 and +1.8% at mu = 0.1.
  - It also never damps content concentrated on element-interface nodes.
  - Applied every step, this can feed a growing mode on element-interface nodes. It was seen in forced 3D BOMEX columns: a horizontally uniform mode that grew faster for larger `mu`, and that was absent with the filter off or with `"legendre"`.
  - Selecting it prints a warning.

These properties cannot all be had at once. No element-local filter can do all three of the following, other than the identity:
1. leave all polynomials of degree < N unchanged;
2. leave element-end values unchanged;
3. never increase the solution's mass-weighted L2 norm.

The reason: a filter with property 3 that fixes the degree < N space must map the remaining mode `L_N` to a multiple of itself. `L_N` is nonzero at the element ends, so damping it changes the end values. `"hierarchical"` keeps 1 and 2 and gives up 3. `"legendre"` keeps 1 and 3 and gives up 2, which the DSS averaging then restores.

The basis option applies to the Legendre–Gauss–Lobatto elements; the Laguerre (semi-infinite) filter is unchanged.

## Files

- `src/kernel/operators/filter.jl`: `init_filter` (filter matrices), `filter!` (application), `filter_step!` / `filter_step_gpu!` (once-per-step entry points)
- `src/kernel/solvers/TimeIntegrators.jl`: the filter callback
- `src/io/mod_inputs.jl`: defaults

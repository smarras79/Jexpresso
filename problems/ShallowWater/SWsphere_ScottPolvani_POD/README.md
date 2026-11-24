# SWsphere_ScottPolvani_POD — Scott & Polvani on the Galewsky grid, for POD

> **This case needs the forcing, and the forcing is not on this branch yet.**
> It runs end to end, the POD chain works, and the fluid stays at rest because
> nothing applies a forcing. See *Status* below. It is committed in this state
> so that the grid, the timing and the POD block are ready the moment
> `sphere_forcing.jl` lands — not because it produces a result today.

```julia
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_ScottPolvani_POD")
```

## What it is

[`SWsphere_ScottPolvani`](../SWsphere_ScottPolvani/README.md) — forced-
dissipative shallow-water turbulence with giant-planet parameters, after Scott
& Polvani (2007) — moved onto **the grid the Galewsky case runs on** and
configured for Proper Orthogonal Decomposition.

| | `SWsphere_ScottPolvani` | here |
|---|---|---|
| grid | `cubed_sphere_64x64.msh`, `npoin ≈ 6.1e5` | **`cubed_sphere.msh`, 10×10 panels, `npoin = 15 002`** |
| `SP_NROT` default | 500 | **20** |
| POD fields | `[:vorticity]` | **`[:vorticity, "h", :velocity]`** |
| snapshot memory | 300 MB (one field) | **29 MB (all three)** |
| raster | 1440 × 720 | 720 × 360 |

The three fields are the point. `:vorticity` is a **derivative**; `h` and
`velocity` are primitive. `ζ = ∇ₛ×u` amplifies grid-scale velocity noise by the
wavenumber, so its temporal fluctuation can be dominated by grid-scale content
while the primitive fields are smooth — and the report's `sub-element content`
line makes that visible per field. On the Galewsky case, where this comparison
can be run today ([`SWsphere_POD`](../SWsphere_POD/README.md)), it reads:

| field | sub-element content (fluctuation / mode 1) |
|---|---|
| `vorticity` | 0.581 / 0.594 |
| `h` | 0.366 / 0.376 |
| `velocity` | 0.426 / 0.404 |

## Status: what it does today, measured

`user_source.jl` of this case states where the forcing lives:

```
# drawn once per step, so they live in src/kernel/operators/sphere_forcing.jl
# and enter the RHS after assembly (switched by :lsphere_forcing in the deck).
```

That file is not in `src/`, and nothing under `src/` reads `:lsphere_forcing`,
`:forcing_epsilon`, `:rayleigh_friction` or `:radiative_relaxation`. The deck
sets them; nothing consumes them. Run on this grid, `SP_NROT=20`, 1012 steps to
16.54 days:

```
 # state: at rest, h = H everywhere; all motion comes from the forcing
 # step 1012  t = 1429214.7 s (16.542 d)  max|ζ| = 2.617e-17  |(φu)·x̂| = 7.95e-22

 #   vorticity: Σλ = 8.403652e-20 , max|q| = 2.689e-17
 #   h        : max|q| = 1.594e+04 , max|q_k - q_1| = 5.150e-09 , ratio = 3.23e-13
 #              <-- THE FIELD BARELY MOVED: the modes below are of the differences
 #   velocity : max|q| = 2.620e-11
```

Machine zero after sixteen days. Nothing here is broken — the equations, the
grid, the metrics, the time loop and the whole POD chain run, and the POD says
plainly that it was handed a field that does not move. The physics is simply
absent from the branch.

## What it should produce once the forcing lands

Scott & Polvani's section 7 runs organise small-scale forcing into banded zonal
jets with coherent vortices between them (their Figs. 13 and 14). For the
decomposition that means:

* the **temporal mean** is the jet structure of Fig. 13;
* the leading **modes** are the vortices and jet meanderings of Fig. 14 — and a
  structure that travels (a Rossby wave riding a jet) appears as a **pair** of
  modes of nearly equal energy in quadrature, a circle in the `(a₁,a₂)` phase
  portrait;
* the **spectrum** measures how many degrees of freedom the equilibrated flow
  really has, which is the size a reduced-order model of this regime would need.

The check that the modes mean something is the `sub-element content` line: the
mode should match the fluctuation it came from, and the primitive fields should
come out markedly smoother than the vorticity.

Full documentation: [`docs/POD.md`](../../../docs/POD.md).

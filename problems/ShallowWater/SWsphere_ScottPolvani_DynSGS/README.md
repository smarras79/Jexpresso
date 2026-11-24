# SWsphere_ScottPolvani_DynSGS — the Scott & Polvani giant-planet turbulence with DynSGS

The same problem as [`SWsphere_ScottPolvani`](../SWsphere_ScottPolvani/README.md)
— Scott & Polvani (2007), forced-dissipative shallow-water turbulence on the
sphere with Jupiter / Saturn / Neptune parameters, on the same cubed-sphere
grid, with the same equations, forcing, large-scale dissipation, planet
selector, environment overrides (`SP_PLANET`, `SP_NROT`, `SP_MESH`) and
outputs — with **one change: the small-scale stabilisation**.

| | `SWsphere_ScottPolvani` | `SWsphere_ScottPolvani_DynSGS` |
|---|---|---|
| small-scale dissipation | constant `ν∇ₛ²(φu)`, `ν = 2e-6 a²Ω`, **plus** the modal filter every step | **DynSGS**: a residual-based `ν_e∇ₛ²(φu)` with a coefficient computed per element and per RK stage; no filter |
| acts | everywhere, always | where and when the discrete equations are not satisfied |
| free numbers | `ν`, filter `α`, order, cutoff | `C₁ = 1`, `C₂ = 0.5` (the paper's), a normalisation floor |

Read the sibling's README for everything about the physics; this one only
covers the stabilisation.

## Run it

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_ScottPolvani_DynSGS")
```

## The model

Marras, Nazarov & Giraldo (2015), *Stabilized high-order Galerkin methods
based on a parameter-free dynamic SGS model for LES*, J. Comput. Phys. **301**,
77–101 — the model the flat cases run with `:visc_model => DSGS()` — carried to
the shell's own right-hand side in `src/kernel/operators/sphere_dsgs.jl`. Per
element `e`:

```
ν_res|e = C₁ Δ² max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω
ν_max|e = C₂ Δ  ‖ |u| + √φ ‖∞,e
ν|e     = max(0, min(ν_max|e, ν_res|e))
```

* `R_i = ∂q_i/∂t − RHS_i(q)` is the strong-form residual of equation `i` of
  `q = [φ, φu, φv, φw]`: the BDF2 derivative over the last two **completed**
  steps (the history is rolled once per step, never per stage) minus the
  assembled, `M⁻¹`-scaled inviscid right-hand side of the current stage,
  forcing and Rayleigh friction included.
* `Δ` is the element's shortest corner edge over `ngl`, the LGL spacing it
  resolves, from the 3-D coordinates.
* `⟨q_i⟩` and `‖q_i − ⟨q_i⟩‖∞,Ω` are the area mean and the deviation norm over
  the domain, rank-local unless `:ldsgs_global_norms`.
* the cap `ν_max` is what keeps the explicit step stable by construction,
  `Δt ν/Δ² ≤ C₂·CFL`; the time loop also feeds the cap on the initial state
  into the diffusive branch of its CFL condition.

`ν_e` multiplies the per-equation `:μ = [0, 1, 1, 1]` (momentum only, Marras
Eq. 10) in the same Laplace–Beltrami weak form the constant-`ν` path uses, so
the only difference in the operator is the coefficient. On a DynSGS run the
RHS is assembled in two passes — inviscid first, then `ν_e` from its residual,
then the viscous term — while the constant-`ν` path keeps its single pass and
is bit-for-bit what it was.

**The floor.** The paper divides by `‖q − ⟨q⟩‖∞,Ω`, which is *exactly zero* for
the resting uniform layer this problem starts from; the residual over zero
would pin `ν` at the cap everywhere until the flow developed. Each denominator
therefore has a floor at `:dsgs_floor` (default `1e-3`) times its natural
scale, `φ̄` for the geopotential and `φ̄√φ̄` for the momentum — the same idea as
the momentum floor of the flat θ path — and plays no role once real deviations
exceed it. A resting start has zero residual and begins inviscid.

## What to look at

Every `:ndiagnostics_prints` steps, one more line:

```
DynSGS: ν_e max = … m²/s  mean = …  at the cap: …% of the elements
```

Compare `max` against the sibling's constant `ν = 2e-6 a²Ω` (`1.8e6 m²/s` on
Jupiter) and watch the cap fraction: a run sitting at the cap in most elements
is telling you the grid is too coarse for what it is being asked to carry. The
VTK files carry `dsgs_nu`, the largest `ν_e` among the elements at each node,
next to `vorticity`, `pv` and `vorticity_forcing`, so *where* the model acts
is on the same file as the flow.

## Measured, Jupiter, the shipped grid

20 rotations, serial, same seed as the sibling's runs:

| | constant `ν` + filter | DynSGS |
|---|---|---|
| KE/mass at 2 rotations [m²/s²] | 0.91 | 0.69 |
| KE/mass at 20 rotations [m²/s²] | 7.19 | 4.16 |
| coefficient [m²/s] | `1.8e6` everywhere | `ν_e` max `2.6e7`, mean `0.9–1.3e7`, 0% of the elements at the cap |
| mass, drift | `1e-13`, `1e-9` | `6e-14`, `1e-9` |

So on this grid, in the spin-up, DynSGS is the more dissipative of the two:
its mean coefficient is 5–7× the sibling's constant `ν`, and the kinetic
energy at 20 rotations is 42% lower. That is the residual at work — the
forcing is redrawn every step, so the flow at the forced scale is never a
solution of the equations the residual is measured against, and the model
treats it as under-resolved. Nothing sits at the cap, i.e. the number is the
residual's own. `C₁` scales it linearly (`:dsgs_C1`) if a run is to be pushed
towards the sibling's level; a finer grid (`SP_MESH`) lowers `Δ²` and with it
`ν_res` quadratically, which is where the model's selectivity pays.

## Switches

| key | meaning |
|---|---|
| `:lvisc => true`, `:visc_model => DSGS()` | the model |
| `:μ` | per-equation multiplier of `ν_e`; `[0,1,1,1]` = momentum only |
| `:dsgs_C1`, `:dsgs_C2` | the paper's 1.0 and 0.5 |
| `:dsgs_floor` | normalisation floor, fraction of the natural scales (1e-3) |
| `:ldsgs_global_norms` | domain norms instead of rank-local ones |
| `:lfilter` | off here; switch on if `ζ` shows grid-scale noise |

## Tests

`julia --project=. test/test_sphere_dsgs.jl` checks the switch, that a resting
layer gets `ν_e = 0`, that `ν_e ≤ ν_max` always and is positive and finite on a
random state, that the two-pass RHS equals inviscid + separately assembled
viscous term, that the DynSGS viscous term only removes energy, and that the
constant-`ν` path is unchanged.

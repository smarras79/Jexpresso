# `CompEuler/ffs_step_M7` — Mach-7 forward-facing step

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "ffs_step_M7")
```

This is `CompEuler/ffs_step` with **M∞ = 7 instead of 3**, and nothing else
changed that is not forced by that one number.

## Why it exists

`rampCaoEtAl2021` on an unstretched grid goes non-finite at `t = 4.8e-7`
(step 481 at its `Δt = 1e-9`), in slot 1 — ρ — on every rank at once. Two
things are new in that deck at the same time, the Mach number *and* the ramp
geometry/grid, so the failure names neither of them.

`ffs_step` runs. Raising **only** the Mach number on `ffs_step` therefore
splits the question in two:

- **this deck reaches `tend`** → Mach 7 is not by itself the problem, and the
  ramp's grid, its boundary conditions and its leading-edge/corner treatment
  are what to look at next;
- **this deck dies** → the same failure is reproduced on a configuration whose
  mesh, boundary conditions, fluxes and DynSGS settings are *all* validated at
  Mach 3, which is a far smaller thing to debug — and the levers listed in the
  `:μ` block of `user_inputs.jl` get tried here rather than on the ramp.

Mach 7 is the rung just below the ramp's 7.7: close enough that a deck
surviving here makes the Mach number an unlikely sole culprit.

## What differs from `ffs_step`

Four values. Everything else is identical — mesh, boundary conditions, fluxes,
primitives, `:nop`, `:init_refine_lvl`, `:μ`, `:Pr`, `:dsgs_sensor`,
`:dsgs_hold_steps => 0`, `:dsgs_norms => "domain"`, `:energy_equation`.

| setting | `ffs_step` | here | why |
|---|---|---|---|
| `M∞` (`initialize.jl`) | 3.0 | **7.0** | the point of the case |
| `:Δt` | 1.0e-7 | **5.0e-8** | `\|u\|+c` doubles, 1371.6 → 2743.3 m/s, so the advective limit halves |
| `:tend` | 8.0e-3 | **3.5e-3** | the same tunnel flow-through count, 2.74 → 2.80 |
| `:diagnostics_at_times` step | 5.0e-5 | **2.5e-5** | the same output cadence in flow-through units |

`user_flux.jl`, `user_source.jl`, `user_primitives.jl` are byte-identical
copies; `user_bc.jl` and `initialize.jl` are identical in code and differ only
in comments and in the one `M∞` literal. `ffs_step_transfinite.{geo,msh}` are
copies too (the `.geo` differs only in its flow-condition comment), kept here
so the case is self-contained.

Keeping `:dsgs_hold_steps => 0` is not incidental. `ffs_step` needs it because
its initial condition is *not* smooth — a supersonic stream started impulsively
against the step — so holding ν at zero while the BDF2 history fills integrates
the most violent steps of the run with no dissipation at all. That argument is
stronger at Mach 7, not weaker.

## Free stream

Built in `ffs_freestream()` from `PhysConst` (γ = cp/cv = 1.398), so the stream
is exactly Mach 7 for the gas the solver integrates:

| | M = 3 | M = 7 |
|---|---|---|
| ρ∞ [kg/m³] | 1.20494 | 1.20494 |
| c∞ [m/s] | 342.91 | 342.91 |
| u∞ [m/s] | 1028.7 | **2400.4** |
| ρE∞ [J/m³] | 8.92e5 | **3.73e6** |
| stagnation T [K] | 818 | **3152** |
| normal-shock ρ₂/ρ₁ | 3.86 | **5.46** |
| flow-through, 3 m [s] | 2.92e-3 | 1.25e-3 |

The gas stays calorically perfect: this is the ideal-gas Euler system, no
dissociation and no vibrational excitation, so the run is a **numerical** test
at Mach 7 and not a physical model of Mach-7 air. (The ramp deck, by contrast,
runs `:lsutherland => true` because it has a real boundary layer to resolve;
this one is inviscid apart from the DynSGS shock capturing, exactly as
`ffs_step` is.)

## Why `:μ` was *not* raised

`ffs_step`'s own sweep — reproduced in `user_inputs.jl`, and **every row of it
is at Mach 3** — shows that more DynSGS dissipation helps only when `Δt` is cut
to match, because what binds on this configuration is the *viscous* step limit,
not the advective one. At Mach 7 both sides of `μΔt/(ρΔx²)` move on their own
and cancel:

- `μ_max = C_max·Δ·(|u|+c)`, the bound DynSGS saturates at the step corner,
  **doubles** with the wave speed;
- `Δt` **halves**.

So the viscous number lands where `ffs_step` measured it, with
`:μ => [1.0, 4.0, 4.0, 4.0]` carried over untouched. Raising `:μ` on top of the
`Δt` cut is the failure mode that sweep documents.

## Cost, and a cheap first look

`3.5e-3 / 5.0e-8` = **70 000 steps** — the same order as `ffs_step`'s 80 000,
at the same per-step cost.

The question this deck exists to answer does not need the full run: the step
corner is where `ffs_step` failed in *every* row of its sweep, and it fails
early or not at all. **`:tend => 2.0e-4` (4000 steps) already answers it.**

`JEXPRESSO_STEP_HEARTBEAT=1` turns on a per-step trace without editing the deck.

## What to look at

`:lschlieren => true` writes `schlieren` and `schlieren_grad_rho` into the VTU;
colour `schlieren` with a **reversed** greyscale in ParaView. Against the
familiar Mach-3 picture, expect the bow shock to stand closer to the step, the
shock layer to be thinner, and the roof reflection to strike further upstream.

## Debugging sweep

Every knob is an environment variable (defaults reproduce the deck, so an unset
environment is the baseline):

```bash
J=julia --project=. src/Jexpresso.jl CompEuler ffs_step_M7

JEXPRESSO_M7_TEND=1.0e-3 $J                               # baseline, short
JEXPRESSO_M7_SENSOR=residual JEXPRESSO_M7_TEND=1.0e-3 $J  # (a) element residual, not |∂ₜq|
JEXPRESSO_M7_FILTER=0.005   JEXPRESSO_M7_TEND=1.0e-3 $J   # (a') kill the top mode instead
JEXPRESSO_M7_NORMS=rank     JEXPRESSO_M7_TEND=1.0e-3 $J   # (b) shrink the normalising Ω
JEXPRESSO_M7_CMAX=2.0       JEXPRESSO_M7_TEND=1.0e-3 $J   # (c) raise the cap
JEXPRESSO_M7_MU1=4.0        JEXPRESSO_M7_TEND=1.0e-3 $J   # more β∇ρ density diffusion
```

| variable | default | notes |
|---|---|---|
| `JEXPRESSO_M7_SENSOR` | `legacy` | `residual` = stage-consistent element residual |
| `JEXPRESSO_M7_NORMS` | `domain` | only `domain`/`rank` exist for this kernel; `element` is DSGS_MHD only |
| `JEXPRESSO_M7_MU1` | 1.0 | slot 1, the β∇ρ density diffusion |
| `JEXPRESSO_M7_MU` | 4.0 | slots 2-4, applied *after* the cap |
| `JEXPRESSO_M7_CMAX` | 0.5 | `μ_cap = Cmax·Δ·ρ_max·(|u|+c)` |
| `JEXPRESSO_M7_CMIN` | 0.0 | unconditional floor `Cmin·Δ·ρ_max·(|u|+c)`; cell Re = 1/Cmin |
| `JEXPRESSO_M7_FILTER` | 0.0 | Boyd-Vandeven blend μ_x; see below |
| `JEXPRESSO_M7_DT` / `_TEND` / `_REF` | 5.0e-8 / 3.5e-3 / 0 | |

### What the filter actually is at `:nop => 4`

**A top-mode killer and nothing else.** Boyd-Vandeven only acts on modes
`k > 2n/3`, so at `n = 4` the transfer weights are `[1, 1, 1, 0.9957, 0]`:
modes 0-2 exactly untouched, mode 3 cut by 0.4% at *full* strength, mode 4
annihilated. `filter_type` `exp` and `quad` collapse to the same operator at
this order.

So `μ_x` is a **rate**, not an amplitude. `filter!` runs inside `rhs!`
(`rhs.jl:799`), i.e. once per RK stage, five times a step, and the top mode
decays as `(1-μ_x)` per call — e-folding in `1/(5·μ_x)` steps:

| μ_x | top mode e-folds in | left after 15 000 steps |
|---|---|---|
| 0.01 | 20 steps | 0 |
| 0.005 | 40 steps | 1e-164 |
| 0.001 | 200 steps | 3e-33 |
| 1e-4 | 2000 steps | 5e-4 |
| 8e-6 | 25 000 steps (one flow-through) | 0.55 |

There is no setting that is both gentle and useful: anything fast enough to
catch a Gibbs mode is a complete P4 → P3 truncation over the run, and anything
slow enough to leave the mode alive cannot catch it. `0.005` is the value to
try, and if it does not help, `0.05` will not either.

## Measured so far

| run | died at | vs baseline |
|---|---|---|
| baseline, coarse grid (`:linitial_refine => false`) | 7.5e-4 | — |
| `JEXPRESSO_M7_FILTER=0.005` | 1.1e-3 | 1.47× |

**The filter result is the informative one.** At `:nop => 4` that setting is a
*complete* removal of the top mode (e-fold 40 steps, nothing left after a few
hundred) — the sharpest available test of "this is undamped grid-scale
ringing". It bought 1.47× and the spurious structures downstream of the corner
survived it. So the defect does **not** live in the top mode: it is in modes
0–3, which the P4 Boyd–Vandeven filter leaves untouched (weights
`[1, 1, 1, 0.9957, 0]`). That is an element-scale structure, not a 2Δx mode.

### Why the DynSGS knobs each buy a factor and none of them fixes it

The trouble is in the **expansion**, and a residual sensor is structurally
quiet there *by design*: a Prandtl–Meyer fan is a smooth solution, so
`μ_res ≈ 0` in it — correctly. What the corner `(0.6, 0.2)` sheds is an
entropy layer, a contact-type feature that convects downstream and never
self-heals, and nothing in a shock-capturing sensor is aimed at it.

This is a documented property of this exact problem, not a Jexpresso defect.
Woodward & Colella (1984) §IV call the step corner *a singular point of the
flow* and apply a special fix to the cells next to it — resetting their
entropy and stagnation enthalpy to the upstream value — precisely because the
spurious entropy layer otherwise streams downstream and corrupts the Mach
stem. That is at Mach 3. At Mach 7 the expansion is far stronger and the
density it expands into far lower.

Three instruments actually target this, in increasing order of work:

1. **`:dsgs_Cmin`** — a background viscosity floor, unconditional, acting on
   every mode, scaling as `Δ·(|u|+c)`. The direct answer to "nothing damps the
   structures downwind". One line.
2. **De-singularize the geometry** — round the corner over 1–2 elements in
   `ffs_step_transfinite.geo`. A point singularity is the worst case for a
   high-order method; a small radius makes the fan non-centered. The ramp's
   leading edge is the same kind of point, so this transfers.
3. **A Woodward–Colella corner fix** — reset entropy and stagnation enthalpy in
   the nodes adjacent to the corner. `user_source!` receives `x, y`, so this is
   implementable in the case files with no kernel change.

## Status

**Not yet run.** The deck was written and checked by inspection against the
current `ffs_step`; no Julia was available in the environment it was written
in, so it has not been executed or even parsed. First run is the test.

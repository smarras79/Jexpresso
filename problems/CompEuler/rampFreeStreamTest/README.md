# rampFreeStreamTest

**This is not a flow. It is a one-property unit test of the discretisation.**

A uniform free stream is an exact solution of the Euler equations — every
flux is constant, so `∂F/∂x + ∂G/∂y = 0` *pointwise*. A correct
discretisation must return an RHS of **exactly zero** at every node, for
ever. This case starts from the uniform Mach-7.7 stream of
`rampCaoEtAl2021` on that case's own mesh, with everything else switched
off, and measures whether it stays there.

## Why it exists

Eight Mach-7 runs across three geometries (`ffs_step_M7`,
`ffs_step_M7_round`, `shock_circle_M7`, `rampCaoEtAl2021_M7`) have now
failed the same way, and **none of them ever measured this property.**

The trigger was the positivity report on `rampCaoEtAl2021_M7`, which named
the first negative pressure in the run:

```
first at (x, y) = (0.19659258262890678, 0.08585606812651451)
```

- `x = 0.1 + 0.1·cos15° = 0.19659258262890683` — the **outflow plane**, to
  sixteen digits.
- `y` is `2.5836e-5` below the top-right corner `(0.19659, 0.08588)` —
  **exactly one wall-normal spacing** of `ramp15_uniform.msh`.

So the first negative pressure was one node below the corner where the
free-stream `top` Dirichlet meets the `outflow` boundary that imposes
*nothing* — in undisturbed free stream, at step ~200, after 1.7 mm of flow
travel, with no shock within 20 cm of it. **Nothing physical happens
there.**

## Why it matters at Mach 7 and not at Mach 3

`p = (γ−1)(ρE − ρ|u|²/2)` is a difference of two nearly equal numbers. At
M = 7.7 the internal energy is only **5.7%** of the total energy
(`ρE = 33466`, `ke = 31558` in SI for this stream), so a relative error in
`ρE` is amplified into the pressure by

```
(γ−1)ρE/p = 1 + γ(γ−1)M²/2   =  3.5 at M = 3,  17.5 at M = 7.7
```

A free-stream-preservation error small enough for `ffs_step` at Mach 3 to
absorb is five times larger in the pressure here — and the pressure only
has 5.7% of headroom before it goes negative.

## The experiment

Two constants in `user_bc.jl`:

| | `FSP_OUTFLOW` | `FSP_WALL` |
|---|---|---|
| **default (run 1)** | `:nothing` — the production condition | `:freestream` — wall removed |
| control (run 2) | `:freestream` — every boundary node overwritten every stage | `:freestream` |

Run 1 is the hypothesis. Run 2 is the control.

## Reading the result

The output carries `dp = p − p∞` and `dp_rel = (p − p∞)/p∞`. Colour by
`dp_rel` in ParaView and **read its range off the Information tab.** That
range at the last frame is the entire result.

| `dp_rel` range | conclusion |
|---|---|
| ~1e-14, flat in time | Free stream is preserved. This line of enquiry is closed and the Mach-7 failures have another cause. |
| larger, growing, **confined to the outflow plane** | The scheme manufactures a source where nothing is imposed — a missing or wrong boundary term in the inviscid assembly. That is a **bug**, not a stabilisation shortfall, and it explains every failure in this campaign. |
| larger, growing, **spread over the interior** | The volume operator or the metrics are at fault and the outflow is innocent. |

If run 1 shows growth and run 2 (control) is machine zero, the outflow is
confirmed. If both grow, it is the volume operator.

## What is off, and why each one is load-bearing

- `:lvisc => false` — an artificial viscosity would *damp the error being
  measured*.
- `:lpositivity => false` — **must** stay false. The repair exists to hide
  negative pressure; here negative pressure is the signal.
- `:lfilter => false`, `:lsource => false` — same argument.
- `:lkep => false` — the baseline is the production inviscid operator. Set
  it `true` with `:volume_flux => ranocha()` for a third run: flux
  differencing carries its own free-stream-preservation requirement on the
  metric terms, and `user_flux.jl` already has the methods.

The mesh, `:nop`, the state and `:Δt` are the production ones unchanged, so
a result here transfers directly to `rampCaoEtAl2021_M7`.

## Running

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "rampFreeStreamTest")
```

500 steps of Δt = 1e-9, ten output frames. Seconds of wall clock, not hours.
